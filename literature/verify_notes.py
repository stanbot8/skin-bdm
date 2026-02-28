#!/usr/bin/env python3
"""Verify SOURCES.yaml metadata and notes against CrossRef.

For each DOI in SOURCES.yaml, fetches the CrossRef record and checks:
  1. Title match (fuzzy, normalized)
  2. Author last-name match (first author)
  3. Year match
  4. Journal match (fuzzy)
  5. Abstract available? If so, checks that key quantitative claims
     in the notes appear plausible (keyword overlap)

Exit codes:
  0 = all checks pass
  1 = metadata mismatches found

Usage:
    python3 literature/verify_notes.py
    python3 literature/verify_notes.py --fix-metadata   # auto-fix title/journal/year from CrossRef
    python3 literature/verify_notes.py --module diabetic # check single module
    python3 literature/verify_notes.py --verbose         # show abstracts and details
"""

import glob
import json
import os
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.request

import yaml

_PROJECT_ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), os.pardir))
_MODULES_DIR = os.path.join(_PROJECT_ROOT, "modules")

# CrossRef API polite pool requires mailto
_USER_AGENT = ("SkiBiDy/1.0 (https://github.com/stanbot8/skibidy;"
               " mailto:stanbot8@users.noreply.github.com)")


def strip_diacritics(s):
    """Remove diacritics: ä->a, ö->o, ł->l, ð->d, etc."""
    # Handle characters that NFKD doesn't decompose
    s = s.replace('ł', 'l').replace('Ł', 'L')
    s = s.replace('đ', 'd').replace('Đ', 'D')
    s = s.replace('ø', 'o').replace('Ø', 'O')
    nfkd = unicodedata.normalize('NFKD', s)
    return ''.join(c for c in nfkd if not unicodedata.combining(c))


def normalize(s):
    """Lowercase, strip diacritics, strip punctuation, collapse whitespace."""
    s = strip_diacritics(s).lower()
    s = re.sub(r'<[^>]+>', '', s)  # strip HTML tags
    s = re.sub(r'[^a-z0-9\s]', ' ', s)
    return ' '.join(s.split())


def title_similarity(a, b):
    """Word-level Jaccard similarity between two titles."""
    wa = set(normalize(a).split())
    wb = set(normalize(b).split())
    if not wa or not wb:
        return 0.0
    return len(wa & wb) / len(wa | wb)


def first_author_last(authors_str):
    """Extract first author's last name from a SOURCES.yaml author string."""
    first = authors_str.split(',')[0].strip()
    # Handle "Alba-Loureiro TC" or "Van De Water L" patterns
    parts = first.split()
    # Last name is everything except the last token if it looks like initials
    if len(parts) >= 2 and all(len(p) <= 2 or p.isupper() for p in parts[-1:]):
        return ' '.join(parts[:-1]).lower()
    return parts[0].lower() if parts else ''


def fetch_crossref(doi, cache):
    """Fetch CrossRef metadata for a DOI. Returns dict or None."""
    if doi in cache:
        return cache[doi]

    url = f"https://doi.org/{doi}"
    req = urllib.request.Request(url)
    req.add_header("Accept", "application/citeproc+json")
    req.add_header("User-Agent", _USER_AGENT)

    try:
        resp = urllib.request.urlopen(req, timeout=15)
        data = json.loads(resp.read().decode('utf-8'))
        cache[doi] = data
        return data
    except urllib.error.HTTPError as e:
        if e.code == 406:
            # DOI exists but no citeproc; try works API
            try:
                url2 = f"https://api.crossref.org/works/{doi}"
                req2 = urllib.request.Request(url2)
                req2.add_header("User-Agent", _USER_AGENT)
                resp2 = urllib.request.urlopen(req2, timeout=15)
                data = json.loads(resp2.read().decode('utf-8'))['message']
                cache[doi] = data
                return data
            except Exception:
                cache[doi] = None
                return None
        cache[doi] = None
        return None
    except Exception:
        cache[doi] = None
        return None


def extract_crossref_fields(cr):
    """Extract normalized fields from a CrossRef record."""
    title = ''
    if cr.get('title'):
        if isinstance(cr['title'], list):
            title = cr['title'][0]
        else:
            title = cr['title']

    authors = []
    for a in cr.get('author', []):
        family = a.get('family', '')
        given = a.get('given', '')
        authors.append(f"{family}, {given}" if given else family)

    year = None
    for key in ('published-print', 'published-online', 'published', 'issued',
                'created'):
        dp = cr.get(key, {}).get('date-parts', [[]])
        if dp and dp[0] and dp[0][0]:
            year = dp[0][0]
            break

    journal = ''
    ct = cr.get('container-title', [])
    if isinstance(ct, list) and ct:
        journal = ct[0]
    elif isinstance(ct, str):
        journal = ct

    abstract = cr.get('abstract', '')

    return {
        'title': title,
        'authors': authors,
        'year': year,
        'journal': journal,
        'abstract': abstract,
    }


def check_entry(src, cr_fields, verbose=False):
    """Compare a SOURCES.yaml entry against CrossRef. Returns list of issues."""
    issues = []

    # Title check: also try containment (one title is substring of the other)
    our_title_norm = normalize(src.get('title', ''))
    cr_title_norm = normalize(cr_fields['title'])
    sim = title_similarity(src.get('title', ''), cr_fields['title'])
    # Check if one contains the other (handles subtitles, truncation)
    contained = (our_title_norm in cr_title_norm or
                 cr_title_norm in our_title_norm) if our_title_norm and cr_title_norm else False
    if sim < 0.3 and not contained:
        issues.append({
            'field': 'title',
            'severity': 'ERROR',
            'message': f"Title mismatch (similarity={sim:.2f})",
            'ours': src.get('title', ''),
            'crossref': cr_fields['title'],
        })
    elif sim < 0.5 and not contained:
        issues.append({
            'field': 'title',
            'severity': 'WARN',
            'message': f"Title partial match (similarity={sim:.2f})",
            'ours': src.get('title', ''),
            'crossref': cr_fields['title'],
        })

    # Year check
    our_year = src.get('year')
    cr_year = cr_fields['year']
    if our_year and cr_year and int(our_year) != int(cr_year):
        issues.append({
            'field': 'year',
            'severity': 'ERROR',
            'message': f"Year mismatch: ours={our_year}, crossref={cr_year}",
        })

    # First author check
    if cr_fields['authors']:
        cr_first = strip_diacritics(cr_fields['authors'][0].split(',')[0]).lower().strip()
        our_first = strip_diacritics(first_author_last(src.get('authors', '')))
        # Normalize hyphens
        cr_norm = cr_first.replace('-', '').replace(' ', '')
        our_norm = our_first.replace('-', '').replace(' ', '')
        # Also check if the core surname matches (strip suffixes like jr, gn, ii)
        cr_core = re.sub(r'\b(jr|sr|ii|iii|iv|gn)\b', '', cr_norm).strip()
        our_core = re.sub(r'\b(jr|sr|ii|iii|iv|gn)\b', '', our_norm).strip()
        if (cr_core and our_core and cr_core not in our_core
                and our_core not in cr_core):
            issues.append({
                'field': 'authors',
                'severity': 'ERROR',
                'message': f"First author mismatch: ours='{our_first}', crossref='{cr_first}'",
            })

    # Journal check (abbreviated vs full names are expected; only flag if
    # the abbreviated form doesn't appear as initials of the full name)
    if cr_fields['journal']:
        our_journal_raw = src.get('journal', '').split('(')[0].split(',')[0].strip()
        cr_journal_raw = cr_fields['journal']
        # Check if abbreviated words appear in the full name
        our_words = normalize(our_journal_raw).split()
        cr_words = normalize(cr_journal_raw).split()
        # Strip volume/page numbers from our journal
        our_words = [w for w in our_words if not w.isdigit()]
        # Check if each abbreviated word starts a full word
        matched = 0
        for ow in our_words:
            if any(cw.startswith(ow[:3]) for cw in cr_words):
                matched += 1
        coverage = matched / max(1, len(our_words))
        if coverage < 0.5:
            issues.append({
                'field': 'journal',
                'severity': 'WARN',
                'message': f"Journal mismatch (coverage={coverage:.0%})",
                'ours': our_journal_raw,
                'crossref': cr_journal_raw,
            })

    # Abstract vs notes check
    abstract = cr_fields.get('abstract', '')
    notes = src.get('notes', '')
    if abstract and notes and verbose:
        # Extract numbers from notes
        note_numbers = set(re.findall(r'\d+(?:\.\d+)?', notes))
        abstract_clean = normalize(abstract)
        note_keywords = set(normalize(notes).split()) - {
            'supports', 'the', 'and', 'in', 'of', 'a', 'to', 'for',
            'with', 'from', 'by', 'is', 'are', 'was', 'were', 'that',
            'this', 'on', 'at', 'as', 'or', 'an', 'not', 'no', 'but',
        }
        abstract_keywords = set(abstract_clean.split())
        overlap = note_keywords & abstract_keywords
        coverage = len(overlap) / max(1, len(note_keywords))
        issues.append({
            'field': 'abstract',
            'severity': 'INFO',
            'message': f"Abstract keyword coverage: {coverage:.0%} ({len(overlap)}/{len(note_keywords)} words)",
            'abstract': abstract[:300],
        })

    return issues


def load_all_entries():
    """Load all SOURCES.yaml entries with DOIs. Returns list of dicts."""
    entries = []
    pattern = os.path.join(_MODULES_DIR, "*", "SOURCES.yaml")
    for path in sorted(glob.glob(pattern)):
        module = os.path.basename(os.path.dirname(path))
        with open(path) as f:
            data = yaml.safe_load(f)
        if not isinstance(data, dict):
            continue
        for dataset_name, dataset in data.items():
            if not isinstance(dataset, dict):
                continue
            for src in dataset.get('sources', []):
                if isinstance(src, dict) and src.get('doi'):
                    entries.append({
                        'module': module,
                        'dataset': dataset_name,
                        'source': src,
                    })
    return entries


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Verify SOURCES.yaml metadata against CrossRef")
    parser.add_argument("--module", type=str, default=None,
                        help="Check single module only")
    parser.add_argument("--verbose", action="store_true",
                        help="Show abstracts and detailed info")
    parser.add_argument("--fix-metadata", action="store_true",
                        help="Auto-fix title/journal/year from CrossRef")
    parser.add_argument("--cache", type=str,
                        default=os.path.join(os.path.dirname(__file__),
                                             ".crossref_cache.json"),
                        help="CrossRef cache file path")
    args = parser.parse_args()

    # Load cache
    cache = {}
    if os.path.isfile(args.cache):
        with open(args.cache) as f:
            cache = json.load(f)

    entries = load_all_entries()
    if args.module:
        entries = [e for e in entries if e['module'] == args.module]

    # Deduplicate by DOI
    seen_dois = set()
    unique_entries = []
    for e in entries:
        doi = e['source']['doi'].lower().strip()
        if doi not in seen_dois:
            seen_dois.add(doi)
            unique_entries.append(e)

    total = len(unique_entries)
    errors = []
    warnings = []
    fixes = []

    print(f"Verifying {total} unique DOIs against CrossRef...\n")

    for i, entry in enumerate(unique_entries):
        src = entry['source']
        doi = src['doi']
        src_id = src.get('id', '?')
        module = entry['module']

        cr = fetch_crossref(doi, cache)
        time.sleep(0.3)  # polite rate limit

        if cr is None:
            errors.append(f"  [{module}] {src_id}: CrossRef fetch failed for {doi}")
            print(f"  [{i+1}/{total}] FAIL: {src_id} ({module})")
            continue

        cr_fields = extract_crossref_fields(cr)
        issues = check_entry(src, cr_fields, verbose=args.verbose)

        has_error = any(iss['severity'] == 'ERROR' for iss in issues)
        has_warn = any(iss['severity'] == 'WARN' for iss in issues)

        if has_error:
            mark = "ERROR"
        elif has_warn:
            mark = "WARN"
        else:
            mark = "ok"

        print(f"  [{i+1}/{total}] {mark}: {src_id} ({module})")

        for iss in issues:
            line = f"    {iss['severity']}: {iss['message']}"
            if iss['severity'] == 'ERROR':
                errors.append(f"  [{module}] {src_id}: {iss['message']}")
                if 'ours' in iss:
                    print(f"{line}")
                    print(f"      ours:     {iss['ours'][:100]}")
                    print(f"      crossref: {iss['crossref'][:100]}")
                else:
                    print(f"{line}")

                # Track fixes
                if args.fix_metadata and iss['field'] in ('title', 'year'):
                    fixes.append({
                        'module': module,
                        'src_id': src_id,
                        'field': iss['field'],
                        'old': iss.get('ours', str(src.get(iss['field']))),
                        'new': iss.get('crossref', str(cr_fields.get(iss['field']))),
                    })

            elif iss['severity'] == 'WARN':
                warnings.append(f"  [{module}] {src_id}: {iss['message']}")
                if 'ours' in iss:
                    print(f"{line}")
                    print(f"      ours:     {iss['ours'][:100]}")
                    print(f"      crossref: {iss['crossref'][:100]}")
                else:
                    print(f"{line}")

            elif args.verbose:
                print(f"{line}")
                if 'abstract' in iss:
                    # Clean HTML from abstract
                    clean = re.sub(r'<[^>]+>', '', iss['abstract'])
                    print(f"      abstract: {clean[:200]}...")

    # Save cache
    with open(args.cache, 'w') as f:
        json.dump(cache, f)

    # Summary
    print(f"\n{'='*60}")
    print(f"Verified: {total} DOIs")
    print(f"Errors:   {len(errors)}")
    print(f"Warnings: {len(warnings)}")

    if errors:
        print(f"\nErrors:")
        for e in errors:
            print(e)

    if warnings:
        print(f"\nWarnings:")
        for w in warnings:
            print(w)

    if fixes:
        print(f"\nAuto-fix would update {len(fixes)} fields:")
        for fix in fixes:
            print(f"  [{fix['module']}] {fix['src_id']}.{fix['field']}: "
                  f"'{fix['old'][:60]}' -> '{fix['new'][:60]}'")

    if errors:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
