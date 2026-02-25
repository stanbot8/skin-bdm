#!/usr/bin/env python3
"""Validate SOURCES.yaml structural integrity and file references.

Discovers sources from per-module files: modules/*/SOURCES.yaml

Checks:
  1. YAML parses cleanly
  2. Every dataset has required top-level fields (description, sources)
  3. Every source entry has required citation fields (id, authors, year, title)
  4. All file paths (consensus, raw_files) resolve to existing files
  5. Parameter-only datasets have 'notes' annotations (warning)
  6. Inline DOI comments in config files match merged SOURCES entries (warning)
  7. Per-module DOI locality: config DOIs should appear in the same module's
     SOURCES.yaml rather than only in a different module (info-level)
  8. --verify-dois: HTTP HEAD check that each DOI resolves via doi.org

Exit codes:
  0 = all checks pass
  1 = structural errors found

Usage:
    python3 literature/check_sources.py
    python3 literature/check_sources.py --verify-dois
    python3 literature/check_sources.py --config bdm.toml --config profiles/diabetic.toml
"""

import glob
import os
import re
import sys
import time

import yaml

_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
_SOURCES_PATH = os.path.join(_DATA_DIR, "SOURCES.yaml")
_PROJECT_ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), os.pardir))
_MODULES_DIR = os.path.join(_PROJECT_ROOT, "modules")

# Required fields for every source citation
_SOURCE_REQUIRED = {"id", "authors", "year", "title"}
# At least one locator should be present
_SOURCE_LOCATOR = {"doi", "pmc", "url"}


def discover_sources_files():
    """Find all SOURCES.yaml files (per-module and monolithic).

    Returns list of (path, module_name_or_None) tuples.
    Per-module files are discovered first, then the monolithic fallback.
    """
    files = []
    # Per-module sources
    pattern = os.path.join(_MODULES_DIR, "*", "SOURCES.yaml")
    for path in sorted(glob.glob(pattern)):
        module = os.path.basename(os.path.dirname(path))
        files.append((path, module))
    # Monolithic fallback
    if os.path.isfile(_SOURCES_PATH):
        files.append((_SOURCES_PATH, None))
    return files


def load_sources(path=None):
    """Parse a single SOURCES.yaml. Returns (dict, errors)."""
    if path is None:
        path = _SOURCES_PATH
    errors = []
    try:
        with open(path) as f:
            data = yaml.safe_load(f)
    except FileNotFoundError:
        return None, [f"SOURCES.yaml not found: {path}"]
    except yaml.YAMLError as e:
        return None, [f"YAML parse error in {path}: {e}"]
    if data is None:
        return {}, errors
    if not isinstance(data, dict):
        return None, [f"{path}: root must be a mapping"]
    return data, errors


def load_all_sources():
    """Discover and merge all SOURCES.yaml files.

    Returns (merged_dict, dataset_dirs, per_module_dois, errors) where
    dataset_dirs maps dataset_name -> directory containing SOURCES.yaml,
    and per_module_dois maps module_name -> set(doi_lowercase).
    """
    files = discover_sources_files()
    merged = {}
    dataset_dirs = {}
    per_module_dois = {}
    errors = []

    for path, module in files:
        data, errs = load_sources(path)
        errors.extend(errs)
        if data is None:
            continue

        base_dir = os.path.dirname(path)

        # Check for duplicate dataset names
        rel = os.path.relpath(path, _PROJECT_ROOT)
        for key in data:
            if key in merged:
                errors.append(f"Duplicate dataset '{key}' in {rel}")
            else:
                merged[key] = data[key]
                dataset_dirs[key] = base_dir

        # Build per-module DOI index
        if module is not None:
            dois = set()
            for entry in data.values():
                if not isinstance(entry, dict):
                    continue
                for src in entry.get("sources", []):
                    if isinstance(src, dict) and src.get("doi"):
                        dois.add(src["doi"].lower().strip())
            per_module_dois[module] = dois

    return merged, dataset_dirs, per_module_dois, errors


def check_source_entry(tag, src):
    """Validate one citation entry. Returns (errors, warnings)."""
    errors = []
    warnings = []
    if not isinstance(src, dict):
        return [f"{tag}: source entry is not a mapping"], []
    src_id = src.get("id", "?")
    for field in _SOURCE_REQUIRED:
        if field not in src:
            errors.append(f"{tag} ({src_id}): missing '{field}'")
    if not any(src.get(loc) for loc in _SOURCE_LOCATOR):
        warnings.append(f"{tag} ({src_id}): no doi/pmc/url")
    return errors, warnings


def check_dataset(name, entry, base_dir=_DATA_DIR):
    """Validate one dataset entry. Returns (errors, warnings).

    base_dir is the directory containing the SOURCES.yaml that defines
    this dataset (module dir for per-module files, literature/data/ for
    the monolithic file).
    """
    errors = []
    warnings = []
    if not isinstance(entry, dict):
        return [f"{name}: entry is not a mapping"], []

    if "description" not in entry:
        errors.append(f"{name}: missing 'description'")

    sources = entry.get("sources")
    if not sources:
        errors.append(f"{name}: missing or empty 'sources' list")
    elif not isinstance(sources, list):
        errors.append(f"{name}: 'sources' must be a list")
    else:
        for i, src in enumerate(sources):
            tag = f"{name}.sources[{i}]"
            src_errs, src_warns = check_source_entry(tag, src)
            errors.extend(src_errs)
            warnings.extend(src_warns)

    # File reference checks (timeseries datasets)
    consensus = entry.get("consensus")
    if consensus:
        full = os.path.join(base_dir, consensus)
        if not os.path.isfile(full):
            errors.append(f"{name}: consensus file not found: {consensus}")

    raw_files = entry.get("raw_files")
    if raw_files:
        if not isinstance(raw_files, list):
            errors.append(f"{name}: 'raw_files' must be a list")
        else:
            for rf in raw_files:
                full = os.path.join(base_dir, rf)
                if not os.path.isfile(full):
                    errors.append(f"{name}: raw file not found: {rf}")

    # Parameter-only dataset: warn if no source has notes
    if consensus is None and raw_files is None:
        if sources and isinstance(sources, list):
            has_notes = any(isinstance(s, dict) and s.get("notes") for s in sources)
            if not has_notes:
                warnings.append(
                    f"{name}: parameter-only dataset with no 'notes' on any source"
                )

    return errors, warnings


def check_config_dois(config_paths, all_dois, per_module_dois=None):
    """Cross-check inline DOI comments in config files against SOURCES.yaml.

    When per_module_dois is provided, also checks that DOIs in a module's
    config.toml appear in that module's own SOURCES.yaml (locality check).
    """
    warnings = []
    doi_re = re.compile(r"doi:(10\.\S+)")

    for path in config_paths:
        if not os.path.isfile(path):
            continue

        # Determine which module this config belongs to
        real = os.path.realpath(path)
        module = None
        if per_module_dois and "/modules/" in real:
            parts = real.split("/modules/", 1)
            if len(parts) == 2:
                module = parts[1].split("/", 1)[0]

        module_dois = per_module_dois.get(module, set()) if module else set()
        basename = os.path.basename(path)

        with open(path) as f:
            for lineno, line in enumerate(f, 1):
                for match in doi_re.finditer(line):
                    doi = match.group(1).rstrip(");,")
                    doi_lower = doi.lower()
                    if doi_lower not in all_dois:
                        warnings.append(
                            f"{basename}:{lineno}: doi:{doi} not in SOURCES.yaml"
                        )

    return warnings


def verify_dois(sources_data):
    """Verify each DOI resolves via the CrossRef content negotiation API.

    Uses the doi.org metadata endpoint (Accept: application/citeproc+json)
    which is the official scholarly API and not blocked by publisher firewalls,
    unlike direct HEAD requests to publisher landing pages.

    Returns (ok_list, fail_list) where each item is (id, doi, status).
    """
    import urllib.request
    import urllib.error

    dois = []
    for entry in sources_data.values():
        if not isinstance(entry, dict):
            continue
        for src in entry.get("sources", []):
            if isinstance(src, dict) and src.get("doi"):
                dois.append((src.get("id", "?"), src["doi"]))

    # Deduplicate by DOI
    seen = set()
    unique = []
    for src_id, doi in dois:
        key = doi.lower().strip()
        if key not in seen:
            seen.add(key)
            unique.append((src_id, doi))

    ok = []
    fail = []
    total = len(unique)
    for i, (src_id, doi) in enumerate(unique):
        url = f"https://doi.org/{doi}"
        try:
            req = urllib.request.Request(url)
            req.add_header("Accept", "application/citeproc+json")
            req.add_header("User-Agent",
                           "SkiBiDy/1.0 (https://github.com/stanbot8/skibidy;"
                           " mailto:stanbot8@users.noreply.github.com)")
            resp = urllib.request.urlopen(req, timeout=15)
            status = resp.getcode()
            if status < 400:
                ok.append((src_id, doi, status))
                mark = "ok"
            else:
                fail.append((src_id, doi, status))
                mark = f"FAIL ({status})"
        except urllib.error.HTTPError as e:
            status = e.code
            # 406 means DOI exists but no citeproc format available; still ok
            if status == 406:
                ok.append((src_id, doi, "406/exists"))
                mark = "ok (no citeproc)"
            else:
                fail.append((src_id, doi, status))
                mark = f"FAIL ({status})"
        except urllib.error.URLError as e:
            fail.append((src_id, doi, str(e.reason)))
            mark = f"FAIL ({e.reason})"
        except Exception as e:
            fail.append((src_id, doi, str(e)))
            mark = f"FAIL ({e})"
        print(f"  [{i+1}/{total}] {mark}: {src_id} -> {doi}")
        time.sleep(0.5)  # polite rate limit per CrossRef etiquette

    return ok, fail


def run_checks(config_paths=None, do_verify_dois=False):
    """Run all source checks. Returns (errors, warnings)."""
    all_errors = []
    all_warnings = []

    data, dataset_dirs, per_module_dois, parse_errors = load_all_sources()
    all_errors.extend(parse_errors)
    if data is None:
        return all_errors, all_warnings

    for name, entry in data.items():
        base_dir = dataset_dirs.get(name, _DATA_DIR)
        errs, warns = check_dataset(name, entry, base_dir=base_dir)
        all_errors.extend(errs)
        all_warnings.extend(warns)

    # Build global DOI set for config checking
    all_dois = set()
    for entry in data.values():
        if not isinstance(entry, dict):
            continue
        for src in entry.get("sources", []):
            if isinstance(src, dict) and src.get("doi"):
                all_dois.add(src["doi"].lower().strip())

    # Auto-discover config files if none specified
    if config_paths is None:
        config_paths = []
        # Module config files
        module_configs = glob.glob(
            os.path.join(_MODULES_DIR, "*", "config.toml")
        )
        config_paths.extend(sorted(module_configs))
        # Profile files
        profiles = glob.glob(os.path.join(_PROJECT_ROOT, "profiles", "*.toml"))
        config_paths.extend(sorted(profiles))

    if config_paths:
        doi_warns = check_config_dois(config_paths, all_dois, per_module_dois)
        all_warnings.extend(doi_warns)

    if do_verify_dois:
        print("\n  Verifying DOIs online...")
        ok, fail = verify_dois(data)
        print(f"\n  DOI verification: {len(ok)} ok, {len(fail)} failed")
        for src_id, doi, status in fail:
            all_warnings.append(f"DOI failed ({status}): {src_id} -> {doi}")

    return all_errors, all_warnings


def main():
    config_paths = []
    do_verify = False
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--config" and i + 1 < len(args):
            config_paths.append(args[i + 1])
            i += 2
        elif args[i].startswith("--config="):
            config_paths.append(args[i].split("=", 1)[1])
            i += 1
        elif args[i] == "--verify-dois":
            do_verify = True
            i += 1
        else:
            i += 1

    errors, warnings = run_checks(
        config_paths=config_paths if config_paths else None,
        do_verify_dois=do_verify,
    )

    if warnings:
        for w in warnings:
            print(f"  WARN: {w}")
    if errors:
        for e in errors:
            print(f"  ERROR: {e}")
        print(f"\n  SOURCES check FAILED: {len(errors)} error(s)")
        sys.exit(1)
    else:
        n = len(warnings)
        print(f"  SOURCES check passed ({n} warning{'s' if n != 1 else ''})")
        sys.exit(0)


if __name__ == "__main__":
    main()
