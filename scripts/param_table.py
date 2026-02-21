#!/usr/bin/env python3
"""Generate supplementary parameter table from module config files.

Parses all modules/*/config.toml files to extract parameter names, values,
descriptions, and literature references. Outputs Markdown and CSV formats.

Usage:
    python3 scripts/param_table.py                    # Markdown to stdout
    python3 scripts/param_table.py --csv              # CSV to stdout
    python3 scripts/param_table.py --out table.md     # Write to file
"""

import glob
import os
import re
import sys


def parse_toml_params(path):
    """Parse a module config.toml file into parameter records.

    Returns list of dicts: {module, section, param, value, description, ref}
    """
    records = []
    module = os.path.basename(os.path.dirname(path))

    section = ""
    section_comment = ""
    with open(path) as f:
        for line in f:
            stripped = line.strip()

            # Skip empty lines
            if not stripped:
                section_comment = ""
                continue

            # Section comment (# heading above parameters)
            if stripped.startswith("# ") and "=" not in stripped:
                section_comment = stripped[2:].strip()
                continue

            # Table header
            if stripped.startswith("["):
                section = stripped.strip("[]")
                continue

            # Skip pure comment lines that are continuations of previous
            if stripped.startswith("#"):
                continue

            # Parameter line: name = value  # description (ref, doi:...)
            match = re.match(r'^(\w+)\s*=\s*(.+?)(?:\s*#\s*(.*))?$', stripped)
            if not match:
                continue

            param = match.group(1)
            value = match.group(2).strip()
            comment = match.group(3) or ""

            # Extract DOI from comment
            doi_match = re.search(r'doi:(\S+)', comment)
            doi = doi_match.group(1).rstrip(")") if doi_match else ""

            # Extract short reference from comment (Author Year pattern)
            ref_match = re.search(
                r'([A-Z][a-z]+(?:\s+(?:&|et al\.?)\s+[A-Z][a-z]+)?\s+\d{4})',
                comment
            )
            ref = ref_match.group(1) if ref_match else ""

            # Clean description: strip ref/DOI citation from end
            desc = comment
            if doi_match:
                desc = desc[:doi_match.start()]
            # Remove trailing "(Author et al. Year)" or "(Author Year)"
            desc = re.sub(r'\s*\([A-Z][^)]*\d{4}[^)]*\)\s*$', '', desc)
            # Remove trailing "(Author et al. Year" (unclosed paren)
            desc = re.sub(r'\s*\([A-Z][^)]*\d{4}[^)]*$', '', desc)
            # Remove trailing "; Author Year" chains
            desc = re.sub(r';\s*[A-Z][a-z].*$', '', desc)
            desc = desc.strip().rstrip("(; ,)")
            # Remove unbalanced trailing open parens
            if desc.count("(") > desc.count(")"):
                desc = desc[:desc.rfind("(")].strip().rstrip("; ,")

            records.append({
                "module": module,
                "section": section_comment or section,
                "param": param,
                "value": value,
                "description": desc,
                "ref": ref,
                "doi": doi,
            })

    return records


def format_markdown(all_records):
    """Format records as a Markdown table grouped by module."""
    lines = []
    lines.append("# Supplementary Table: Model Parameters")
    lines.append("")
    lines.append("| Module | Parameter | Value | Description | Reference |")
    lines.append("|--------|-----------|-------|-------------|-----------|")

    current_module = ""
    for r in all_records:
        mod = r["module"] if r["module"] != current_module else ""
        current_module = r["module"]

        ref = r["ref"]
        if r["doi"]:
            ref = f"[{ref}](https://doi.org/{r['doi']})" if ref else ""

        lines.append(
            f"| {mod} | `{r['param']}` | {r['value']} "
            f"| {r['description']} | {ref} |"
        )

    return "\n".join(lines) + "\n"


def format_csv(all_records):
    """Format records as CSV."""
    lines = []
    lines.append("module,parameter,value,description,reference,doi")
    for r in all_records:
        desc = r["description"].replace('"', '""')
        ref = r["ref"].replace('"', '""')
        lines.append(
            f'"{r["module"]}","{r["param"]}","{r["value"]}",'
            f'"{desc}","{ref}","{r["doi"]}"'
        )
    return "\n".join(lines) + "\n"


def main():
    project_root = os.path.join(os.path.dirname(__file__), os.pardir)
    config_pattern = os.path.join(project_root, "modules", "*", "config.toml")
    config_files = sorted(glob.glob(config_pattern))

    if not config_files:
        print("No module config files found.", file=sys.stderr)
        sys.exit(1)

    all_records = []
    for path in config_files:
        records = parse_toml_params(path)
        all_records.extend(records)

    # Parse CLI args
    fmt = "markdown"
    out_path = None
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--csv":
            fmt = "csv"
        elif args[i] == "--out" and i + 1 < len(args):
            out_path = args[i + 1]
            i += 1
        i += 1

    if fmt == "csv":
        output = format_csv(all_records)
    else:
        output = format_markdown(all_records)

    if out_path:
        with open(out_path, "w") as f:
            f.write(output)
        print(f"Wrote {len(all_records)} parameters to {out_path}")
    else:
        print(output)

    print(f"Total: {len(all_records)} parameters from {len(config_files)} modules",
          file=sys.stderr)


if __name__ == "__main__":
    main()
