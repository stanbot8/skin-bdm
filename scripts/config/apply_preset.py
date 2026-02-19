#!/usr/bin/env python3
"""Apply a preset/profile overlay to bdm.toml.

Reads the overlay TOML file and replaces matching key = value lines in
bdm.toml.  Only keys present in the overlay are touched; everything else
(comments, formatting, other sections) is preserved.

Supports dotted section names like [skin.wound], [skin.diabetic], etc.
Keys are matched by (section, key) tuple so 'enabled' under [skin.wound]
does not collide with 'enabled' under [skin.tumor].

Usage:
    python3 scripts/config/apply_preset.py presets/wound.toml
    python3 scripts/config/apply_preset.py profiles/diabetic.toml bdm.toml
"""

import re
import sys


def parse_overrides(preset_path):
    """Extract (section, key) -> value pairs from a preset file."""
    overrides = {}
    current_section = ""
    with open(preset_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Skip array-of-tables
            if line.startswith("[["):
                continue
            # Match section headers: [skin], [skin.wound], etc.
            section_match = re.match(r"^\[([\w.]+)\]$", line)
            if section_match:
                current_section = section_match.group(1)
                continue
            m = re.match(r"^(\w+)\s*=\s*(.+)$", line)
            if m and current_section.startswith("skin"):
                overrides[(current_section, m.group(1))] = m.group(2)
    return overrides


def apply_overrides(toml_path, overrides):
    """Replace matching lines in skin sections of the TOML file."""
    with open(toml_path) as f:
        lines = f.readlines()

    applied = set()
    current_section = ""
    for i, line in enumerate(lines):
        stripped = line.lstrip()
        # Track section headers (skip array-of-tables like [[visualize_diffusion]])
        if stripped.startswith("[["):
            continue
        section_match = re.match(r"^\[([\w.]+)\]", stripped)
        if section_match:
            current_section = section_match.group(1)
            continue
        # Only replace keys inside skin sections
        if not current_section.startswith("skin"):
            continue
        if stripped.startswith("#") or "=" not in stripped:
            continue
        m = re.match(r"^(\s*)(\w+)\s*=\s*(.+)$", line)
        if m and (current_section, m.group(2)) in overrides:
            key = m.group(2)
            section_key = (current_section, key)
            indent = m.group(1)
            old_val = m.group(3)
            new_val = overrides[section_key]
            # Preserve old inline comment only if new value has none
            comment = ""
            new_has_comment = re.search(r"\s+#\s+", new_val)
            if not new_has_comment:
                comment_match = re.search(r"\s+#\s+", old_val)
                if comment_match:
                    comment = old_val[comment_match.start():]
            lines[i] = f"{indent}{key} = {new_val}{comment}\n"
            applied.add(section_key)

    # Append missing keys to the correct section (or create section)
    missing = set(overrides.keys()) - applied
    if missing:
        # Find the last line index for each section
        section_ends = {}
        current = ""
        for i, line in enumerate(lines):
            stripped = line.lstrip()
            if stripped.startswith("[["):
                continue
            sm = re.match(r"^\[([\w.]+)\]", stripped)
            if sm:
                current = sm.group(1)
            if current:
                section_ends[current] = i

        # Group missing keys by section
        by_section = {}
        for section, key in sorted(missing):
            by_section.setdefault(section, []).append(key)

        appended = set()
        for section, keys in by_section.items():
            insert_lines = []
            if section not in section_ends:
                # Section doesn't exist; create it
                insert_lines.append(f"\n[{section}]\n")
            for key in keys:
                val = overrides[(section, key)]
                insert_lines.append(f"{key} = {val}\n")
                appended.add((section, key))

            if section in section_ends:
                idx = section_ends[section] + 1
            else:
                idx = len(lines)
            for j, il in enumerate(insert_lines):
                lines.insert(idx + j, il)
            # Shift section_ends for any sections after the insertion point
            shift = len(insert_lines)
            for s in section_ends:
                if section_ends[s] >= idx:
                    section_ends[s] += shift

        applied |= appended

    with open(toml_path, "w") as f:
        f.writelines(lines)

    return applied


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/config/apply_preset.py <overlay.toml> [bdm.toml]")
        sys.exit(1)

    preset_path = sys.argv[1]
    toml_path = sys.argv[2] if len(sys.argv) > 2 else "bdm.toml"

    overrides = parse_overrides(preset_path)
    if not overrides:
        print(f"No overrides found in {preset_path}")
        sys.exit(1)

    applied = apply_overrides(toml_path, overrides)
    applied_strs = [f"{s}.{k}" for s, k in sorted(applied)]
    print(f"Applied {preset_path}: {', '.join(applied_strs)}")


if __name__ == "__main__":
    main()
