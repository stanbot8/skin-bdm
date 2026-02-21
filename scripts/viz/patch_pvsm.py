#!/usr/bin/env python3
"""Patch a ParaView state file with settings defined in paraview.toml.

Reads [colormap.*], [opacity], [camera], [glyph], and [rendering]
from paraview.toml and applies them to the .pvsm XML file generated
by BioDynaMo.  Project-agnostic: all visual configuration lives in
the TOML file.

Uses only stdlib (no tomllib/tomli dependency).

Usage:
    python3 scripts/viz/patch_pvsm.py [path/to/file.pvsm] [path/to/paraview.toml]
"""

import ast
import re
import sys
import xml.etree.ElementTree as ET

# ParaView ColorSpace enum values
COLOR_SPACE_MAP = {
    "RGB": "0",
    "HSV": "1",
    "Lab": "2",
    "Diverging": "3",
    "Step": "5",
}


# ---------------------------------------------------------------------------
# Minimal TOML reader (stdlib only, handles the subset we need)
# ---------------------------------------------------------------------------

def _parse_value(raw):
    """Parse a TOML value string into a Python object."""
    raw = raw.strip()
    if raw.lower() == "true":
        return True
    if raw.lower() == "false":
        return False
    # Strip inline comment (not inside quotes or brackets)
    depth = 0
    in_str = False
    for i, ch in enumerate(raw):
        if ch == '"' and (i == 0 or raw[i - 1] != '\\'):
            in_str = not in_str
        elif not in_str:
            if ch in '([':
                depth += 1
            elif ch in ')]':
                depth -= 1
            elif ch == '#' and depth == 0:
                raw = raw[:i].rstrip()
                break
    # Try ast.literal_eval (handles numbers, strings, lists)
    try:
        return ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        # Bare string (unquoted TOML string -- shouldn't happen in our config)
        return raw


def _set_nested(d, keys, value):
    """Set d[k1][k2][...] = value, creating dicts as needed."""
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


def parse_toml_sections(path):
    """Parse a TOML file into nested dicts.

    Handles [dotted.sections], key = value, and comments.
    Skips [[array_of_tables]] (BioDynaMo's own config sections).
    """
    data = {}
    section_keys = []

    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue

            # [[array_of_tables]] -- skip (BioDynaMo handles these)
            if stripped.startswith('[['):
                section_keys = []
                continue

            # [section.path]
            m = re.match(r'^\[([^\]]+)\]', stripped)
            if m:
                section_keys = m.group(1).strip().split('.')
                # Ensure section exists
                d = data
                for k in section_keys:
                    d = d.setdefault(k, {})
                continue

            # key = value
            m = re.match(r'^(\w+)\s*=\s*(.+)$', stripped)
            if m:
                key = m.group(1)
                val = _parse_value(m.group(2))
                full_keys = section_keys + [key]
                _set_nested(data, full_keys, val)

    return data


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_colormaps(vis_cfg):
    """Extract named colormaps from visualization.colormap.* sections.

    Returns list of (name, color_space, points) where points is
    [(value, r, g, b), ...] sorted by value.
    """
    cm_section = vis_cfg.get("colormap", {})
    if not isinstance(cm_section, dict):
        return []

    result = []
    for name, entries in cm_section.items():
        if not isinstance(entries, dict):
            continue
        color_space = entries.get("color_space", "Step")
        # Collect cN = [value, R, G, B] entries
        points = []
        for key, val in sorted(entries.items()):
            if key.startswith('c') and key[1:].isdigit() and isinstance(val, list):
                points.append(tuple(val))
        if points:
            result.append((name, color_space, points))

    return result


# ---------------------------------------------------------------------------
# PVSM patching (XML manipulation)
# ---------------------------------------------------------------------------

def flatten_rgb_points(points):
    """Flatten [(val, r, g, b), ...] into [val, r, g, b, val, r, g, b, ...]."""
    flat = []
    for v, r, g, b in points:
        flat.extend([v, r, g, b])
    return flat


def patch_lut(proxy_elem, points, color_space="Step"):
    """Replace RGBPoints and ColorSpace in a PVLookupTable proxy element."""
    flat = flatten_rgb_points(points)
    cs_val = COLOR_SPACE_MAP.get(color_space, "5")

    for prop in proxy_elem.findall("Property"):
        name = prop.get("name")

        if name == "RGBPoints":
            for child in list(prop):
                if child.tag == "Element":
                    prop.remove(child)
            prop.set("number_of_elements", str(len(flat)))
            for i, val in enumerate(flat):
                el = ET.SubElement(prop, "Element")
                el.set("index", str(i))
                el.set("value", str(val))

        elif name == "ColorSpace":
            for child in prop.findall("Element"):
                child.set("value", cs_val)

        elif name == "NumberOfTableValues":
            for child in prop.findall("Element"):
                child.set("value", str(len(points)))

        elif name == "Discretize":
            for child in prop.findall("Element"):
                child.set("value", "1")

        elif name == "AutomaticRescaleRangeMode":
            for child in prop.findall("Element"):
                child.set("value", "-1")

        elif name == "ScalarRangeInitialized":
            for child in prop.findall("Element"):
                child.set("value", "1")


def patch_opacity(root, opacity_points):
    """Set opacity ramp on all PiecewiseFunction proxies."""
    flat = []
    for val, opac in opacity_points:
        flat.extend([val, opac, 0.5, 0.0])

    patched = 0
    for proxy in root.iter("Proxy"):
        if proxy.get("type") != "PiecewiseFunction":
            continue
        for prop in proxy.findall("Property"):
            if prop.get("name") != "Points":
                continue
            for child in list(prop):
                if child.tag == "Element":
                    prop.remove(child)
            prop.set("number_of_elements", str(len(flat)))
            for i, val in enumerate(flat):
                el = ET.SubElement(prop, "Element")
                el.set("index", str(i))
                el.set("value", str(val))
            patched += 1
    return patched


def patch_agent_representation(root, lut_id):
    """Set agent representations that have stratum_ to use the given LUT."""
    patched = 0
    for proxy in root.iter("Proxy"):
        for prop in proxy.findall("Property"):
            if prop.get("name") != "ColorArrayName":
                continue
            domain = prop.find("Domain")
            if domain is None:
                continue
            has_stratum = any(
                s.get("text") == "stratum_" for s in domain.findall("String")
            )
            if not has_stratum:
                continue

            for el in prop.findall("Element"):
                if int(el.get("index", "-1")) == 4:
                    el.set("value", "stratum_")

            for lut_prop in proxy.findall("Property"):
                if lut_prop.get("name") != "LookupTable":
                    continue
                lut_prop.set("number_of_elements", "1")
                for child in list(lut_prop):
                    if child.tag == "Proxy":
                        lut_prop.remove(child)
                proxy_ref = ET.Element("Proxy")
                proxy_ref.set("value", str(lut_id))
                lut_prop.insert(0, proxy_ref)
                break

            patched += 1
            break
    return patched


def patch_agent_solid_colors(root, agent_colors):
    """Set solid (uniform) color on agent representation proxies.

    agent_colors: dict mapping agent type substring to [R, G, B].
    Matches against RegistrationName of GeometryRepresentation proxies.
    """
    patched = 0
    for proxy in root.iter("Proxy"):
        if proxy.get("type") != "GeometryRepresentation":
            continue
        reg_name = ""
        for prop in proxy.findall("Property"):
            if prop.get("name") == "RegistrationName":
                for el in prop.findall("Element"):
                    reg_name = el.get("value", "")
        if not reg_name:
            continue

        for agent_name, rgb in agent_colors.items():
            if agent_name not in reg_name:
                continue
            for prop in proxy.findall("Property"):
                if prop.get("name") == "DiffuseColor":
                    for el in prop.findall("Element"):
                        idx = int(el.get("index", "-1"))
                        if 0 <= idx < 3:
                            el.set("value", str(rgb[idx]))
                    patched += 1
                    break
            break
    return patched


def patch_glyph_resolution(root, glyph_cfg):
    """Set sphere glyph resolution on all SphereSource proxies."""
    theta = glyph_cfg.get("theta_resolution", 8)
    phi = glyph_cfg.get("phi_resolution", 8)

    patched = 0
    for proxy in root.iter("Proxy"):
        if proxy.get("type") != "SphereSource":
            continue
        for prop in proxy.findall("Property"):
            name = prop.get("name")
            if name == "ThetaResolution":
                for el in prop.findall("Element"):
                    el.set("value", str(theta))
            elif name == "PhiResolution":
                for el in prop.findall("Element"):
                    el.set("value", str(phi))
        patched += 1
    return patched


def patch_diffusion_representation(root, rep_value):
    """Set representation mode on all UniformGridRepresentation proxies."""
    patched = 0
    for proxy in root.iter("Proxy"):
        if proxy.get("type") != "UniformGridRepresentation":
            continue
        for prop in proxy.findall("Property"):
            if prop.get("name") == "Representation":
                for el in prop.findall("Element"):
                    el.set("value", rep_value)
                patched += 1
                break
    return patched


def patch_camera(root, cam_cfg):
    """Set camera position in the RenderView proxy."""
    position = cam_cfg.get("position", [-40, -80, 30])
    focal_point = cam_cfg.get("focal_point", [15, 15, 15])
    view_up = cam_cfg.get("view_up", [0, 0, 1])
    view_angle = cam_cfg.get("view_angle", 30)

    for proxy in root.iter("Proxy"):
        if proxy.get("type") != "RenderView":
            continue

        camera_props = {
            "CameraPosition": position,
            "CameraPositionInfo": position,
            "CameraFocalPoint": focal_point,
            "CameraFocalPointInfo": focal_point,
            "CameraViewUp": view_up,
            "CameraViewUpInfo": view_up,
        }

        for prop in proxy.findall("Property"):
            name = prop.get("name")
            if name in camera_props:
                values = camera_props[name]
                for el in prop.findall("Element"):
                    idx = int(el.get("index", "-1"))
                    if 0 <= idx < len(values):
                        el.set("value", str(values[idx]))
            elif name == "CameraViewAngle":
                for el in prop.findall("Element"):
                    el.set("value", str(view_angle))
            elif name == "Camera3DManipulators":
                for el in prop.findall("Element"):
                    if el.get("value") == "4":
                        el.set("value", "5")
        break


def get_lut_registration_name(proxy):
    """Extract RegistrationName from a PVLookupTable proxy."""
    for prop in proxy.findall("Property"):
        if prop.get("name") == "RegistrationName":
            for el in prop.findall("Element"):
                return el.get("value", "")
    return ""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    pvsm_path = sys.argv[1] if len(sys.argv) > 1 else "output/skibidy/skibidy.pvsm"
    toml_path = sys.argv[2] if len(sys.argv) > 2 else "paraview.toml"

    cfg = parse_toml_sections(toml_path)
    colormaps = load_colormaps(cfg)

    if not colormaps:
        print(f"No [colormap.*] defined in {toml_path}, skipping color patch.")
        return

    tree = ET.parse(pvsm_path)
    root = tree.getroot()

    # First colormap is the default
    default_name, default_cs, default_pts = colormaps[0]

    # Patch lookup tables
    patched_luts = 0
    first_lut_id = None
    for proxy in root.iter("Proxy"):
        if proxy.get("type") != "PVLookupTable":
            continue

        reg_name = get_lut_registration_name(proxy)

        # Find matching colormap by name substring
        matched = None
        for cm_name, cm_cs, cm_pts in colormaps:
            stem = cm_name.rstrip("s")
            if stem in reg_name:
                matched = (cm_cs, cm_pts)
                break

        cs, pts = matched or (default_cs, default_pts)
        patch_lut(proxy, pts, cs)

        if first_lut_id is None:
            first_lut_id = proxy.get("id")
        patched_luts += 1

    if patched_luts == 0:
        print("WARNING: No PVLookupTable proxies found in PVSM file.")
        sys.exit(1)

    # Opacity
    opacity_cfg = cfg.get("opacity", {})
    if opacity_cfg.get("points"):
        patch_opacity(root, opacity_cfg["points"])

    # Agent representation (stratum-based LUT)
    if first_lut_id:
        patch_agent_representation(root, first_lut_id)

    # Agent solid colors
    agent_color_cfg = cfg.get("agent_color", {})
    if agent_color_cfg:
        ac = {name: vals["color"] for name, vals in agent_color_cfg.items()
              if isinstance(vals, dict) and "color" in vals}
        if ac:
            n_ac = patch_agent_solid_colors(root, ac)
            if n_ac:
                print(f"  Agent colors: {n_ac} patched")

    # Glyph resolution
    glyph_cfg = cfg.get("glyph")
    if glyph_cfg:
        patch_glyph_resolution(root, glyph_cfg)

    # Diffusion grid representation
    render_cfg = cfg.get("rendering", {})
    diff_rep = render_cfg.get("diffusion_representation")
    if diff_rep:
        patch_diffusion_representation(root, diff_rep)

    # Camera
    camera_cfg = cfg.get("camera")
    if camera_cfg:
        patch_camera(root, camera_cfg)

    tree.write(pvsm_path, xml_declaration=False)

    # Summary
    print(f"Patched {patched_luts} lookup table(s) in {pvsm_path}")
    for cm_name, cm_cs, cm_pts in colormaps:
        print(f"  {cm_name}: {len(cm_pts)} color(s), space={cm_cs}")
    if glyph_cfg:
        print(f"  Glyph: theta={glyph_cfg.get('theta_resolution', 8)} phi={glyph_cfg.get('phi_resolution', 8)}")
    if diff_rep:
        print(f"  Diffusion: {diff_rep}")
    if camera_cfg:
        print(f"  Camera: pos={camera_cfg.get('position')} angle={camera_cfg.get('view_angle')}")


if __name__ == "__main__":
    main()
