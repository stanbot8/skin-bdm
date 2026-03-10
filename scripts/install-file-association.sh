#!/usr/bin/env bash
# Register .skibidy file association so double-clicking opens the dashboard.
# Run once after cloning the repo.

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DASHBOARD="$SCRIPT_DIR/dashboard.py"

# 1. Register MIME type for .skibidy files
MIME_DIR="$HOME/.local/share/mime"
mkdir -p "$MIME_DIR/packages"
cat > "$MIME_DIR/packages/skibidy.xml" << 'MIMEEOF'
<?xml version="1.0" encoding="UTF-8"?>
<mime-info xmlns="http://www.freedesktop.org/standards/shared-mime-info">
  <mime-type type="application/x-skibidy">
    <comment>SkiBiDy Study Project</comment>
    <glob pattern="*.skibidy"/>
  </mime-type>
</mime-info>
MIMEEOF
update-mime-database "$MIME_DIR" 2>/dev/null || true

# 2. Install desktop entry with correct path
APPS_DIR="$HOME/.local/share/applications"
mkdir -p "$APPS_DIR"
sed "s|SCRIPT_PATH|$SCRIPT_DIR|g" "$SCRIPT_DIR/skibidy-open.desktop" > "$APPS_DIR/skibidy-open.desktop"

# 3. Associate .skibidy files with the desktop entry
xdg-mime default skibidy-open.desktop application/x-skibidy 2>/dev/null || true

echo "Done. Double-click any .skibidy file to open it in the dashboard."
