#!/bin/bash
# Assemble PNG frames from output/animation/ into an animated GIF.
#
# Workflow:
#   1. Run simulation (./run.sh)
#   2. In ParaView: File > Save Animation > output/animation/frame.png
#   3. Run this script: ./scripts/viz/render_gif.sh
#
# Usage:
#     ./scripts/viz/render_gif.sh [fps]
#
# Requires: ffmpeg
# Output:   output/skibidy.gif

set -eo pipefail

FPS="${1:-6}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
FRAME_DIR="$PROJECT_DIR/output/animation"
GIF_OUT="$PROJECT_DIR/output/skibidy.gif"

command -v ffmpeg >/dev/null || { echo "Error: ffmpeg not found (apt install ffmpeg)"; exit 1; }

# Find PNGs (ParaView names them frame.0000.png or frame_0000.png, etc.)
N_FRAMES=$(find "$FRAME_DIR" -maxdepth 1 -name "*.png" 2>/dev/null | wc -l)
if [ "$N_FRAMES" -eq 0 ]; then
    echo "Error: No PNGs found in $FRAME_DIR/"
    echo ""
    echo "Export frames from ParaView first:"
    echo "  File > Save Animation > $FRAME_DIR/frame.png"
    exit 1
fi

echo "=== Assembling GIF ($N_FRAMES frames, ${FPS} fps) ==="

# Build sorted file list (handles any naming convention)
ffmpeg -y -framerate "$FPS" \
    -pattern_type glob -i "$FRAME_DIR/*.png" \
    -vf "fps=$FPS,scale=640:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=128:stats_mode=diff[p];[s1][p]paletteuse=dither=bayer:bayer_scale=3" \
    "$GIF_OUT" 2>&1 | tail -3

echo ""
echo "Done: $GIF_OUT ($(du -h "$GIF_OUT" | cut -f1))"
