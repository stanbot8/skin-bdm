#!/bin/bash
if [ -z "$BDMSYS" ]; then
  echo "ERROR: BioDynaMo not sourced. Run: source <path>/bin/thisbdm.sh"
  exit 1
fi
cd "$(dirname "$0")/.."
rm -rf output/N3bdm*
biodynamo build && ./build/skibidy-test
