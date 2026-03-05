#!/bin/bash
# Build BioDynaMo from the stanbot8 fork (dev/stanbot8 branch).
# After running this, source the fork build before running skibidy:
#   source ~/biodynamo/build/bin/thisbdm.sh
set -e

BDM_DIR="$HOME/biodynamo"
FORK_REMOTE="_stanbot8"
FORK_BRANCH="dev/stanbot8"

cd "$BDM_DIR"

# Ensure fork remote exists
if ! git remote get-url "$FORK_REMOTE" &>/dev/null; then
  echo "Adding fork remote..."
  git remote add "$FORK_REMOTE" https://github.com/stanbot8/biodynamo.git
fi

git fetch "$FORK_REMOTE" "$FORK_BRANCH"

# Checkout the fork branch
git checkout "$FORK_BRANCH" 2>/dev/null || git checkout -b "$FORK_BRANCH" "$FORK_REMOTE/$FORK_BRANCH"
git reset --hard "$FORK_REMOTE/$FORK_BRANCH"

echo "Building BioDynaMo from $FORK_REMOTE/$FORK_BRANCH..."
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
make -j"$(nproc)" 2>&1 | tail -5
if [ ${PIPESTATUS[0]} -ne 0 ]; then
  echo "BDM build failed."
  exit 1
fi

echo ""
echo "Done. Source the build before running skibidy:"
echo "  source ~/biodynamo/build/bin/thisbdm.sh"
