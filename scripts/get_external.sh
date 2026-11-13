#!/bin/bash
# Fetch the external Reburn modules at the pinned hashes (merge plan section 2.1; Reburn discussion #789).
# Usage: scripts/get_external.sh   (from the repo root). HPC nodes without GitHub: run on a login node or unpack a vendored tarball.
set -e
cd "$(dirname "$0")/.."
FW_HASH=31051b0; BU_HASH=6d56e29
for spec in "Fireweed $FW_HASH" "Burnup $BU_HASH"; do
  set -- $spec
  if [ ! -d external/$1/.git ]; then git clone -q https://github.com/JoshuaRady/$1.git external/$1; fi
  git -C external/$1 fetch -q origin
  git -C external/$1 checkout -q $2
  echo "external/$1 at $(git -C external/$1 rev-parse --short HEAD)"
done
