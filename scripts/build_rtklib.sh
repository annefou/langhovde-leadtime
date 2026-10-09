#!/usr/bin/env bash
# Builds RTKLIB command-line tools (convbin, rnx2rtkp) from a pinned release of the
# rtklibexplorer fork into .tools/rtklib/bin. RTKLIB is not packaged on conda-forge.
set -euo pipefail
TAG=v2.5.1
SHA256=ad8e56f64dd2b71e2eb8098c27925a65b027b161eae366cf5afb71402348d2e1
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DEST=$ROOT/.tools/rtklib
mkdir -p "$DEST/bin" "$DEST/src"
TGZ=$DEST/src/RTKLIB-$TAG.tar.gz
[ -f "$TGZ" ] || curl -sSL -o "$TGZ" "https://github.com/rtklibexplorer/RTKLIB/archive/refs/tags/$TAG.tar.gz"
GOT=$(sha256sum "$TGZ" | cut -d' ' -f1)
echo "RTKLIB $TAG sha256 $GOT"
if [ "$GOT" != "$SHA256" ]; then echo "checksum mismatch" >&2; exit 1; fi
tar -xzf "$TGZ" -C "$DEST/src"
SRC=$(ls -d "$DEST"/src/RTKLIB-*/ | head -1)
for app in convbin rnx2rtkp; do
  make -s -C "$SRC/app/consapp/$app/gcc" LDLIBS=-lm >/dev/null 2>&1
  cp "$SRC/app/consapp/$app/gcc/$app" "$DEST/bin/"
done
"$DEST/bin/convbin" 2>&1 | head -2 || true
