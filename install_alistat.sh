#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# install_alistat.sh — Download and build AliStat from source
# https://github.com/thomaskf/AliStat
#
# Platform: macOS / Linux only
# Usage:
#   ./install_alistat.sh              # build & install to ./tools/alistat
#   ./install_alistat.sh /usr/local/bin   # install to custom directory
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

INSTALL_DIR="${1:-$(cd "$(dirname "$0")" && pwd)/tools}"
BUILD_DIR="$(mktemp -d)"

echo "=== AliStat Installer ==="
echo "  Build directory : $BUILD_DIR"
echo "  Install target  : $INSTALL_DIR/alistat"
echo ""

# ── Check prerequisites ─────────────────────────────────────────────────────
for cmd in wget tar make g++ ; do
    if ! command -v "$cmd" &>/dev/null; then
        # On macOS, g++ is provided by Xcode command line tools (clang++)
        if [[ "$cmd" == "g++" ]] && command -v c++ &>/dev/null; then
            continue
        fi
        echo "ERROR: '$cmd' not found. Please install it first."
        exit 1
    fi
done

# ── Download ─────────────────────────────────────────────────────────────────
echo "Downloading AliStat source..."
wget -q --show-progress \
    "https://github.com/thomaskf/AliStat/archive/refs/heads/master.tar.gz" \
    -O "$BUILD_DIR/AliStat.tar.gz"

# ── Extract & build ─────────────────────────────────────────────────────────
echo "Extracting..."
tar -zxf "$BUILD_DIR/AliStat.tar.gz" -C "$BUILD_DIR"

echo "Building..."
cd "$BUILD_DIR/AliStat-master"
make -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)"

# ── Install ──────────────────────────────────────────────────────────────────
mkdir -p "$INSTALL_DIR"
cp alistat "$INSTALL_DIR/alistat"
chmod +x "$INSTALL_DIR/alistat"

# ── Clean up ─────────────────────────────────────────────────────────────────
rm -rf "$BUILD_DIR"

echo ""
echo "✔ AliStat installed to: $INSTALL_DIR/alistat"
echo ""
echo "  To use with the pipeline, either:"
echo "    1. Add '$INSTALL_DIR' to your \$PATH, or"
echo "    2. Set the path in config/config.yaml:"
echo "       tools:"
echo "         alistat: \"$INSTALL_DIR/alistat\""
echo ""
