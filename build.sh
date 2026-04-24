#!/usr/bin/env bash
# build.sh  –  builds seqmatcher Docker image for the native host architecture.
#
# Usage:
#   ./build.sh                        # build only
#   ./build.sh --push --registry=ghcr.io/your-org   # build + push
#
# The mambaorg/micromamba base image is multi-arch; bioconda resolves the
# correct package binaries (arm64 / amd64) at solve time automatically.
# No BUILDPLATFORM / TARGETPLATFORM args are needed.

set -euo pipefail

IMAGE="seqmatcher:latest"
REGISTRY=""
PUSH=false

for arg in "$@"; do
    case $arg in
        --push)         PUSH=true ;;
        --image=*)      IMAGE="${arg#*=}" ;;
        --registry=*)   REGISTRY="${arg#*=}" ;;
    esac
done

ARCH="$(uname -m)"
case "$ARCH" in
    arm64|aarch64) PLATFORM="linux/arm64" ;;
    x86_64|amd64)  PLATFORM="linux/amd64" ;;
    *)
        echo "Unsupported architecture: $ARCH" >&2; exit 1 ;;
esac

echo ">>> Host arch : $ARCH  →  $PLATFORM"
echo ">>> Image tag : $IMAGE"
echo ""

docker build \
    --platform "$PLATFORM" \
    -t "$IMAGE" \
    .

echo ""
echo "✓ Build complete → $IMAGE ($PLATFORM)"

# ── Smoke-test ────────────────────────────────────────────────────────────────
echo ""
echo ">>> Smoke-test …"
docker run --rm --platform "$PLATFORM" "$IMAGE" bash -c "
    python - << 'PYEOF'
import sys, pyBigWig, pyfaidx, logomaker, Bio, numpy, pandas
print('Python   :', sys.version.split()[0])
print('biopython:', Bio.__version__)
print('pyBigWig :', pyBigWig.__version__)
print('pyfaidx  :', pyfaidx.__version__)
print('logomaker:', logomaker.__version__)
print('numpy    :', numpy.__version__)
print('pandas   :', pandas.__version__)
print('All Python imports OK.')
PYEOF
    echo -n 'bigWigInfo: '; bigWigInfo 2>&1 | head -1 || true
"

# ── Optional push ─────────────────────────────────────────────────────────────
if $PUSH && [ -n "$REGISTRY" ]; then
    FULL_TAG="${REGISTRY}/${IMAGE}"
    docker tag  "$IMAGE" "$FULL_TAG"
    docker push "$FULL_TAG"
    echo "✓ Pushed → $FULL_TAG"
fi
