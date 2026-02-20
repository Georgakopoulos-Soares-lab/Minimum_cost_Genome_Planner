#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# install_sdsl.sh — clone and build sdsl-lite v2.1.1
# ─────────────────────────────────────────────────────────────────────────────
# Usage:
#   bash install_sdsl.sh [PREFIX]
#
# PREFIX defaults to ./deps/sdsl (relative to this script).
# After installation, compile the planners with:
#   make               (if PREFIX is ./deps/sdsl  — auto-detected)
#   make SDSL_PREFIX=<PREFIX>   (if you chose a custom PREFIX)
#
# Notes:
#   - Requires git, cmake, gcc/g++.
#   - Pins v2.1.1 which fixes the louds_tree.hpp member-name bug present
#     in some older checkouts (m_select1 → m_bv_select1).
#   - Uses -DBUILD_PORTABLE=ON so the resulting library is not tied to
#     the CPU architecture used during compilation (safe for ARM and x86).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PREFIX="${1:-${SCRIPT_DIR}/deps/sdsl}"
VENDOR_DIR="${SCRIPT_DIR}/vendor"
BUILD_DIR="${SCRIPT_DIR}/.build_sdsl_tmp"

SDSL_VER="v2.1.1"
SDSL_GIT_URL="https://github.com/simongog/sdsl-lite.git"
SDSL_SRC="${VENDOR_DIR}/sdsl-lite"

echo "Installing sdsl-lite ${SDSL_VER} → ${PREFIX}"
mkdir -p "${PREFIX}" "${VENDOR_DIR}"

# Fast path: already installed
if [[ -f "${PREFIX}/include/sdsl/csa_wt.hpp" && -f "${PREFIX}/lib/libsdsl.a" ]]; then
    echo "SDSL already installed in ${PREFIX} — skipping."
    exit 0
fi

# Clone if not already present
if [[ ! -d "${SDSL_SRC}/.git" ]]; then
    echo "Cloning sdsl-lite ${SDSL_VER}..."
    rm -rf "${SDSL_SRC}"
    git clone --recursive "${SDSL_GIT_URL}" "${SDSL_SRC}"
fi

cd "${SDSL_SRC}"
git fetch --tags --force
git checkout "${SDSL_VER}"
git submodule update --init --recursive

# Force GCC/G++ to avoid NVHPC ABI / CMake path issues
CC_BIN="$(command -v gcc  || true)"
CXX_BIN="$(command -v g++ || true)"
if [[ -z "${CC_BIN}" || -z "${CXX_BIN}" ]]; then
    echo "ERROR: gcc/g++ not found in PATH; cannot build sdsl-lite." >&2
    exit 2
fi
export CC="${CC_BIN}"
export CXX="${CXX_BIN}"

echo "Configuring CMake..."
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"
cmake -S "${SDSL_SRC}" -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="${CC}" \
    -DCMAKE_CXX_COMPILER="${CXX}" \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DBUILD_PORTABLE=ON

echo "Building..."
cmake --build "${BUILD_DIR}" --target sdsl
cmake --build "${BUILD_DIR}" --target install

rm -rf "${BUILD_DIR}"

echo ""
echo "Done.  SDSL installed to: ${PREFIX}"
if [[ "${PREFIX}" == "${SCRIPT_DIR}/deps/sdsl" ]]; then
    echo "Now run:  make"
else
    echo "Now run:  make SDSL_PREFIX=${PREFIX}"
fi
