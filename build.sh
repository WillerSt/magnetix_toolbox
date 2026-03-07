#!/bin/bash
# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

# Build script for magnetix_toolbox
# Usage: ./build.sh [options]
# Options:
#   -j <num>          Number of parallel jobs (default: auto-detect)
#   -c, --clean       Clean build directory before building
#   -v, --verbose     Verbose output
#   -h, --help        Show this help message

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
JOBS=$(nproc)
CLEAN=0
VERBOSE=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -j)
            JOBS="$2"
            shift 2
            ;;
        -c|--clean)
            CLEAN=1
            shift
            ;;
        -v|--verbose)
            VERBOSE="--verbose"
            shift
            ;;
        -h|--help)
            sed -n '2,/^$/p' "$0" | tail -n +1 | sed 's/# //'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  magnetix_toolbox Build Script${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check prerequisites
echo -e "${YELLOW}Checking prerequisites...${NC}"
command -v cmake &> /dev/null || { echo -e "${RED}Error: cmake not found${NC}"; exit 1; }
command -v make &> /dev/null || { echo -e "${RED}Error: make not found${NC}"; exit 1; }
echo -e "${GREEN}✓ Prerequisites OK${NC}"
echo ""

# Ensure dependencies are available
echo -e "${YELLOW}Checking for FEniCSx and dependencies...${NC}"
echo -e "${YELLOW}(These should be in your environment or CMAKE_PREFIX_PATH)${NC}"
if command -v module &> /dev/null; then
    echo -e "${YELLOW}Tip: If not available, try: module load dolfinx${NC}"
fi
echo ""

# Create or clean build directory
BUILD_DIR="${SCRIPT_DIR}/build"
if [ $CLEAN -eq 1 ]; then
    echo -e "${YELLOW}Cleaning build directories...${NC}"
    rm -rf "$BUILD_DIR"
    rm -rf "${SCRIPT_DIR}/install"
    echo -e "${GREEN}✓ Build directories cleaned${NC}"
fi

if [ ! -d "$BUILD_DIR" ]; then
    echo -e "${YELLOW}Creating build directory...${NC}"
    mkdir -p "$BUILD_DIR"
    echo -e "${GREEN}✓ Build directory created${NC}"
fi

# Build hysteresis_model
echo ""
echo -e "${YELLOW}Step 1/3: Building hysteresis_model...${NC}"
HYST_BUILD_DIR="${SCRIPT_DIR}/build/hysteresis_model_build"
mkdir -p "$HYST_BUILD_DIR"
cd "$HYST_BUILD_DIR"
cmake ${VERBOSE} -DCMAKE_INSTALL_PREFIX="${SCRIPT_DIR}/install" "${SCRIPT_DIR}/hysteresis_model" >/dev/null
cmake --build . -j "${JOBS}" ${VERBOSE}
cmake --build . --target install ${VERBOSE} >/dev/null
echo -e "${GREEN}✓ hysteresis_model built successfully${NC}"

# Prepare CMAKE_PREFIX_PATH for next steps
export CMAKE_PREFIX_PATH="${SCRIPT_DIR}/install:${CMAKE_PREFIX_PATH}"

# Build fenicsx_tools
echo ""
echo -e "${YELLOW}Step 2/3: Building fenicsx_tools...${NC}"
FENICS_BUILD_DIR="${SCRIPT_DIR}/build/fenicsx_tools_build"
mkdir -p "$FENICS_BUILD_DIR"
cd "$FENICS_BUILD_DIR"
cmake ${VERBOSE} -DCMAKE_INSTALL_PREFIX="${SCRIPT_DIR}/install" -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" "${SCRIPT_DIR}/fenicsx_tools/library" >/dev/null
cmake --build . -j "${JOBS}" ${VERBOSE}
cmake --build . --target install ${VERBOSE} >/dev/null
echo -e "${GREEN}✓ fenicsx_tools built successfully${NC}"

# Build examples
echo ""
echo -e "${YELLOW}Step 3/3: Building examples...${NC}"
EXAMPLES_BUILD_DIR="${SCRIPT_DIR}/build/examples_build"
mkdir -p "$EXAMPLES_BUILD_DIR"
cd "$EXAMPLES_BUILD_DIR"
cmake ${VERBOSE} -DCMAKE_INSTALL_PREFIX="${SCRIPT_DIR}/install" -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" "${SCRIPT_DIR}/examples/magnetostatic_2D" >/dev/null
cmake --build . -j "${JOBS}" ${VERBOSE}
echo -e "${GREEN}✓ examples built successfully${NC}"

# Move examples executable
echo ""
echo -e "${YELLOW}Organizing build artifacts...${NC}"
if [ -f "${EXAMPLES_BUILD_DIR}/magnetostatic_2D_exec" ]; then
    cp "${EXAMPLES_BUILD_DIR}/magnetostatic_2D_exec" "${SCRIPT_DIR}/examples/magnetostatic_2D/"
    echo -e "${GREEN}✓ Executable placed at examples/magnetostatic_2D/magnetostatic_2D_exec${NC}"
fi

echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}Build successful!${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Setup environment instructions
INSTALL_PREFIX="${SCRIPT_DIR}/install"
echo -e "${YELLOW}To setup your environment, add to ~/.bashrc or ~/.zshrc:${NC}"
echo ""
echo -e "${GREEN}export CMAKE_PREFIX_PATH=\"${INSTALL_PREFIX}/dpc_hysteresis-0.1:${INSTALL_PREFIX}/fenicsx_magnetics_toolbox-0.10:\$CMAKE_PREFIX_PATH\"${NC}"
echo -e "${GREEN}export PYTHONPATH=\"${INSTALL_PREFIX}/fenicsx_magnetics_toolbox-0.10/python:\$PYTHONPATH\"${NC}"
echo ""

echo -e "${YELLOW}Or run to update your current shell:${NC}"
echo ""
echo -e "${GREEN}export CMAKE_PREFIX_PATH=\"${INSTALL_PREFIX}/dpc_hysteresis-0.1:${INSTALL_PREFIX}/fenicsx_magnetics_toolbox-0.10:\$CMAKE_PREFIX_PATH\"${NC}"
echo -e "${GREEN}export PYTHONPATH=\"${INSTALL_PREFIX}/fenicsx_magnetics_toolbox-0.10/python:\$PYTHONPATH\"${NC}"
echo ""

echo -e "${YELLOW}Run TEAM Problem 32 example:${NC}"
echo ""
echo -e "${GREEN}cd examples/magnetostatic_2D/TEAM_Problem_32${NC}"
echo -e "${GREEN}mpirun -np 4 ../magnetostatic_2D_exec --scen TeamProblem32_case3_default.xml${NC}"
echo ""
