#!/usr/bin/env bash

set -e

clean=0
defs=()

for arg in "$@"; do
    case "$arg" in
        clean)
            clean=1
            ;;
        --real=float)
            defs+=("-DREAL_FLOAT")
            ;;
        --real=double)
            defs+=("-DREAL_DOUBLE")
            ;;
        --real=long-double)
            defs+=("-DREAL_LONG_DOUBLE")
            ;;
        --bit=small)
            defs+=("-DSMALL_BIT")
            ;;
        --bit=large)
            defs+=("-DLARGE_BIT")
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: ./scripts/build.sh [clean] [--real=float|double|long-double] [--bit=small|large]"
            exit 1
            ;;
    esac
done

if [[ "$clean" -eq 1 ]]; then
    rm -rf build
fi

cmake -S . -B build \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DCMAKE_CXX_FLAGS="${defs[*]}"

cmake --build build

ln -sf build/compile_commands.json compile_commands.json

# usage example:
#
# chmod +x scripts/build.sh
# 
# ./scripts/build.sh clean --real=double --bit=large
# ./scripts/build.sh --real=float --bit=small
