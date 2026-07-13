#!/usr/bin/env bash

set -e

if [[ "$1" == "clean" ]]; then
    rm -rf build
fi

cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build

ln -sf build/compile_commands.json compile_commands.json

# chmod +x scripts/build.sh
# 
# ./scripts/build.sh
# ./scripts/build.sh clean
