#!/bin/bash -e

PROJ_DIR="$(dirname "$0")/.."

if [ $# -eq 0 ] ; then
  find "${PROJ_DIR}/src" -type f \
    \( -name \*.c -o -name \*.cc -o -name \*.cpp \) \
    -print0 \
    | xargs -0 \
    clang-tidy -format-style=google -p "$PROJ_DIR" -header-filter=src/.*
else
  clang-tidy -format-style=google -p "$PROJ_DIR" -header-filter=src/.* \
    "$@"
fi

