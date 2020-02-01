#!/bin/bash -e

PROJ_DIR="$(dirname "$0")/.."

DIFF_CMD=diff
DIFF_PARAMS=-u

if [ -t 1 ] ; then 
  if [ -x /usr/bin/colordiff ] ; then
    DIFF_CMD=colordiff
    PAGER=${PAGER:=less -r}
  else
    PAGER=${PAGER:=less}
  fi
else
  PAGER=cat
fi


format_diff() {
  local f="$1"
  "$DIFF_CMD" $DIFF_PARAMS "$f" <(clang-format -style=google "$f")
}

export DIFF_CMD
export DIFF_PARAMS
export -f format_diff
find "${PROJ_DIR}/src" -type f \
  \( -name \*.c -o -name \*.cc -o -name \*.cpp -o -name \*.h \) \
  -print0 \
  | xargs -0 -n1 bash -c 'format_diff "$1"' {} \
  | $PAGER

