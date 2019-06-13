#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: sh c [sim subdir]"
    exit 1
fi
echo $#

cat main.cpp | grep '// RUN ' | cut -f 3- -d ' ' | sh /dev/stdin "$@"
