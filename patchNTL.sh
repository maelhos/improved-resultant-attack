#!/bin/bash

file="ntl/include/NTL/FFT.h"
tempfile="$file.tmp"

if [[ ! -f "$file" ]]; then
    echo "Error: File '$file' not found."
    exit 1
fi

sed -e 's/#if (25 <= NTL_FFTMaxRootBnd)/#if (36 <= NTL_FFTMaxRootBnd)/' \
    -e 's/#define NTL_FFTMaxRoot (25)/#define NTL_FFTMaxRoot (36)/' "$file" > "$tempfile"

mv "$tempfile" "$file"
echo "File modified successfully."