#!/bin/sh
curl -fsSL https://zenodo.org/records/14674981/files/SourceCode_2024_11_07.zip | bsdtar xv
./build.sh SourceCode
