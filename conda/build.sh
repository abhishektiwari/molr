#!/bin/bash

# Create version file with the build version
VERSION="${PKG_VERSION:-0.1.0}"
echo "version = \"$VERSION\"" > molr/_version.py

# Install the package
$PYTHON -m pip install . -vv