#!/usr/bin/env bash

set -e  # Exit immediately if a command exits with a non-zero status.
set -u  # Treat unset variables as an error.

# Set the Cplex installer name and the install path for Cplex.
CPLEX_INSTALLER="cplex_studio2211.linux_x86_64.bin"
CPLEX_INSTALL_DIR="cplex_studio2211"
CPLEX_INSTALL_PATH=$(realpath "${CPLEX_INSTALL_DIR}")

# Set up virtual environment.
python3 -m venv .env
# shellcheck disable=SC1091
source .env/bin/activate

# Install Cplex and its Python API.
chmod u+x "${CPLEX_INSTALLER}"
"./${CPLEX_INSTALLER}" -i silent -DINSTALLER_UI=silent -DLICENSE_ACCEPTED=TRUE -DUSER_INSTALL_DIR="${CPLEX_INSTALL_PATH}"
python3 "${CPLEX_INSTALL_PATH}/python/setup.py" install

# Install Compass.
pip3 install --editable .
