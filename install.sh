#!/bin/bash

if [ "$(uname)" == "Darwin" ]; then
    # Download latest Miniconda distribution for MacOS        
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ./lib/miniconda.sh
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Download latest Miniconda distribution for Linux
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./lib/miniconda.sh
fi
# Silent install Miniconda into HOME directory
bash ./lib/miniconda.sh -b -p $HOME/miniconda
# Add conda activation script to bashrc
echo -e "\n#Adding conda activation script to .bashrc \n. $HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
# Refresh terminal
source ~/.bashrc
# Update conda
~/miniconda/bin/conda update -y conda
# Install Anaconda dependencies and packages
~/miniconda/bin/conda install -y anaconda-client anaconda-build conda-build
# Create NeST virtual environment
~/miniconda/bin/conda env create -n nest -f lib/nest_env.yaml
# Remove Miniconda installation script
rm lib/miniconda.sh
