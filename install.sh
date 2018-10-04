#!/bin/bash

if [ "$(uname)" == "Darwin" ]; then
    # Do something under Mac OS X platform        
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Do something under GNU/Linux platform
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    # Do something under 32 bits Windows NT platform
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    # Do something under 64 bits Windows NT platform
fi
# Download latest Miniconda distribution
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./lib/miniconda.sh
# Silent install Miniconda into HOME directory
bash ./lib/miniconda.sh -b -p $HOME/miniconda
# Add conda activation script to bashrc
echo -e "\n#Adding conda activation script to .bashrc \n. $HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
# Refresh terminal
source ~/.bashrc
# Update conda
conda update -y conda
# Install Anaconda dependencies and packages
conda install -y anaconda-client anaconda-build conda-build
# Create NeST virtual environment
conda env create -n nest -f lib/nest_env.yaml
# Remove Miniconda installation script
rm lib/miniconda.sh
