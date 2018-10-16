#!/bin/bash

if [ "$(uname)" == "Darwin" ]; then
    # Download latest Miniconda distribution for MacOS        
    curl -o lib/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh 


elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Download latest Miniconda distribution for Linux
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./lib/miniconda.sh

fi

    # Silent install Miniconda into HOME directory
    bash ./lib/miniconda.sh -b -p $HOME/miniconda

# Add conda activation script to bashrc
if [ -f ~/.bash_profile ]; then
    echo "Appending conda path to bash_profile"
    echo -e "\n#Adding conda activation script to .bashr_profile \n. $HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bash_profile
elif [ -f ~/.bashrc ]; then
    echo "Appending conda path to bashrc"
    echo -e "\n#Adding conda activation script to .bashrc \n. $HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
else
    echo "Could not find bashrc or bash_profile, please append this to your bashrc or profile files"
    echo ". $HOME/miniconda/etc/profile.d/conda.sh"
fi

# Update conda
~/miniconda/bin/conda update -y conda
# Install Anaconda dependencies and packages
~/miniconda/bin/conda install -y anaconda-client anaconda-build conda-build
# Create NeST virtual environment
~/miniconda/bin/conda env create -n nest -f lib/nest_env.yaml
# Remove Miniconda installation script
rm lib/miniconda.sh

