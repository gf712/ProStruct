#!/usr/bin/env bash

cd $HOME

sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -y
sudo apt-get install gcc-7 g++-7 autotools-dev libicu-dev build-essential libbz2-dev libblas-dev liblapack-dev -y
# increase priority of gcc7 and g++7
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 150 --slave /usr/bin/g++ g++ /usr/bin/g++-7

# install armadillo 7.8 manually (trusty only has older versions)
git clone -b 7.800.x https://gitlab.com/conradsnicta/armadillo-code.git
cd armadillo-code/
./configure && make
sudo make install