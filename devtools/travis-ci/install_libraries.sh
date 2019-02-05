#!/usr/bin/env bash

cd $HOME

sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -y
sudo apt-get install gcc-5 g++-5 autotools-dev libicu-dev build-essential libbz2-dev libblas-dev liblapack-dev -y
# increase priority of gcc5 and g++5
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 150 --slave /usr/bin/g++ g++ /usr/bin/g++-5

# install armadillo 7.8 manually (trusty only has older versions)
git clone -b 7.800.x https://github.com/conradsnicta/armadillo-code.git
cd armadillo-code/
./configure && make
sudo make install