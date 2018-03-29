#!/usr/bin/env bash

cd $HOME

sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -y
sudo apt-get install gcc-5 g++-5
# increase priority of gcc5 and g++5
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 150 --slave /usr/bin/g++ g++ /usr/bin/g++-5

# install boost 1.59 manually (this is the minimum requirement)
wget -O boost_1_59_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.gz/download
tar -xzf boost_1_59_0.tar.gz
cd boost_1_66_0/
sudo apt-get install autotools-dev libicu-dev build-essential libbz2-dev -y
./bootstrap.sh --with-libraries=test
./b2 -j4
sudo ./b2 install
