sudo apt-get update
sudo apt-get install build-essential git automake libtool
git clone https://github.com/chokkan/liblbfgs.git
cd liblbfgs
./autogen.sh
./configure
make
sudo make install
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
