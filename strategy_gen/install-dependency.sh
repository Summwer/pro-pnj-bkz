
apt-get update
apt-get install -y autoconf automake libtool libsysfs-dev build-essential dh-autoreconf  pkg-config
apt-get install libmpfr-dev #ubuntu
apt-get install m4 #should install m4 first
apt-get install libgmp-dev
# apt-get install libboost-all-dev
apt-get install libalglib-dev
apt-get install libcurl4-openssl-dev libcurl4-doc libidn11-dev libkrb5-dev libldap2-dev librtmp-dev libssh2-1-dev libssl-dev zlib1g-dev



git clone https://github.com/fplll/fplll.git
cd fplll
./autogen.sh
./configure
make
sudo make install
make check
cd ..


./compile-framework.sh
