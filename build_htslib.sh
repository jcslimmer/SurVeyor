tar -xjf htslib-1.19.1.tar.bz2
cd htslib-1.19.1
autoheader
autoconf
./configure --prefix=`pwd`
make
make install
