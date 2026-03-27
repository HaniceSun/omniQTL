set -ex

export CPPFLAGS="-I$PREFIX/include"
export LDFLAGS="-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib"

./configure --prefix=$PREFIX --libdir=$PREFIX/lib --with-x=no
make
make install

cd $SRC_DIR/src/nmath/standalone
make
cp libRmath* $PREFIX/lib/R/lib
