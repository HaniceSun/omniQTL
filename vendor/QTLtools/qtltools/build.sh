set -ex

export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig
export CPPFLAGS="-I$PREFIX/include"
export LDFLAGS="-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib"

export BOOST_INC=$PREFIX/include
export BOOST_LIB=$PREFIX/lib
export RMATH_INC=$PREFIX/lib/R/include
export RMATH_LIB=$PREFIX/lib/R/lib
export HTSLD_INC=$PREFIX/include
export HTSLD_LIB=$PREFIX/lib

sed -i '/^BOOST_INC=$/d' Makefile
sed -i '/^BOOST_LIB=$/d' Makefile
sed -i '/^RMATH_INC=$/d' Makefile
sed -i '/^RMATH_LIB=$/d' Makefile
sed -i '/^HTSLD_INC=$/d' Makefile
sed -i '/^HTSLD_LIB=$/d' Makefile
sed -i 's|\$(BOOST_LIB)/libboost_iostreams.a \$(BOOST_LIB)/libboost_program_options.a|-lboost_iostreams -lboost_program_options -ldeflate -ldl|g' Makefile

make CXX="$CXX $CPPFLAGS $LDFLAGS"
mkdir -p $PREFIX/bin
cp bin/QTLtools $PREFIX/bin
cp -r scripts $PREFIX/bin
