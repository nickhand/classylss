#!/bin/sh -e

CLASS_VERSION=$1; shift
PREFIX="$1"; shift

TMP="tmp-class-v$CLASS_VERSION"
LOGFILE="build.log"

START=$(pwd)
mkdir -p $PREFIX;
ROOT=`dirname $0`/../
cd $ROOT/depends; mkdir -p $TMP 
if ! [ -f $ROOT/depends/class-v$CLASS_VERSION.tar.gz ]; then
wget https://github.com/lesgourg/class_public/archive/v$CLASS_VERSION.tar.gz \
    -O $ROOT/depends/class-v$CLASS_VERSION.tar.gz 
fi

gzip -dc $ROOT/depends/class-v$CLASS_VERSION.tar.gz | tar xf - -C $TMP
cd $TMP/class_public-$CLASS_VERSION

make libclass.a
cp -r include $START/$PREFIX/
mkdir -p $START/$PREFIX/lib
cp libclass.a $START/$PREFIX/lib

# (
#
# ../class-${PFFT_VERSION}/configure --prefix=$PREFIX --disable-shared --enable-static  \
# --disable-fortran --disable-doc --enable-mpi ${OPTIMIZE} &&
# make -j 4   &&
# make install && echo "PFFT_DONE"
# ) 2>&1 > ${LOGFILE}.double
#
# if ! grep PFFT_DONE ${LOGFILE}.double > /dev/null; then
#     tail ${LOGFILE}.double
#     exit 1
# fi
# (
# mkdir -p single;cd single
# ../pfft-${PFFT_VERSION}/configure --prefix=$PREFIX --enable-single --disable-shared --enable-static  \
# --disable-fortran --disable-doc --enable-mpi $2 ${OPTIMIZE1} &&
# make -j 4  &&
# make install && echo "PFFT_DONE"
# ) 2>&1 > ${LOGFILE}.single
#
# if ! grep PFFT_DONE ${LOGFILE}.single > /dev/null; then
#     tail ${LOGFILE}.single
#     exit 1
# fi
