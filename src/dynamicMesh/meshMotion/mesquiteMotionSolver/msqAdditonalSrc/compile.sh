#!/bin/bash

rm Makefile
rm Makefile.in
autoreconf -fi

./configure \
    --enable-shared             \
    --disable-imesh             \
    --disable-igeom             \
    --disable-irel              \
    --without-cppunit           \
    --enable-trap-fpe           \
    --disable-function-timers   \
    --enable-api-doc=HTML       \
    --enable-user-guide=PDF

make mostlyclean-generic
make

# make DESTDIR=/fu/bar install
