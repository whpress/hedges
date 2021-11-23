#!/bin/sh

g++ -fPIC -fpermissive -w -c NRpyDNAcode.cpp -o NRpyDNAcode.o -I/usr/include/python2.7 \
 -I/usr/local/lib/python2.7/dist-packages/numpy/core/include
g++ -shared NRpyDNAcode.o -o NRpyDNAcode.so

g++ -fPIC -fpermissive -w -c NRpyRS.cpp -o NRpyRS.o -I/usr/include/python2.7 \
 -I/usr/local/lib/python2.7/dist-packages/numpy/core/include
g++ -shared NRpyRS.o -o NRpyRS.so

echo "done"
