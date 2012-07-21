#!/bin/sh

rm lorenzShared
touch lorenzShared

make rodes

gdb rodes

#./rodes 1 lorenzShared 1255 757 8 > log_1.txt
