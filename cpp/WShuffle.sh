#!/bin/bash -x

g++ ./WShuffle.cpp -o WShuffle.out

./WShuffle.out ../python/data/Gplus/edges.csv 10000 1 0.1 0.9 3.22 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv 20000 1 0.1 0.9 3.84 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv 30000 1 0.1 0.9 4.22 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv 40000 1 0.1 0.9 4.49 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv 50000 1 0.1 0.9 4.71 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv 60000 1 0.1 0.9 4.88 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv 70000 1 0.1 0.9 5.03 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv 80000 1 0.1 0.9 5.16 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv 90000 1 0.1 0.9 5.28 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv 100000 1 0.1 0.9 5.38 6 1 -1 5-1 
./WShuffle.out ../python/data/Gplus/edges.csv -1 1 0.1 0.9 5.20 8 1 -1 5-1 
