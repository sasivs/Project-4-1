#!/bin/bash -x
g++ ./PermuteTriangleCountingFisher.cpp -o PermuteTriangleCountingFisher.out
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 10000 1 3.5 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 20000 1 4.14 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 30000 1 4.52 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 40000 1 4.8 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 50000 1 5.02 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 60000 1 5.19 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 70000 1 5.34 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 80000 1 5.48 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 90000 1 5.59 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 100000 1 5.56 0 1 1 1.5 1-1 1

# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 100000 1 5.56 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 100000 1 5.56 0 1 1 1 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 100000 1 5.56 0 1 1 0.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 100000 1 5.56 0 1 1 0 1-1 1

# Deg Clipping
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 10000 1 3.22 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 10000 1 3.22 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 10000 1 3.22 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 10000 1 3.22 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 10000 1 3.22 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 10000 1 3.22 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 10000 1 3.22 0.1 1 1 3 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 20000 1 3.84 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 20000 1 3.84 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 20000 1 3.84 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 20000 1 3.84 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 20000 1 3.84 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 20000 1 3.84 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 20000 1 3.84 0.1 1 1 3 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 30000 1 4.22 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 30000 1 4.22 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 30000 1 4.22 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 30000 1 4.22 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 30000 1 4.22 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 30000 1 4.22 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 30000 1 4.22 0.1 1 1 3 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 40000 1 4.49 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 40000 1 4.49 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 40000 1 4.49 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 40000 1 4.49 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 40000 1 4.49 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 40000 1 4.49 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 40000 1 4.49 0.1 1 1 3 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 50000 1 4.71 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 50000 1 4.71 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 50000 1 4.71 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 50000 1 4.71 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 50000 1 4.71 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 50000 1 4.71 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 50000 1 4.71 0.1 1 1 3 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 60000 1 4.88 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 60000 1 4.88 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 60000 1 4.88 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 60000 1 4.88 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 60000 1 4.88 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 60000 1 4.88 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 60000 1 4.88 0.1 1 1 3 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 70000 1 5.03 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 70000 1 5.03 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 70000 1 5.03 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 70000 1 5.03 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 70000 1 5.03 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 70000 1 5.03 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 70000 1 5.03 0.1 1 1 3 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 80000 1 5.16 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 80000 1 5.16 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 80000 1 5.16 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 80000 1 5.16 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 80000 1 5.16 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 80000 1 5.16 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 80000 1 5.16 0.1 1 1 3 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 90000 1 5.28 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 90000 1 5.28 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 90000 1 5.28 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 90000 1 5.28 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 90000 1 5.28 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 90000 1 5.28 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 90000 1 5.28 0.1 1 1 3 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 100000 1 5.38 0.1 1 1 0 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 100000 1 5.38 0.1 1 1 0.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 100000 1 5.38 0.1 1 1 1 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 100000 1 5.38 0.1 1 1 1.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 100000 1 5.38 0.1 1 1 2 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 100000 1 5.38 0.1 1 1 2.5 5-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 100000 1 5.38 0.1 1 1 3 5-1 1


# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 200000 1 6.25 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 200000 1 6.25 0 1 1 1 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 200000 1 6.25 0 1 1 0.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 200000 1 6.25 0 1 1 0 1-1 1
# USed delta = 10^-8 for n=3*10^6
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 300000 1 6.53 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 300000 1 6.53 0 1 1 1 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 300000 1 6.53 0 1 1 0.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 300000 1 6.53 0 1 1 0 1-1 1

# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 500000 1 7.04 0 1 1 1.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 500000 1 7.04 0 1 1 1 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 500000 1 6.53 0 1 1 0.5 1-1 1
# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 500000 1 6.53 0 1 1 0 1-1 1

# ./PermuteTriangleCountingFisher.out ../python/data/IMDB/edges.csv 100000 1 5.56 0 1 1 1.5 1-1 1

# :: for running shuffler triangle counting
# :: <File Name> <Edges file> <NodeNum> <Eps_2)> <Eps_l> <Eps_1> <c> <t> <Sampling prob> <ItrNum-FixPerm(Only for NodeNum!=-1)> <Algorithm>
# :: .\PermuteTriangleCounting.exe ..\python\data\Gplus\edges.csv -1 1 5.2 0.1 1 -1 0 1-1 2

# :: for running new algorithm with empirical estimation
# :: <File Name> <Edges file> <NodeNum> <Eps)> <Eps_l> <Eps_1> <c> <t> <Sampling prob> <ItrNum-FixPerm(Only for NodeNum!=-1)> <Algorithm>
# :: <Eps_1> <c> <t> are not used for alg 1
# :: Note that delta is not for calculation, used for calculating only Eps_l
# ::.\PermuteTriangleCounting.exe ..\python\data\Gplus\edges.csv 15000 1 3.87 0 1 1 3 1-1 1
