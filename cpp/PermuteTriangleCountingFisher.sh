#!/bin/bash -x
g++ ./PermuteTriangleCountingFisher.cpp -o PermuteTriangleCountingFisher.out
# ./PermuteTriangleCountingFisher.out ../python/data/Gplus/edges.csv 10000 1 3.5 0 1 1 1.5 1-1 1
# ./PermuteTriangleCounting.out ../python/data/Gplus/edges.csv 20000 1 4.14 0 1 1 1.5 1-1 1
# ./PermuteTriangleCounting.out ../python/data/Gplus/edges.csv 30000 1 4.52 0 1 1 1.5 1-1 1
# ./PermuteTriangleCounting.out ../python/data/Gplus/edges.csv 40000 1 4.8 0 1 1 1.5 1-1 1
./PermuteTriangleCounting.out ../python/data/Gplus/edges.csv 50000 1 5.02 0 1 1 1.5 1-1 1
# .\PermuteTriangleCounting.exe ..\python\data\Gplus\edges.csv 15000 3.87 2 1-1 1
# .\PermuteTriangleCounting.exe ..\python\data\Gplus\edges.csv 15000 3.87 1 1-1 1
# .\PermuteTriangleCounting.exe ..\python\data\Gplus\edges.csv 15000 3.87 1.5 1-1 1
# .\PermuteTriangleCounting.exe ..\python\data\Gplus\edges.csv 15000 3.87 0.5 1-1 1

# :: for running shuffler triangle counting
# :: <File Name> <Edges file> <NodeNum> <Eps_2)> <Eps_l> <Eps_1> <c> <t> <Sampling prob> <ItrNum-FixPerm(Only for NodeNum!=-1)> <Algorithm>
# :: .\PermuteTriangleCounting.exe ..\python\data\Gplus\edges.csv -1 1 5.2 0.1 1 -1 0 1-1 2

# :: for running new algorithm with empirical estimation
# :: <File Name> <Edges file> <NodeNum> <Eps)> <Eps_l> <Eps_1> <c> <t> <Sampling prob> <ItrNum-FixPerm(Only for NodeNum!=-1)> <Algorithm>
# :: <Eps_1> <c> <t> are not used for alg 1
# :: Note that delta is not for calculation, used for calculating only Eps_l
# ::.\PermuteTriangleCounting.exe ..\python\data\Gplus\edges.csv 15000 1 3.87 0 1 1 3 1-1 1
