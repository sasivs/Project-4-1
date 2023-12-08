g++ .\Clique4Counting.cpp -o Clique4Counting
:: for running Shuffler 4-clique counting
:: <File Name> <Edges file> <NodeNum> <Eps_2> <Eps_l> <Eps_1> <c> <t> <Sampling prob> <ItrNum-FixPerm(Only for NodeNum!=-1)> <Algorithm>
.\Clique4Counting.exe ..\python\data\Gplus\edges.csv -1 0.9 5.2 0.1 1 -1 1 1-1 2

:: for running 4-clique non interactive without arr (RR)
:: <File Name> <Edges file> <NodeNum> <Eps> <Eps_l> <Eps_1> <c> <t> <Sampling prob> <ItrNum-FixPerm(Only for NodeNum!=-1)> <Algorithm>
:: <Eps_l> <Eps_1> <c> <t> <Sampling prob> are not used for this alg

:: for running 4-clique non interactive with arr
:: <File Name> <Edges file> <NodeNum> <Eps> <Eps_l> <Eps_1> <c> <t> <Sampling prob> <ItrNum-FixPerm(Only for NodeNum!=-1)> <Algorithm>
:: <Eps_l> <Eps_1> <c> <t> are not used for this alg

