[ Info: Problem Parameters Successfully Loaded
[ Info: Matrices for robust oprimizations successfully constructed. 
[ Info: CCGA Modeling tools has been loaded. 
Set parameter Username
Academic license - for non-commercial use only - expires 2023-07-14
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 3273 rows, 4569 columns and 16473 nonzeros
Model fingerprint: 0x7e0b8c79
Variable types: 4065 continuous, 504 integer (504 binary)
Coefficient statistics:
  Matrix range     [3e-02, 8e+02]
  Objective range  [0e+00, 0e+00]
  Bounds range     [5e+01, 2e+02]
  RHS range        [1e+00, 1e+08]
Presolve removed 1107 rows and 3695 columns
Presolve time: 0.01s
Presolved: 2166 rows, 874 columns, 7191 nonzeros
Variable types: 384 continuous, 490 integer (490 binary)
Found heuristic solution: objective 0.0000000

Explored 0 nodes (0 simplex iterations) in 0.02 seconds (0.04 work units)
Thread count was 8 (of 8 available processors)

Solution count 1: 0 

Optimal solution found (tolerance 1.00e-03)
Best objective 0.000000000000e+00, best bound 0.000000000000e+00, gap 0.0000%

User-callback calls 299, time in user-callback 0.00 sec
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 840 rows, 577 columns and 1962 nonzeros
Model fingerprint: 0x6be4e664
Variable types: 145 continuous, 432 integer (432 binary)
Coefficient statistics:
  Matrix range     [1e+00, 1e+00]
  Objective range  [1e+00, 1e+00]
  Bounds range     [5e+01, 5e+01]
  RHS range        [1e+00, 1e+00]
Found heuristic solution: objective -0.0000000
Presolve removed 840 rows and 577 columns
Presolve time: 0.00s
Presolve: All rows and columns removed

Explored 0 nodes (0 simplex iterations) in 0.00 seconds (0.00 work units)
Thread count was 1 (of 8 available processors)

Solution count 2: 7200 -0 

Optimal solution found (tolerance 1.00e-03)
Best objective 7.200000000000e+03, best bound 7.200000000000e+03, gap 0.0000%

User-callback calls 1772, time in user-callback 0.00 sec
[ Info: Outter Forloop itr=1
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
┌ Warning: Sparse Vee Feature has been DEPRECATED, this option is now useless. 
└ @ Main ~/Desktop/repos/Robust-Optimization/src/ccga_modeling_fsp_fmp.jl:133
[ Info: [2022-10-28 14:07:19.724] Inner loop is initialized with fmp, and we are solving the initial fmp. 
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 19669 rows, 8003 columns and 67181 nonzeros
Model fingerprint: 0x5fc1acb2
Variable types: 7715 continuous, 288 integer (288 binary)
Coefficient statistics:
  Matrix range     [3e-02, 1e+08]
  Objective range  [1e+00, 1e+00]
  Bounds range     [1e+00, 1e+00]
  RHS range        [1e+00, 1e+08]
Presolve removed 5881 rows and 2037 columns
Presolve time: 0.06s
Presolved: 13788 rows, 5966 columns, 47176 nonzeros
Variable types: 5822 continuous, 144 integer (144 binary)
Found heuristic solution: objective -0.0000000

Root relaxation: objective 7.389107e+03, 16256 iterations, 0.96 seconds (2.62 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0 7389.10665    0  144   -0.00000 7389.10665      -     -    1s
H    0     0                    2287.0651454 5693.83162   149%     -    1s
H    0     0                    2833.0063627 5693.83162   101%     -    1s
     0     0 5693.83162    0  144 2833.00636 5693.83162   101%     -    1s
H    0     0                    3493.0383070 5693.83162  63.0%     -    1s
H    0     0                    4110.6362741 5681.38306  38.2%     -    1s
     0     0 5681.38306    0  144 4110.63627 5681.38306  38.2%     -    1s
     0     0 5680.65840    0  144 4110.63627 5680.65840  38.2%     -    1s
     0     0 5680.65840    0  144 4110.63627 5680.65840  38.2%     -    1s
     0     0 5513.47529    0  144 4110.63627 5513.47529  34.1%     -    1s
     0     0 5509.35030    0  144 4110.63627 5509.35030  34.0%     -    2s
     0     0 5507.48129    0  144 4110.63627 5507.48129  34.0%     -    2s
     0     0 5506.67851    0  144 4110.63627 5506.67851  34.0%     -    2s
     0     0 5495.61602    0  144 4110.63627 5495.61602  33.7%     -    2s
H    0     0                    4370.4554214 5493.81015  25.7%     -    2s
     0     0 5493.81015    0  144 4370.45542 5493.81015  25.7%     -    2s
     0     0 5493.21043    0  144 4370.45542 5493.21043  25.7%     -    2s
     0     0 5485.42886    0  144 4370.45542 5485.42886  25.5%     -    2s
     0     0 5482.83568    0  144 4370.45542 5482.83568  25.5%     -    2s
     0     0 5482.77234    0  144 4370.45542 5482.77234  25.5%     -    2s
     0     0 5482.03400    0  144 4370.45542 5482.03400  25.4%     -    2s
     0     0 5481.90550    0  144 4370.45542 5481.90550  25.4%     -    2s
     0     0 5481.88108    0  144 4370.45542 5481.88108  25.4%     -    3s
     0     0 5481.88108    0  144 4370.45542 5481.88108  25.4%     -    3s
     0     2 5481.88108    0  144 4370.45542 5481.88108  25.4%     -    3s
H   15    24                    5295.4752426 5455.58502  3.02%   333    3s
   100   101 5345.69728   13  116 5295.47524 5455.58502  3.02%   349    5s
   528   336 5401.08204    6  118 5295.47524 5434.59358  2.63%   266   10s
H  709   407                    5295.4752436 5424.53229  2.44%   260   12s
H  711   407                    5342.8100124 5424.53229  1.53%   260   12s
   877   429 infeasible   11      5342.81001 5423.31956  1.51%   257   18s
  1055   422 infeasible   10      5342.81001 5420.03927  1.45%   258   20s
  1643   434 5358.64441    8  120 5342.81001 5390.67045  0.90%   252   25s
H 1725   420                    5342.8100134 5390.42966  0.89%   251   26s
H 1726   420                    5342.8100162 5390.42966  0.89%   251   26s
H 1746   420                    5342.8100171 5390.42966  0.89%   250   26s
  1856   423 infeasible   13      5342.81002 5390.42966  0.89%   250   30s
  2247   435 infeasible   12      5342.81002 5383.95993  0.77%   250   35s
  2883   480 infeasible   11      5342.81002 5379.70918  0.69%   247   40s
  3458   456 5348.28639   13   99 5342.81002 5375.36271  0.61%   245   45s
  3892   422 5351.12594    8  132 5342.81002 5372.90635  0.56%   242   51s
  4550   348 infeasible   13      5342.81002 5369.42423  0.50%   239   55s
  5341   196 infeasible   14      5342.81002 5364.10572  0.40%   235   60s

Cutting planes:
  Gomory: 14
  MIR: 50
  Flow cover: 40
  RLT: 695

Explored 6084 nodes (1411793 simplex iterations) in 63.17 seconds (155.44 work units)
Thread count was 8 (of 8 available processors)

Solution count 10: 5342.81 5342.81 5295.48 ... -0

Optimal solution found (tolerance 1.00e-03)
Best objective 5.342810017097e+03, best bound 5.342810017097e+03, gap 0.0000%

User-callback calls 19633, time in user-callback 0.02 sec
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
[ Info: [2022-10-28 14:08:23.139] FSP is made and we are solving it. 
┌ Warning: Sparse_vee option for FMP has been deprecated for the FSP, the option is now useless. 
└ @ Main ~/Desktop/repos/Robust-Optimization/src/ccga_modeling_fsp_fmp.jl:49
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 2697 rows, 1153 columns and 10935 nonzeros
Model fingerprint: 0xe45944ca
Variable types: 1081 continuous, 72 integer (72 binary)
Coefficient statistics:
  Matrix range     [3e-02, 4e+02]
  Objective range  [1e+00, 1e+00]
  Bounds range     [0e+00, 0e+00]
  RHS range        [1e+00, 1e+08]
Found heuristic solution: objective 1500.0000000
Presolve removed 841 rows and 672 columns
Presolve time: 0.01s
Presolved: 1856 rows, 481 columns, 8148 nonzeros
Variable types: 409 continuous, 72 integer (72 binary)

Root relaxation: objective 2.994430e+01, 1049 iterations, 0.02 seconds (0.04 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

*    0     0               0      29.9442955   29.94430  0.00%     -    0s

Explored 1 nodes (1049 simplex iterations) in 0.03 seconds (0.06 work units)
Thread count was 8 (of 8 available processors)

Solution count 2: 29.9443 1500 

Optimal solution found (tolerance 1.00e-03)
Best objective 2.994429551956e+01, best bound 2.994429551956e+01, gap 0.0000%

User-callback calls 236, time in user-callback 0.00 sec
[ Info: (FSP Lower, FMP Upper) = (29.94429551955827, 5342.810017097392) at itr = 1
[ Info: Inner CCGA forloop terminated due to a positive lower of bound of: 29.94429551955827 from FSP which is higher than: ϵ=0.1.
[ Info: [2022-10-28 14:08:27.873] Introduced cut to the msp and we are solving it. 
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 3538 rows, 1729 columns and 14448 nonzeros
Model fingerprint: 0x1666be43
Variable types: 1225 continuous, 504 integer (504 binary)
Coefficient statistics:
  Matrix range     [3e-02, 8e+02]
  Objective range  [1e+00, 1e+00]
  Bounds range     [5e+01, 5e+01]
  RHS range        [1e+00, 1e+08]

MIP start from previous solve did not produce a new incumbent solution

Presolve removed 1348 rows and 831 columns
Presolve time: 0.02s
Presolved: 2190 rows, 898 columns, 8208 nonzeros
Variable types: 408 continuous, 490 integer (490 binary)
Found heuristic solution: objective 3153.0127505

Root relaxation: objective 7.200000e+03, 676 iterations, 0.01 seconds (0.01 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0 7200.00000    0    6 3153.01275 7200.00000   128%     -    0s
H    0     0                    7200.0000000 7200.00000  0.00%     -    0s

Cutting planes:
  MIR: 1
  Flow cover: 4

Explored 1 nodes (1049 simplex iterations) in 0.07 seconds (0.10 work units)
Thread count was 8 (of 8 available processors)

Solution count 2: 7200 3153.01 

Optimal solution found (tolerance 1.00e-03)
Best objective 7.200000000000e+03, best bound 7.200000000000e+03, gap 0.0000%

User-callback calls 416, time in user-callback 0.00 sec
[ Info: Objective value of msp settled at: 7200.0. 
[ Info: Outter Forloop itr=2
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
┌ Warning: Sparse Vee Feature has been DEPRECATED, this option is now useless. 
└ @ Main ~/Desktop/repos/Robust-Optimization/src/ccga_modeling_fsp_fmp.jl:133
[ Info: [2022-10-28 14:08:28.045] Inner loop is initialized with fmp, and we are solving the initial fmp. 
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 19669 rows, 8003 columns and 67233 nonzeros
Model fingerprint: 0x476de500
Variable types: 7715 continuous, 288 integer (288 binary)
Coefficient statistics:
  Matrix range     [3e-02, 1e+08]
  Objective range  [1e+00, 1e+00]
  Bounds range     [1e+00, 1e+00]
  RHS range        [1e+00, 1e+08]
Presolve removed 5915 rows and 2009 columns
Presolve time: 0.06s
Presolved: 13754 rows, 5994 columns, 47174 nonzeros
Variable types: 5850 continuous, 144 integer (144 binary)
Found heuristic solution: objective 295.0000000

Root relaxation: objective 9.450970e+03, 14671 iterations, 0.84 seconds (2.52 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0 9450.97003    0  126  295.00000 9450.97003  3104%     -    1s
H    0     0                    6646.8282817 8335.35767  25.4%     -    1s
     0     0 8335.35767    0  126 6646.82828 8335.35767  25.4%     -    1s
H    0     0                    6985.4245629 8335.35767  19.3%     -    1s
     0     0 8315.61199    0  126 6985.42456 8315.61199  19.0%     -    1s
     0     0 8312.91666    0  126 6985.42456 8312.91666  19.0%     -    1s
     0     0 8312.66363    0  126 6985.42456 8312.66363  19.0%     -    1s
     0     0 8304.20265    0  126 6985.42456 8304.20265  18.9%     -    1s
H    0     0                    6992.9944184 8303.47008  18.7%     -    1s
     0     0 8303.47008    0  126 6992.99442 8303.47008  18.7%     -    1s
     0     0 7799.19146    0  107 6992.99442 7799.19146  11.5%     -    1s
     0     0 7666.24802    0  114 6992.99442 7666.24802  9.63%     -    1s
     0     0 7645.39345    0  102 6992.99442 7645.39345  9.33%     -    1s
     0     0 7640.16472    0  108 6992.99442 7640.16472  9.25%     -    1s
     0     0 7639.32828    0  108 6992.99442 7639.32828  9.24%     -    1s
     0     0 7639.32828    0  108 6992.99442 7639.32828  9.24%     -    1s
     0     0 7464.08551    0  114 6992.99442 7464.08551  6.74%     -    2s
     0     0 7452.31254    0  114 6992.99442 7452.31254  6.57%     -    2s
     0     0 7444.32224    0  114 6992.99442 7444.32224  6.45%     -    2s
     0     0 7441.26135    0  114 6992.99442 7441.26135  6.41%     -    2s
     0     0 7440.20708    0  108 6992.99442 7440.20708  6.40%     -    2s
     0     0 7439.58740    0  108 6992.99442 7439.58740  6.39%     -    2s
     0     0 7369.47805    0  115 6992.99442 7369.47805  5.38%     -    2s
H    0     0                    7016.5166203 7369.47805  5.03%     -    2s
     0     0 7356.74019    0  115 7016.51662 7356.74019  4.85%     -    2s
     0     0 7355.64325    0  115 7016.51662 7355.64325  4.83%     -    2s
     0     0 7355.38544    0  115 7016.51662 7355.38544  4.83%     -    2s
     0     0 7328.08984    0  102 7016.51662 7328.08984  4.44%     -    3s
H    0     0                    7016.8798074 7328.08984  4.44%     -    3s
     0     0 7323.93479    0   90 7016.87981 7323.93479  4.38%     -    3s
     0     0 7322.65315    0   90 7016.87981 7322.65315  4.36%     -    3s
     0     0 7301.76849    0   96 7016.87981 7301.76849  4.06%     -    3s
     0     0 7300.28220    0   90 7016.87981 7300.28220  4.04%     -    3s
     0     0 7299.25021    0   96 7016.87981 7299.25021  4.02%     -    3s
     0     0 7298.92904    0   96 7016.87981 7298.92904  4.02%     -    3s
     0     0 7289.14230    0   84 7016.87981 7289.14230  3.88%     -    3s
     0     0 7287.40564    0   90 7016.87981 7287.40564  3.86%     -    3s
     0     0 7287.24836    0   84 7016.87981 7287.24836  3.85%     -    3s
     0     0 7283.88351    0   78 7016.87981 7283.88351  3.81%     -    3s
     0     0 7283.19203    0   72 7016.87981 7283.19203  3.80%     -    3s
     0     0 7282.13402    0   72 7016.87981 7282.13402  3.78%     -    3s
     0     0 7280.91231    0   78 7016.87981 7280.91231  3.76%     -    3s
     0     0 7280.46939    0   78 7016.87981 7280.46939  3.76%     -    3s
     0     0 7271.48611    0   60 7016.87981 7271.48611  3.63%     -    3s
     0     0 7270.60850    0   54 7016.87981 7270.60850  3.62%     -    4s
     0     0 7263.24000    0   42 7016.87981 7263.24000  3.51%     -    4s
     0     0 7263.00172    0   42 7016.87981 7263.00172  3.51%     -    4s
     0     0 7251.86754    0   60 7016.87981 7251.86754  3.35%     -    4s
     0     0 7250.64879    0   53 7016.87981 7250.64879  3.33%     -    4s
     0     0 7250.64338    0   48 7016.87981 7250.64338  3.33%     -    4s
     0     0 7249.32127    0   54 7016.87981 7249.32127  3.31%     -    4s
     0     0 7249.32127    0   48 7016.87981 7249.32127  3.31%     -    4s
     0     2 7249.32127    0   48 7016.87981 7249.32127  3.31%     -    4s
    47    57 7194.24120    7   39 7016.87981 7206.32204  2.70%   137    5s
  1459   498 7043.63615   24   21 7016.87981 7101.97983  1.21%  74.0   10s
  2611   407     cutoff   15      7016.87981 7067.39094  0.72%  73.0   16s

Cutting planes:
  Gomory: 9
  Lift-and-project: 8
  MIR: 99
  Flow cover: 21
  RLT: 629

Explored 3380 nodes (251901 simplex iterations) in 17.43 seconds (28.01 work units)
Thread count was 8 (of 8 available processors)

Solution count 6: 7016.88 7016.52 6992.99 ... 295

Optimal solution found (tolerance 1.00e-03)
Best objective 7.016879807350e+03, best bound 7.016879807350e+03, gap 0.0000%

User-callback calls 10515, time in user-callback 0.01 sec
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
[ Info: [2022-10-28 14:08:45.514] FSP is made and we are solving it. 
┌ Warning: Sparse_vee option for FMP has been deprecated for the FSP, the option is now useless. 
└ @ Main ~/Desktop/repos/Robust-Optimization/src/ccga_modeling_fsp_fmp.jl:49
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 2697 rows, 1153 columns and 10935 nonzeros
Model fingerprint: 0x62b25563
Variable types: 1081 continuous, 72 integer (72 binary)
Coefficient statistics:
  Matrix range     [3e-02, 4e+02]
  Objective range  [1e+00, 1e+00]
  Bounds range     [0e+00, 0e+00]
  RHS range        [1e+00, 1e+08]
Found heuristic solution: objective 1500.0000000
Presolve removed 837 rows and 672 columns
Presolve time: 0.01s
Presolved: 1860 rows, 481 columns, 8159 nonzeros
Variable types: 409 continuous, 72 integer (72 binary)

Root relaxation: objective 7.810297e+01, 575 iterations, 0.01 seconds (0.01 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0   78.10297    0    1 1500.00000   78.10297  94.8%     -    0s
H    0     0                      78.1029729   78.10297  0.00%     -    0s
     0     0   78.10297    0    1   78.10297   78.10297  0.00%     -    0s

Explored 1 nodes (575 simplex iterations) in 0.02 seconds (0.03 work units)
Thread count was 8 (of 8 available processors)

Solution count 2: 78.103 1500 

Optimal solution found (tolerance 1.00e-03)
Best objective 7.810297285962e+01, best bound 7.810297285962e+01, gap 0.0000%

User-callback calls 232, time in user-callback 0.00 sec
[ Info: (FSP Lower, FMP Upper) = (78.1029728596188, 7016.87980735026) at itr = 1
[ Info: Inner CCGA forloop terminated due to a positive lower of bound of: 78.1029728596188 from FSP which is higher than: ϵ=0.1.
[ Info: [2022-10-28 14:08:45.555] Introduced cut to the msp and we are solving it. 
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 6235 rows, 2881 columns and 26790 nonzeros
Model fingerprint: 0x565ddc92
Variable types: 2305 continuous, 576 integer (576 binary)
Coefficient statistics:
  Matrix range     [3e-02, 8e+02]
  Objective range  [1e+00, 1e+00]
  Bounds range     [5e+01, 5e+01]
  RHS range        [1e+00, 1e+08]

MIP start from previous solve did not produce a new incumbent solution

Presolve removed 2407 rows and 1529 columns
Presolve time: 0.06s
Presolved: 3828 rows, 1352 columns, 15638 nonzeros
Variable types: 792 continuous, 560 integer (560 binary)
Found heuristic solution: objective 2667.6796734

Root relaxation: objective 5.757649e+03, 2654 iterations, 0.09 seconds (0.14 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0 5757.64874    0  210 2667.67967 5757.64874   116%     -    0s
H    0     0                    4099.0967121 5337.65818  30.2%     -    0s
     0     0 5337.65818    0  230 4099.09671 5337.65818  30.2%     -    0s
H    0     0                    4411.1111111 5337.65818  21.0%     -    0s
     0     0 5335.41236    0  230 4411.11111 5335.41236  21.0%     -    0s
     0     0 5335.39930    0  233 4411.11111 5335.39930  21.0%     -    0s
     0     0 5287.72451    0  267 4411.11111 5287.72451  19.9%     -    0s
     0     0 5286.42164    0  264 4411.11111 5286.42164  19.8%     -    0s
     0     0 5286.40073    0  267 4411.11111 5286.40073  19.8%     -    0s
     0     0 5253.90968    0  257 4411.11111 5253.90968  19.1%     -    0s
     0     0 5246.19581    0  223 4411.11111 5246.19581  18.9%     -    0s
     0     0 5244.27756    0  236 4411.11111 5244.27756  18.9%     -    0s
     0     0 5243.94410    0  234 4411.11111 5243.94410  18.9%     -    0s
     0     0 5243.85347    0  244 4411.11111 5243.85347  18.9%     -    0s
     0     0 5233.46989    0  267 4411.11111 5233.46989  18.6%     -    0s
     0     0 5228.94665    0  270 4411.11111 5228.94665  18.5%     -    0s
     0     0 5228.36951    0  273 4411.11111 5228.36951  18.5%     -    0s
     0     0 5228.30817    0  273 4411.11111 5228.30817  18.5%     -    0s
     0     0 5222.61927    0  273 4411.11111 5222.61927  18.4%     -    0s
     0     0 5222.43702    0  276 4411.11111 5222.43702  18.4%     -    0s
     0     0 5222.09016    0  273 4411.11111 5222.09016  18.4%     -    1s
     0     0 5221.48546    0  273 4411.11111 5221.48546  18.4%     -    1s
     0     0 5221.48546    0  273 4411.11111 5221.48546  18.4%     -    1s
     0     0 5219.44468    0  269 4411.11111 5219.44468  18.3%     -    1s
     0     0 5219.27874    0  274 4411.11111 5219.27874  18.3%     -    1s
     0     0 5216.81823    0  273 4411.11111 5216.81823  18.3%     -    1s
     0     0 5216.20816    0  273 4411.11111 5216.20816  18.3%     -    1s
     0     0 5216.12562    0  275 4411.11111 5216.12562  18.2%     -    1s
     0     0 5215.94927    0  273 4411.11111 5215.94927  18.2%     -    1s
     0     0 5215.94690    0  275 4411.11111 5215.94690  18.2%     -    1s
     0     0 5215.65345    0  275 4411.11111 5215.65345  18.2%     -    1s
     0     0 5215.65345    0  274 4411.11111 5215.65345  18.2%     -    1s
     0     2 5215.65345    0  274 4411.11111 5215.65345  18.2%     -    1s
   947   811     cutoff   72      4411.11111 5083.59301  15.2%   119    5s
  1495  1157 5044.90638   10  222 4411.11111 5044.90638  14.4%   124   11s
  1759  1336 4717.65571   29  219 4411.11111 5041.24451  14.3%   132   15s
  3009  1930 4705.30330   35  192 4411.11111 5031.12670  14.1%   135   20s
  4738  2556     cutoff   66      4411.11111 4941.66952  12.0%   132   25s
  6396  3707 4448.13842   37  130 4411.11111 4899.70843  11.1%   133   30s
  8480  5101 4692.66703   25  216 4411.11111 4867.55866  10.3%   131   35s
 10128  6201     cutoff   38      4411.11111 4847.82858  9.90%   131   40s
 11879  7400 4615.18618   24  169 4411.11111 4827.50912  9.44%   130   45s
 13505  8298 4450.51683   35  170 4411.11111 4812.53294  9.10%   130   50s
 15810  9610 4460.93867   27  206 4411.11111 4797.14687  8.75%   128   55s
 18424 10852 4437.80851   37  221 4411.11111 4776.91395  8.29%   128   61s
 20422 12060 4535.80903   41  207 4411.11111 4766.12392  8.05%   128   65s
 22964 13295 4604.25959   28  157 4411.11111 4755.18013  7.80%   128   71s
 25459 14315 4472.21303   28  193 4411.11111 4740.67561  7.47%   128   76s
 27909 15279 4506.43580   27  207 4411.11111 4729.26576  7.21%   128   81s
 29535 16113 4506.60536   39  218 4411.11111 4722.55602  7.06%   128  100s
 32092 17002 4446.38862   33  194 4411.11111 4713.49354  6.86%   128  105s
 34490 17976 4522.73782   25  207 4411.11111 4704.92340  6.66%   128  110s
 37030 18875 4511.44928   38  176 4411.11111 4694.55213  6.43%   128  115s
 39362 19527 4429.67263   27  163 4411.11111 4686.99664  6.25%   128  120s
 40862 19718 4437.74032   35  214 4411.11111 4682.87193  6.16%   129  125s
 42447 20259     cutoff   23      4411.11111 4677.45349  6.04%   129  130s
 44742 20890 4444.38587   32  168 4411.11111 4671.10086  5.89%   129  136s
 47260 21389 4448.95245   26  186 4411.11111 4663.66304  5.73%   129  141s
 48879 21725 4439.68016   23  175 4411.11111 4659.38323  5.63%   130  145s
 51223 22158     cutoff   33      4411.11111 4652.74510  5.48%   130  151s
 52615 22350 4519.23344   48  157 4411.11111 4649.41121  5.40%   130  155s
 54946 22675 4411.50382   29  176 4411.11111 4642.30972  5.24%   131  161s
 56505 22842 4509.56813   39  197 4411.11111 4637.87602  5.14%   131  165s
 59021 23106     cutoff   63      4411.11111 4631.66377  5.00%   131  171s
 60622 23309 4434.67682   19  186 4411.11111 4627.14171  4.90%   132  175s
 63079 23394     cutoff   24      4411.11111 4620.96638  4.76%   132  180s

Cutting planes:
  Cover: 1
  Implied bound: 14
  Projected implied bound: 1
  Clique: 3
  MIR: 440
  Flow cover: 465
  Relax-and-lift: 32

Explored 63441 nodes (8382474 simplex iterations) in 180.01 seconds (353.30 work units)
Thread count was 8 (of 8 available processors)

Solution count 3: 4411.11 4099.1 2667.68 

Time limit reached
Best objective 4.411111111111e+03, best bound 4.620136940331e+03, gap 4.7386%

User-callback calls 144386, time in user-callback 0.04 sec
[ Info: Objective value of msp settled at: 4411.111111111111. 
[ Info: Outter Forloop itr=3
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
┌ Warning: Sparse Vee Feature has been DEPRECATED, this option is now useless. 
└ @ Main ~/Desktop/repos/Robust-Optimization/src/ccga_modeling_fsp_fmp.jl:133
[ Info: [2022-10-28 14:11:45.715] Inner loop is initialized with fmp, and we are solving the initial fmp. 
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 19669 rows, 8003 columns and 67419 nonzeros
Model fingerprint: 0x918546b6
Variable types: 7715 continuous, 288 integer (288 binary)
Coefficient statistics:
  Matrix range     [3e-02, 1e+08]
  Objective range  [1e+00, 1e+00]
  Bounds range     [1e+00, 1e+00]
  RHS range        [1e+00, 1e+08]
Presolve removed 5762 rows and 1870 columns
Presolve time: 0.06s
Presolved: 13907 rows, 6133 columns, 47913 nonzeros
Variable types: 5989 continuous, 144 integer (144 binary)
Found heuristic solution: objective -0.0000000

Root relaxation: objective 4.411111e+03, 12812 iterations, 0.73 seconds (2.11 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0 4411.11111    0  144   -0.00000 4411.11111      -     -    0s
     0     0 3200.73748    0  143   -0.00000 3200.73748      -     -    1s
H    0     0                     121.4940897 3200.73748  2534%     -    1s
     0     0 3164.25591    0  144  121.49409 3164.25591  2504%     -    1s
     0     0 3153.16351    0  144  121.49409 3153.16351  2495%     -    1s
     0     0 3152.78342    0  144  121.49409 3152.78342  2495%     -    1s
     0     0 3148.30285    0  144  121.49409 3148.30285  2491%     -    1s
     0     0 3143.21700    0  144  121.49409 3143.21700  2487%     -    1s
H    0     0                     146.9150455 3136.81804  2035%     -    1s
     0     0 3136.81804    0  144  146.91505 3136.81804  2035%     -    1s
     0     0 2655.30957    0  143  146.91505 2655.30957  1707%     -    1s
     0     0 2587.80230    0  140  146.91505 2587.80230  1661%     -    1s
     0     0 2575.02146    0  141  146.91505 2575.02146  1653%     -    1s
     0     0 2567.91319    0  141  146.91505 2567.91319  1648%     -    1s
     0     0 2562.93477    0  141  146.91505 2562.93477  1645%     -    1s
     0     0 2557.26012    0  141  146.91505 2557.26012  1641%     -    2s
     0     0 2557.25089    0  141  146.91505 2557.25089  1641%     -    2s
     0     0 2248.88828    0  144  146.91505 2248.88828  1431%     -    2s
     0     0 2226.77200    0  144  146.91505 2226.77200  1416%     -    3s
     0     0 2213.20200    0  144  146.91505 2213.20200  1406%     -    3s
     0     0 2203.04897    0  144  146.91505 2203.04897  1400%     -    3s
     0     0 2198.56721    0  144  146.91505 2198.56721  1396%     -    3s
     0     0 2195.18486    0  144  146.91505 2195.18486  1394%     -    3s
     0     0 2195.02860    0  144  146.91505 2195.02860  1394%     -    3s
     0     0 2140.04327    0  144  146.91505 2140.04327  1357%     -    4s
     0     0 2129.82525    0  144  146.91505 2129.82525  1350%     -    4s
     0     0 2127.79892    0  144  146.91505 2127.79892  1348%     -    4s
     0     0 2126.76792    0  144  146.91505 2126.76792  1348%     -    4s
     0     0 2126.41853    0  144  146.91505 2126.41853  1347%     -    4s
H    0     0                     170.1838615 2126.41853  1149%     -    4s
     0     0 2112.46210    0  144  170.18386 2112.46210  1141%     -    4s
     0     0 2109.14170    0  144  170.18386 2109.14170  1139%     -    4s
     0     0 2107.48654    0  144  170.18386 2107.48654  1138%     -    4s
     0     0 2107.02949    0  144  170.18386 2107.02949  1138%     -    5s
     0     0 2102.68851    0  144  170.18386 2102.68851  1136%     -    5s
     0     0 2101.97720    0  144  170.18386 2101.97720  1135%     -    5s
     0     0 2101.33906    0  144  170.18386 2101.33906  1135%     -    5s
     0     0 2101.33906    0  144  170.18386 2101.33906  1135%     -    5s
     0     2 2101.33906    0  144  170.18386 2101.33906  1135%     -    5s
   175   191 1652.28877   14  130  170.18386 2000.93519  1076%   389   10s
   302   319 1375.71230   24  120  170.18386 2000.93519  1076%   345   15s
*  591   575             120     216.0200991 2000.93519   826%   285   17s
   710   725  713.71312   53   83  216.02010 2000.93519   826%   299   23s
H  749   704                     260.7907578 2000.93519   667%   299   23s
H  772   734                     312.4713840 2000.93519   540%   297   24s
   834   808  479.66437   66   62  312.47138 2000.93519   540%   292   25s
  1311  1263 1141.47536   30  114  312.47138 1994.64823   538%   294   30s
H 1597  1457                     355.9832370 1994.64823   460%   291   32s
  1693  1498  469.19045   80  144  355.98324 1968.41464   453%   292   38s
  1696  1500 1112.94975   36  143  355.98324 1968.41464   453%   291   40s
  1706  1507 1710.31254   13  144  355.98324 1918.13964   439%   289   45s
  1720  1516  935.13481   39  144  355.98324 1657.95592   366%   287   50s
H 1727  1444                     362.9262524 1598.86077   341%   286   53s
H 1728  1371                     386.6142369 1539.00745   298%   286   54s
H 1728  1302                     404.6142369 1539.00745   280%   286   54s
  1732  1305  856.87932   46  144  404.61424 1487.08075   268%   285   55s
  1744  1313 1043.99339   34  144  404.61424 1378.44842   241%   283   60s
H 1750  1250                     453.0761389 1306.30656   188%   282   64s
  1753  1252 1277.72981    7  144  453.07614 1277.72981   182%   282   65s
H 1756  1190                     532.8038133 1264.36162   137%   281   67s
  1766  1197 1246.97991    7  144  532.80381 1246.97991   134%   280   70s
  1777  1205 1244.47575   27  143  532.80381 1244.47575   134%   314   75s
  1784  1210 1240.60751   17  144  532.80381 1240.60751   133%   313   80s
H 1784  1149                     555.9993581 1239.07538   123%   313   80s
  1793  1158 1227.92573   19  143  555.99936 1237.40250   123%   329   91s
  1795  1161 1201.14620   20  142  555.99936 1226.01937   121%   330   98s
  1807  1173 1175.58451   22  140  555.99936 1199.46225   116%   337  100s
H 1823  1126                     592.6132831 1197.23437   102%   345  102s
H 1826  1070                     677.7713418 1197.23437  76.6%   346  102s
  1848  1089 1141.80451   24  138  677.77134 1197.23437  76.6%   358  105s
  1922  1144 1090.84265   28  134  677.77134 1197.23437  76.6%   377  110s
  2055  1247 1044.11815   35  126  677.77134 1197.23437  76.6%   391  115s
  2213  1303  965.67002   44  104  677.77134 1197.23437  76.6%   396  120s
  2349  1402  842.07516   56   61  677.77134 1197.23437  76.6%   406  125s
  2445  1448  811.35230   66   55  677.77134 1197.23437  76.6%   411  136s
H 2471  1387                     684.8222575 1197.23437  74.8%   412  136s
  2563  1485  720.12516   82   33  684.82226 1197.23437  74.8%   416  140s
  2769  1617 1084.00401   25  137  684.82226 1179.99787  72.3%   418  146s
  2870  1616 1023.93810   27  135  684.82226 1179.99787  72.3%   416  151s
  2954  1658  989.41631   29  133  684.82226 1179.99787  72.3%   421  176s
H 2972  1533                     716.3689223 1179.99787  64.7%   422  176s
  3042  1570  965.80631   30  132  716.36892 1179.99787  64.7%   425  180s

Cutting planes:
  Gomory: 28
  Lift-and-project: 176
  MIR: 24
  Flow cover: 96
  RLT: 631

Explored 3071 nodes (1337451 simplex iterations) in 180.01 seconds (439.34 work units)
Thread count was 8 (of 8 available processors)

Solution count 10: 716.369 684.822 677.771 ... 362.926

Time limit reached
Best objective 7.163689223135e+02, best bound 1.179997871901e+03, gap 64.7193%

User-callback calls 51543, time in user-callback 0.07 sec
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
[ Info: [2022-10-28 14:14:45.885] FSP is made and we are solving it. 
┌ Warning: Sparse_vee option for FMP has been deprecated for the FSP, the option is now useless. 
└ @ Main ~/Desktop/repos/Robust-Optimization/src/ccga_modeling_fsp_fmp.jl:49
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 2697 rows, 1153 columns and 10935 nonzeros
Model fingerprint: 0x8a67b17b
Variable types: 1081 continuous, 72 integer (72 binary)
Coefficient statistics:
  Matrix range     [3e-02, 4e+02]
  Objective range  [1e+00, 1e+00]
  Bounds range     [0e+00, 0e+00]
  RHS range        [1e+00, 1e+08]
Found heuristic solution: objective 1456.2129493
Found heuristic solution: objective 1456.2129492
Presolve removed 745 rows and 672 columns
Presolve time: 0.01s
Presolved: 1952 rows, 481 columns, 8412 nonzeros
Variable types: 409 continuous, 72 integer (72 binary)

Root relaxation: objective 6.274510e+00, 790 iterations, 0.01 seconds (0.02 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0    6.27451    0    2 1456.21295    6.27451   100%     -    0s
H    0     0                       6.2745098    6.27451  0.00%     -    0s
     0     0    6.27451    0    2    6.27451    6.27451  0.00%     -    0s

Explored 1 nodes (790 simplex iterations) in 0.03 seconds (0.04 work units)
Thread count was 8 (of 8 available processors)

Solution count 2: 6.27451 1456.21 

Optimal solution found (tolerance 1.00e-03)
Best objective 6.274509803922e+00, best bound 6.274509803922e+00, gap 0.0000%

User-callback calls 240, time in user-callback 0.00 sec
[ Info: (FSP Lower, FMP Upper) = (6.274509803921578, 716.3689223135495) at itr = 1
[ Info: Inner CCGA forloop terminated due to a positive lower of bound of: 6.274509803921578 from FSP which is higher than: ϵ=0.1.
[ Info: [2022-10-28 14:14:45.934] Introduced cut to the msp and we are solving it. 
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 8932 rows, 4033 columns and 39132 nonzeros
Model fingerprint: 0x1793abbb
Variable types: 3385 continuous, 648 integer (648 binary)
Coefficient statistics:
  Matrix range     [3e-02, 8e+02]
  Objective range  [1e+00, 1e+00]
  Bounds range     [5e+01, 5e+01]
  RHS range        [1e+00, 1e+08]

MIP start from previous solve produced solution with objective 3611.11 (0.10s)
Loaded MIP start from previous solve with objective 3611.11

Presolve removed 3465 rows and 2227 columns
Presolve time: 0.09s
Presolved: 5467 rows, 1806 columns, 22791 nonzeros
Variable types: 1176 continuous, 630 integer (630 binary)

Root relaxation: objective 4.411111e+03, 3671 iterations, 0.14 seconds (0.20 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0 4411.11111    0   41 3611.11156 4411.11111  22.2%     -    0s
     0     0 4411.11111    0   48 3611.11156 4411.11111  22.2%     -    0s
     0     0 4411.11111    0   48 3611.11156 4411.11111  22.2%     -    0s
     0     0 4411.11111    0   30 3611.11156 4411.11111  22.2%     -    1s
     0     0 4411.11111    0   42 3611.11156 4411.11111  22.2%     -    1s
     0     0 4411.11111    0   41 3611.11156 4411.11111  22.2%     -    1s
     0     0 4411.11111    0   50 3611.11156 4411.11111  22.2%     -    1s
     0     0 4411.11111    0   45 3611.11156 4411.11111  22.2%     -    2s
     0     0 4411.11111    0   56 3611.11156 4411.11111  22.2%     -    2s
     0     0 4411.11111    0   42 3611.11156 4411.11111  22.2%     -    2s
     0     0 4411.11111    0   69 3611.11156 4411.11111  22.2%     -    2s
     0     0 4411.11111    0   50 3611.11156 4411.11111  22.2%     -    2s
     0     0 4411.11111    0   47 3611.11156 4411.11111  22.2%     -    3s
     0     2 4411.11111    0   45 3611.11156 4411.11111  22.2%     -    3s
    81    90 4411.11111   11   79 3611.11156 4411.11111  22.2%   428    5s
   578   999 4411.11111   41  131 3611.11156 4411.11111  22.2%   312   10s
H  913   999                    3611.1115589 4411.11111  22.2%   268   10s
H 1641  1375                    3663.1947035 4411.11111  20.4%   237   13s
  1645  1378 4411.11111   31   42 3663.19470 4411.11111  20.4%   237   15s
  1664  1396 4411.11111   17   66 3663.19470 4411.11111  20.4%   256   20s
  1763  1513 4411.11111   23  106 3663.19470 4411.11111  20.4%   279   26s
  2029  2494 4322.62953   35  155 3663.19470 4411.11111  20.4%   301   31s
  3412  2471 4411.11111   27  105 3663.19470 4411.11111  20.4%   256   35s
H 3667  2619                    3664.7277237 4411.11111  20.4%   271   37s
H 3703  2540                    3665.6346076 4411.11111  20.3%   273   37s
H 3777  2453                    3666.0766502 4411.11111  20.3%   272   37s
H 3930  2342                    3667.6096705 4411.11111  20.3%   271   37s
  4072  2411 4411.11111   36  152 3667.60967 4411.11111  20.3%   277   44s
H 4099  2345                    3667.9182354 4411.11111  20.3%   278   44s
H 4100  2204                    3724.3517448 4411.11111  18.4%   278   44s
H 4188  1981                    3774.8737965 4411.11111  16.9%   282   51s
  4940  2797 4020.62612   62  126 3774.87380 4411.11111  16.9%   282   56s
H 5512  2796                    3798.4720174 4411.11111  16.1%   285   59s
  5604  3361 4411.11111   36  119 3798.47202 4411.11111  16.1%   290   62s
  6376  3366 4411.11111   43  146 3798.47202 4411.11111  16.1%   287   68s
H 6377  3363                    3801.0390769 4411.11111  16.1%   287   68s
  6387  3741 4411.11111   44  174 3801.03908 4411.11111  16.1%   288   71s
  7663  4348 4411.11111   37  167 3801.03908 4411.11111  16.1%   283   82s
H 7669  4337                    3804.4049838 4411.11111  15.9%   283   82s
  7678  4937 4411.11111   38  175 3804.40498 4411.11111  15.9%   285   85s
  9082  5689 4341.52742   34  131 3804.40498 4411.11111  15.9%   289   91s
  9429  6219 4411.11111   37  123 3804.40498 4411.11111  15.9%   291   95s
 10640  7123 4411.11111   48  138 3804.40498 4411.11111  15.9%   301  102s
 11225  7752 4411.11111   37   82 3804.40498 4411.11111  15.9%   302  107s
 12045  8318 4411.11111   56  187 3804.40498 4411.11111  15.9%   304  112s
 12748  8325 4383.82591   45  160 3804.40498 4411.11111  15.9%   308  128s
 12756  9271 4371.02164   46  137 3804.40498 4411.11111  15.9%   308  133s
 13947 10173 4411.11111   42  150 3804.40498 4411.11111  15.9%   301  138s
 15158 11010 3948.32874  101   75 3804.40498 4411.11111  15.9%   298  147s
H15191 11010                    3804.4051594 4411.11111  15.9%   298  147s
 16250 12119 4146.00765   52  135 3804.40516 4411.11111  15.9%   301  153s
 17620 13444 4312.93242   40  175 3804.40516 4411.11111  15.9%   302  160s
 19359 14697     cutoff   65      3804.40516 4411.11111  15.9%   297  165s
 21092 15466 3891.68637   53  176 3804.40516 4411.11111  15.9%   290  171s
 22172 16409 3911.07659   83  122 3804.40516 4411.11111  15.9%   291  176s
 23484 16805 4411.11111   33  121 3804.40516 4411.11111  15.9%   289  180s

Cutting planes:
  Gomory: 3
  Implied bound: 1
  Projected implied bound: 1
  MIR: 13
  Flow cover: 64
  Relax-and-lift: 1

Explored 24002 nodes (6950243 simplex iterations) in 180.01 seconds (332.87 work units)
Thread count was 8 (of 8 available processors)

Solution count 10: 3804.41 3804.4 3801.04 ... 3665.63

Time limit reached
Best objective 3.804405159393e+03, best bound 4.411111111111e+03, gap 15.9475%

User-callback calls 74914, time in user-callback 0.05 sec
[ Info: Objective value of msp settled at: 3804.405159392837. 
[ Info: Outter Forloop itr=4
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
┌ Warning: Sparse Vee Feature has been DEPRECATED, this option is now useless. 
└ @ Main ~/Desktop/repos/Robust-Optimization/src/ccga_modeling_fsp_fmp.jl:133
[ Info: [2022-10-28 14:17:46.043] Inner loop is initialized with fmp, and we are solving the initial fmp. 
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 19669 rows, 8003 columns and 67365 nonzeros
Model fingerprint: 0x3ae07378
Variable types: 7715 continuous, 288 integer (288 binary)
Coefficient statistics:
  Matrix range     [3e-06, 1e+08]
  Objective range  [1e+00, 1e+00]
  Bounds range     [1e+00, 1e+00]
  RHS range        [3e-06, 1e+08]
Warning: Model contains large matrix coefficient range
         Consider reformulating model or setting NumericFocus parameter
         to avoid numerical issues.
Presolve removed 5809 rows and 1912 columns
Presolve time: 0.09s
Presolved: 13860 rows, 6091 columns, 47684 nonzeros
Variable types: 5947 continuous, 144 integer (144 binary)
Found heuristic solution: objective 0.0000000

Root relaxation: objective 3.804405e+03, 13005 iterations, 0.74 seconds (2.07 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0 3804.40516    0  144    0.00000 3804.40516      -     -    0s
H    0     0                     146.2721794 2709.40968  1752%     -    1s
     0     0 2709.40968    0  144  146.27218 2709.40968  1752%     -    1s
     0     0 2689.24143    0  144  146.27218 2689.24143  1739%     -    1s
     0     0 2670.25629    0  144  146.27218 2670.25629  1726%     -    1s
     0     0 2663.14238    0  144  146.27218 2663.14238  1721%     -    1s
     0     0 2640.49553    0  144  146.27218 2640.49553  1705%     -    1s
     0     0 2640.45280    0  144  146.27218 2640.45280  1705%     -    1s
H    0     0                     263.8962863 2208.68125   737%     -    1s
     0     0 2208.68125    0  142  263.89629 2208.68125   737%     -    1s
     0     0 2086.51041    0  143  263.89629 2086.51041   691%     -    1s
     0     0 2071.92610    0  143  263.89629 2071.92610   685%     -    1s
     0     0 2070.76754    0  143  263.89629 2070.76754   685%     -    1s
     0     0 2070.49537    0  143  263.89629 2070.49537   685%     -    1s
     0     0 1776.53759    0  141  263.89629 1776.53759   573%     -    2s
     0     0 1699.31091    0  142  263.89629 1699.31091   544%     -    2s
     0     0 1674.23662    0  142  263.89629 1674.23662   534%     -    2s
     0     0 1672.99180    0  141  263.89629 1672.99180   534%     -    2s
     0     0 1672.32117    0  141  263.89629 1672.32117   534%     -    2s
     0     0 1566.50169    0  136  263.89629 1566.50169   494%     -    2s
     0     0 1555.79517    0  135  263.89629 1555.79517   490%     -    2s
     0     0 1553.93349    0  134  263.89629 1553.93349   489%     -    2s
     0     0 1553.21855    0  135  263.89629 1553.21855   489%     -    2s
     0     0 1540.49388    0  128  263.89629 1540.49388   484%     -    3s
     0     0 1537.63026    0  129  263.89629 1537.63026   483%     -    3s
     0     0 1536.46927    0  129  263.89629 1536.46927   482%     -    3s
     0     0 1536.22635    0  129  263.89629 1536.22635   482%     -    3s
     0     0 1531.80156    0  128  263.89629 1531.80156   480%     -    3s
     0     0 1531.33990    0  127  263.89629 1531.33990   480%     -    3s
     0     0 1530.72712    0  122  263.89629 1530.72712   480%     -    3s
     0     0 1530.72712    0  117  263.89629 1530.72712   480%     -    3s
H    0     0                     267.3000928 1530.72712   473%     -    3s
     0     2 1530.72712    0  117  267.30009 1530.72712   473%     -    3s
H   10    16                     348.3257500 1467.69078   321%   197    4s
H   32    40                     369.2836255 1420.28452   285%   238    4s
   120   132 1235.82627   11  111  369.28363 1420.28452   285%   181    5s
H  446   464                     393.0119163 1420.28452   261%   192    7s
H  450   464                     412.6188707 1420.28452   244%   191    7s
H  838   785                     434.9850281 1418.28033   226%   181   16s
H  847   777                     442.3164529 1418.28033   221%   181   16s
H 1200  1118                     443.1101197 1418.28033   220%   178   19s
  1339  1235  623.41532   32   95  443.11012 1418.28033   220%   179   20s
  1446  1287  771.11414   32  144  443.11012 1386.17377   213%   179   25s
  1471  1304  779.04711   31  119  443.11012 1242.58011   180%   176   30s
  1487  1318 1208.28861   10   98  443.11012 1231.29087   178%   195   37s
  1489  1321 1194.10112   11  120  443.11012 1196.52286   170%   196   40s
  1617  1417 1044.96596   20   91  443.11012 1165.35501   163%   220   45s
H 1732  1424                     472.6064913 1165.35501   147%   233   47s
H 1737  1359                     491.2871160 1165.35501   137%   233   47s
  1941  1458  839.11570   37   71  491.28712 1165.35501   137%   244   50s
H 1962  1391                     505.2006467 1165.35501   131%   245   50s
  2003  1429  796.47739   41   69  505.20065 1165.35501   131%   248   58s
  2163  1588  691.17777   53   63  505.20065 1165.35501   131%   252   60s
  2658  1795 1091.88340   15  102  505.20065 1146.30025   127%   251   65s
  2995  1924  947.28407   23   88  505.20065 1146.30025   127%   259   78s
H 2996  1867                     505.2006752 1146.30025   127%   259   78s
  3003  1930  950.55577   24   90  505.20068 1146.30025   127%   260   81s
  3324  2095  753.61241   37   75  505.20068 1146.30025   127%   267   85s
* 3740  2219              97     506.3905451 1111.70193   120%   271   89s
  3814  2255  952.42494   21   89  506.39055 1111.70193   120%   273   90s
  4205  2486  759.74472   33   75  506.39055 1108.17711   119%   284   95s
  4535  2677  808.60593   27   84  506.39055 1105.29153   118%   290  100s
H 4680  2586                     528.6525054 1105.29153   109%   292  104s
H 4682  2338                     573.0483214 1105.29153  92.9%   292  104s
  4739  2414  919.73094   27   80  573.04832 1105.29153  92.9%   292  106s
  4987  2587 1052.89534   17   96  573.04832 1099.36631  91.8%   296  110s
  5420  2921  817.01567   30   81  573.04832 1094.21756  90.9%   306  117s
  5737  3182  887.57556   27   84  573.04832 1085.24637  89.4%   311  122s
  6048  3404  829.50158   24   91  573.04832 1085.24637  89.4%   316  139s
H 6065  3404                     573.0483612 1085.24637  89.4%   316  139s
H 6136  3404                     573.0483784 1085.24637  89.4%   318  139s
  6187  3501  628.02579   33   83  573.04838 1084.63698  89.3%   319  141s
  6573  3819  651.85917   44   69  573.04838 1078.53632  88.2%   323  147s
  6781  3960  791.37120   30   87  573.04838 1073.40317  87.3%   325  150s
  6990  3968  979.72432   21   89  573.04838 1072.85562  87.2%   327  163s
  6998  4081  956.76811   22   88  573.04838 1072.85562  87.2%   328  168s
  7129  4256  738.87579   30   80  573.04838 1072.85562  87.2%   332  171s
  7314  4340  585.03552   37   75  573.04838 1070.34720  86.8%   337  175s
  7700  4535  968.90697   20   92  573.04838 1065.64434  86.0%   342  180s

Cutting planes:
  Gomory: 90
  Lift-and-project: 9
  MIR: 62
  Flow cover: 204
  RLT: 1066

Explored 7761 nodes (2681658 simplex iterations) in 180.01 seconds (433.26 work units)
Thread count was 8 (of 8 available processors)

Solution count 10: 573.048 573.048 528.653 ... 442.316

Time limit reached
Best objective 5.730483783850e+02, best bound 1.065644343251e+03, gap 85.9606%

User-callback calls 47883, time in user-callback 0.07 sec
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
[ Info: [2022-10-28 14:20:46.115] FSP is made and we are solving it. 
┌ Warning: Sparse_vee option for FMP has been deprecated for the FSP, the option is now useless. 
└ @ Main ~/Desktop/repos/Robust-Optimization/src/ccga_modeling_fsp_fmp.jl:49
Set parameter MIPGap to value 0.001
Set parameter TimeLimit to value 180
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 2697 rows, 1153 columns and 10935 nonzeros
Model fingerprint: 0x9a9029f6
Variable types: 1081 continuous, 72 integer (72 binary)
Coefficient statistics:
  Matrix range     [3e-02, 4e+02]
  Objective range  [1e+00, 1e+00]
  Bounds range     [0e+00, 0e+00]
  RHS range        [3e-06, 1e+08]
Found heuristic solution: objective 1465.6155639
Found heuristic solution: objective 1465.6155639
Presolve removed 773 rows and 672 columns
Presolve time: 0.01s
Presolved: 1924 rows, 481 columns, 8335 nonzeros
Variable types: 409 continuous, 72 integer (72 binary)

Root relaxation: objective 6.917523e+00, 832 iterations, 0.01 seconds (0.02 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

*    0     0               0       6.9175228    6.91752  0.00%     -    0s

Explored 1 nodes (832 simplex iterations) in 0.03 seconds (0.04 work units)
Thread count was 8 (of 8 available processors)

Solution count 2: 6.91752 1465.62 

Optimal solution found (tolerance 1.00e-03)
Best objective 6.917522788014e+00, best bound 6.917522788014e+00, gap 0.0000%

User-callback calls 236, time in user-callback 0.00 sec
[ Info: (FSP Lower, FMP Upper) = (6.917522788014266, 573.0483783850335) at itr = 1
[ Info: Inner CCGA forloop terminated due to a positive lower of bound of: 6.917522788014266 from FSP which is higher than: ϵ=0.1.
[ Info: [2022-10-28 14:20:46.168] Introduced cut to the msp and we are solving it. 
Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[arm])
Thread count: 8 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 11629 rows, 5185 columns and 51474 nonzeros
Model fingerprint: 0xbfd7b486
Variable types: 4465 continuous, 720 integer (720 binary)
Coefficient statistics:
  Matrix range     [3e-02, 8e+02]
  Objective range  [1e+00, 1e+00]
  Bounds range     [5e+01, 5e+01]
  RHS range        [1e+00, 1e+08]

MIP start from previous solve produced solution with objective 3378.14 (0.19s)
Loaded MIP start from previous solve with objective 3378.14

Presolve removed 4527 rows and 2922 columns
Presolve time: 0.11s
Presolved: 7102 rows, 2263 columns, 27862 nonzeros
Variable types: 1560 continuous, 703 integer (700 binary)

Root relaxation: objective 3.804405e+03, 5213 iterations, 0.29 seconds (0.50 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0 3804.40516    0   19 3378.13811 3804.40516  12.6%     -    0s
     0     0 3804.40516    0   46 3378.13811 3804.40516  12.6%     -    1s
     0     0 3804.40516    0   46 3378.13811 3804.40516  12.6%     -    1s
     0     0 3804.40516    0   22 3378.13811 3804.40516  12.6%     -    1s
     0     0 3804.40516    0   22 3378.13811 3804.40516  12.6%     -    1s
     0     0 3804.40516    0   23 3378.13811 3804.40516  12.6%     -    1s
     0     0 3804.40516    0   40 3378.13811 3804.40516  12.6%     -    1s
     0     0 3804.40516    0   24 3378.13811 3804.40516  12.6%     -    2s
     0     0 3804.40516    0   31 3378.13811 3804.40516  12.6%     -    2s
     0     0 3804.40516    0   22 3378.13811 3804.40516  12.6%     -    2s
     0     0 3804.40516    0   22 3378.13811 3804.40516  12.6%     -    2s
     0     2 3804.40516    0   22 3378.13811 3804.40516  12.6%     -    2s
    76   109 3804.40516   10   34 3378.13811 3804.40516  12.6%   816    5s
   353   704 3489.99387   46   53 3378.13811 3804.40516  12.6%   725   14s
H  752   704                    3433.1677142 3804.40516  10.8%   698   14s
