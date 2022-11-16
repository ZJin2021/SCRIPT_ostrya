//Parameters for the coalescence simulation program : simcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NPOP0
NPOP1
//Haploid samples sizes and samples age
26
52
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
6
//Migration matrix 0
0 MR0M01
MR0M10 0
//Migration matrix 1
0 MR1M01
MR1M10 0
//Migration matrix 2
0 MR2M01
MR2M10 0
//Migration matrix 3
0 MR3M01
MR3M10 0
//Migration matrix 4
0 MR4M01
MR4M10 0
//Migration matrix 5
0 0
0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
5 historical event
TH1$ 0 0 0 RS1R00 0 5
TH2$ 1 1 0 RS2R11 0 2
TH3$ 1 1 0 RS3R11 0 5
TH4$ 0 0 0 RS4R00 0 4
TH5$ 0 1 1 RS5R01 0 5
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, n
FREQ 1 0 2.182e-8
