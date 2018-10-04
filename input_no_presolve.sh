module load cplex
#ampl -x 10000 init.mod constraints.mod data_files/data_SD13_C8_H4_J16_K8_W2_S2_L2.dat input_solve.run
#ampl -x 10000 init.mod constraints.mod data_files/data_SD0_C4_H2_J8_K4_W2_S2_L2.dat input_results.run | tee results/SD0_C4_H2_J8_K4_W2_S2_L2.txt

DATA_FILE=SD0_C4_H2_J8_K4_W2_S2_L2
TYPE=_no_presolve
ampl -x 10000 init.mod no_presolve_constraints.mod data_files/data_$DATA_FILE.dat input_no_presolve.run | tee results/$DATA_FILE$TYPE.txt
