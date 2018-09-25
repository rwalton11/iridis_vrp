module load cplex
#DATE_TIME=$(date +%Y%m%d-%H%M%S)
# echo >> result_log.txt
# echo -n  $DATE_TIME >> result_log.txt
# ampl -x 10000 init.mod constraints.mod data_files/data_SD2_C4_H2_J8_K4_W2_S2_L2.dat input_solve.run | tee  results/multiple_$DATE_TIME.txt
ampl -x 10000 init.mod constraints.mod data_files/data_SD0_C4_H2_J8_K4_W2_S2_L2.dat input_solve.run
# ampl -x 10000 rand_init.mod constraints.mod rand.dat multiple_solve.run>/home/rw9g11/VRP/results/multiple_$DATE_TIME.txt
#echo $DATE_TIME
