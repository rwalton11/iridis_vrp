module load cplex
DATE_TIME=$(date +%Y%m%d-%H%M%S)
echo >> result_log.txt
echo -n  $DATE_TIME >> result_log.txt
 ampl -x 10000 rand_init.mod constraints.mod rand.dat multiple_solve.run | tee  results/multiple_$DATE_TIME.txt
# ampl -x 10000 rand_init.mod constraints.mod rand.dat multiple_solve.run>/home/rw9g11/VRP/results/multiple_$DATE_TIME.txt
echo -n  '	improved' >> result_log.txt
echo $DATE_TIME
