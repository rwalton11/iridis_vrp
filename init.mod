#  Define sets in the problem

# --------------------------------------------------------------------------------------
# first we define all the sets and paramters that are read in from the .dat file
# except for those sets required to read in others sets/paramters
# --------------------------------------------------------------------------------------

/* 
----------------------------------------------------------------------------------------
 The sets and parameters in the following section are for use with the rand.dat file
 to randomly generate data with minimal input.
----------------------------------------------------------------------------------------
*/

param log_file = 'output.txt';
#print "  " > log_file;
#print "---------------------------------" > ;
#print "beginning of .mod file";
#print "---------------------------------";
#print "  ";

print "start of log file: "&ctime() > log_file;
printf "%s\n%s\n%s\n%s\n%s", "  ","------------------------------","beginning of .mod file","-------------------------------","  ">log_file;


#option randseed 0 > log_file;
#option randseed 0;
#option randseed 306004749497;
#option randseed;
option presolve 0;
#print "randseed 1";

param jobs;
param depots;
param techs;
param skldom;
param skllvl;

param max_travel;

set JOBS = 1..jobs;
set DEPOTS = jobs+1..jobs+depots;
set TECHS = 1..techs;
set SKL_DOM = 1..skldom;
set SKL_LVL = 1..skllvl;
set ALL_NODES = JOBS union DEPOTS;


param weeks_count >= 0 integer;		# number of weeks in the full planning period
param days_week >= 0 integer;		# number of working days per week

set DAYS = 1..days_week ordered; 		# days in each week {1, ..., 5}
					# ***not read in but required for others***

param hours_day{DAYS} >= 0;# number of working hours for each day of the week

param hours_week = sum{d in DAYS} hours_day[d];

# sets/parameters relating to job duration

set TECH_COUNT = 1..2;			# possible number of techs that cna be assigned to jobs
					# ***not read in but required for others***
set WEEKS = 
	1..weeks_count ordered; 			# all weeks in planning period

set HOURS = 
	1..hours_week*weeks_count ordered;		# all working hours across full planning period

# set of hours in week w
set HOURS_W{w in WEEKS} within HOURS = 
	(w-1)*hours_week+1
	..
	w*hours_week ordered;

# set of hours in week w, day d
# if statement to handle sum to d-1 when d=1

set HOURS_WD{w in WEEKS, d in DAYS} within HOURS = 
	(w-1)*hours_week
	+(if d = 1 
		then 0 
		else sum{j in 1..d-1}hours_day[j]
	 )
	 +1
	 ..
	 (w-1)*hours_week
	 +sum{j in 1..d}hours_day[j] ordered;


# params relating to the feasibility of techs performing jobs
# values of these are claculated in .dat file
param job_tech{JOBS,TECHS} default 0; 		# parameter takes value 1 if tech k can do job j alone, 0 otherwise
param job_tech2{JOBS,TECHS,TECHS} default 0;	# parameter takes value 1 if techs k and k' can do job j together, 0 otherwise 


# params relating to the number of tech allowed on each job

param tech_min{JOBS} >=1 integer;	# the minimum number of techs allowed to work on job j
param tech_max{JOBS} <=2 integer;	# the maximum number of techs allowed to work on job j


param duration{JOBS,TECH_COUNT} ;	# the duration of job j when n techs are assigned to it


# parameters relating to time windows

param start{JOBS} >= 0;					# earliest allowable start time of job j
param finish{j in JOBS} >= start[j]+duration[j,1];	# latest allowable finish time of job j 

param earliest_leave{j in ALL_NODES} = 
	if j in JOBS then
		start[j]+(if duration[j,2]>=0 then duration [j,2] else duration [j,1])
	else 1;

param latest_arrive{j in JOBS} = 
	finish[j]-(if duration[j,2]>=0 then duration [j,2] else duration [j,1]);

# parameters relating to number of technicians allowed on each job

# parameters relating to travel time

param travel{i in ALL_NODES,j in ALL_NODES} >=0; 

# parameters relating to skill levels and requirements

param skl_tech{TECHS,SKL_DOM,SKL_LVL} >=0, <=1 integer;	# skill level of tech k in domain s
param skl_req{JOBS,SKL_DOM,SKL_LVL} >=0, <=2 integer;	# skill requierement of job j in domain s
 

# parameters relating to work hours

param work_hours{TECHS}>=0;

#param work_hours{TECHS} = 
#	hours_week;		# the weekly working hours of each technician

# paramters relating to base/depot nodes


param base{TECHS} >= 0; 

param tech_base{k in TECHS,h in DEPOTS} = 
	if base[k]=h 
		then 1 
		else 0;	# tech k has base node h

param night_away_travel;
param night_away{j in JOBS, k in TECHS} = 
	if travel[j,base[k]]>=night_away_travel and (job_tech[j,k]=1 or sum{kp in TECHS}job_tech2[j,k,kp]>=1)
		then 1
		else 0;

param lateness{j in JOBS} >= 0;


param test;				# used in the calculation of other parameters for comparison tests


# --------------------------------------------------------------------------------------
# parameters relating to the objective function
# --------------------------------------------------------------------------------------

param weight{1..3}; 		# weighting of each term in the objective function
				# 1: travel time/distance, 2: nights away, 3: lateness

param cost_travel;		# cost per distance/time pf travel

param cost_nights;		# cost per night away

param cost_late; 		# 'cost' of job lateness (per time unit)

# --------------------------------------------------------------------------------------
# now we perform calculations of other sets and parameters
# --------------------------------------------------------------------------------------



# sets/parameters relating to days/weeks in the planning period
/*
set WEEKS = 
	1..weeks_count ordered; 			# all weeks in planning period

set HOURS = 
	1..hours_week*weeks_count ordered;		# all working hours across full planning period

# set of hours in week w
set HOURS_W{w in WEEKS} within HOURS = 
	(w-1)*hours_week+1
	..
	w*hours_week ordered;

# set of hours in week w, day d
# if statement to handle sum to d-1 when d=1

set HOURS_WD{w in WEEKS, d in DAYS} within HOURS = 
	(w-1)*hours_week
	+(if d = 1 
		then 0 
		else sum{j in 1..d-1}hours_day[j]
	 )
	 +1
	 ..
	 (w-1)*hours_week
	 +sum{j in 1..d}hours_day[j] ordered;
*/
/*
param night_split{j in JOBS, n in TECH_COUNT, t in  HOURS_W[1]} = 
	if tech_min[j]<=n<=tech_max[j] then 
		floor(((if t mod hours_day[1]>0 then t mod hours_day[1] else hours_day[1])+duration[j,n]-1)/hours_day[1])
		else 0;
*/

set JOB_TECH within {JOBS,TECHS} default {};		# all sets of permissible tech and job combinations (j,k)

set JOB_TECH2 within {JOBS,k in TECHS,kp in TECHS:k!=kp} default {};	# all sets of permissible tech pair and job combinations (j,k,k')

param max_hours = max{j in JOBS}finish[j];

param week_t{t in 1..max_hours}=ceil(t/hours_week);

param day_t{t in 1..(weeks_count+1)*hours_week}; # has to be calculated in .dat file due to requirement of let command

param night_split{j in JOBS, n in TECH_COUNT, t in HOURS} = 
	if tech_min[j]<=n<=tech_max[j] then 
		floor(((if (t-(week_t[t]-1)*hours_week) mod hours_day[1]>0 then (t-(week_t[t]-1)*hours_week) mod hours_day[1] else hours_day[1])+duration[j,n]-1)/hours_day[1])
		else 0;

printf "%s\n%s\n%s\n%s\n%s", "  ","------------------------------","end of .mod file","-------------------------------","  ">log_file;




