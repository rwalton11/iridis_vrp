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
option randseed 0;
#option randseed 306004749497;
option randseed;
option presolve 0;
#print "randseed 1";

param rand_locs;
param rand_jobs;
param rand_depots;
param rand_techs;
param rand_skldom;
param rand_skllvl;

param max_travel;

set JOBS = 1..rand_jobs;
set DEPOTS = rand_jobs+1..rand_jobs+rand_depots;
set TECHS = 1..rand_techs;
set SKL_DOM = 1..rand_skldom;
set SKL_LVL = 1..rand_skllvl;
set ALL_NODES = JOBS union DEPOTS;
set LOCS = 1..rand_locs; 

/*

# sets relating to nodes in the network

set JOBS;				# nodes relating to jobs to be completed
set DEPOTS;				# nodes relaitng to base location of technicians
#set ALL_NODES = JOBS union DEPOTS;

# sets relating to technicians

set TECHS;				# all techinicans

# sets relating to skill levels

set SKL_DOM;				# all skill domains
set SKL_LVL;				# all skill levels

*/

# sets relating to time

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


# params realting to the number of tech allowed on each job

#param tech_min{JOBS} >=1 integer;	# the minimum number of techs allowed to work on job j
#param tech_max{j in JOBS} >= tech_min{j} <=2 integer;	# the maximum number of techs allowed to work on job j

param tech_min{JOBS} =
	if Uniform01()<=0.8
		then 1
		else 2;


param tech_max{j in JOBS} = 
	max(tech_min[j],
	    if Uniform01()<=0.6 
	    	then 1 
		else 2
	   );

#param duration{JOBS,TECH_COUNT} >= 0;	# the duration of job j when n techs are assigned to it
# this is being split up for ease of generating random data

# duration of jobs with 1 tech
# only defined for jobs that can be completed by one tech
#param duration1{j in JOBS:tech_min[j]=1} = 
#	floor(Uniform(1,hours_week));

param dur_mean = 6;
param dur_stdv = 6;

param duration1{j in JOBS} = 
	if tech_min[j]=1
		then min(max(floor(Normal(dur_mean,dur_stdv)),1),hours_week-2*max_travel)
#		then floor(Uniform(1,hours_week-2*max_travel))
	else -1;

# duration of jobs with 2 techs
# only defined for jobs that can be completed by two techs
# half of duration with 1 if tech_min=1 and tech_max=2
param duration2{j in JOBS} =
	if tech_max[j]=1
		then -1
		else 	if tech_min[j]=2 
#				then floor(Uniform(1,hours_week-2*max_travel))
				then min(max(floor(Normal(dur_mean,dur_stdv)),1),hours_week-2*max_travel)
	     			else ceil(duration1[j]/2);

param duration{j in JOBS, n in 1..2} = 
	if n=1 
		then duration1[j]
		else duration2[j];


# parameters relating to time windows

param window{j in JOBS} = if weeks_count=1 then 1 else (if Uniform01()<0.7 then 2 else 1);

param start{j in JOBS} = 
	floor((if window[j]=1 then weeks_count else weeks_count-1)*Uniform01())*hours_week+1;

param finish{j in JOBS} = 
	start[j]+window[j]*hours_week-1;
# --- full week time windows




/*
--- irregular random time windows

param start{j in JOBS} =
	min(floor(Uniform(1,card(HOURS))),card(HOURS)-(if duration[j,1]>0 then duration[j,1] else duration[j,2]));

param finish{j in JOBS} = 
	min(floor(Uniform(1,card(HOURS)))
	+start[j]
	+ if tech_min[j]=1 
		then duration[j,1] 
		else duration[j,2],
	card(HOURS));

*/

param earliest_leave{j in ALL_NODES} = 
	if j in JOBS then
		start[j]+(if duration[j,2]>=0 then duration [j,2] else duration [j,1])
	else 0;

param latest_arrive{j in JOBS} = 
	finish[j]-(if duration[j,2]>=0 then duration [j,2] else duration [j,1]);
#param start{JOBS} >= 0;					# earliest allowable start time of job j
#param finish{j in JOBS} >= start[j]+duration[j,1];	# latest allowable finish time of job j 

# parameters relating to number of technicians allowed on each job

# parameters relating to travel time

param coords{j in LOCS, x in 1..2}=
	Uniform(0,1);

param job_loc{j in JOBS}=
	floor(Uniform(1,rand_locs+1));

param job_coord{j in JOBS, x in 1..2}=
	coords[job_loc[j],x];

param dpt_coord{d in DEPOTS, x in 1..2}=
	Uniform(0,1);

param all_coord{j in ALL_NODES, x in 1..2}=
	if j in JOBS 
		then job_coord[j,x]
		else dpt_coord[j,x]; 

#param travel_tri{i in ALL_NODES, j in ALL_NODES:i<=j}=
#	if i!=j 
#		then floor(Uniform(1,max_travel+1)) 
#		else 0;# travel time (or distance?) between all nodes in the problem


#param travel{i in ALL_NODES,j in ALL_NODES} = 
#	if i<=j
#		then travel_tri[i,j]
#		else travel_tri[j,i];

param travel{i in ALL_NODES, j in ALL_NODES}=
	ceil(3*sqrt((all_coord[i,1]-all_coord[j,1])**2+(all_coord[i,2]-all_coord[j,2])**2));



# parameters relating to skill levels and requirements

#param skl_tech{TECHS,SKL_DOM,SKL_LVL} >=0, <=1 integer;	# skill level of tech k in domain s
#param skl_req{JOBS,SKL_DOM,SKL_LVL} >=0, <=1 integer;	# skill requierement of job j in domain s
 
# random generation of skill levels, translated in to skill matrix on line below

# assign skills based on Uniform Distribution
#param skills{TECHS,SKL_DOM}=
#	floor(Uniform(0,rand_skllvl+1));

# assign skills based on Normal Distribution (with mean skllvl+1)
param skills{TECHS,SKL_DOM}=
	min(rand_skllvl,floor(Normal(rand_skllvl+1,(rand_skllvl/2)^2)));


param skl_tech{k in TECHS,s in SKL_DOM,l in SKL_LVL}= 
	if skills[k,s]>=l 
		then 1 
		else 0;

param reqs{JOBS,SKL_DOM}=
	floor(Uniform(0,rand_skllvl+1));

param skl_req{j in JOBS,s in SKL_DOM,l in SKL_LVL} = 
	if reqs[j,s]>= l 
		then (if tech_min[j]=2 
			    then 2 
			    else 1)
		else 0;


# parameters relating to work hours

param work_hours{TECHS} = 
	hours_week;		# the wweekly working hours of each technician

# paramters relating to base/depot nodes

# random genertion of depot nodes, translated in to matrix on line below
param base{TECHS} = 
	floor(Uniform(rand_jobs+1,rand_jobs+rand_depots+1));

param tech_base{k in TECHS,h in DEPOTS} = 
	if base[k]=h 
		then 1 
		else 0;	# tech k has base node h

param night_away_travel;
param night_away{j in JOBS, k in TECHS} = 
	if travel[j,base[k]]>=night_away_travel and (job_tech[j,k]=1 or sum{kp in TECHS}job_tech2[j,k,kp]>=1)
		then 1
		else 0;

param lateness{j in JOBS} = 
	max(0,floor(Normal(0,3)));


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




