/*
file containing parameters for the VRP  with time windows and skill levels

Author: Ruth Walton
*/

# ------------------------------------------------------------
# declaration of variables
# ------------------------------------------------------------


printf "%s\n%s\n%s\n%s\n%s", "  ","------------------------------","beginning of constraints file","-------------------------------","  ">log_file;

/*
set set_LEAVE = {i in ALL_NODES,j in ALL_NODES,k in TECHS,t in HOURS:
	(if i in JOBS then job_tech[i,k]=1 or if i in JOBS then sum{kp in TECHS}job_tech2[i,k,kp]>=1 or if i in DEPOTS then base[k]=i) and
	(if j in JOBS then job_tech[j,k]=1 or if j in JOBS then sum{kp in TECHS}job_tech2[j,k,kp]>=1 or if j in DEPOTS then base[k]=j) and
	if (i in JOBS and j in JOBS) then start[j]<=t+travel[i,j] and
#	if (i in JOBS and j in JOBS) then t+travel[i,j] - finish [j] + (if duration[j,2]>=0 then duration[j,2] else duration[j,1]) <= 0
	if (i in JOBS and j in JOBS) then t+travel[i,j] - finish [j] <=0 and 
	if (j in JOBS and duration[j,2]>=0) then t+travel[i,j]+duration[j,2]<=hours_week else t+travel[i,j]+duration[j,1]<=hours_week
	};
*/

set set_LEAVE = {i in ALL_NODES,j in ALL_NODES,k in TECHS,t in HOURS:
	i!=j and
	(if i in JOBS then job_tech[i,k]=1 or if i in JOBS then sum{kp in TECHS}job_tech2[i,k,kp]>=1 or if i in DEPOTS then base[k]=i) and
	(if j in JOBS then job_tech[j,k]=1 or if j in JOBS then sum{kp in TECHS}job_tech2[j,k,kp]>=1 or if j in DEPOTS then base[k]=j) and
	earliest_leave[i]-t<=0 and
	if j in JOBS then start[j]-t-travel[i,j]<=0 and
	if j in JOBS then t+travel[i,j]-latest_arrive[j]<=0 and
#	if j in JOBS then day_t[t+travel[i,j]]+max{n in TECH_COUNT}night_split[j,n,t+travel[i,j]]<=5 and
	if j in JOBS then day_t[t+travel[i,j]]+(if tech_max[j]=2 then night_split[j,2,t+travel[i,j]] else night_split[j,1,t+travel[i,j]])<=5 and
	day_t[t]-day_t[t+travel[i,j]]<=0
	};

/*
for {t in 74..75, i in 6..6, j in 4..4}{
	printf "testing set_LEAVE for t=%d", t;
	printf "earliest_leave[%d] <= t: %d <= %d",i,earliest_leave[6],t;
	printf "start[%d] <= t+travel[%d,%d]: %d <= %d + %d",j,i,j, start[4],t,travel[6,4];
	printf "t +";
	print "  ";
};
*/
set set_ARRIVE = {i in ALL_NODES,j in ALL_NODES,k in TECHS,t in HOURS:
	i!=j and
	(if i in JOBS then job_tech[i,k]=1 or if i in JOBS then sum{kp in TECHS}job_tech2[i,k,kp]>=1 or if i in DEPOTS then base[k]=i) and
	(if j in JOBS then job_tech[j,k]=1 or if j in JOBS then sum{kp in TECHS}job_tech2[j,k,kp]>=1 or if j in DEPOTS then base[k]=j) and
	earliest_leave[i]-t+travel[i,j]<=0 and
	if j in JOBS then start[j]-t<=0 and
	if j in JOBS then t-latest_arrive[j]<=0 and
#	if j in JOBS then day_t[t]+max{n in TECH_COUNT}night_split[j,n,t]<=5 and
	if j in JOBS then day_t[t]+(if tech_max[j]=2 then night_split[j,2,t] else night_split[j,1,t])<=5 and
	day_t[max(1,t-travel[i,j])]-day_t[t]<=0
	};

	
set set_FEASIBLE_HOURS = {j in JOBS, k in TECHS, t in HOURS:
	(job_tech[j,k]=1 or sum{kp in TECHS:k!=kp}job_tech2[j,k,kp]>=1) and
	start[j]<=t<=finish[j]	
	};

set set_FEASIBLE_WEEKS = {j in JOBS, k in TECHS, w in WEEKS:
	(job_tech[j,k]=1 or sum{kp in TECHS:k!=kp}job_tech2[j,k,kp]>=1) and
	week_t[start[j]]<=w<=week_t[finish[j]]	
	};

set set_FEASIBLE_DAYS = {j in JOBS, k in TECHS, w in WEEKS, d in DAYS:
	(job_tech[j,k]=1 or sum{kp in TECHS:k!=kp}job_tech2[j,k,kp]>=1) and
	week_t[start[j]]<=w<=week_t[finish[j]] and
	#(j,k,w) in set_FEASIBLE_WEEKS and 
	#(if w = week_t[start[j]] then d>=day_t[start[j]]) and
	(if w=week_t[finish[j]] then d<=day_t[finish[j]]) and
	(if w=week_t[start[j]] then day_t[start[j]]-d<=0)
	};


# variable x in formulation
var var_leave{(i,j,k,t) in set_LEAVE} binary;

# variable bar(x) in formulation
var var_arrive{(i,j,k,t) in set_ARRIVE} binary;

# variable y in formulation
var var_tracking{(j,k,t) in set_FEASIBLE_HOURS} binary;

# variable u in formulation
var var_start{j in JOBS} >= 0;

# variable z in formulation
var var_twotech{j in JOBS} binary;

# variable lambda in formulation
#var var_night_away{(j,k,w) in set_FEASIBLE_WEEKS, d in DAYS:
var var_night_away{(j,k,w,d) in set_FEASIBLE_DAYS:
	night_away[j,k]=1
	} binary;

# variable sigma in formultion
#var var_hours_tech{(j,k,w) in set_FEASIBLE_WEEKS,d in DAYS};
var var_hours_tech{(j,k,w,d) in set_FEASIBLE_DAYS} >= 0;

# variable rho in formultion
#var var_job_split{(j,k,w) in set_FEASIBLE_WEEKS, d in DAYS
var var_job_split{(j,k,w,d) in set_FEASIBLE_DAYS
#	:d!=last(DAYS)
	} binary;

# variable p in formulation
var var_tech_pair{j in JOBS, k in TECHS, kp in TECHS:
	job_tech2[j,k,kp]=1
	} binary;


# variable epsilon in formulation
var var_late{j in JOBS} binary;


# ------------------------------------------------------------
# objective function
# ------------------------------------------------------------

minimize objective:
	weight[1]*cost_travel*(sum{(i,j,k,t) in set_LEAVE}travel[i,j]*var_leave[i,j,k,t])+
	weight[2]*cost_nights*(sum{(j,k,w,d) in set_FEASIBLE_DAYS: night_away[j,k]=1}var_night_away[j,k,w,d])+
	weight[3]*cost_late*(sum{j in JOBS}var_late[j])+
#	weight[3]*cost_late*(sum{j in JOBS}var_late[j]*(var_start[j]+(if duration[j,1]>=0 then (1-var_twotech[j])*duration[j,1])+(if duration[j,2]>=0 then (var_twotech[j]*duration[j,2]))-finish[j]))
#		  	    (sum{(i,j,k,t) in set_ARRIVE}t*var_arrive[i,j,k,t])+
0;

# ------------------------------------------------------------
# declaration of constraints
# ------------------------------------------------------------


# ------------------------------------------------------------
# constraints relating to number of technicians (lb, ub, pairs etc)
# ------------------------------------------------------------

# each job is scheduled, with at least the minimum number of required techs
subject to c_tech_min{j in JOBS}:
	sum {i in ALL_NODES, k in TECHS, t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t] >= 
		tech_min[j];

# each job is scheduled, with at most the maximum number of allowed techs
subject to c_tech_max{j in JOBS}:
	sum {i in ALL_NODES, k in TECHS, t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t] <= 
		tech_max[j];

# number of techs arriving and leaving each job is equal
subject to c_tech_inout{j in ALL_NODES,k in TECHS}:
	sum {i in ALL_NODES, t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t] =
		sum {h in ALL_NODES, t in HOURS:(j,h,k,t) in set_LEAVE}var_leave[j,h,k,t];

# this is a relaxation of the constraint above, summed across k
#subject to c_tech_inout{j in ALL_NODES}:
#	sum {i in ALL_NODES, k in TECHS, t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t] =
#		sum {h in ALL_NODES, k in TECHS, t in HOURS:(j,h,k,t) in set_LEAVE}var_leave[j,h,k,t];

# number of techs arriving at each job is equal to the number completing it
subject to c_twotech{j in JOBS}:
	sum {i in ALL_NODES, k in TECHS, t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t] = 
		var_twotech[j]+1;

# tech can only be assigned in a pair to a job if it is assigned to that job
#subject to c_pair_ub{(j,k) in JOB_TECH:tech_max[j]=2}:
subject to c_pair_ub{j in JOBS, k in TECHS:tech_max[j]=2 and (job_tech[j,k]=1 or sum{kp in TECHS}job_tech2[j,k,kp]>=1)}:
	sum{kp in TECHS:k!=kp and job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp] <=
		sum {i in ALL_NODES,t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t];

# for jobs completed by two techs, exactly one pair is assigned to this job
subject to c_pair_balance{j in JOBS}:
	2*var_twotech[j] = 
		sum{k in TECHS, kp in TECHS:k!=kp and job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp];

subject to c_pair_sym{j in JOBS, k in TECHS, kp in TECHS:job_tech2[j,k,kp]=1 and k<kp}:
	var_tech_pair[j,k,kp]=var_tech_pair[j,kp,k];


# ------------------------------------------------------------
# constraints relating to time windows and duration of jobs 
# ------------------------------------------------------------


# start time of each job is after earliest start
subject to c_start_time{j in JOBS}:
	var_start[j]>=
		start[j];

# start time of each job allows for completion of job before deadline
subject to c_finish_time{j in JOBS}:
	var_start[j] + (if duration[j,2]>=0 then duration[j,2]*(var_twotech[j]) else duration[j,1]*(1-var_twotech[j])) <= 
		finish[j] + lateness[j]*var_late[j];

# ------------------------------------------------------------

#start time of each job is at least after the tech has arrived there
#correct constraint
subject to c_tw_arrive{(i,j,k,t) in set_ARRIVE:j in JOBS}:
	var_start[j]>=
		t*var_arrive[i,j,k,t];

# relaxation without z
subject to c_tw_arrive_r{j in JOBS, t in HOURS}:
	var_start[j]>=
		t*(sum{i in ALL_NODES, k in TECHS:(i,j,k,t) in set_ARRIVE}var_arrive[i,j,k,t]);

# relaxation with z
subject to c_tw_arrive_rz{j in JOBS, t in HOURS}:
	var_start[j]>=
		t*(sum{i in ALL_NODES, k in TECHS:(i,j,k,t) in set_ARRIVE}var_arrive[i,j,k,t]-var_twotech[j]);

# ------------------------------------------------------------
#tech canonly leave a job once it has been completed
#correct constraint
subject to c_tw_leave{(j,h,k,t) in set_LEAVE:j in JOBS}:
	var_start[j]+(if duration[j,1]>=0 then duration[j,1]*(1-var_twotech[j]))+(if duration[j,2]>=0 then duration[j,2]*var_twotech[j])<=
		card(HOURS)+(t-card(HOURS))*var_leave[j,h,k,t];

# relaxation without z
subject to c_tw_leave_r{j in JOBS, t in HOURS}:
	var_start[j]+(if duration[j,1]>=0 then duration[j,1]*(1-var_twotech[j]))+(if duration[j,2]>=0 then duration[j,2]*var_twotech[j])<=
		card(HOURS)+(t-card(HOURS))*(sum{h in ALL_NODES, k in TECHS:(j,h,k,t) in set_LEAVE}var_leave[j,h,k,t]);

# relaxation with z
subject to c_tw_leave_rz{j in JOBS, t in HOURS}:
	var_start[j]+(if duration[j,1]>=0 then duration[j,1]*(1-var_twotech[j]))+(if duration[j,2]>=0 then duration[j,2]*var_twotech[j])<=
		card(HOURS)+(t-card(HOURS))*(sum{h in ALL_NODES, k in TECHS:(j,h,k,t) in set_LEAVE}var_leave[j,h,k,t]-var_twotech[j]);


# ------------------------------------------------------------

subject to c_return_home{k in TECHS,w in WEEKS}:
	sum{j in JOBS, t in HOURS_W[w]:(base[k],j,k,t) in set_LEAVE}var_leave[base[k],j,k,t]=
		1;


subject to c_return_home_new{k in TECHS}:
	sum{j in JOBS, t in HOURS:(base[k],j,k,t) in set_LEAVE}var_leave[base[k],j,k,t]=
		1;

# tech leaving job i at time t will arrive at j at time t+(travel time)
subject to c_arcs{(i,j,k,t) in set_LEAVE:t+travel[i,j]<= last(HOURS)}:
	var_leave[i,j,k,t] = 
		var_arrive[i,j,k,t+travel[i,j]];

# the total time spent in each job is equal to the duration of that job (by tracking variables)
subject to c_tracking_balance{j in JOBS}:
	sum{k in TECHS, t in HOURS:(j,k,t) in set_FEASIBLE_HOURS} var_tracking[j,k,t] = 
		duration[j,1]*(1-var_twotech[j]) + 2*duration[j,2]*var_twotech[j];

# the total time spent in each job by each technician is equal to the duration of that job
subject to c_tracking_tech_balance{j in JOBS, k in TECHS:(job_tech[j,k]=1 or sum{kp in TECHS:k!=kp}job_tech2[j,k,kp]>=1)}:
	sum{t in HOURS:(j,k,t) in set_FEASIBLE_HOURS} var_tracking[j,k,t]=
		duration[j,1]*(sum{i in JOBS,t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t]- sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp])+
		duration[j,2]*sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp];

# the total time spen} in each job is equal to the duration of that job (by hours_tech variables)
subject to c_hours_balance{j in JOBS}:
	#sum{k in TECHS, w in WEEKS, d in DAYS:(j,k,w) in set_FEASIBLE_WEEKS} 
	sum{k in TECHS, w in WEEKS, d in DAYS:(j,k,w,d) in set_FEASIBLE_DAYS} 
	var_hours_tech[j,k,w,d] = 
		duration[j,1]*(1-var_twotech[j]) + 2*duration[j,2]*var_twotech[j];

# the total time spent on each job by each technician is equal to the duration of that job
subject to c_hours_tech_balance{j in JOBS, k in TECHS:(job_tech[j,k]=1 or sum{kp in TECHS:k!=kp}job_tech2[j,k,kp]>=1)}:
	sum{w in WEEKS, d in DAYS:(j,k,w,d) in set_FEASIBLE_DAYS} var_hours_tech[j,k,w,d]=
		duration[j,1]*(sum{i in ALL_NODES,t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t]- sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp])+
		duration[j,2]*sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp];

# if a job is completed by a pair, the respective hours/tech variables must be the same for each day
subject to c_hours_tech_pair_1{j in JOBS,k in TECHS, kp in TECHS, w in WEEKS, d in DAYS:k!=kp and (j,k,w,d) in set_FEASIBLE_DAYS and (j,kp,w,d) in set_FEASIBLE_DAYS and job_tech2[j,k,kp]=1}:
	var_hours_tech[j,k,w,d]-var_hours_tech[j,kp,w,d]<=
		hours_week*(1-var_tech_pair[j,k,kp]);

subject to c_hours_tech_pair_2{j in JOBS,k in TECHS, kp in TECHS, w in WEEKS, d in DAYS:k!=kp and (j,k,w,d) in set_FEASIBLE_DAYS and (j,kp,w,d) in set_FEASIBLE_DAYS and job_tech2[j,k,kp]=1}:
	var_hours_tech[j,kp,w,d]-var_hours_tech[j,k,w,d]<=
		hours_week*(1-var_tech_pair[j,k,kp]);


# the start time of consecutive jobs must be in the correct order
# for this constraint, M_1 (as in the formulation) is the total number of hours in the planning period
subject to c_job_order{i in JOBS, j in JOBS}:
	var_start[i]+
#	duration[i,1]*(1-var_twotech[i]) + 2*duration[i,2]*var_twotech[i]+
	duration[i,1]*(1-var_twotech[i]) + duration[i,2]*var_twotech[i]+
	travel[i,j]<=
		var_start[j]+
#		card(HOURS)*(1-sum{k in TECHS,t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t]+var_twotech[j]);
		2*card(HOURS)*(1-sum{k in TECHS,t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t]+var_twotech[j]);

# the daily working capacity of each technician does not exceed their contracted hours
subject to c_day_capacity{k in TECHS, w in WEEKS, d in DAYS}:
	#sum{j in JOBS:(j,k,w) in set_FEASIBLE_WEEKS}
	sum{j in JOBS:(j,k,w,d) in set_FEASIBLE_DAYS}
	(sum{t in HOURS_WD[w,d], i in ALL_NODES:(i,j,k,t) in set_LEAVE}travel[i,j]*var_leave[i,j,k,t] + 
	var_hours_tech[j,k,w,d] + 
	sum{t in HOURS_WD[w,d]:(j,base[k],k,t) in set_LEAVE} travel[j,base[k]]*var_leave[j,base[k],k,t]) <=
		hours_day[d];

# ------------------------------------------------------------
# constraints relating to nights away and job splits 
# ------------------------------------------------------------

# lb on the number of nights away per job and tech
subject to c_night_away_lb{j in JOBS, k in TECHS:night_away[j,k]=1}:
	#sum{w in WEEKS, d in DAYS:(j,k,w) in set_FEASIBLE_WEEKS}var_night_away[j,k,w,d] >=
	sum{w in WEEKS, d in DAYS:(j,k,w,d) in set_FEASIBLE_DAYS}var_night_away[j,k,w,d] >=
		night_split[j,1,first(HOURS_W[1])]*
			(sum{i in ALL_NODES, t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t]-
			sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp])+
		night_split[j,2,first(HOURS_W[1])]*
			sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp];

# ub on the number of nights away per job and tech
subject to c_night_away_ub{j in JOBS, k in TECHS:night_away[j,k]=1}:
	#sum{w in WEEKS, d in DAYS:(j,k,w) in set_FEASIBLE_WEEKS}var_night_away[j,k,w,d] <=
	sum{w in WEEKS, d in DAYS:(j,k,w,d) in set_FEASIBLE_DAYS}var_night_away[j,k,w,d] <=
		night_split[j,1,last(HOURS_W[1])]*
			(sum{i in ALL_NODES, t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t]-
			sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp])+
		night_split[j,2,last(HOURS_W[1])]*
			sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp];

# lb on the number of night splits per job and tech
subject to c_job_split_lb{j in JOBS, k in TECHS}:
	#sum{w in WEEKS, d in DAYS:(j,k,w) in set_FEASIBLE_WEEKS and d!=last(DAYS)}var_job_split[j,k,w,d] >=
	sum{w in WEEKS, d in DAYS:(j,k,w,d) in set_FEASIBLE_DAYS and d!=last(DAYS)}var_job_split[j,k,w,d] >=
		night_split[j,1,first(HOURS_W[1])]*
			(sum{i in ALL_NODES, t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t]-
			sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp])+
		night_split[j,2,first(HOURS_W[1])]*
			sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp];

# ub on the number of night splits per job and tech
subject to c_job_split_ub{j in JOBS, k in TECHS}:
	#sum{w in WEEKS, d in DAYS:(j,k,w) in set_FEASIBLE_WEEKS and d!=last(DAYS)}var_job_split[j,k,w,d] <=
	sum{w in WEEKS, d in DAYS:(j,k,w,d) in set_FEASIBLE_DAYS and d!=last(DAYS)}var_job_split[j,k,w,d] <=
		night_split[j,1,last(HOURS_W[1])]*
			(sum{i in ALL_NODES, t in HOURS:(i,j,k,t) in set_LEAVE}var_leave[i,j,k,t]-
			sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp])+
		night_split[j,2,last(HOURS_W[1])]*
			sum{kp in TECHS:job_tech2[j,k,kp]=1}var_tech_pair[j,k,kp];

# if a job is started at time t, tracking variables must be one for the duration of the job
# this version uses the arrive variable
subject to c_consec_tracking_a{j in JOBS, k in TECHS, t in HOURS}:
	sum{tp in t..t+max(duration[j,1],duration[j,2]):(j,k,tp) in set_FEASIBLE_HOURS}var_tracking[j,k,tp]>=
		max(duration[j,1],duration[j,2])*sum{i in ALL_NODES:(i,j,k,t) in set_ARRIVE}var_arrive[i,j,k,t]-
		(max(duration[j,1],duration[j,2])-duration[j,2])*var_twotech[j];

# if a job is started at time t, tracking variables must be one for the duration of the job
# this version uses the leave variable
subject to c_consec_tracking_l{j in JOBS, k in TECHS, t in HOURS}:
	sum{tp in t..t+max(duration[j,1],duration[j,2]):(j,k,tp) in set_FEASIBLE_HOURS}var_tracking[j,k,tp]>=
		max(duration[j,1],duration[j,2])*sum{i in ALL_NODES:(i,j,k,t-travel[i,j]) in set_LEAVE}var_leave[i,j,k,t-travel[i,j]]-
		(max(duration[j,1],duration[j,2])-duration[j,2])*var_twotech[j];

# if a job is started at time t, job split variables must be one for the duration of the job
# this version uses the arrive variable
subject to c_consec_split_a{j in JOBS, k in TECHS, t in HOURS}:
	sum{d in day_t[t]..day_t[t]+max{n in TECH_COUNT}night_split[j,n,t]-1:
				(j,k,t) in set_FEASIBLE_HOURS and
				day_t[t]+max{n in TECH_COUNT}night_split[j,n,t]-1<=4 and
				(j,k,week_t[t],d) in set_FEASIBLE_DAYS
				}
	var_job_split[j,k,week_t[t],d]>=
		(max{n in TECH_COUNT}night_split[j,n,(if t mod hours_week>0 then t mod hours_week else hours_week)])*
		    	sum{i in ALL_NODES:(i,j,k,t) in set_ARRIVE}var_arrive[i,j,k,t]-
		((max{n in TECH_COUNT}night_split[j,n,(if t mod hours_week>0 then t mod hours_week else hours_week)])-
		     night_split[j,2,(if t mod hours_week >0 then t mod hours_week else hours_week)])*
		     	var_twotech[j];


# if a job is started at time t, job split variables must be one for the duration of the job
# this version uses the leave variable
subject to c_consec_split_l{j in JOBS, k in TECHS, t in HOURS}:
	sum{d in day_t[t]..day_t[t]+max{n in TECH_COUNT}night_split[j,n,t]-1:
				(j,k,t) in set_FEASIBLE_HOURS and
				day_t[t]+max{n in TECH_COUNT}night_split[j,n,t]-1<=4 and
				(j,k,week_t[t],d) in set_FEASIBLE_DAYS
				}
	var_job_split[j,k,week_t[t],d]>=
		(max{n in TECH_COUNT}night_split[j,n,(if t mod hours_week>0 then t mod hours_week else hours_week)])*
		    	sum{i in ALL_NODES:(i,j,k,t-travel[i,j]) in set_LEAVE}var_leave[i,j,k,t-travel[i,j]]-
		((max{n in TECH_COUNT}night_split[j,n,(if t mod hours_week>0 then t mod hours_week else hours_week)])-
		     night_split[j,2,(if t mod hours_week >0 then t mod hours_week else hours_week)])*
		     	var_twotech[j];

# if a job is started at time t, tracking variables must be one for the duration of the job
# this version uses the arrive variable
subject to c_consec_hours_a{j in JOBS, k in TECHS, t in HOURS}:
	sum{d in day_t[t]..day_t[t]+max{n in TECH_COUNT}night_split[j,n,t]-1:
				(j,k,t) in set_FEASIBLE_HOURS and
				day_t[t]+max{n in TECH_COUNT}night_split[j,n,t]-1<=4 and
				(j,k,week_t[t],d) in set_FEASIBLE_DAYS
				}
	var_hours_tech[j,k,week_t[t],d]>=
		max(duration[j,1],duration[j,2])*sum{i in ALL_NODES:(i,j,k,t) in set_ARRIVE}var_arrive[i,j,k,t]-
		(max(duration[j,1],duration[j,2])-duration[j,2])*var_twotech[j];

# if a job is started at time t, tracking variables must be one for the duration of the job
# this version uses the leave variable
subject to c_consec_hours_l{j in JOBS, k in TECHS, t in HOURS}:
	sum{d in day_t[t]..day_t[t]+max{n in TECH_COUNT}night_split[j,n,t]-1:
				(j,k,t) in set_FEASIBLE_HOURS and
				day_t[t]+max{n in TECH_COUNT}night_split[j,n,t]-1<=4 and
				(j,k,week_t[t],d) in set_FEASIBLE_DAYS
				}
	var_hours_tech[j,k,week_t[t],d]>=
		max(duration[j,1],duration[j,2])*sum{i in ALL_NODES:(i,j,k,t-travel[i,j]) in set_LEAVE}var_leave[i,j,k,t-travel[i,j]]-
		(max(duration[j,1],duration[j,2])-duration[j,2])*var_twotech[j];

#just to make the listing of constraints easier
subject to c_list_end:
	1=1;


printf "%s\n%s\n%s\n%s\n%s", "  ","------------------------------","end of constraints file","-------------------------------","  ">log_file;
