//********relax indices information (control)********
//relax variable collection. 0= relax; 1=don't relax
double relax_slack_variable1=0.0; //0 regard as no slack variable
double relax_slack_variable2=0.0; //regard as no slack variable
double relax_gen_variable=1.0;    //0 regard gen as known variables; 1 regard as unknown
double relax_candidateline_variable = 1.0;   //empty matrix means no candidate line; 1 means candidate lines exist

double load_change=1.0;           //(control)

//********given parameter********
// indices information
double line_Cfrom=1.0 * relax_candidateline_variable; // Candidate line info indices
double line_Cto=2*relax_candidateline_variable;
double line_Ccost=3*relax_candidateline_variable;
double line_Creactance=4*relax_candidateline_variable;
double line_Climit=5*relax_candidateline_variable;
double gen_busnum=1; // gen info indices (generation at bus #)
double gen_min=2;
double gen_max=3;
double gen_fixed=4;
double gen_cost=5;
double load_busnum=1;// load info indices
double load_fixed=2;

// basic number info
double nbus = BUS_NUM; //        ****(control)
double nCline = row_cline;
double nGen = row_gen;
double nload = row_load;
double nMLC = 1; //nMaxLineConnect ****(control)
double nCline_total = nCline*nMLC;
double nplan_year=10;// numbers of planning years ****(control)
//     Fkl_1            Xkl+Gk_ts           Gk       S1 S2 thata  I
double nAV=(nCline*nMLC)+(nCline*nMLC)*2+(nGen)+(nbus*3)+nCline_total; //numbers of  annual variable

//basic gen & load info
load_P=load_info(:,load_fixed)*load_change;
load_increase_factor=0.01;%                                       (control)
discount_rate=0.1; % discount rate
M=7*10^4;%                                                         (control)
M_kl=7*10^4;% candidate branch flow different M meanning same value
loss_penalty=10^4;
ref_position=1; % reference bus position
duration_time=8760; %year duration time (hour) 365*24*60
capacity_factor=0.6;% capacity factor of generator 50%
million_transfer=10^6;% transfer unit from dollar to million dollar
candidate_line_pool=[1:8]';
max_open_num=0;%nCline: all lines are controlable. 1: 1 line is controable (control)
