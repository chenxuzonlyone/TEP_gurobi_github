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
double *load_P;
Load_struct_read (load_P, Load_info, 'c' ,row_load, row_load, 0, load_fixed);// load at each bus
Array_coef_multiply(load_P, load_change, row_load); // load with coeff.

double load_increase_factor=0.01; // (control)
double discount_rate=0.1;       // discount rate
double M=7*10^4; //(control)
double M_kl=7*10^4;    //candidate branch flow different M meanning same value
double loss_penalty=10^4;
double ref_position=1; //reference bus position
double duration_time=8760; //year duration time (hour) 365*24*60
double capacity_factor=0.6; // capacity factor of generator 50%
double million_transfer=10^6; // transfer unit from dollar to million dollar

double *candidate_line_pool;
Array_ascend(candidate_line_pool, nCline);

