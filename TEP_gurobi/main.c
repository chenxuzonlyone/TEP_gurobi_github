#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gurobi_c.h"
#include "read.h"
#include <string.h>

#define COL_LINE 5
#define COL_GEN 5
#define COL_LOAD 2
#define COL_BUS 1
#define BUS_NUM 6
#define MAXSTR 128
#define DELIMINATE_NUM  10e5

int main()
{
    //********This part related to cline info********
    int row_cline, col_cline;
    
    //FILE *f_cline_stream = fopen("/Users/zhangcaihua/Desktop/TEP_gurobi/data/cline_info.csv", "r");/////////
    FILE *f_cline_stream = fopen("/Users/zhangcaihua/Desktop/TEP_gurobi_github/data/cline_info.csv", "r");
    file_size(f_cline_stream, &row_cline, &col_cline);// get rows and cols from file
    
    Cline_struct Cline_info;
    double*Cline_data_space;
    int space_delimiter_line[col_cline];// this is Mac version
    //int space_delimiter_line[COL_LINE]; // this is Windows version
    delimiter_array(space_delimiter_line, row_cline, col_cline);
    
    Cline_data_space = calloc(row_cline*col_cline, sizeof(double));
    Cline_info.line_Cfrom = &Cline_data_space[space_delimiter_line[0]];
    Cline_info.line_Cto = &Cline_data_space[space_delimiter_line[1]];
    Cline_info.line_Ccost = &Cline_data_space[space_delimiter_line[2]];
    Cline_info.line_Creactance = &Cline_data_space[space_delimiter_line[3]];
    Cline_info.line_Climit = &Cline_data_space[space_delimiter_line[4]];
    
    Cline_data_read(f_cline_stream, row_cline, col_cline, &Cline_info);//read data from file
    
//    for (int i = 0; i< row_cline; i++)
//    {
//        
//        printf("%f\t", Cline_info.line_Cfrom[i]);
//        printf("%f\t", Cline_info.line_Cto[i]);
//        printf("%f\t", Cline_info.line_Ccost[i]);
//        printf("%f\t", Cline_info.line_Creactance[i]);
//        printf("%f\t", Cline_info.line_Climit[i]);
//        printf("\n");
//    }
    //********Cline info ending********
    
    
    //********This part related to gen info********
    int row_gen, col_gen;
    FILE *f_gen_stream = fopen("/Users/zhangcaihua/Desktop/TEP_gurobi_github/data/gen_info.csv", "r");
    file_size(f_gen_stream, &row_gen, &col_gen);// get rows and cols from file
    
    Gen_struct Gen_info;
    double*Gen_data_space;
    int space_delimiter_gen[col_gen];//this is Mac version
    //int space_delimiter_gen[COL_GEN]; // this is Windows version
    delimiter_array(space_delimiter_gen, row_gen, col_gen);
    
    Gen_data_space = calloc(row_gen*col_gen, sizeof(double));
    Gen_info.gen_busnum = &Gen_data_space[space_delimiter_gen[0]];
    Gen_info.gen_min = &Gen_data_space[space_delimiter_gen[1]];
    Gen_info.gen_max = &Gen_data_space[space_delimiter_gen[2]];
    Gen_info.gen_fixed = &Gen_data_space[space_delimiter_gen[3]];
    Gen_info.gen_cost = &Gen_data_space[space_delimiter_gen[4]];
    
    Gen_data_read(f_gen_stream, row_gen, col_gen, &Gen_info);//read data from file
    
//    for (int i = 0; i< row_gen; i++)
//    {
//        
//        printf("%f\t", Gen_info.gen_busnum[i]);
//        printf("%f\t", Gen_info.gen_min[i]);
//        printf("%f\t", Gen_info.gen_max[i]);
//        printf("%f\t", Gen_info.gen_fixed[i]);
//        printf("%f\t", Gen_info.gen_cost[i]);
//        printf("\n");
//    }
    //********gen info ending********
    
    //********This part related to load info********
    int row_load, col_load;
    FILE *f_load_stream = fopen("/Users/zhangcaihua/Desktop/TEP_gurobi_github/data/load_info.csv", "r");
    file_size(f_load_stream, &row_load, &col_load);// get rows and cols from file
    
    Load_struct Load_info;
    double*Load_data_space;
    int space_delimiter_load[col_load]; // this is Mac version
    //int space_delimiter_load[COL_LOAD];// this is Windows version
    delimiter_array(space_delimiter_load, row_load, col_load);
    
    Load_data_space = calloc(row_load*col_load, sizeof(double));
    Load_info.load_busnum= &Load_data_space[space_delimiter_load[0]];
    Load_info.load_fixed = &Load_data_space[space_delimiter_load[1]];
    
    Load_data_read(f_load_stream, row_load, col_load, &Load_info);//read data from file
    
//    for (int i = 0; i< row_load; i++)
//    {
//        
//        printf("%f\t", Load_info.load_busnum[i]);
//        printf("%f\t", Load_info.load_fixed[i]);
//        printf("\n");
//    }
    //********load info ending********
    
    //********This part related to number setting********
    
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
    double nVariable_total = nplan_year*nAV;
    
    //basic gen & load info
    double load_P[(int)nload];
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
    double million_transfer=pow(10, 6); // transfer unit from dollar to million dollar
    
    double candidate_line_pool[(int)nCline];
    Array_ascend(candidate_line_pool, nCline);
    /********This is end related to number setting********/
    
    //********parameter calculation********
    //bus related incidence matrix
    
    // bus branch incidence matrix
    double Kl_C[(int)(nbus*nCline)];
    Array_initial(Kl_C, nbus*nCline); //initialize array
    //candidate
    Kl_C_set(Kl_C, Cline_info, nbus, nCline);
    
    
   
    // bus generation incidence matrix
    double Kp[(int)(nbus*nGen)];
    Array_initial(Kp, nbus*nGen); //initialize array
    Kp_set(Kp, Gen_info, nbus, nGen);
    
    // bus pseudo generator incidence matrix (ts: transmission switching)
    double Kp_ts[(int)(nbus*nCline)];
    Array_initial(Kp_ts, nbus*nCline); //initialize array
    Kp_ts_set(Kp_ts, Cline_info, nbus, nCline);
    
    // bus load incidence matrix
    double Kd[(int)(nbus*nload)];
    Array_initial(Kp_ts, nbus*nload); //initialize array
    Kd_set(Kd, Load_info, nbus, nload);

    // bus slack variable incidence matrix
    double Kr1[(int)(nbus*nbus)];
    Array_initial(Kr1, nbus*nbus);
    double Kr2[(int)(nbus*nbus)];
    Array_initial(Kr2, nbus*nbus);
    eye(Kr1, nbus);
    eye(Kr2, nbus);
    
    // bus angle(theta) incidence matrix
    double Ka[(int)(nbus*nbus)];
    Array_initial(Ka, nbus*nbus);
    eye(Ka, nbus);
    
    // branch reactance matrix
    double X_C[(int)(nCline*nCline)];
    Array_initial(X_C, nCline*nCline);
    
    double resist_array[(int)nCline];
    Array_initial(resist_array, nCline);
    
    Cline_struct_read (resist_array, Cline_info, 'c', row_cline, col_cline, 0, line_Creactance);
    Array_coef_multiply(load_P, load_change, row_load); // load with coeff.
    
    diag(X_C, nCline, resist_array);
    //********end of parameter calculation********
    
    
    //********This part deal with Gurobi ********
    
    //****This part deal with objective function setting****
    double f[(int)nVariable_total];
    Array_initial(f, nVariable_total);
    double lb[(int)nVariable_total];
    Array_initial(lb, nVariable_total);
    double ub[(int)nVariable_total];
    Array_initial(ub, nVariable_total);
    
    //set the candidate line status
    double line_initial[(int)nCline_total];
    ones(line_initial, nCline_total); // regard all lines are built already
    for (int i = 0; i < nCline_total; i++) { // choose the candidate lines
        int exist_line = (int)candidate_line_pool[i];
        line_initial[exist_line-1]=0; // Due to the array strats from 0 to (Max-1)
    }
//    for (int i =0 ; i < nCline_total; i++) {
//        printf("line_init_status %f\n",line_initial[i]);
//    }
    
    
    //line_initial=ones(nCline_total,1);% In fact, this constraint is used to set the limit of binary variable
    //line_initial(candidate_line_pool)=0;% the meaning is to let the candidate line binary variable limit from 0 to 1
    
    //Array_initial(f, nVariable_total);
    Array_initial(f, nVariable_total);
    int var_pt = 0;
    
    for (int i_year = 1; i_year <= nplan_year; i_year++) {
        for (int i = 0; i<nCline_total;i++){ //flowC_F_t
            f[var_pt] = 0;//Obj. func.
            lb[var_pt] = -GRB_INFINITY; //lower bound
            ub[var_pt] =  GRB_INFINITY; //upper bound
            var_pt=var_pt+1;
        }
        
        
        for (int i = 0; i<nCline_total; i++) { //x_F_t
            f[var_pt] = 0;
            lb[var_pt] = line_initial[i]; //lower bound (existing lines' lower bound set to 1)
            ub[var_pt] =  1; //upper bound
            var_pt=var_pt+1;
        }
        
        double gen_discount_coeff=(capacity_factor*duration_time)/(pow(1+discount_rate, (double)(i_year-1))*million_transfer);
        for (int i = 0; i<nGen; i++) {//gen_F_t
            f[var_pt] = Gen_info.gen_cost[i] * gen_discount_coeff;
            lb[var_pt] = 0; //lower bound
            ub[var_pt] = GRB_INFINITY; //upper bound
            var_pt=var_pt+1;
        }
        
        double s_discount_coeff =(loss_penalty*duration_time)/(pow(1+discount_rate, (double)(i_year-1))*million_transfer);
        for (int i = 0; i<nbus; i++) {//s1_F_t
            f[var_pt] =  s_discount_coeff * relax_slack_variable1;
            lb[var_pt] = 0; //lower bound
            ub[var_pt] = GRB_INFINITY; //upper bound
            var_pt=var_pt+1;
        }
        for (int i = 0; i<nbus; i++) {//s2_F_t
            f[var_pt] =  s_discount_coeff * relax_slack_variable1;
            lb[var_pt] = 0; //lower bound
            ub[var_pt] = GRB_INFINITY; //upper bound
            var_pt=var_pt+1;
        }
        
        for (int i = 0; i<nbus; i++) {//theta_F_t
            f[var_pt] =  0;
        
            if (i == (int)(ref_position-1)) { //select the reference bus
                lb[var_pt] = 0; //lower bound
                ub[var_pt] = 0; //upper bound
            }
            lb[var_pt] = -GRB_INFINITY; //lower bound
            ub[var_pt] =  GRB_INFINITY; //upper bound
            var_pt=var_pt+1;
        }
        
        for (int i = 0; i<nCline_total; i++) {//Kp_ts_F_t (pseudo generator)
            f[var_pt] =  0;
            lb[var_pt] = -GRB_INFINITY; //lower bound
            ub[var_pt] =  GRB_INFINITY; //upper bound
            var_pt=var_pt+1;
        }
        
        for (int i = 0; i<nCline_total; i++) {//I_F_t
            f[var_pt] =  Cline_info.line_Ccost[i]/(pow(1+discount_rate, (double)(i_year-1)));
            lb[var_pt] = 0; //lower bound
            ub[var_pt] = 1; //upper bound
            var_pt=var_pt+1;
        }
    }

    //testing
    //    for(int i =0;i<f_pt;i++){
    //        printf("i %d \t",i+1);
    //        printf("%f \n", f[i]);
    //    }
    

    //****This part deal with Gurobi setting****
    GRBenv *env = NULL;
    GRBmodel *model = NULL;
    int error = 0;
    char vname[MAXSTR];
    
    // Create environment
    error = GRBloadenv(&env, "TEP_GUROBI.log");
    //if (error) goto QUIT;
    
    // Create initial model
    error = GRBnewmodel(env, &model, "TEP_GUROBI", (int)nVariable_total,
                        NULL, NULL, NULL, NULL, NULL);
    //if (error) goto QUIT;
    
    //Update model
    error = GRBupdatemodel(model);
    //if (error) goto QUIT;
    
    // ******** This will set the TYPE, COEFFICIENT, and UPPER,LOWER BOUND of all variable ********
    for (int i_total = 0; i_total < nVariable_total; )
    {
        //Fkl_1
        for (int i_C_Cline = 0; i_C_Cline<nCline_total; i_C_Cline++) {
            error = GRBsetcharattrelement(model, "VType", i_total, GRB_CONTINUOUS);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "Obj", i_total, f[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "LB", i_total, lb[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "UB", i_total, ub[i_total]);
            //if (error) goto QUIT;
            
            sprintf(vname, "Fkl_%i", i_total);
            error = GRBsetstrattrelement(model, "VarName", i_total, vname);
            //if (error) goto QUIT;
            i_total++;
        }
        //printf("i_total %d \n", i_total);
        //Xkl
        for (int i_B_Xkl = 0; i_B_Xkl<nCline_total; i_B_Xkl++) {
            error = GRBsetcharattrelement(model, "VType", i_total, GRB_BINARY);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "Obj", i_total, f[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "LB", i_total, lb[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "UB", i_total, ub[i_total]);
            //if (error) goto QUIT;
            
            sprintf(vname, "Xkl_%i", i_total);
            error = GRBsetstrattrelement(model, "VarName", i_total, vname);
            //if (error) goto QUIT;
            i_total++;
        }
        //printf("i_total %d \n", i_total);
        //Gk
        for (int i_C_Gk = 0; i_C_Gk<nGen; i_C_Gk++) {
            error = GRBsetcharattrelement(model, "VType", i_total, GRB_CONTINUOUS);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "Obj", i_total, f[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "LB", i_total, lb[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "UB", i_total, ub[i_total]);
            //if (error) goto QUIT;
            
            sprintf(vname, "Gkl_%i", i_total);
            error = GRBsetstrattrelement(model, "VarName", i_total, vname);
            //if (error) goto QUIT;
            i_total++;
        }
        //printf("i_total %d \n", i_total);
        //Rk1
        for (int i_C_Rk1 = 0; i_C_Rk1<nbus; i_C_Rk1++) {
            error = GRBsetcharattrelement(model, "VType", i_total, GRB_CONTINUOUS);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "Obj", i_total, f[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "LB", i_total, lb[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "UB", i_total, ub[i_total]);
            //if (error) goto QUIT;
            
            sprintf(vname, "Rk1_%i", i_total);
            error = GRBsetstrattrelement(model, "VarName", i_total, vname);
            //if (error) goto QUIT;
            i_total++;
        }
        //printf("i_total %d \n", i_total);
        //Rk2
        for (int i_C_Rk2 = 0; i_C_Rk2<nbus; i_C_Rk2++) {
            error = GRBsetcharattrelement(model, "VType", i_total, GRB_CONTINUOUS);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "Obj", i_total, f[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "LB", i_total, lb[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "UB", i_total, ub[i_total]);
            //if (error) goto QUIT;
            
            sprintf(vname, "Rk1_%i", i_total);
            error = GRBsetstrattrelement(model, "VarName", i_total, vname);
            //if (error) goto QUIT;
            i_total++;
        }
        //printf("i_total %d \n", i_total);
        //theta
        for (int i_C_theta = 0; i_C_theta<nbus; i_C_theta++) {
            error = GRBsetcharattrelement(model, "VType", i_total, GRB_CONTINUOUS);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "Obj", i_total, f[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "LB", i_total, lb[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "UB", i_total, ub[i_total]);
            //if (error) goto QUIT;
            
            sprintf(vname, "Rk2_%i", i_total);
            error = GRBsetstrattrelement(model, "VarName", i_total, vname);
            //if (error) goto QUIT;
            i_total++;
        }
        //printf("i_total %d \n", i_total);
        //Gk_ts
        for (int i_C_Gk_ts = 0; i_C_Gk_ts<nCline_total; i_C_Gk_ts++) {
            error = GRBsetcharattrelement(model, "VType", i_total, GRB_CONTINUOUS);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "Obj", i_total, f[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "LB", i_total, lb[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "UB", i_total, ub[i_total]);
            //if (error) goto QUIT;
            
            sprintf(vname, "Gk_ts_%i", i_total);
            error = GRBsetstrattrelement(model, "VarName", i_total, vname);
            //if (error) goto QUIT;
            i_total++;
        }
        //printf("i_total %d \n", i_total);
        //I
        for (int i_B_I = 0; i_B_I<nCline_total; i_B_I++) {
            error = GRBsetcharattrelement(model, "VType", i_total, GRB_BINARY);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "Obj", i_total, f[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "LB", i_total, lb[i_total]);
            //if (error) goto QUIT;
            error = GRBsetdblattrelement(model, "UB", i_total, ub[i_total]);
            //if (error) goto QUIT;
            
            sprintf(vname, "I_%i", i_total);
            error = GRBsetstrattrelement(model, "VarName", i_total, vname);
            //if (error) goto QUIT;
            i_total++;
        }
        //printf("i_total %d \n", i_total);

    }
    
    //Update model due to lazy model update strategy
    error = GRBupdatemodel(model);
    //if (error) goto QUIT;
    GRBwrite (model, "groubi_obj.lp" );
    GRBwrite (model, "groubi_obj.rlp" );
    
//    for (int i=0; i<nVariable_total; i++) {
//         error = GRBgetdblattr(model, "", &objcon);
//    }
    
    
    // ******** The end of set the TYPE, COEFFICIENT, and UPPER,LOWER BOUND of all variable ********
    
    // ******** The begining of set Nodal Balance constraints ********
    //For the A part of Ax<=b
    double non_zero_num = 0;
    double p_position[BUFFER_MAX];
    double n_position[BUFFER_MAX];
    double p_value[BUFFER_MAX];
    double n_value[BUFFER_MAX];
    ones(p_position, BUFFER_MAX);
    ones(n_position, BUFFER_MAX);
    Array_initial(p_value, BUFFER_MAX);
    Array_initial(n_value, BUFFER_MAX);
    Array_coef_multiply(p_position, DELIMINATE_NUM, BUFFER_MAX);//a huge number to tell the different
    Array_coef_multiply(n_position, DELIMINATE_NUM, BUFFER_MAX);
    int p_pt = 0;
    int n_pt = 0;
    int position_pt = 0;
    
    //Fkl_1
    for (int i = 0; i < nCline; i++) {
        //printf("KL_c %f\n",Kl_C[(int)(0+i*nCline)]);
        if (Kl_C[(int)(i + 0)] > 0) { //0+(i-1): the 0 can be used as variable
            ++non_zero_num;
            p_position[p_pt] = position_pt;
            p_value[p_pt] = Kl_C[(int)(i + 0)];
            ++p_pt;
        }
        if (Kl_C[(int)(i + 0)] < 0) { //0+(i-1): the 0 can be used as variable
            ++non_zero_num;
            n_position[n_pt] = position_pt;
            n_value[n_pt] = Kl_C[(int)(i + 0)];
            ++n_pt;
        }
        ++position_pt;
    }
    
    //Xkl
    position_pt = position_pt + nCline;
    
    //Gk
    for (int i = 0; i < nGen; i++) {
        
        if (-(Kp[(int)(i + 0)]) > 0) { //0+(i-1): the 0 can be used as variable
            ++non_zero_num;
            p_position[p_pt] = position_pt;
            p_value[p_pt] = -(Kp[(int)(i + 0)]);
            ++p_pt;
        }
        if (-(Kp[(int)(i + 0)]) < 0) { //0+(i-1): the 0 can be used as variable
            ++non_zero_num;
            n_position[n_pt] = position_pt;
            n_value[n_pt] = -(Kp[(int)(i + 0)]);
            ++n_pt;
        }
        ++position_pt;
    }
    
    //Rk1
    position_pt = position_pt + nbus;
    
    //Rk2
    position_pt = position_pt + nbus;
    
    //theta
    position_pt = position_pt + nbus;
    
    //Gk_ts (it is the reverse of candidate line flow)
    for (int i = 0; i < nCline; i++) {
        if (-(Kl_C[(int)(i + 0)]) > 0) { //0+(i-1): the 0 can be used as variable
            ++non_zero_num;
            p_position[p_pt] = position_pt;
            p_value[p_pt] = -(Kl_C[(int)(i + 0)]);
            ++p_pt;
        }
        if (-(Kl_C[(int)(i + 0)]) < 0) { //0+(i-1): the 0 can be used as variable
            ++non_zero_num;
            n_position[n_pt] = position_pt;
            n_value[n_pt] = -(Kl_C[(int)(i + 0)]);
            ++n_pt;
        }
        ++position_pt;
    }
    
    
    //I
    position_pt = position_pt + nCline;
    
    
    int ind_t[(int)non_zero_num];
    double val_t[(int)non_zero_num];
    Array_initial_int(ind_t, non_zero_num);
    Array_initial(val_t, non_zero_num);
    int ind_pt = 0;
    
    for (int i = 0; i < (int) p_pt; i++) { //positive coefficient
        ind_t[ind_pt] = p_position[i];
        val_t[ind_pt] = p_value[i];
        ++ind_pt;
    }
    for (int i = 0; i < (int) n_pt; i++) { //negative coefficient
        ind_t[ind_pt] = n_position[i];
        val_t[ind_pt] = n_value[i];
        ++ind_pt;
    }
    for (int i=0; i<7; i++) {
        printf("val %f\t",val_t[i]);
        printf("ind %d\n",ind_t[i]);
    }
    
    //For the b part of Ax<=b
    double load_increase_coef = pow((1+load_increase_factor), 0); // This will be changed due to the years
    double load_Pt[(int)nload];
    double load_NB_t[(int)nload];
    Array_initial(load_NB_t, nload);
    
    Array_coef_multiply(load_Pt, load_increase_coef, nload);
    

    
    
    
    char constraint_name[MAXSTR];
    sprintf(constraint_name, "%i_constraint", 1);
    // Add a constraint
    error = GRBaddconstr(model, (int)non_zero_num, ind_t, val_t, GRB_LESS_EQUAL, load_Pt[0], constraint_name);
    //Update model due to lazy model update strategy
    error = GRBupdatemodel(model);
    //if (error) goto QUIT;
    GRBwrite (model, "groubi_obj.lp" );
    GRBwrite (model, "groubi_obj.rlp" );
    
    
    /* Solve */
    error = GRBoptimize(model);
    if (error) goto QUIT;
    
//********This part will deal with code Finish or Crash********
QUIT:
    
    // Error reporting
    
    if (error)
    {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }
    
    // Free model
    
    GRBfreemodel(model);
    
    // Free environment
    
    GRBfreeenv(env);
    
    //********This part will close file and clean memory********
    fclose(f_cline_stream);
    free(Cline_data_space);
    
    fclose(f_gen_stream);
    free(Gen_data_space);
    
    fclose(f_load_stream);
    free(Load_data_space);
    return 0;
}
