//#include<stdio.h>
//#include<math.h>
//#include "read.h"
//int inverse_matrix_size_global ; // To avoid multiple linker definitions, just one declaration of your global symbol must be present across your compilation units
//
//int main()
//{
//    //double a[9]={2,3,5,2,1,5,8,3,9};
//    double a[16]={9,3,4,6,4,2,3,5,9,8,6,5,3,4,6,9};
//    double b[16];
//    double k, d;
//
//    //k = 3.0;
//    inverse_matrix_size_global = 4; // should be equal to the value of k (but after setting, it should not be changed)
//    k = 4.0; // size of inverse matrix (not global and can be change)
//
//    d = determinant(a, k);
//    if (d == 0)
//        printf("\nInverse of Entered Matrix is not possible\n");
//    else
//        cofactor(a, k, b);
//
//    Matrix_display(b, k, k);
//}

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
int inverse_matrix_size_global ; // To avoid multiple linker definitions, just one declaration of your global symbol must be present across your compilation units

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
    double nplan_year=5;// numbers of planning years ****(control)
    //     Fkl_1            Xkl+Gk_ts           Gk       S1 S2 thata  I
    double nAV=(nCline*nMLC)+(nCline*nMLC)*2+(nGen)+(nbus*3)+nCline_total; //numbers of  annual variable
    double nVariable_total = nplan_year*nAV;
    
    //basic gen & load info
    double load_P[(int)nload];
    Load_struct_read (load_P, Load_info, 'c' ,row_load, row_load, 0, load_fixed);// load at each bus
    Array_coef_multiply(load_P, load_change, row_load); // load with coeff.
    
    double load_increase_factor=0.01; // (control)
    double discount_rate=0.1;       // discount rate
    double M=7*pow(10, 4); //(control)
    double M_kl=7*pow(10, 4);    //candidate branch flow different M meanning same value
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
    Array_initial(Kd, nbus*nload); //initialize array
    Kd_set(Kd, Load_info, nbus, nload);
    
    //    for (int i = 0; i<6; i++) {
    //        for (int j = 0; j<6;j++) {
    //            printf("%f \t", Kd[i*6+j]);
    //        }
    //        printf("\n");
    //    }
    
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
    for (int NB_year = 0; NB_year < (int)nplan_year; NB_year++) {
        
        //**** One year constraints ****
        for (int NB_bus_t = 0; NB_bus_t < nbus; NB_bus_t++) { // this will add number of bus constraints
            
            //For the A part of Ax<=b
            double non_zero_num = 0;
            double p_NB_position[BUFFER_MAX];
            double n_NB_position[BUFFER_MAX];
            double p_NB_value[BUFFER_MAX];
            double n_NB_value[BUFFER_MAX];
            
            Array_initial(p_NB_value, BUFFER_MAX);
            Array_initial(n_NB_value, BUFFER_MAX);
            Array_initial(p_NB_position, BUFFER_MAX);
            Array_initial(n_NB_position, BUFFER_MAX);
            
            int p_pt = 0;
            int n_pt = 0;
            int position_pt = 0;
            //Fkl_1
            for (int i = 0; i < nCline; i++) {
                //printf("KL_c %f\n",Kl_C[(int)(0+i*nCline)]);
                if (Kl_C[(int)(i + NB_bus_t*nCline)] > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    p_NB_position[p_pt] = position_pt;
                    p_NB_value[p_pt] = Kl_C[(int)(i + NB_bus_t*nCline)];
                    ++p_pt;
                }
                if (Kl_C[(int)(i + + NB_bus_t*nCline)] < 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    n_NB_position[n_pt] = position_pt;
                    n_NB_value[n_pt] = Kl_C[(int)(i + NB_bus_t*nCline)];
                    ++n_pt;
                }
                ++position_pt;
            }
            
            //Xkl
            position_pt = position_pt + nCline;
            
            //Gk
            for (int i = 0; i < nGen; i++) {
                
                if (-(Kp[(int)(i + NB_bus_t*nGen)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    p_NB_position[p_pt] = position_pt;
                    p_NB_value[p_pt] = -(Kp[(int)(i + NB_bus_t*nGen)]);
                    ++p_pt;
                }
                if (-(Kp[(int)(i + NB_bus_t*nGen)]) < 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    n_NB_position[n_pt] = position_pt;
                    n_NB_value[n_pt] = -(Kp[(int)(i + NB_bus_t*nGen)]);
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
                if (-(Kl_C[(int)(i + NB_bus_t*nCline)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    p_NB_position[p_pt] = position_pt;
                    p_NB_value[p_pt] = -(Kl_C[(int)(i + NB_bus_t*nCline)]);
                    ++p_pt;
                }
                if (-(Kl_C[(int)(i + NB_bus_t*nCline)]) < 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    n_NB_position[n_pt] = position_pt;
                    n_NB_value[n_pt] = -(Kl_C[(int)(i + NB_bus_t*nCline)]);
                    ++n_pt;
                }
                ++position_pt;
            }
            
            //I
            position_pt = position_pt + nCline_total;
            
            //Set positive & negative coefficient (ind_t) and positive & negative position(val_t)
            int ind_t[(int)non_zero_num];
            double val_t[(int)non_zero_num];
            Array_initial_int(ind_t, non_zero_num);
            Array_initial(val_t, non_zero_num);
            int ind_pt = 0;
            
            for (int i = 0; i < (int) p_pt; i++) { //positive coefficient
                ind_t[ind_pt] = p_NB_position[i] + (int)nAV*NB_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                //printf("%d\n",ind_pt + (int)nAV*NB_year);
                val_t[ind_pt] = p_NB_value[i];// But the value of cofficients are not change
                ++ind_pt;
            }
            for (int i = 0; i < (int) n_pt; i++) { //negative coefficient
                ind_t[ind_pt] = n_NB_position[i]  + (int)nAV*NB_year;// The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                val_t[ind_pt] = n_NB_value[i];// But the value of cofficients are not change
                ++ind_pt;
            }
            //        printf("ind_pt number: %d\n", ind_pt);
            //        for (int i = 0; i< ind_pt; i++) {
            //            printf("%d\n",ind_t[i]);
            //        }
            
            //For the b part of Ax<=b
            double load_increase_coef = pow((1+load_increase_factor), NB_year); //Loads are changing according to the years change.
            //printf("load_increase_coef %f \n",load_increase_coef);
            double load_Pt[(int)nload];
            double load_NB_t[(int)nbus];
            Array_initial(load_Pt, nload);
            Array_initial(load_NB_t, nload);
            
            Matrix_multiply(load_Pt, load_P, nload, 1.0, &load_increase_coef, 1.0, 1.0);//load_Pt=load_P*load_increase_coef
            
            Matrix_multiply(load_NB_t, Kd, nbus, nload, load_Pt, nload, 1.0);//load_NB_t=Kd*load_Pt
            
            
            //char constraint_name[MAXSTR];
            //sprintf(constraint_name, "%i_constraint", NB_bus_t+1);
            // Add a constraint
            //printf("load_NB_t %f\n", load_NB_t[NB_bus_t]);
            //printf("non_zero %f\n", non_zero_num);
            
            //*****************************
            error = GRBaddconstr(model, (int)non_zero_num, ind_t, val_t, GRB_EQUAL, -load_NB_t[NB_bus_t], NULL);//load is negative value at this moment
            //*****************************
            //Update model due to lazy model update strategy
            
            error = GRBupdatemodel(model);
            //if (error) goto QUIT;
            GRBwrite (model, "groubi_obj.lp" );
            GRBwrite (model, "groubi_obj.rlp" );
        }// NODAL BALANCE CONSTRAINT: Each year constraints are added
        //printf("\n");
        
        //****The ending of One year constraints ****
    }// NODAL BALANCE CONSTRAINT: Total planning years constraints are added
    // ******** The ending of set Nodal Balance constraints ********
    
    
    // ******** The begining of set STATUS CHANGE constraints ********
    for (int SC_year = 0; SC_year < (int)nplan_year; SC_year++) {
        
        //**** One year constraints ****
        for (int SC_nCline_t = 0; SC_nCline_t < nCline_total; SC_nCline_t++) { // this will add number of bus constraints
            
            //For the A part of Ax<=b
            double non_zero_num = 0;
            double p_SC_position[BUFFER_MAX];
            double p_SC_value[BUFFER_MAX];
            
            Array_initial(p_SC_value, BUFFER_MAX);
            Array_initial(p_SC_position, BUFFER_MAX);
            
            int p_pt = 0;
            int position_pt = 0;
            //Fkl_1
            position_pt = position_pt + nCline_total;
            
            //Xkl
            double x_SC_C_t[(int)(nCline_total * nCline_total)];
            double x_SC_P_t[(int)(nCline_total * nCline_total)];
            Array_initial(x_SC_C_t,nCline_total * nCline_total);
            Array_initial(x_SC_P_t,nCline_total * nCline_total);
            eye(x_SC_C_t, nCline_total);
            eye(x_SC_P_t, nCline_total);
            
            for (int i = 0; i < nCline_total; i++) {
                if ((x_SC_C_t[(int)(i + SC_nCline_t*nCline_total)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    p_SC_position[p_pt] = position_pt;
                    p_SC_value[p_pt] = (x_SC_C_t[(int)(i + SC_nCline_t*nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++p_pt;
                    if (SC_year > 0) { // Add the PREVIOUS status back into the constraints
                        ++non_zero_num;
                        p_SC_position[p_pt] = position_pt - nAV;
                        p_SC_value[p_pt] = -(x_SC_C_t[(int)(i + SC_nCline_t*nCline_total)]);
                        //printf("%f\n", p_SC_position[p_pt]);
                        ++p_pt;
                    }
                    
                }
                ++position_pt;
            }
            
            //Gk
            position_pt = position_pt + nGen;
            
            //Rk1
            position_pt = position_pt + nbus;
            
            //Rk2
            position_pt = position_pt + nbus;
            
            //theta
            position_pt = position_pt + nbus;
            
            //Gk_ts (it is the reverse of candidate line flow)
            position_pt = position_pt + nCline_total;
            
            //I
            for (int i = 0; i < nCline_total; i++) {
                if ((x_SC_C_t[(int)(i + SC_nCline_t*nCline_total)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    p_SC_position[p_pt] = position_pt;
                    p_SC_value[p_pt] = -(x_SC_C_t[(int)(i + SC_nCline_t*nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++p_pt;
                }
                ++position_pt;
            }
            
            //Set positive & negative coefficient (ind_t) and positive & negative position(val_t)
            int ind_t[(int)non_zero_num];
            double val_t[(int)non_zero_num];
            Array_initial_int(ind_t, non_zero_num);
            Array_initial(val_t, non_zero_num);
            int ind_pt = 0;
            
            for (int i = 0; i < (int) p_pt; i++) { //positive coefficient
                ind_t[ind_pt] = p_SC_position[i] + (int)nAV*SC_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                //printf("%d\n",ind_pt + (int)nAV*NB_year);
                val_t[ind_pt] = p_SC_value[i];// But the value of cofficients are not change
                ++ind_pt;
            }
            //                    printf("ind_pt number: %d\n", ind_pt);
            //                    for (int i = 0; i< ind_pt; i++) {
            //                        printf("%d\n",ind_t[i]);
            //                    }
            
            //For the b part of Ax<=b
            // The elements in b are always "0"
            
            //Add constraint
            //*****************************
            error = GRBaddconstr(model, (int)non_zero_num, ind_t, val_t, GRB_EQUAL, 0.0, NULL);//load is negative value at this moment
            //*****************************
            
            //Update model due to lazy model update strategy
            error = GRBupdatemodel(model);
            //if (error) goto QUIT;
            GRBwrite (model, "groubi_obj.lp" );
            GRBwrite (model, "groubi_obj.rlp" );
        }// STATUS CHANGE CONSTRAINTS: Each year constraints are added
        //printf("\n");
        //****The ending of One year constraints ****
    }// STATUS CHANGE CONSTRAINTS: Total planning years constraints are added
    // ******** The ending of set STATUS CHANGE constraints ********
    
    
    // ******** The begining of set INSTALLATION STATUS MAINTAIN constraints ********
    for (int SM_year = 0; SM_year < (int)nplan_year; SM_year++) {
        
        //**** One year constraints ****
        for (int SM_nCline_t = 0; SM_nCline_t < nCline_total; SM_nCline_t++) { // this will add number of bus constraints
            
            //For the A part of Ax<=b
            double non_zero_num = 0;
            double p_SM_position[BUFFER_MAX];
            double p_SM_value[BUFFER_MAX];
            
            Array_initial(p_SM_value, BUFFER_MAX);
            Array_initial(p_SM_position, BUFFER_MAX);
            
            int p_pt = 0;
            int position_pt = 0;
            //Fkl_1
            position_pt = position_pt + nCline_total;
            
            //Xkl
            double x_SM_C_t[(int)(nCline_total * nCline_total)];
            double x_SM_P_t[(int)(nCline_total * nCline_total)];
            Array_initial(x_SM_C_t,nCline_total * nCline_total);
            Array_initial(x_SM_P_t,nCline_total * nCline_total);
            eye(x_SM_C_t, nCline_total);
            eye(x_SM_P_t, nCline_total);
            
            for (int i = 0; i < nCline_total; i++) {
                if ((x_SM_C_t[(int)(i + SM_nCline_t*nCline_total)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    p_SM_position[p_pt] = position_pt;
                    p_SM_value[p_pt] = -(x_SM_C_t[(int)(i + SM_nCline_t*nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++p_pt;
                    if (SM_year > 0) { // Add the PREVIOUS status back into the constraints
                        ++non_zero_num;
                        p_SM_position[p_pt] = position_pt - nAV;
                        p_SM_value[p_pt] =  (x_SM_C_t[(int)(i + SM_nCline_t*nCline_total)]);
                        //printf("%f\n", p_SC_position[p_pt]);
                        ++p_pt;
                    }
                    
                }
                ++position_pt;
            }
            
            //Gk
            position_pt = position_pt + nGen;
            
            //Rk1
            position_pt = position_pt + nbus;
            
            //Rk2
            position_pt = position_pt + nbus;
            
            //theta
            position_pt = position_pt + nbus;
            
            //Gk_ts (it is the reverse of candidate line flow)
            position_pt = position_pt + nCline_total;
            
            //I
            position_pt = position_pt + nCline_total;
            
            
            //Set positive & negative coefficient (ind_t) and positive & negative position(val_t)
            int ind_t[(int)non_zero_num];
            double val_t[(int)non_zero_num];
            Array_initial_int(ind_t, non_zero_num);
            Array_initial(val_t, non_zero_num);
            int ind_pt = 0;
            
            for (int i = 0; i < (int) p_pt; i++) { //positive coefficient
                ind_t[ind_pt] = p_SM_position[i] + (int)nAV*SM_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                //printf("%d\n",ind_pt + (int)nAV*NB_year);
                val_t[ind_pt] = p_SM_value[i];// But the value of cofficients are not change
                ++ind_pt;
            }
            //                    printf("ind_pt number: %d\n", ind_pt);
            //                    for (int i = 0; i< ind_pt; i++) {
            //                        printf("%d\n",ind_t[i]);
            //                    }
            
            //For the b part of Ax<=b
            // The elements in b are always "0"
            
            //Add constraint
            //*****************************
            error = GRBaddconstr(model, (int)non_zero_num, ind_t, val_t, GRB_LESS_EQUAL, 0.0, NULL);//load is negative value at this moment
            //*****************************
            
            //Update model due to lazy model update strategy
            error = GRBupdatemodel(model);
            //if (error) goto QUIT;
            GRBwrite (model, "groubi_obj.lp" );
            GRBwrite (model, "groubi_obj.rlp" );
        }// STATUS CHANGE CONSTRAINTS: Each year constraints are added
        //printf("\n");
        
        //****The ending of One year constraints ****
    }// STATUS MAINTAIN CONSTRAINTS: Total planning years constraints are added
    // ******** The ending of set STATUS MAINTAIN constraints ********
    
    
    
    // ******** The begining of set POWER FLOW constraints ********
    for (int PF_year = 0; PF_year < (int)nplan_year; PF_year++) {
        
        //**** One year constraints ****
        for (int PF_nCline_t = 0; PF_nCline_t < nCline_total; PF_nCline_t++) { // this will add number of bus constraints
            
            //For the A part of Ax<=b
            double non_zero_num = 0;
            double nonzero_PF_position[BUFFER_MAX];
            double nonzero_PF_Pvalue[BUFFER_MAX];
            double nonzero_PF_Nvalue[BUFFER_MAX];
            
            Array_initial(nonzero_PF_position, BUFFER_MAX);
            Array_initial(nonzero_PF_Pvalue, BUFFER_MAX);
            Array_initial(nonzero_PF_Nvalue, BUFFER_MAX);
            
            int nonzero_pt = 0;
            int position_pt = 0;
            
            //Fkl_1
            double flowC_CPF_t[(int)(nCline_total * nCline_total)];
            Array_initial(flowC_CPF_t,nCline_total * nCline_total);
            eye(flowC_CPF_t, nCline_total);
            
            for (int i = 0; i < nCline_total; i++) {
                if ((flowC_CPF_t[(int)(i + PF_nCline_t*nCline_total)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_PF_position[nonzero_pt] = position_pt;
                    nonzero_PF_Pvalue[nonzero_pt] = (flowC_CPF_t[(int)(i + PF_nCline_t*nCline_total)]);
                    nonzero_PF_Nvalue[nonzero_pt] =-(flowC_CPF_t[(int)(i + PF_nCline_t*nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //Xkl
            Array_coef_multiply(flowC_CPF_t, M, (nCline_total * nCline_total));

            
            for (int i = 0; i < nCline_total; i++) {
                if ((flowC_CPF_t[(int)(i + PF_nCline_t*nCline_total)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_PF_position[nonzero_pt] = position_pt;
                    nonzero_PF_Pvalue[nonzero_pt] = (flowC_CPF_t[(int)(i + PF_nCline_t*nCline_total)]);
                    nonzero_PF_Nvalue[nonzero_pt] = (flowC_CPF_t[(int)(i + PF_nCline_t*nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //Gk
            position_pt = position_pt + nGen;
            
            //Rk1
            position_pt = position_pt + nbus;
            
            //Rk2
            position_pt = position_pt + nbus;
            
            //theta
            double X_C_inverse[(int)(nCline_total*nCline_total)];// Inverse X_C
            Array_initial(X_C_inverse, nCline_total*nCline_total);
            Matrix_diagonal_inverse(X_C, nCline_total, nCline_total*nCline_total, X_C_inverse);
            
            double Kl_C_transpose[(int)(nCline*nbus)]; // Transpose Kl_c
            Array_initial(Kl_C_transpose, nCline*nbus);
            Matrix_transpose(Kl_C, nbus, nCline, Kl_C_transpose);
            
            double theta_CPF_t[(int)(nCline_total*nbus)]; // define theta variable's coefficients
            Matrix_multiply(theta_CPF_t, X_C_inverse, nCline_total, nCline_total, Kl_C_transpose, nCline_total, nbus);//theta cal.
            //Matrix_display(theta_CPF_t, nCline_total, nbus);
            
            
            for (int i = 0; i < nbus; i++) { // search for each column of a certain row
                if ((theta_CPF_t[(int)(i + PF_nCline_t*nbus)]) != 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_PF_position[nonzero_pt] = position_pt;
                    nonzero_PF_Pvalue[nonzero_pt] = -(theta_CPF_t[(int)(i + PF_nCline_t*nbus)]);
                    nonzero_PF_Nvalue[nonzero_pt] =  (theta_CPF_t[(int)(i + PF_nCline_t*nbus)]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //Gk_ts (it is the reverse of candidate line flow)
            position_pt = position_pt + nCline_total;
            
            //I
            position_pt = position_pt + nCline_total;
            
            //Set positive coefficient (ind_t) and positive position(val_t)
            int p_ind_t[(int)non_zero_num];
            double p_val_t[(int)non_zero_num];
            Array_initial_int(p_ind_t, non_zero_num);
            Array_initial(p_val_t, non_zero_num);
            int p_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                p_ind_t[p_ind_pt] = nonzero_PF_position[i] + (int)nAV*PF_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                p_val_t[p_ind_pt] = nonzero_PF_Pvalue[i];// But the value of cofficients are not change
                ++p_ind_pt;
            }
            
            //Set negative coefficient (ind_t) and negative position(val_t)
            int n_ind_t[(int)non_zero_num];
            double n_val_t[(int)non_zero_num];
            Array_initial_int(n_ind_t, non_zero_num);
            Array_initial(n_val_t, non_zero_num);
            int n_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                n_ind_t[n_ind_pt] = nonzero_PF_position[i] + (int)nAV*PF_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                n_val_t[n_ind_pt] = nonzero_PF_Nvalue[i];// But the value of cofficients are not change
                ++n_ind_pt;
            }

            //For the b part of Ax<=b
            // The elements in b are always "M"
            
            //Add constraint
            //*****************************
            error = GRBaddconstr(model, (int)non_zero_num, p_ind_t, p_val_t, GRB_LESS_EQUAL, M, NULL);//load is negative value at this moment
            error = GRBaddconstr(model, (int)non_zero_num, n_ind_t, n_val_t, GRB_LESS_EQUAL, M, NULL);//load is negative value at this moment
            //*****************************
            
            //Update model due to lazy model update strategy
            error = GRBupdatemodel(model);
            //if (error) goto QUIT;
            GRBwrite (model, "groubi_obj.lp" );
            GRBwrite (model, "groubi_obj.rlp" );
        }// STATUS CHANGE CONSTRAINTS: Each year constraints are added
        //printf("\n");
        //printf("\n");
        //****The ending of One year constraints ****
    }// POWER FLOW CONSTRAINTS: Total planning years constraints are added
    // ******** The ending of set POWER FLOW constraints ********
   
    
    // ******** The begining of BRANCH LIMIT constraints ********
    for (int BL_year = 0; BL_year < (int)nplan_year; BL_year++) {
        
        //**** One year constraints ****
        for (int BL_nCline_t = 0; BL_nCline_t < nCline_total; BL_nCline_t++) { // this will add number of bus constraints
            
            //For the A part of Ax<=b
            double non_zero_num = 0;
            double nonzero_BL_position[BUFFER_MAX];
            double nonzero_BL_Pvalue[BUFFER_MAX];
            double nonzero_BL_Nvalue[BUFFER_MAX];
            
            Array_initial(nonzero_BL_position, BUFFER_MAX);
            Array_initial(nonzero_BL_Pvalue, BUFFER_MAX);
            Array_initial(nonzero_BL_Nvalue, BUFFER_MAX);
            
            int nonzero_pt = 0;
            int position_pt = 0;
            
            //Fkl_1
            double flowC_CBL_t[(int)(nCline_total * nCline_total)];
            Array_initial(flowC_CBL_t,nCline_total * nCline_total);
            eye(flowC_CBL_t, nCline_total);
            
            for (int i = 0; i < nCline_total; i++) {
                if ((flowC_CBL_t[(int)(i + BL_nCline_t*nCline_total)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_BL_position[nonzero_pt] = position_pt;
                    nonzero_BL_Pvalue[nonzero_pt] = (flowC_CBL_t[(int)(i + BL_nCline_t*nCline_total)]);
                    nonzero_BL_Nvalue[nonzero_pt] =-(flowC_CBL_t[(int)(i + BL_nCline_t*nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //Xkl
            double x_CBL_t[(int)(nCline_total * nCline_total)];
            Array_initial(x_CBL_t,nCline_total * nCline_total);
            eye(x_CBL_t, nCline_total);
            Array_coef_multiply(x_CBL_t, M, (nCline_total * nCline_total));
            
            for (int i = 0; i<nCline_total; i++) {
                x_CBL_t[(int)(i*nCline_total+i)] =  Cline_info.line_Climit[i] - x_CBL_t[(int)(i*nCline_total+i)];
            }
            
            //Matrix_display(x_CBL_t, nCline_total, nCline_total);
            
            for (int i = 0; i < nCline_total; i++) { //nCline_total number of column
                if ((x_CBL_t[(int)(i + BL_nCline_t*nCline_total)]) != 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_BL_position[nonzero_pt] = position_pt;
                    nonzero_BL_Pvalue[nonzero_pt] = -(x_CBL_t[(int)(i + BL_nCline_t*nCline_total)]);
                    nonzero_BL_Nvalue[nonzero_pt] = -(x_CBL_t[(int)(i + BL_nCline_t*nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //Gk
            position_pt = position_pt + nGen;
            
            //Rk1
            position_pt = position_pt + nbus;
            
            //Rk2
            position_pt = position_pt + nbus;
            
            //theta
            position_pt = position_pt + nbus;

            
            //Gk_ts (it is the reverse of candidate line flow)
            position_pt = position_pt + nCline_total;
            
            //I
            position_pt = position_pt + nCline_total;
            
            
            //Set positive coefficient (ind_t) and positive position(val_t)
            int p_ind_t[(int)non_zero_num];
            double p_val_t[(int)non_zero_num];
            Array_initial_int(p_ind_t, non_zero_num);
            Array_initial(p_val_t, non_zero_num);
            int p_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                p_ind_t[p_ind_pt] = nonzero_BL_position[i] + (int)nAV*BL_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                p_val_t[p_ind_pt] = nonzero_BL_Pvalue[i];// But the value of cofficients are not change
                ++p_ind_pt;
            }
            
            //Set negative coefficient (ind_t) and negative position(val_t)
            int n_ind_t[(int)non_zero_num];
            double n_val_t[(int)non_zero_num];
            Array_initial_int(n_ind_t, non_zero_num);
            Array_initial(n_val_t, non_zero_num);
            int n_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                n_ind_t[n_ind_pt] = nonzero_BL_position[i] + (int)nAV*BL_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                n_val_t[n_ind_pt] = nonzero_BL_Nvalue[i];// But the value of cofficients are not change
                ++n_ind_pt;
            }
            
            
            
            //For the b part of Ax<=b
            // The elements in b are always "M"
            
            //Add constraint
            //*****************************
            error = GRBaddconstr(model, (int)non_zero_num, p_ind_t, p_val_t, GRB_LESS_EQUAL, M, NULL);//load is negative value at this moment
            error = GRBaddconstr(model, (int)non_zero_num, n_ind_t, n_val_t, GRB_LESS_EQUAL, M, NULL);//load is negative value at this moment
            //*****************************
            
            //Update model due to lazy model update strategy
            error = GRBupdatemodel(model);
            //if (error) goto QUIT;
            GRBwrite (model, "groubi_obj.lp" );
            GRBwrite (model, "groubi_obj.rlp" );
        }// STATUS CHANGE CONSTRAINTS: Each year constraints are added
        //printf("\n");
        //printf("\n");
        //****The ending of One year constraints ****
    }//BRANCH LIMIT: Total planning years constraints are added
    // ******** The ending of set BRANCH LIMIT constraints ********

    
    // ******** The begining of Generation Limit constraints ********
    for (int GL_year = 0; GL_year < (int)nplan_year; GL_year++) {
        
        //**** One year constraints ****
        for (int GL_nGen_t = 0; GL_nGen_t < nGen; GL_nGen_t++) { // this will add number of bus constraints
            
            //For the A part of Ax<=b
            double non_zero_num = 0;
            double nonzero_GL_position[BUFFER_MAX];
            double nonzero_GL_Pvalue[BUFFER_MAX];
            double nonzero_GL_Nvalue[BUFFER_MAX];
            
            Array_initial(nonzero_GL_position, BUFFER_MAX);
            Array_initial(nonzero_GL_Pvalue, BUFFER_MAX);
            Array_initial(nonzero_GL_Nvalue, BUFFER_MAX);
            
            int nonzero_pt = 0;
            int position_pt = 0;
            
            //Fkl_1
            position_pt = position_pt + nCline_total;
            
            //Xkl
            position_pt = position_pt + nCline_total;
            
            //Gk
            double gen_GL_t[(int)(nGen * nGen)];
            Array_initial(gen_GL_t,nGen * nGen);
            eye(gen_GL_t, nGen);
            //Matrix_display(gen_GL_t, nGen, nGen);
            
            for (int i = 0; i < nGen; i++) {
                if ((gen_GL_t[(int)(i + GL_nGen_t*nGen)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_GL_position[nonzero_pt] = position_pt;
                    nonzero_GL_Pvalue[nonzero_pt] = (gen_GL_t[(int)(i + GL_nGen_t*nGen)]);
                    nonzero_GL_Nvalue[nonzero_pt] =-(gen_GL_t[(int)(i + GL_nGen_t*nGen)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
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
            position_pt = position_pt + nCline_total;
            
            //I
            position_pt = position_pt + nCline_total;
            
            
            //Set positive coefficient (ind_t) and positive position(val_t)
            int p_ind_t[(int)non_zero_num];
            double p_val_t[(int)non_zero_num];
            Array_initial_int(p_ind_t, non_zero_num);
            Array_initial(p_val_t, non_zero_num);
            int p_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                p_ind_t[p_ind_pt] = nonzero_GL_position[i] + (int)nAV*GL_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                p_val_t[p_ind_pt] = nonzero_GL_Pvalue[i];// But the value of cofficients are not change
                ++p_ind_pt;
            }
            
            //Set negative coefficient (ind_t) and negative position(val_t)
            int n_ind_t[(int)non_zero_num];
            double n_val_t[(int)non_zero_num];
            Array_initial_int(n_ind_t, non_zero_num);
            Array_initial(n_val_t, non_zero_num);
            int n_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                n_ind_t[n_ind_pt] = nonzero_GL_position[i] + (int)nAV*GL_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                n_val_t[n_ind_pt] = nonzero_GL_Nvalue[i];// But the value of cofficients are not change
                ++n_ind_pt;
            }
            
            
            //For the b part of Ax<=b
            // The b is the Gmin and Gmax
            
            //Add constraint
            //*****************************
            error = GRBaddconstr(model, (int)non_zero_num, p_ind_t, p_val_t, GRB_LESS_EQUAL, Gen_info.gen_max[GL_nGen_t], NULL);//load is negative value at this moment
            error = GRBaddconstr(model, (int)non_zero_num, n_ind_t, n_val_t, GRB_LESS_EQUAL,-Gen_info.gen_min[GL_nGen_t], NULL);//load is negative value at this moment
            //*****************************
            
            //Update model due to lazy model update strategy
            error = GRBupdatemodel(model);
            //if (error) goto QUIT;
            GRBwrite (model, "groubi_obj.lp" );
            GRBwrite (model, "groubi_obj.rlp" );
        }// STATUS CHANGE CONSTRAINTS: Each year constraints are added
        //printf("\n");
        //printf("\n");
        //****The ending of One year constraints ****
    }//GENERATION LIMIT: Total planning years constraints are added
    // ******** The ending of set GENERATION LIMIT constraints ********
    
    
    // ******** The begining of PSEUDO GENERATOR LINE RELATIONSHIP constraints ********
    for (int PS_year = 0; PS_year < (int)nplan_year; PS_year++) {
        
        //**** One year constraints ****
        for (int PS_nCline_t = 0; PS_nCline_t < nCline_total; PS_nCline_t++) { // this will add number of bus constraints
            
            //For the A part of Ax<=b
            double non_zero_num = 0;
            double nonzero_PS_position[BUFFER_MAX];
            double nonzero_PS_Pvalue[BUFFER_MAX];
            double nonzero_PS_Nvalue[BUFFER_MAX];
            
            Array_initial(nonzero_PS_position, BUFFER_MAX);
            Array_initial(nonzero_PS_Pvalue, BUFFER_MAX);
            Array_initial(nonzero_PS_Nvalue, BUFFER_MAX);
            
            int nonzero_pt = 0;
            int position_pt = 0;
            
            //Fkl_1
            double flowC_PS_t[(int)(nCline_total * nCline_total)];
            Array_initial(flowC_PS_t,nCline_total * nCline_total);
            eye(flowC_PS_t, nCline_total);
            
            for (int i = 0; i < nCline_total; i++) {
                if ((flowC_PS_t[(int)(i + PS_nCline_t*nCline_total)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_PS_position[nonzero_pt] = position_pt;
                    nonzero_PS_Pvalue[nonzero_pt] = (flowC_PS_t[(int)(i + PS_nCline_t * nCline_total)]);
                    nonzero_PS_Nvalue[nonzero_pt] =-(flowC_PS_t[(int)(i + PS_nCline_t * nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //Xkl
            double x_PS_t[(int)(nCline_total * nCline_total)];
            Array_initial(x_PS_t,nCline_total * nCline_total);
            eye(x_PS_t, nCline_total);
            Array_coef_multiply(x_PS_t, M, (nCline_total * nCline_total));
            
            for (int i = 0; i < nCline_total; i++) { //nCline_total number of column
                if ((x_PS_t[(int)(i + PS_nCline_t*nCline_total)]) != 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_PS_position[nonzero_pt] = position_pt;
                    nonzero_PS_Pvalue[nonzero_pt] = -(x_PS_t[(int)(i + PS_nCline_t*nCline_total)]);
                    nonzero_PS_Nvalue[nonzero_pt] = -(x_PS_t[(int)(i + PS_nCline_t*nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //Gk
            position_pt = position_pt + nGen;
            
            //Rk1
            position_pt = position_pt + nbus;
            
            //Rk2
            position_pt = position_pt + nbus;
            
            //theta
            position_pt = position_pt + nbus;
            
            
            //Gk_ts (it is the reverse of candidate line flow)
            double Kp_ts_PS_t[(int)(nCline_total * nCline_total)];
            Array_initial(Kp_ts_PS_t,nCline_total * nCline_total);
            eye(Kp_ts_PS_t, nCline_total);
            
            for (int i = 0; i < nCline_total; i++) {
                if ((Kp_ts_PS_t[(int)(i + PS_nCline_t*nCline_total)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_PS_position[nonzero_pt] = position_pt;
                    nonzero_PS_Pvalue[nonzero_pt] =-(Kp_ts_PS_t[(int)(i + PS_nCline_t * nCline_total)]);
                    nonzero_PS_Nvalue[nonzero_pt] = (Kp_ts_PS_t[(int)(i + PS_nCline_t * nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //I
            position_pt = position_pt + nCline_total;
            
            
            //Set positive coefficient (ind_t) and positive position(val_t)
            int p_ind_t[(int)non_zero_num];
            double p_val_t[(int)non_zero_num];
            Array_initial_int(p_ind_t, non_zero_num);
            Array_initial(p_val_t, non_zero_num);
            int p_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                p_ind_t[p_ind_pt] = nonzero_PS_position[i] + (int)nAV*PS_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                p_val_t[p_ind_pt] = nonzero_PS_Pvalue[i];// But the value of cofficients are not change
                ++p_ind_pt;
            }
            
            //Set negative coefficient (ind_t) and negative position(val_t)
            int n_ind_t[(int)non_zero_num];
            double n_val_t[(int)non_zero_num];
            Array_initial_int(n_ind_t, non_zero_num);
            Array_initial(n_val_t, non_zero_num);
            int n_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                n_ind_t[n_ind_pt] = nonzero_PS_position[i] + (int)nAV*PS_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                n_val_t[n_ind_pt] = nonzero_PS_Nvalue[i];// But the value of cofficients are not change
                ++n_ind_pt;
            }
            
            
            //For the b part of Ax<=b
            // The elements in b are always 0
            
            //Add constraint
            //*****************************
            error = GRBaddconstr(model, (int)non_zero_num, p_ind_t, p_val_t, GRB_LESS_EQUAL, 0, NULL);//load is negative value at this moment
            error = GRBaddconstr(model, (int)non_zero_num, n_ind_t, n_val_t, GRB_LESS_EQUAL, 0, NULL);//load is negative value at this moment
            //*****************************
            
            //Update model due to lazy model update strategy
            error = GRBupdatemodel(model);
            //if (error) goto QUIT;
            GRBwrite (model, "groubi_obj.lp" );
            GRBwrite (model, "groubi_obj.rlp" );
        }// STATUS CHANGE CONSTRAINTS: Each year constraints are added
        //printf("\n");
        //printf("\n");
        //****The ending of One year constraints ****
    }//PSEUDO GENERATOR LINE RELATIONSHIP: Total planning years constraints are added
    // ******** The ending of set PSEUDO GENERATOR LINE RELATIONSHIP constraints ********
    
    
    
    // ******** The begining of PSEUDO GENERATOR LIMIT relationship constraints ********
    for (int PSL_year = 0; PSL_year < (int)nplan_year; PSL_year++) {
        
        //**** One year constraints ****
        for (int PSL_nCline_t = 0; PSL_nCline_t < nCline_total; PSL_nCline_t++) { // this will add number of bus constraints
            
            //For the A part of Ax<=b
            double non_zero_num = 0;
            double nonzero_PSL_position[BUFFER_MAX];
            double nonzero_PSL_Pvalue[BUFFER_MAX];
            double nonzero_PSL_Nvalue[BUFFER_MAX];
            
            Array_initial(nonzero_PSL_position, BUFFER_MAX);
            Array_initial(nonzero_PSL_Pvalue, BUFFER_MAX);
            Array_initial(nonzero_PSL_Nvalue, BUFFER_MAX);
            
            int nonzero_pt = 0;
            int position_pt = 0;
            
            //Fkl_1
            position_pt = position_pt + nCline_total;
            
            //Xkl
            double x_PSL_t[(int)(nCline_total * nCline_total)];
            Array_initial(x_PSL_t,nCline_total * nCline_total);
            eye(x_PSL_t, nCline_total);
            Array_coef_multiply(x_PSL_t, M, (nCline_total * nCline_total));
            
            for (int i = 0; i < nCline_total; i++) { //nCline_total number of column
                if ((x_PSL_t[(int)(i + PSL_nCline_t*nCline_total)]) != 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_PSL_position[nonzero_pt] = position_pt;
                    nonzero_PSL_Pvalue[nonzero_pt] =  (x_PSL_t[(int)(i + PSL_nCline_t*nCline_total)]);
                    nonzero_PSL_Nvalue[nonzero_pt] =  (x_PSL_t[(int)(i + PSL_nCline_t*nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //Gk
            position_pt = position_pt + nGen;
            
            //Rk1
            position_pt = position_pt + nbus;
            
            //Rk2
            position_pt = position_pt + nbus;
            
            //theta
            position_pt = position_pt + nbus;
            
            
            //Gk_ts (it is the reverse of candidate line flow)
            double Kp_ts_PSL_t[(int)(nCline_total * nCline_total)];
            Array_initial(Kp_ts_PSL_t,nCline_total * nCline_total);
            eye(Kp_ts_PSL_t, nCline_total);
            
            for (int i = 0; i < nCline_total; i++) {
                if ((Kp_ts_PSL_t[(int)(i + PSL_nCline_t*nCline_total)]) > 0) { //0+(i-1): the 0 can be used as variable
                    ++non_zero_num;
                    nonzero_PSL_position[nonzero_pt] = position_pt;
                    nonzero_PSL_Pvalue[nonzero_pt] = (Kp_ts_PSL_t[(int)(i + PSL_nCline_t * nCline_total)]);
                    nonzero_PSL_Nvalue[nonzero_pt] =-(Kp_ts_PSL_t[(int)(i + PSL_nCline_t * nCline_total)]);
                    //printf("%f\n", p_SC_position[p_pt]);
                    ++nonzero_pt;
                }
                ++position_pt;
            }
            
            //I
            position_pt = position_pt + nCline_total;
            
            
            //Set positive coefficient (ind_t) and positive position(val_t)
            int p_ind_t[(int)non_zero_num];
            double p_val_t[(int)non_zero_num];
            Array_initial_int(p_ind_t, non_zero_num);
            Array_initial(p_val_t, non_zero_num);
            int p_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                p_ind_t[p_ind_pt] = nonzero_PSL_position[i] + (int)nAV*PSL_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                p_val_t[p_ind_pt] = nonzero_PSL_Pvalue[i];// But the value of cofficients are not change
                ++p_ind_pt;
            }
            
            //Set negative coefficient (ind_t) and negative position(val_t)
            int n_ind_t[(int)non_zero_num];
            double n_val_t[(int)non_zero_num];
            Array_initial_int(n_ind_t, non_zero_num);
            Array_initial(n_val_t, non_zero_num);
            int n_ind_pt = 0;
            
            for (int i = 0; i < (int) nonzero_pt; i++) { //positive coefficient
                n_ind_t[n_ind_pt] = nonzero_PSL_position[i] + (int)nAV*PSL_year; // The term "(int)nAV*NB_year" is offset of the year. When year increase, constraints just move in the matrix diagonal
                n_val_t[n_ind_pt] = nonzero_PSL_Nvalue[i];// But the value of cofficients are not change
                ++n_ind_pt;
            }
            
            
            //For the b part of Ax<=b
            // The elements in b are always M
            
            //Add constraint
            //*****************************
            error = GRBaddconstr(model, (int)non_zero_num, p_ind_t, p_val_t, GRB_LESS_EQUAL, M, NULL);//load is negative value at this moment
            error = GRBaddconstr(model, (int)non_zero_num, n_ind_t, n_val_t, GRB_LESS_EQUAL, M, NULL);//load is negative value at this moment
            //*****************************
            
            //Update model due to lazy model update strategy
            error = GRBupdatemodel(model);
            //if (error) goto QUIT;
            GRBwrite (model, "groubi_obj.lp" );
            GRBwrite (model, "groubi_obj.rlp" );
        }// STATUS CHANGE CONSTRAINTS: Each year constraints are added
        //printf("\n");
        //printf("\n");
        //****The ending of One year constraints ****
    }//PSEUDO GENERATOR LINE LIMIT: Total planning years constraints are added
    // ******** The ending of set PSEUDO GENERATOR LINE LIMIT constraints ********
    
    
    
    
    
    
    
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
