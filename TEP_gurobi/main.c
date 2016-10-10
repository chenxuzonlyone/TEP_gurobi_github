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
    
    for (int i = 0; i< row_cline; i++)
    {
        
        printf("%f\t", Cline_info.line_Cfrom[i]);
        printf("%f\t", Cline_info.line_Cto[i]);
        printf("%f\t", Cline_info.line_Ccost[i]);
        printf("%f\t", Cline_info.line_Creactance[i]);
        printf("%f\t", Cline_info.line_Climit[i]);
        printf("\n");
    }
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
    
    for (int i = 0; i< row_gen; i++)
    {
        
        printf("%f\t", Gen_info.gen_busnum[i]);
        printf("%f\t", Gen_info.gen_min[i]);
        printf("%f\t", Gen_info.gen_max[i]);
        printf("%f\t", Gen_info.gen_fixed[i]);
        printf("%f\t", Gen_info.gen_cost[i]);
        printf("\n");
    }
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
    
    for (int i = 0; i< row_load; i++)
    {
        
        printf("%f\t", Load_info.load_busnum[i]);
        printf("%f\t", Load_info.load_fixed[i]);
        printf("\n");
    }
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
    double million_transfer=10^6; // transfer unit from dollar to million dollar
    
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
    
    
    
    /********This part will close file and clean memory********/
    fclose(f_cline_stream);
    free(Cline_data_space);
    
    fclose(f_gen_stream);
    free(Gen_data_space);
    
    fclose(f_load_stream);
    free(Load_data_space);
    return 0;
}

