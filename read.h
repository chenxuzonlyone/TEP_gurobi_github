//
//  read.h
//  TEP_gurobi
//
//  Created by zhangcaihua on 9/19/16.
//  Copyright © 2016 zhangcaihua. All rights reserved.
//

#ifndef read_h
#define read_h
#endif /* read_h */
#define BUFFER_MAX 1024

//For structure like this way, the memory location is continuous.

#include <stdio.h>
typedef struct Cline{
    double *line_Cfrom;
    double *line_Cto;
    double *line_Ccost;
    double *line_Creactance;
    double *line_Climit;
} Cline_struct;

typedef struct Gen{
    double *gen_busnum;
    double *gen_min;
    double *gen_max;
    double *gen_fixed;
    double *gen_cost;
} Gen_struct;


typedef struct Load{
    double *load_busnum;
    double *load_fixed;
} Load_struct;


int abs_value(int value);

/*Read data from FILE*/
int file_size(FILE * fstream, int *row, int *col); //decide the dimension of matrix
void delimiter_array(int* space_delimiter, int row, int col); // generate space array
int Cline_data_read(FILE * fstream, int row, int col, Cline_struct* Cline_info); //read cline data from file into code
int Gen_data_read(FILE * fstream, int row, int col, Gen_struct* info);//read gen data from file into code
int Load_data_read(FILE * fstream, int row, int col, Load_struct* info);//read load data from file into code


/*Read data from STRUCTURE*/
int Cline_struct_read (double *Cline_array, Cline_struct line_info, char mode,double row_struc, double col_struc, double row_need, double col_need);     //read candidate line structure. Row or Col is set to "0" when not intend to read
int Gen_struct_read (double *Gen_array, Gen_struct gen_info, char mode,double row_struc, double col_struc, double row_need, double col_need);           // read candidate gen structure. Row or Col is set to "0" when not intend to read
int Load_struct_read (double *Load_array, Load_struct load_info, char mode,double row_struc, double col_struc, double row_need, double col_need);     // read candidate load structure. Row or Col is set to "0" when not intend to read


/*Array operation*/
int Array_coef_multiply(double *Array_need_multiply, double coeff, double length);//This will let every element in the array multiply the coefficient
int Array_ascend(double *Array_asscending, double length); //This will make an array with contents from 1 to length+1
int Array_initial(double *Array_initialition, double length);
int eye(double *eye_array, double length);
int diag(double *diag_array, double length, double *coeff_array);
int ones(double *ones_array, double length); // This will give you an array with all one

/*Kl_C, Kp parameter calculation*/
int Kl_C_set(double *Kl_C, Cline_struct Cline_info, double nbus, double nCline); // Kl_C
int Kp_set(double *Kp, Gen_struct Gen_info, double nbus, double nGen);// Kp
int Kp_ts_set(double *Kp_ts, Cline_struct Cline_info, double nbus, double nCline);//Kp_ts
int Kd_set(double *Kd, Load_struct Load_info, double nbus, double nload); // kd

