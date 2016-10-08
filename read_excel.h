//
//  read_excel.h
//  TEP_gurobi
//
//  Created by zhangcaihua on 9/19/16.
//  Copyright © 2016 zhangcaihua. All rights reserved.
//

#ifndef read_excel_h
#define read_excel_h
#endif /* read_excel_h */
#define BUFFER_MAX 1024

#include <stdio.h>
typedef struct Cline{
    double *line_Cfrom;
    double *line_Cto;
    double *line_Ccost;
    double *line_Creactance;
    double *line_Climit;
} Cline_struct;


typedef struct Gen{
    double *
    
 } Gen_struct;

int abs_value (int value);
int file_size(FILE * fstream, int *row, int *col); //decide the dimension of matrix
void delimiter_array(int* space_delimiter,int row, int col); // generate space array 
int Cline_data_read(FILE * fstream, int row, int col, Cline_struct* Cline_info); //
