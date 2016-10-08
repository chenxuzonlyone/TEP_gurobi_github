#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gurobi_c.h"
#include "read_excel.h"
#include <string.h>



int main()
{
    int row, col;

    FILE *f_cline_stream = fopen("/Users/zhangcaihua/Desktop/TEP_gurobi/data/cline_info.csv","r");
    file_size(f_cline_stream, &row, &col);// get rows and cols from file
    
    Cline_struct Cline_info;
    double*Cline_data_space;
    int space_delimiter[col];
    delimiter_array(space_delimiter, row, col);

    Cline_data_space = calloc(row*col, sizeof(double));
    Cline_info.line_Cfrom = &Cline_data_space[space_delimiter[0]];
    Cline_info.line_Cto = &Cline_data_space[space_delimiter[1]];
    Cline_info.line_Ccost = &Cline_data_space[space_delimiter[2]];
    Cline_info.line_Creactance = &Cline_data_space[space_delimiter[3]];
    Cline_info.line_Climit = &Cline_data_space[space_delimiter[4]];
    
    
    
    Cline_data_read(f_cline_stream, row, col, &Cline_info);//read data from file
    
    
    for (int i=0; i< row; i++)
    {
        
        printf("%f\t", Cline_info.line_Cfrom[i]);
        printf("%f\t", Cline_info.line_Cto[i]);
        printf("%f\t", Cline_info.line_Ccost[i]);
        printf("%f\t", Cline_info.line_Creactance[i]);
        printf("%f\t", Cline_info.line_Climit[i]);
        printf("\n");
    }
    
    
    
    
    /*This part will close file and clean memory*/
    fclose(f_cline_stream);
    free(Cline_data_space);
    

    return 0;
}
