#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gurobi_c.h"
#include "read_excel.h"
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
    
    
    
    
    
    
    
    /*This part will close file and clean memory*/
    fclose(f_cline_stream);
    free(Cline_data_space);
    
    fclose(f_gen_stream);
    free(Gen_data_space);
    
    fclose(f_load_stream);
    free(Load_data_space);
    return 0;
}

