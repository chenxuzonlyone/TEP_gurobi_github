//
//  read.c
//  TEP_gurobi
//
//  Created by zhangcaihua on 9/19/16.
//  Copyright © 2016 zhangcaihua. All rights reserved.
//

#include "read.h" // This include a global variable used for Inverse of Matrix
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gurobi_c.h"
#include <string.h>

int file_size(FILE * fstream, int *row, int *col)
{
    FILE * fstream_inside = fstream;
    char buffer[BUFFER_MAX];
    char *record, *line;
    int i = 0, j = 0;
    
    if (fstream_inside == NULL)   {
        printf("\n file opening failed ");
        return -1;
    }
    while ((line = fgets(buffer, sizeof(buffer), fstream_inside)) != NULL)
    {
        record = strtok(line, ",");
        //printf("file position %ld\n", ftell(fstream));
        while (record != NULL && i == 0) //just need to test first line
        {
            j++;
            record = strtok(NULL, ",");
            
        }
        i++;
    }
    *row = i;
    *col = j;
    
    rewind(fstream_inside); // It is vary important that to restore file pointer back to start point when finish
    // Or it can assign it's value to a new variable and keep origional pointer safe.
    // Besides, fseek(fp, 0, SEEK_SET) and rewind(fp) have same effect.
    return 0;
}

void delimiter_array(int* space_delimiter, int row, int col)
{
    for (int i = 0; i < col; i++)
    {
        space_delimiter[i] = i*row;
    }
}

// Cline_data_read
int Cline_data_read(FILE * fstream, int row, int col, Cline_struct* info)
{
    FILE * fstream_inside = fstream;
    char buffer[BUFFER_MAX];
    char *record, *line;
    int i = 0, j = 0;
    
    if (fstream_inside == NULL)   {
        printf("\n file opening failed ");
        return -1;
    }
    
    while ((line = fgets(buffer, sizeof(buffer), fstream_inside)) != NULL)
    {
        record = strtok(line, ",");
        j = 0;
        while (record != NULL)
        {
            info->line_Cfrom[i + j*row] = atof(record);
            ++j;
            record = strtok(NULL, ",");
        }
        ++i;
    }
    
    printf("\n");
    rewind(fstream_inside);
    return 0;
}

int Gen_data_read(FILE * fstream, int row, int col, Gen_struct* info)
{
    FILE * fstream_inside = fstream;
    char buffer[BUFFER_MAX];
    char *record, *gen;
    int i = 0, j = 0;
    
    if (fstream_inside == NULL)   {
        printf("\n file opening failed ");
        return -1;
    }
    
    
    while ((gen = fgets(buffer, sizeof(buffer), fstream_inside)) != NULL)
    {
        record = strtok(gen, ",");
        j = 0;
        while (record != NULL)
        {
            info->gen_busnum[i + j*row] = atof(record);
            ++j;
            record = strtok(NULL, ",");
        }
        ++i;
    }
    
    printf("\n");
    rewind(fstream_inside);
    return 0;
}

int Load_data_read(FILE * fstream, int row, int col, Load_struct* info)
{
    FILE * fstream_inside = fstream;
    char buffer[BUFFER_MAX];
    char *record, *load;
    int i = 0, j = 0;
    
    if (fstream_inside == NULL)   {
        printf("\n file opening failed ");
        return -1;
    }
    
    
    while ((load = fgets(buffer, sizeof(buffer), fstream_inside)) != NULL)
    {
        record = strtok(load, ",");
        j = 0;
        while (record != NULL)
        {
            info->load_busnum[i + j*row] = atof(record);
            ++j;
            record = strtok(NULL, ",");
        }
        ++i;
    }
    
    printf("\n");
    rewind(fstream_inside);
    return 0;
}


int Cline_struct_read (double *Cline_array, Cline_struct line_info, char mode,double row_struc, double col_struc, double row_need, double col_need)
{
    // The "Cline_struct" is a structure filled with double type pointer. The mamory address is the continuous
    double *struct_first_member_pt =(line_info.line_Cfrom); //this show the address of first member (pointer)
    int r = (int) row_need;
    int c = (int) col_need;

    //********This will read one ROW data********
    if (mode == 'r') // r = row
    {
        for (int i = 1; i <= (int) col_struc; i++){
            Cline_array[i-1] = ((struct_first_member_pt) + ((i-1)*(int)row_struc))[r-1];
        }
    }
    
    //********This will read one COLUMN data********
    else if (mode == 'c') // c = column
    {
        for (int j = 1; j <= (int)row_struc; j++){
            Cline_array[j-1] = (struct_first_member_pt + ((c-1)*(int)row_struc))[j-1];
        }
    }
    //********This will report ERROR********
    else {printf("Please choose mode r or c\n");}
    
    return 0;
}

int Gen_struct_read (double *Gen_array, Gen_struct gen_info, char mode,double row_struc, double col_struc, double row_need, double col_need)
{
    // The "Gen_struct" is a structure filled with double type pointer. The mamory address is the continuous
    double *struct_first_member_pt =(gen_info.gen_busnum); //this show the address of first member (pointer)
    int r = (int) row_need;
    int c = (int) col_need;
    
    //********This will read one ROW data********
    if (mode == 'r') // r = row
    {
        for (int i = 1; i <= (int) col_struc; i++){
            Gen_array[i-1] = ((struct_first_member_pt) + ((i-1)*(int)row_struc))[r-1];
        }
    }
    
    //********This will read one COLUMN data********
    else if (mode == 'c') // c = column
    {
        for (int j = 1; j <= (int)row_struc; j++){
            Gen_array[j-1] = (struct_first_member_pt + ((c-1)*(int)row_struc))[j-1];
        }
    }
    //********This will report ERROR********
    else {printf("Please choose mode r or c\n");}
    
    return 0;
}

int Load_struct_read (double *Load_array, Load_struct load_info, char mode,double row_struc, double col_struc, double row_need, double col_need)
{
    // The "Load_struct" is a structure filled with double type pointer. The mamory address is the continuous
    double *struct_first_member_pt =(load_info.load_busnum); //this show the address of first member (pointer)
    int r = (int) row_need;
    int c = (int) col_need;
    
    //********This will read one ROW data********
    if (mode == 'r') // r = row
    {
        for (int i = 1; i <= (int) col_struc; i++){
            Load_array[i-1] = ((struct_first_member_pt) + ((i-1)*(int)row_struc))[r-1];
        }
    }
    
    //********This will read one COLUMN data********
    else if (mode == 'c') // c = column
    {
        for (int j = 1; j <= (int)row_struc; j++){
            Load_array[j-1] = (struct_first_member_pt + ((c-1)*(int)row_struc))[j-1];
        }
    }
    //********This will report ERROR********
    else {printf("Please choose mode r or c\n");}
    
    return 0;
}

// This will let every element in the array multiply the coefficient
int Array_coef_multiply(double *Array_need_multiply, double coeff, double length)
{
    int l = (int) length;
    for (int i = 0; i < l; i++){
        Array_need_multiply[i] = Array_need_multiply[i] * coeff;
    }
    return 0;
}

//This will make an array with contents from 1 to length+1
int Array_ascend(double *Array_asscending, double length)
{
    int l = (int) length;
    for (int i = 0; i < l; i++){
        Array_asscending[i] = i+1;
    }
    return 0;
}

// Array initialization
int Array_initial(double *Array_initialition, double length)
{
    int l = (int) length;
    for (int i = 0; i < l; i++) {
        Array_initialition[i] = 0;
    }
    return 0;
}

// Array initialization
int Array_initial_int(int *Array_initialition, double length)
{
    int l = (int) length;
    for (int i = 0; i < l; i++) {
        Array_initialition[i] = (int)0;
    }
    return 0;
}

int Kl_C_set(double *Kl_C, Cline_struct Cline_info, double nbus, double nCline)
{
    for (int i = 0; i < (int)nCline; ++i)
    {
        int l_from =(int)(Cline_info.line_Cfrom[i] -1)*nCline + i; // row * col_num + col
        int l_to =(int)(Cline_info.line_Cto[i] -1)*nCline + i; // row * col_num + col

        Kl_C[l_from] = 1.0;
        Kl_C[l_to] =  -1.0;
        //For testing purpose
        //printf("line from : %d \t", l_from);
        //printf("line to : %d \t", l_to);
        //printf("\n");
    }
    //Output to a file
    FILE *fp;
    fp = fopen("Kl_C.txt", "w+");
    for (int i = 0; i<(int)(nbus*nCline); i++) {
        fprintf(fp, "%f", Kl_C[i]);
        fprintf(fp, "\t");
        if (((i+1)%(int)nCline) == 0){
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
    return 0;
}

int Kp_set(double *Kp, Gen_struct Gen_info, double nbus, double nGen)
{
    for (int i = 0; i < (int)nGen; ++i)
    {
        int g_num =(int)(Gen_info.gen_busnum[i] -1)*nGen + i; // row * row_num + col
        Kp[g_num] = 1.0;
        //For testing purpose
        //printf("gen_bus_number : %d \t", g_num);
        //printf("\n");
    }
    //Output to a file
    FILE *fp;
    fp = fopen("Kp.txt", "w+");
    for (int i = 0; i<(int)(nbus*nGen); i++) {
        fprintf(fp, "%f", Kp[i]);
        fprintf(fp, "\t");
        if (((i+1)%(int)nGen) == 0){
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
    return 0;
}

int Kp_ts_set(double *Kp_ts, Cline_struct Cline_info, double nbus, double nCline)
{
    for (int i = 0; i < (int)nCline; ++i)
    {
        int l_from =(int)(Cline_info.line_Cfrom[i] -1)*nCline + i; // row * row_num + col
        int l_to =(int)(Cline_info.line_Cto[i] -1)*nCline + i; // row * row_num + col
        
        Kp_ts[l_from] = 1.0;
        Kp_ts[l_to] =  -1.0;
        //For testing purpose
        //printf("line from : %d \t", l_from);
        //printf("line to : %d \t", l_to);
        //printf("\n");
    }
    //Output to a file
    FILE *fp;
    fp = fopen("Kp_ts.txt", "w+");
    for (int i = 0; i<(int)(nbus*nCline); i++) {
        fprintf(fp, "%f", Kp_ts[i]);
        fprintf(fp, "\t");
        if (((i+1)%(int)nCline) == 0){
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
    return 0;
}


int Kd_set(double *Kd, Load_struct Load_info, double nbus, double nload)
{
    for (int i = 0; i < (int)nload; ++i)
    {
        int load_num =(int)(Load_info.load_busnum[i] -1)*nload + i; // row * row_num + col
        Kd[load_num] = 1.0;
        //For testing purpose
        //printf("load_bus_number : %d \t", load_num);
        //printf("\n");
    }
    //Output to a file
    FILE *fp;
    fp = fopen("Kd.txt", "w+");
    for (int i = 0; i<(int)(nbus*nload); i++) {
        fprintf(fp, "%f", Kd[i]);
        fprintf(fp, "\t");
        if (((i+1)%(int)nload) == 0){
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
    return 0;
}

int eye(double *eye_array, double length) // only need row or col number due to only work on the diaginal element
{
    for (int i = 0; i < (int)length; ++i)
    {
        int eye_num =(int)(i*length) + i; // row * row_num + col
        eye_array[eye_num] = 1.0;
        //For testing purpose
        //printf("eye : %d \t", eye_num);
        //printf("\n");
    }
    //Output to a file
    /*
    FILE *fp;
    fp = fopen("eye.txt", "w+");
    for (int i = 0; i<(int)(length*length); i++) {
        fprintf(fp, "%f", eye_array[i]);
        fprintf(fp, "\t");
        if (((i+1)%(int)length) == 0){
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
     */
    return  0;
}

int diag(double *diag_array, double length, double *coeff_array)
{
    for (int i = 0; i < (int)length; ++i)
    {
        int diag_num =(int)(i*length) + i; // row * row_num + col
        diag_array[diag_num] = coeff_array[i];
        //For testing purpose
        //printf("eye : %d \t", eye_num);
        //printf("\n");
    }
    /*
    FILE *fp;
    fp = fopen("diag.txt", "w+");
    for (int i = 0; i<(int)(length*length); i++) {
        fprintf(fp, "%f", diag_array[i]);
        fprintf(fp, "\t");
        if (((i+1)%(int)length) == 0){
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
     */
    return 0;
}

int ones(double *ones_array, double length)
{
    for (int i = 0; i <= (int)length; ++i){
        ones_array[i] = 1;
    }
    return 0;
}

int Matrix_multiply(double *outcome, double *matrix_a, double a_row, double a_col, double *matrix_b, double b_row, double b_col)
{
    Array_initial(outcome, a_row*b_col);
    int a_row_int = (int)a_row;
    int a_col_int = (int)a_col;
    //int b_row_int = (int)b_row;
    int b_col_int = (int)b_col;
    
    for (int i_a = 0; i_a < a_row_int; i_a++) {
        for (int j_b = 0; j_b < b_col_int; j_b++) {
            for (int j_a = 0; j_a < a_col_int; j_a++) {
                
                outcome[(i_a)*(b_col_int)+j_b] = (outcome[(i_a)*(b_col_int)+j_b])+matrix_a[j_a+ i_a*(a_col_int)] * matrix_b[(j_a)*b_col_int + j_b];
//                printf("%f\n", matrix_a[j_a+ i_a*(a_col_int)]);
//                printf("%f\n", matrix_b[(j_a)*b_col_int +j_b]);
//                printf("%f\n", matrix_a[j_a+ i_a*(a_col_int)] * matrix_b[(j_a)*b_col_int + j_b]);
            }
        }
    }
    
//    for (int i = 0; i<a_row_int; i++) {
//        for (int j = 0; j<a_col_int;j++) {
//            printf("%f \t", matrix_a[i*a_col_int+j]);
//        }
//        printf("\n");
//    }
//    for (int i = 0; i<b_row_int; i++) {
//        for (int j = 0; j<b_col_int;j++) {
//            printf("%f \t", matrix_b[i*b_col_int+j]);
//        }
//        printf("\n");
//    }
    
    
//    double outcome[9];
//    Array_initial(outcome, 9);
//    double matrix_a[6]={1,2,5,6,9,10};
//    double matrix_b[6]={12,11,10,9,8,7};
//    matrix_multiply(outcome, matrix_a,3,2,matrix_b,2,3);
//    for (int i=0; i<3; i++) {
//        for (int j = 0; j<3; j++) {
//            printf("%f \t",outcome[(i*3) + j]);
//        }
//        printf("\n");
//    }
    return 0;
}

int Matrix_display(double *Matrix, double row, double col)
{
    for (int i =0; i<(int)row; i++) {
        for (int j = 0; j<(int)col; j++) {
            printf("%f\t", Matrix[(int)(i*col+j)]);
        }
        printf("\n");
    }    
    return 0;
}

// This part will calculate the INVERSE OF MATRIX
//http://www.sanfoundry.com/c-program-find-inverse-matrix/
//This C program sorts a given array of integer numbers using Bubble Sort technique. The algorithm gets its name from the way smaller elements “bubble” to the top of the list. Because it only uses comparisons to operate on elements, it is a comparison sort. Time Complexity of this algorithm is O(n2).
/*For calculating Determinant of the Matrix */
//double determinant(double a[size_s][size_s], double k)
double determinant(double *matrix_need_inverse, double sizeof_square_inverse_matrix) // Part of inverse matrix
{
    double s = 1, det = 0;
    double b[inverse_matrix_size_global*inverse_matrix_size_global];
    int i, j, m, n, c;
    if (sizeof_square_inverse_matrix == 1)
    {
        return (matrix_need_inverse[0]);
    }
    else
    {
        det = 0;
        for (c = 0; c < sizeof_square_inverse_matrix; c++)
        {
            m = 0;
            n = 0;
            for (i = 0;i < sizeof_square_inverse_matrix; i++)
            {
                for (j = 0 ;j < sizeof_square_inverse_matrix; j++)
                {
                    //b[i][j] = 0;
                    b[(int)(i*inverse_matrix_size_global+j)] = 0;
                    if (i != 0 && j != c)
                    {
                        //b[m][n] = a[i][j];
                        b[m*inverse_matrix_size_global+n] = matrix_need_inverse[(int)(i*inverse_matrix_size_global+j)];
                        //printf("%f\n",b[m*inverse_matrix_size_global+n]);
                        if (n < (sizeof_square_inverse_matrix - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (matrix_need_inverse[c] * determinant(b, sizeof_square_inverse_matrix - 1));
            s = -1 * s;
        }
    }
    
    return (det);
}

void cofactor(double *matrix_need_inverse, double sizeof_square_inverse_matrix, double *desired_inversed_matrix) // Part of inverse matrix
{
    double b[inverse_matrix_size_global*inverse_matrix_size_global], fac[inverse_matrix_size_global*inverse_matrix_size_global];
    //double b[(int)(f*f)], fac[(int)(f*f)];
    int p, q, m, n, i, j;
    for (q = 0;q < sizeof_square_inverse_matrix; q++)
    {
        for (p = 0;p < sizeof_square_inverse_matrix; p++)
        {
            m = 0;
            n = 0;
            for (i = 0;i < sizeof_square_inverse_matrix; i++)
            {
                for (j = 0;j < sizeof_square_inverse_matrix; j++)
                {
                    if (i != q && j != p)
                    {
                        b[m*inverse_matrix_size_global+n] = matrix_need_inverse[(int)(i*inverse_matrix_size_global+j)];
                        if (n < (sizeof_square_inverse_matrix - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[(int)(q*sizeof_square_inverse_matrix+p)] = pow(-1, q + p) * determinant(b, sizeof_square_inverse_matrix - 1);
        }
    }
    transpose(matrix_need_inverse, fac, sizeof_square_inverse_matrix, desired_inversed_matrix);
}

/*Finding transpose of matrix*/
void transpose(double *matrix_need_inverse, double *fac, double sizeof_square_inverse_matrix, double *desired_inverse_matrix)
{
    int i, j;
    //double b[(int)(r*r)], inverse[(int)(r*r)], d;
    double b[inverse_matrix_size_global*inverse_matrix_size_global], inverse[inverse_matrix_size_global*inverse_matrix_size_global], d;
    
    for (i = 0;i < sizeof_square_inverse_matrix; i++)
    {
        for (j = 0;j < sizeof_square_inverse_matrix; j++)
        {
            b[(int)(i*inverse_matrix_size_global+j)] = fac[(int)(j*inverse_matrix_size_global+i)];
        }
    }
    d = determinant(matrix_need_inverse, sizeof_square_inverse_matrix);
    for (i = 0;i < sizeof_square_inverse_matrix; i++)
    {
        for (j = 0;j < sizeof_square_inverse_matrix; j++)
        {
            inverse[(int)(i*inverse_matrix_size_global+j)] = b[(int)(i*inverse_matrix_size_global+j)] / d;
        }
    }
    printf("\n\n\nThe inverse of matrix is : \n");
    
    for (i = 0;i < sizeof_square_inverse_matrix; i++)
    {
        for (j = 0;j < sizeof_square_inverse_matrix; j++)
        {
            //printf("\t%f", inverse[(int)(i*inverse_matrix_size_global+j)]);
            desired_inverse_matrix[(int)(i*inverse_matrix_size_global+j)]= inverse[(int)(i*inverse_matrix_size_global+j)];
        }
        //printf("\n");
    }
}
// This is the Ending point of calculation the INVERSE OF MATRIX
