//
//  read_excel.c
//  TEP_gurobi
//
//  Created by zhangcaihua on 9/19/16.
//  Copyright © 2016 zhangcaihua. All rights reserved.
//

#include "read_excel.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gurobi_c.h"
#include <string.h>

int file_size(FILE * fstream, int *row, int *col)
{
    FILE * fstream_inside = fstream;
    char buffer[BUFFER_MAX] ;
    char *record,*line;
    int i=0,j=0;
    
    if(fstream_inside == NULL)   {
        printf("\n file opening failed ");
        return -1 ;
    }
    while((line=fgets(buffer,sizeof(buffer),fstream_inside))!=NULL)
    {
        record = strtok(line,",");
        //printf("file position %ld\n", ftell(fstream));
        while(record != NULL && i == 0) //just need to test first line
        {
            j++;
            record = strtok(NULL,",");
            
        }
        i++ ;
    }
    *row = i;
    *col = j;
    
    rewind(fstream_inside); // It is vary important that to restore file pointer back to start point when finish
                            // Or it can assign it's value to a new variable and keep origional pointer safe.
                            // Besides, fseek(fp, 0, SEEK_SET) and rewind(fp) have same effect.
    return 0;
}

void delimiter_array(int* space_delimiter,int row, int col)
{
    for (int i = 0; i < col; i++)
    {
        space_delimiter[i]=i*row;
    }
}


int Cline_data_read(FILE * fstream, int row, int col, Cline_struct* info)
{
    FILE * fstream_inside = fstream;
    char buffer[BUFFER_MAX] ;
    char *record,*line;
    int i=0,j=0;
    
    if(fstream_inside == NULL)   {
        printf("\n file opening failed ");
        return -1 ;
    }
    

    while((line=fgets(buffer,sizeof(buffer),fstream_inside))!=NULL)
    {
        record = strtok(line,",");
        j = 0;
        while(record != NULL)
        {
            info->line_Cfrom[i+j*row] = atof(record);
            ++j;
            record = strtok(NULL,",");
        }
        ++i ;
    }
    
    printf("\n");
    rewind(fstream_inside);
    return 0 ;
}



