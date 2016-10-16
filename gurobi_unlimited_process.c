# include <math.h>
# include <mpi.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <gurobi_c.h>

int main ( int argc, char *argv[] );

void p0_stop_decision(int *stop_decision, int stop_counter);
void p0_send_decision(int process_size,int stop_decision);

void p0_set_input ( double *input, int process_size, int coef_length );
void p0_send_input ( double *input, int process_size, int coef_length );
void p0_receive_output ( double *output, int process_size );

void p1_receive_decision(int *stop_decision_i, int id);
void p1_receive_input (double *input_i, int id, int coef_length);
double p1_compute_output ( double *input_i, int stop_decision_i);
void p1_send_output ( double output_i, int id );

void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
 Purpose:
 
 MAIN is the main program for MPI_MULTITASK.
 
 Discussion:
 
 Message tag 1: P0 sends input to P1
 Message tag 2: P0 sends input to P2
 Message tag 3: P1 sends output to P0.
 Message tag 4: P2 sends output to P0.
 */
{
    int id;
    int ierr;
    int coef_length=3;
    int p;
    double wtime;
    
    /*
     Process 0 is the "monitor".
     It chooses the inputs, and sends them to the workers.
     It waits for the outputs.
     It plots the outputs.
     */
    ierr = MPI_Init ( &argc, &argv );
    
    ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );
    
    ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );
    /*
     Make sure we have enough processes.
     */
    if ( p < 3 )
    {
        printf ( "\n" );
        printf ( "MPI_MULTITASK - Fatal error!\n" );
        printf ( "  Number of available processes must be at least 3!\n" );
        ierr = MPI_Finalize ( );
        exit ( 1 );
    }
    /*
     Run program P0 on process 0, and so on.
     */
    if ( id == 0 )
    {
        timestamp ( );
        
        printf ( "\n" );
        printf ( "MPI_MULTITASK:\n" );
        printf ( "  C / MPI version\n" );
        
        wtime = MPI_Wtime ( );//MPI_Wtime returns a floating-point number of seconds, representing elapsed wall-clock time since some time in the past.
        
        double input[(p-1)*coef_length];
        double output[(p-1)];
        int stop_decision;
        int stop_counter;
        
        for (stop_counter=0;  ; stop_counter++) // this loop will let program p0_stop_decision to decide ending point.
        {
            printf("the stop_counter is %d\n", stop_counter);
            p0_stop_decision( &stop_decision,  stop_counter);
            printf("stop decision is %d ", stop_decision);
            printf("\n");
            p0_send_decision( p, stop_decision);
            
            if (stop_decision == 0){
                break;}
            
            
            
            p0_set_input ( input, p, coef_length ); // p stands for numbers of processes
            
            
            
            
            
            p0_send_input ( input, p, coef_length );
            
            
            
            p0_receive_output ( output, p );
            
        }
        
        
        wtime = MPI_Wtime ( ) - wtime;
        printf ( "  Process 0 time = %g\n", wtime );
        
        
        
        ierr = MPI_Finalize ( );
        
        printf ( "\n" );
        printf ( "MPI_MULTITASK:\n" );
        printf ( "  Normal end of execution.\n" );
        
        
        timestamp ( );
    }
    
    /*
     Process 1 to n works on task 1 to n.
     It receives input from process 0.
     It computes the output.
     It sends the output to process 0.
     */
    else if ( id > 0 )
    {
        
        
        wtime = MPI_Wtime ( );
        double input_i[coef_length];
        double output_i;
        int stop_decision_i;
        
        
        while (1)
        {
            
            p1_receive_decision( &stop_decision_i, id);
            
            printf("stop decision received at subproblem %d is %d \n", id, stop_decision_i);
            
            if (stop_decision_i == 1)
            {
                
                
                printf("the loop is no error at this moment at process %d\n", id);
                
                p1_receive_input (input_i, id, coef_length);
                
                
                output_i = p1_compute_output ( input_i, stop_decision_i);
                
                
                
                p1_send_output ( output_i, id );
                
                wtime = MPI_Wtime ( ) - wtime;
                printf ( "  Process %d time = %g\n",id, wtime );//  For output, you're passing a value, which will be promoted from float to double when passed as a variadic parameter. For input you're passing a pointer, which is not promoted, so you have to tell scanf whether you want to read a float or a double, so for scanf, %f means you want to read a float and %lf means you want to read a double (and, for what it's worth, for a long double, you use %Lf for either printf or scanf).
            }
            else
            {
                printf("MPI communication will terminate \n");
                ierr = MPI_Finalize ( );//All processes must call this routine before exiting. All processes will still exist but may not make any further MPI calls.
                break;
            }
            
        }
        
    }
    
    //EXIT_SUBPROCESS:
    return 0;
    
}

/******************************************************************************/
void p0_stop_decision(int *stop_decision, int stop_counter)
/*
 Purpose: program stop decision
 */
/******************************************************************************/
{
    if (stop_counter <= 3 )
    {
        *stop_decision = 1; //program will continue
    }
    else
    {
        *stop_decision = 0; // program will stop
    }
    //printf("this is %d loop, and the stop decision is %d \n", stop_counter, *stop_decision);
    return;
}
/******************************************************************************/
void p0_send_decision(int process_size,int stop_decision)
/*
 Purpose: Master node send stop decision to slave nodes
 */
/******************************************************************************/
{
    int id;
    int ierr;
    int tag;
    int process_last = process_size - 1; // the last process's number
    
    for (int i=1; i <= process_last; i++)
    {
        id = i;
        tag = i;
        printf("the sent stop_decision%d is %d \n",id ,stop_decision);
        ierr = MPI_Send ( &stop_decision, 1, MPI_INT, id, tag, MPI_COMM_WORLD );
    }
    return;
}

/******************************************************************************/

void p0_set_input ( double *input, int process_size, int coef_length )

/******************************************************************************/
/*
 Purpose:
 
 P0_SET_INPUT sets input.
 
 Licensing:
 
 This code is distributed under the GNU LGPL license.
 
 Modified:
 
 21 October 2011
 
 Author:
 
 John Burkardt
 
 Parameters:
 
 Output, int *INPUT1, *INPUT2, the values of two
 inputs used by tasks 1 and 2.
 */
{
    double obj_coef[coef_length];
    int process_last = process_size - 1; // the last process's number
    
    obj_coef[0]=1; obj_coef[1]=1; obj_coef[2]=2;
    for (int i=1; i<=process_last; i++)
    {
        input[0+(i-1)*coef_length]=obj_coef[0]*i; // change the coef. of obj. function according to rank(rank times original setup)
        input[1+(i-1)*coef_length]=obj_coef[1]*i;
        input[2+(i-1)*coef_length]=obj_coef[2]*i;
    }
    /*
     
     for (int i=1; i<=process_last; i++) // this is used for checking setting parameter
     {
     printf("%2.1f   ",input[0+(i-1)*coef_length]);
     printf("%2.1f   ",input[1+(i-1)*coef_length]);
     printf("%2.1f"  ,input[2+(i-1)*coef_length]);
     printf("\n");
     }
     */
    return;
}
/******************************************************************************/

void p0_send_input ( double *input, int process_size, int coef_length )

/******************************************************************************/
/*
 Purpose:
 
 P0_SEND_INPUT sends input to processes 1 and 2.
 
 Licensing:
 
 This code is distributed under the GNU LGPL license.
 
 Modified:
 
 21 October 2011
 
 Author:
 
 John Burkardt
 
 Parameters:
 
 Input, int INPUT1, INPUT2, the values of two
 inputs used by tasks 1 and 2.
 */
{
    int id;
    int ierr;
    int tag;
    int process_last = process_size - 1; // the last process's number
    
    for (int i=1; i <= process_last; i++)
    {
        id = i;
        tag = i;
        //printf("id is %d, tag is %d \n", id, tag);
        
        ierr = MPI_Send ( &input[0+(i-1)*coef_length], 3, MPI_DOUBLE, id, tag, MPI_COMM_WORLD );
        
    }
    
    
    return;
}
/******************************************************************************/

void p0_receive_output ( double *output, int process_size )

/******************************************************************************/
/*
 Purpose:
 
 P0_RECEIVE_OUTPUT receives output from processes 1 and 2.
 
 Licensing:
 
 This code is distributed under the GNU LGPL license.
 
 Modified:
 
 21 October 2011
 
 Author:
 
 John Burkardt
 
 Parameters:
 
 Output, int OUTPUT1, OUTPUT2, the values of the
 outputs of tasks 1 and 2.
 */
{
    int ierr;
    int process_last = process_size - 1; // the last process's number
    MPI_Status status[process_size-1];//except process 0, which is master process
    
    
    for (int i=1; i<=process_last; i++)
    {
        ierr = MPI_Recv ( &output[i-1], 1, MPI_DOUBLE, i, 100+i,
                         MPI_COMM_WORLD, &status[i-1] ); // slave message id is (100+id)
        printf ( "\n" );
        printf ( "  Process %d returned OUTPUT%d = %g\n", i, i, output[i-1] );
    }
    return;
}

/******************************************************************************/
void p1_receive_decision(int *stop_decision_i, int id)
/*
 Purpose: receive the stop decision from master node
 */
/******************************************************************************/
{
    int recv_source = 0;
    int ierr;
    MPI_Status status;
    int tag;
    
    tag = id;
    
    ierr = MPI_Recv ( stop_decision_i, 1, MPI_INT, recv_source, tag, MPI_COMM_WORLD, &status );
    
    //printf("process %d receive stop decision is %d \n", id, *stop_decision_i);
    return;
}

/******************************************************************************/

void p1_receive_input (double *input_i, int id, int coef_length)

/******************************************************************************/
/*
 Purpose:
 
 P1_RECEIVE_INPUT receives input from process 0.
 
 Licensing:
 
 This code is distributed under the GNU LGPL license.
 
 Modified:
 
 21 October 2011
 
 Author:
 
 John Burkardt
 
 Parameters:
 
 Output, int P1_RECEIVE_INPUT, the value of the parameter.
 */
{
    int recv_source = 0;
    int ierr;
    MPI_Status status;
    int tag;
    
    tag = id;
    
    ierr = MPI_Recv ( input_i, coef_length, MPI_DOUBLE, recv_source, tag, MPI_COMM_WORLD, &status );
    
    /*
     for (int i =0; i<3 ; i++)
     {
     printf(" the coeff.  is %f \n", input_i[i] );
     }
     */
    
    return ;
}
/******************************************************************************/

double p1_compute_output ( double *input_i, int stop_decision_i)

/******************************************************************************/
/*
 Purpose:
 P1_COMPUTE_OUTPUT carries out computation number 1. 
 Output, int P1_COMPUTE_OUTPUT1, the problem output.
 */
{
    static int initial_setup=1; // This will get rid of build model and add variable all the time
    
    
    static double output=NAN;
    static GRBenv   *env   = NULL;
    static GRBmodel *model = NULL;
    static int       error = 0;
    static double    sol[3];
    static int       ind[3];
    static double    val[3];
    static double    obj[3];
    static char      vtype[3];
    static int       optimstatus;
    static double    objval;
    
    
    
    
    if (initial_setup == 1) // start envirnoment setup
    {
        
        
        initial_setup = 0; // close initial setup when it already done
        /* Create environment */
        
        error = GRBloadenv(&env, "mip1.log");
        if (error) goto QUIT;
        
        /*SET output is none after optimization*/
        error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
        if (error) goto QUIT;
        
        /* Create an empty model */
        
        error = GRBnewmodel(env, &model, "mip1", 0, NULL, NULL, NULL, NULL, NULL);
        if (error) goto QUIT;
        
        /* Add variables */
        for(int i=0;i<3;i++)
        {
            obj[i]= *(input_i);
        }
        //obj[0] = 1; obj[1] = 1; obj[2] = 2;
        
        vtype[0] = GRB_BINARY; vtype[1] = GRB_BINARY; vtype[2] = GRB_BINARY;
        error = GRBaddvars(model, 3, 0, NULL, NULL, NULL, obj, NULL, NULL, vtype,
                           NULL);
        
        // program update
        error = GRBupdatemodel(model);
        if (error) goto QUIT;
        
    }
    
    
    else            // after envirnoment setup, no need setup again
    {
        printf("this is the second time that program come in \n");
        int change_index[3];
        // change the coef. of variables
        for(int i=0;i<3;i++)
        {
            obj[i]= *(input_i) ;
            printf("the updated coef. is %f \n", obj[i]);
            change_index[i] = i;
        }
        error = GRBsetdblattrlist(model, GRB_DBL_ATTR_OBJ, 3, change_index, obj);
        if (error) goto QUIT;
        // program update
        error = GRBupdatemodel(model);
        if (error) goto QUIT;
        
        printf("up to this point, program is fine \n");
        
        double obj_updated[3];
        error = GRBgetdblattrlist(model, GRB_DBL_ATTR_OBJ, 3, change_index, obj_updated);
        
        /*
         for(int i=0;i<3;i++)
         {
         printf("the updated coef. is %g", obj_updated[i]);
         }
         */
    }
    
    
    
    
    /* Change objective sense to maximization */
    
    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
    if (error) goto QUIT;
    
    /* Integrate new variables */
    
    error = GRBupdatemodel(model);
    if (error) goto QUIT;
    
    
    /* First constraint: x + 2 y + 3 z <= 4 */
    
    ind[0] = 0; ind[1] = 1; ind[2] = 2;
    val[0] = 1; val[1] = 2; val[2] = 3;
    
    error = GRBaddconstr(model, 3, ind, val, GRB_LESS_EQUAL, 4.0, "c0");
    if (error) goto QUIT;
    
    /* Second constraint: x + y >= 1 */
    
    ind[0] = 0; ind[1] = 1;
    val[0] = 1; val[1] = 1;
    
    error = GRBaddconstr(model, 2, ind, val, GRB_GREATER_EQUAL, 1.0, "c1");
    if (error) goto QUIT;
    
    /* Optimize model */
    
    error = GRBoptimize(model);
    if (error) goto QUIT;
    
    /* Write model to 'mip1.lp' */
    
    error = GRBwrite(model, "mip1.lp");
    if (error) goto QUIT;
    
    /* Capture solution information */
    
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
    if (error) goto QUIT;
    
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if (error) goto QUIT;
    
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 3, sol);
    if (error) goto QUIT;
    
    
    //printf("\nOptimization complete\n");
    if (optimstatus == GRB_OPTIMAL) {
        printf("Optimal objective: %.4e\n", objval);  // silence the output result
        //printf("  x=%.0f, y=%.0f, z=%.0f\n", sol[0], sol[1], sol[2]);
        output = objval;
    } else if (optimstatus == GRB_INF_OR_UNBD) {
        printf("Model is infeasible or unbounded\n");
        output = INFINITY;
    } else {
        printf("Optimization was stopped early\n");
        output = NAN;
    }
    
    
QUIT:
    
    /* Error reporting */
    
    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }
    
    
    if (stop_decision_i==0) // stop_decision = 0 means calculation will stop, and envirnoment, model will be release
    {
        /* Free model */
        
        GRBfreemodel(model);
        
        /* Free environment */
        
        GRBfreeenv(env);
    }
    
    
    return output;
}
/******************************************************************************/

void p1_send_output ( double output_i, int id )

/******************************************************************************/
/*
 Purpose:
 
 P1_SEND_OUTPUT sends output to process 0.
 
 */
{
    int dest;
    int ierr;
    int tag;
    
    dest = 0;
    tag = 100 + id; // the slave send id is (100+id)
    //printf(" process %d will send output %f \n", id, output_i);
    ierr = MPI_Send ( &output_i, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD );
    
    return;
}

void timestamp ( )

/******************************************************************************/
/*
 Purpose:
 TIMESTAMP prints the current YMDHMS date as a time stamp.
 */
{
# define TIME_SIZE 40
    
    static char time_buffer[TIME_SIZE];//it exist through program execution. static storage duration, block scope, and no linkage. it's initialize just onece, when timestamp compiled. And set to 0 if not initialize it.
    const struct tm *tm; //value can't be modified by assignment or by incrementing or decrementing
    //the value of tm can be changed but tm points to a value that must remain constant
    time_t now;//type suitable for storing the calendar time
    
    now = time ( NULL );//calculate the current calender time and encodes it into time_t format
    tm = localtime ( &now );//the value of timer is broken up into the structure tm and expressed in the local time zone
    
    strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );//Formats the time represented in the structure timeptr according to the formatting rules defined in format and stored into str.
    
    printf ( "%s\n", time_buffer );
    
    return;
# undef TIME_SIZE
}


