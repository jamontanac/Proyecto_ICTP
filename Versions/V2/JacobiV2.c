#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>

#define TAG_1 100
#define TAG_2 200
#define MPI_PROC_ROOT 0

void save_gnuplot( double *M, size_t dim );
void evolve( double * matrix, double *matrix_new, int dimension, int loc_size );
double seconds( void );
void Exchange_information(double *matrix,int rank,int nproc,int dimension,int loc_size);
void Set_Boundary_Conditions(double * matrix,double *matrix_new, int loc_size, int dimension,int rank,int nproc,int offset);
void Print_matrix(double *matrix, int loc_size, int dimension);


int main(int argc, char* argv[])
{
  
  int i, j, it;
  
  
  double *matrix, *matrix_new, *tmp_matrix;
  
  
  int dimension = 0, iterations = 0;
  int row_peek = 0, col_peek = 0; //To debug if everything is going correct
  int rank, nproc;
  int loc_size;
  size_t byte_dimension = 0;
  int num_rank;
  int rest, offset;
  double t_start, time_local, Total_time;
  
  dimension = atoi(argv[1]);
  iterations = atoi(argv[2]);
  
  if(argc != 3)
    {
      fprintf(stderr,"\n Wrong number of arguments. Usage: ./a.out dim it n m\n");
      return 1;
    }
  
  
  
  
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &nproc) ;
  
  // Setting a balanced distribution of the matrix
  
  loc_size = dimension/nproc;
  rest = offset = dimension % nproc;
  if(rank < rest){
    loc_size += 1;
    offset = 0;
  }
  
  byte_dimension = sizeof(double) * ( loc_size + 2 ) * ( dimension + 2 );
  matrix = ( double* ) malloc( byte_dimension );
  matrix_new = ( double* ) malloc( byte_dimension );
  memset( matrix, 0, byte_dimension );
  memset( matrix_new, 0, byte_dimension );
  
  //fill initial values
#pragma omp parallel for
  for( i = 1; i <= loc_size; ++i )
    for( j = 1; j <= dimension; ++j )
      matrix[ ( i * ( dimension + 2 ) ) + j ] = 0.5;
  
  
  Set_Boundary_Conditions(matrix,matrix_new,loc_size,dimension,rank,nproc,offset);
  
  MPI_Request request_1;
  MPI_Status  status_1;
  
  MPI_Request request_2;
  MPI_Status  status_2;
  
  if(rank==MPI_PROC_ROOT){
    int Destination_Source = rank+1;
    MPI_Irecv(&matrix[(loc_size+1)*(dimension+2)],dimension+2,MPI_DOUBLE,Destination_Source,TAG_1,MPI_COMM_WORLD,&request_1);
    
    // Evolve the lower frontier
#pragma omp parallel for
    for(i = 1; i< dimension;i++){
      matrix_new[(loc_size)*(dimension+2)+i] = ( 0.25 ) *
	( matrix[ ( ( loc_size - 1 ) * ( dimension + 2 ) ) + i ]
	  + matrix[ ( loc_size * ( dimension + 2 ) ) + ( i + 1 ) ] +
	  matrix[ ( ( loc_size + 1 ) * ( dimension + 2 ) ) + i ] +
	  matrix[ ( loc_size * ( dimension + 2 ) ) + ( i - 1 ) ] );
    }
    
    
    MPI_Isend(&matrix[(loc_size)*(dimension+2)],dimension+2,MPI_DOUBLE,Destination_Source,TAG_2,MPI_COMM_WORLD,&request_1);
    
    
    //Evolve the Bulk
#pragma omp parallel for
    for( i = 1 ; i <= loc_size-1; ++i )
      for( j = 1; j <= dimension; ++j )
	matrix_new[ ( i * ( dimension + 2 ) ) + j ] = ( 0.25 ) * 
	  ( matrix[ ( ( i - 1 ) * ( dimension + 2 ) ) + j ] + 
	    matrix[ ( i * ( dimension + 2 ) ) + ( j + 1 ) ] + 	  
	    matrix[ ( ( i + 1 ) * ( dimension + 2 ) ) + j ] + 
	    matrix[ ( i * ( dimension + 2 ) ) + ( j - 1 ) ] ); 
    
    MPI_Wait(&request_1, &status_1);
  }

  
  else if(rank==nproc-1){
    int Destination_Source = rank-1;
    MPI_Irecv(&matrix[0],dimension+2,MPI_DOUBLE,Destination_Source,TAG_1,MPI_COMM_WORLD,&request_2);
    
    // Evolve the upper frontier
#pragma omp parallel for
    for(i = 1; i< dimension;i++){
      matrix_new[(dimension+2)+i] = ( 0.25 ) *
	( matrix[ i ]
	  + matrix[ (( dimension + 2 ) ) + ( i + 1 ) ] +
	  matrix[ ( ( 1+1 ) * ( dimension + 2 ) ) + i ] +
	  matrix[ (( dimension + 2 ) ) + ( i - 1 ) ] );
    }
    
    MPI_Isend(&matrix[(dimension+2)],dimension+2,MPI_DOUBLE,Destination_Source,TAG_2,MPI_COMM_WORLD,&request_2);
    
    //Evolve the Bulk
#pragma omp parallel for
    for( i = 2 ; i <= loc_size; ++i )
      for( j = 1; j <= dimension; ++j )
	matrix_new[ ( i * ( dimension + 2 ) ) + j ] = ( 0.25 ) * 
	  ( matrix[ ( ( i - 1 ) * ( dimension + 2 ) ) + j ] + 
	    matrix[ ( i * ( dimension + 2 ) ) + ( j + 1 ) ] + 	  
	    matrix[ ( ( i + 1 ) * ( dimension + 2 ) ) + j ] + 
	    matrix[ ( i * ( dimension + 2 ) ) + ( j - 1 ) ] );
    
    MPI_Wait(&request_2, &status_2);
  }


  
  else{
    int Destination_Source_up=rank-1;
    int Destination_Source_down = rank+1;
    
    //Receiving
    MPI_Irecv(&matrix[(loc_size+1)*(dimension+2)],dimension+2,MPI_DOUBLE,Destination_Source_down,TAG_1,MPI_COMM_WORLD,&request_1);
    MPI_Irecv(&matrix[0],dimension+2,MPI_DOUBLE,Destination_Source_down,TAG_1,MPI_COMM_WORLD,&request_2);
    
    // Evolve the lower frontier
#pragma omp parallel for
    for(i = 1; i< dimension;i++){
      matrix_new[(loc_size)*(dimension+2)+i] = ( 0.25 ) *
	( matrix[ ( ( loc_size - 1 ) * ( dimension + 2 ) ) + i ]
	  + matrix[ ( loc_size * ( dimension + 2 ) ) + ( i + 1 ) ] +
	  matrix[ ( ( loc_size + 1 ) * ( dimension + 2 ) ) + i ] +
	  matrix[ ( loc_size * ( dimension + 2 ) ) + ( i - 1 ) ] );
    }
    
    // Evolve the upper frontier
#pragma omp parallel for
    for(i = 1; i< dimension;i++){
      matrix_new[(dimension+2)+i] = ( 0.25 ) *
	( matrix[ i ]
	  + matrix[ (( dimension + 2 ) ) + ( i + 1 ) ] +
	  matrix[ ( ( 1+1 ) * ( dimension + 2 ) ) + i ] +
	  matrix[ (( dimension + 2 ) ) + ( i - 1 ) ] );	
    }
    
    MPI_Isend(&matrix[(dimension+2)],dimension+2,MPI_DOUBLE,Destination_Source_up,TAG_2,MPI_COMM_WORLD,&request_2);
    MPI_Isend(&matrix[(loc_size)*(dimension+2)],dimension+2,MPI_DOUBLE,Destination_Source_down,TAG_2,MPI_COMM_WORLD,&request_1);
    
    //evolve the bulk
#pragma omp parallel for
    for( i = 2 ; i <= loc_size-1; ++i )
      for( j = 1; j <= dimension; ++j )
	matrix_new[ ( i * ( dimension + 2 ) ) + j ] = ( 0.25 ) * 
	  ( matrix[ ( ( i - 1 ) * ( dimension + 2 ) ) + j ] + 
	    matrix[ ( i * ( dimension + 2 ) ) + ( j + 1 ) ] + 	  
	    matrix[ ( ( i + 1 ) * ( dimension + 2 ) ) + j ] + 
	    matrix[ ( i * ( dimension + 2 ) ) + ( j - 1 ) ] );
    
    MPI_Wait(&request_2, &status_2);
    MPI_Wait(&request_1, &status_1);
    
  }
  
  
  
  
  /* t_start = seconds(); */
  /* for( it = 0; it < iterations; ++it ) */
  /*   { */
  /*     Exchange_information(matrix,rank,nproc,dimension,loc_size); */
  /*     evolve( matrix, matrix_new, dimension,loc_size ); */
  /*     tmp_matrix = matrix; */
  /*     matrix = matrix_new; */
  /*     matrix_new = tmp_matrix; */
  /*   } */
  /* time_local = seconds()-t_start; */
  /* Total_time=0; */
  
  /* MPI_Reduce(&time_local,&Total_time,1,MPI_DOUBLE,MPI_SUM,MPI_PROC_ROOT,MPI_COMM_WORLD);   */
  /* if(rank==MPI_PROC_ROOT) */
  /*   printf("%d\t%.15f\n",nproc,Total_time/nproc); */
  
  
  if(rank==MPI_PROC_ROOT){
  
    FILE *file;
    file = fopen("data.txt","w");
    for(i =1;i <= loc_size;i++){
      for(j =1 ;j <= dimension;j++){
  	fprintf(file,"%f\t",matrix[i*(dimension+2)+j]);
      }
      fprintf(file,"\n");
    }
    /* Print_matrix(matrix,loc_size,dimension); */
    
    for(int pe = 1; pe < nproc; pe++){
      MPI_Recv(matrix,(dimension+2)*(loc_size+2),MPI_DOUBLE,pe,pe,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  
      if( rest && pe >= rest){
  	for(i =1;i <= loc_size;i++){
  	  for(j =1 ;j <= dimension;j++){
  	    fprintf(file,"%f\t",matrix[i*(dimension+2)+j]);
  	  }
  	  fprintf(file,"\n");
  	}
      }
  	else{
  	  for(i =1;i <= loc_size;i++){
  	    for(j =1 ;j <= dimension;j++){
  	      fprintf(file,"%f\t",matrix[i*(dimension+2)+j]);
  	    }
  	    fprintf(file,"\n");
  	  }
  	}/* Print_matrix(matrix,loc_size,dimension); */
    }
    fclose( file );
  }
  else
    MPI_Send(matrix,(dimension+2)*(loc_size+2),MPI_DOUBLE,MPI_PROC_ROOT,rank,MPI_COMM_WORLD);
  
  free( matrix );
  free( matrix_new );
	
  MPI_Finalize();
  return 0;
}


/* -------------------------Evolve -------------------- */
void Print_matrix(double *matrix, int loc_size, int dimension)
{
  int i,j;
  for(i =1;i <= loc_size;i++){
  	for(j =1 ;j <= dimension;j++){
  	  printf("%f\t",matrix[i*(dimension+2)+j]);
  	}
  	printf("\n");
      }
}

void evolve( double * matrix, double *matrix_new, int dimension, int loc_size )
{  
  size_t i , j;
  #pragma omp parallel for
  for( i = 1 ; i <= loc_size; ++i )
    for( j = 1; j <= dimension; ++j )
      matrix_new[ ( i * ( dimension + 2 ) ) + j ] = ( 0.25 ) * 
	( matrix[ ( ( i - 1 ) * ( dimension + 2 ) ) + j ] + 
	  matrix[ ( i * ( dimension + 2 ) ) + ( j + 1 ) ] + 	  
	  matrix[ ( ( i + 1 ) * ( dimension + 2 ) ) + j ] + 
	  matrix[ ( i * ( dimension + 2 ) ) + ( j - 1 ) ] ); 
}

/* -------------------------Set Boundaries -------------------- */
void Set_Boundary_Conditions(double * matrix,double *matrix_new, int loc_size, int dimension,int rank,int nproc,int offset)
{
  int i;
  double increment = 100.0 / ( dimension + 1 );
  #pragma omp parallel for
  for(i = 1; i<=loc_size; i++){
    matrix[i*(dimension + 2)] = (rank*loc_size + i + offset)*increment;
    matrix_new[i*(dimension + 2)] = (rank*loc_size + i + offset)*increment;
    }
  
  if(rank==nproc-1){
    #pragma omp parallel for
    for( i=1; i <= dimension+1; ++i ){
      matrix[ ( ( loc_size + 1 ) * ( dimension + 2 ) ) + ( dimension + 1 - i ) ] = i * increment;
      matrix_new[ ( ( loc_size + 1 ) * ( dimension + 2 ) ) + ( dimension + 1 - i ) ] = i * increment;
    }
  }
}

/* -------------------------Prepare for gnuplot -------------------- */

void save_gnuplot( double *M, size_t dimension )
{
  size_t i , j;
  const double h = 0.1;
  FILE *file;
  file = fopen( "solution.dat", "w" );
  for( i = 0; i < dimension + 2; ++i )
    for( j = 0; j < dimension + 2; ++j )
      fprintf(file, "%f\t%f\t%f\n", h * j, -h * i, M[ ( i * ( dimension + 2 ) ) + j ] );
  fclose( file );
}


/* -------------------------Measure the time -------------------- */
double seconds()
{
    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}
/* -------------------------Communicate information  -------------------- */
void Exchange_information(double *matrix,int rank,int nproc,int dimension,int loc_size)
{
  if(rank==MPI_PROC_ROOT){
      int Destination_Source = rank+1;
      MPI_Sendrecv(&matrix[(loc_size)*(dimension+2)],dimension+2,
  		   MPI_DOUBLE,Destination_Source,TAG_1,
  		   &matrix[(loc_size+1)*(dimension+2)],dimension+2,
  		   MPI_DOUBLE,Destination_Source,TAG_2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
  
  else if(rank ==nproc-1){
    int Destination_Source = rank-1;
    MPI_Sendrecv(&matrix[(dimension+2)],dimension+2,
  		 MPI_DOUBLE,Destination_Source,TAG_2,
  		 &matrix[0],dimension+2,
  		 MPI_DOUBLE,Destination_Source,TAG_1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  else{
    int Destination_Source_up=rank-1;
    int Destination_Source_down = rank+1;

    
    //Sending up
    MPI_Sendrecv(&matrix[(dimension+2)],dimension+2,
  		 MPI_DOUBLE,Destination_Source_up,TAG_2,
  		 &matrix[0],dimension+2,
  		 MPI_DOUBLE,Destination_Source_up,TAG_1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    // Sending down
    MPI_Sendrecv(&matrix[(loc_size)*(dimension+2)],dimension+2,
  		   MPI_DOUBLE,Destination_Source_down,TAG_1,
  		   &matrix[(loc_size+1)*(dimension+2)],dimension+2,
  		   MPI_DOUBLE,Destination_Source_down,TAG_2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
  }
}
