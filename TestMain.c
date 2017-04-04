
#include <sys/times.h>
#include <stdlib.h>
#include <stdio.h>
#include <semaphore.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>


double start, stop, used, mf;

double ftime(void);
//void multiply (double **a, double **b, double **c, int n);

sem_t *semaphores;
double **a, **b;
double *c;

double ftime (void)
{
    struct tms t;

    times ( &t );

    return (t.tms_utime + t.tms_stime) / 100.0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function     : printMatrix
// Description  : helper function to print a single matrix
//
// Inputs       : a double pointer to a matrix, the size
// Outputs      : N/A
void printMatrix(double **matrix, int size){
	int i, j;
	for(i = 0; i <size; i++){
		for(j = 0; j < size; j++){
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}
		printf("\n");
}
void printOneDMatrix(double *matrix, int size){
	int i, j;
	for(i = 0; i <size; i++){
		for(j = 0; j < size; j++){
			printf("%f ", matrix[i*size + j]);
		}
		printf("\n");
	}
		printf("\n");

}


////////////////////////////////////////////////////////////////////////////////
//
// Function     : min
// Description  : helper function to determine the smaller of two values
//
// Inputs       : two ints
// Outputs      : the smaller int
int min(int i, int j){
	if(i>j){
		return j;
	}
	else{
		return i;
	}
}


void multiply (double **a, double **b, double *c, int n)
{
   int i, j, k;

   for (i=0; i<n; i++)
   {
     for (j=0; j<n; j++)

         c[i*n+j] = 0;
    }

    for (i=0; i<n; i++)
    {
       for (j=0; j<n; j++)
       {
         for (k=0; k<n; k++)
           c[i*n+j]= c[i*n+j] + a[i][k] * b[k][j];
        }
     }
  }


////////////////////////////////////////////////////////////////////////////////
//
// Function     : transpose
// Description  : reformats the matrix by switching rows and columns
//
// Inputs       : double pointer to matrix
// Outputs      : N/A
void transpose(double **matrix, int size){
	int i, j;
	double temp;

	for(i = 0; i < size; i++){
		for(j = i+1; j<size;j++){
			temp = matrix[i][j];
			matrix[i][j] = matrix[j][i];
			matrix[j][i] = temp;
		}
	}

}


////////////////////////////////////////////////////////////////////////////////
//
// Function     : transposeMult
// Description  : multiplies matrices with locality
//
// Inputs       : a double pointer to a regular matrix, a double pointer to a transposed matrix, a matrix to output to, and the size of the matrices
// Outputs      : N/A
void transposeMult(double **matrix, double **transMatrix, double *output, int size){
	int i, j, k;

		//fill output with zeros
	   for (i=0; i<size; i++)
	   {
	     for (j=0; j<size; j++)

	         output[i*size+j] = 0;
	    }

	   //multiply
	   for (i=0; i<size; i++)
	   {
		   for (j=0; j<size; j++)
	       {
			   for (k=0; k<size; k++)
			   {
				   output[i*size+j]= output[i*size+j] + matrix[i][k] * transMatrix[j][k];
			   }
	       }
	   }
}


////////////////////////////////////////////////////////////////////////////////
//
// Function     : blockMult
// Description  : splits matrices into blocks and multiplies the smaller chunks
//
// Inputs       : two double pointers to matrices, a matrix to output to, the block size, and the size of the matrices
// Outputs      : N/A
void blockMult(double **matrixA, double **matrixB, double *output, int blockSize, int size){
    int i, j, k, l, m, n, imin=0, jmin=0, kmin=0;


    //Fill output matrices with zeros
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            output[i*size+j]=0;
        }
    }

    //Blocking
    for(i=0; i < size; i += blockSize){

        imin = min(i + blockSize, size);

        for(j=0; j < size; j += blockSize){

            jmin = min(j + blockSize, size);

            for(k=0; k < size; k += blockSize){

                kmin = min( k + blockSize, size);
                //for each block
                for (l=i; l < imin; l++){
                    for(m=j; m < jmin; m++){
                        for(n=k; n < kmin; n++){
                            output[l*size+m] = output[l*size+m] + matrixA[l][n] * matrixB[n][m];
                        }
                    }
                }
            }
        }
    }
}



////////////////////////////////////////////////////////////////////////////////
//
// Function     : setupSharedMem
// Description  : creates semaphores and shared memory for a 1D matrix 
//
// Inputs       : int size for the dimensions of the matrix
// Outputs      : N/A
void setupSharedMem(int size){
	int shmfd;

	shmfd = shm_open("/mcdermottsjSemaphores", O_RDWR | O_CREAT, 0666);
	if(shmfd < 0){
		printf("error creating semaphores. Error code %d\n", shmfd);
		exit(1);
	}
	ftruncate(shmfd, size * size * sizeof(sem_t));
	semaphores = (sem_t *)mmap (NULL, size*size*sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED, shmfd,0);
	if(semaphores == NULL){
		printf("semaphores failed");
		exit(1);
	}
	close(shmfd);
	shm_unlink("/mcdermottsjSemaphores");
	//TODO remove this vvv

/*
	shmfd = shm_open("/mcdermottsjSharedA", O_RDWR | O_CREAT, 0666);
	if(shmfd < 0){
		printf("error creating shared mem A. Error code %d\n", shmfd);
		exit(1);
	}
	ftruncate(shmfd, size * size * sizeof(double));
	sharedA = (double (*))mmap (NULL, size * size * sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED, shmfd,0);
	if(sharedA == NULL){
		printf("shared mem A failed");
		exit(1);
	}
	close(shmfd);
	shm_unlink("/mcdermottsjSharedA");



	shmfd = shm_open("/mcdermottsjSharedB", O_RDWR | O_CREAT, 0666);
	if(shmfd < 0){
		printf("error creating shared mem B. Error code %d\n", shmfd);
		exit(1);
	}
	ftruncate(shmfd, size * size * sizeof(double));
	sharedB = (double (*))mmap (NULL, size*size*sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED, shmfd,0);
	if(sharedB == NULL){
		printf("shared mem B failed");
		exit(1);
	}
	close(shmfd);
	shm_unlink("/mcdermottsjSharedB");*/




	shmfd = shm_open("/mcdermottsjSharedC", O_RDWR | O_CREAT, 0666);
	if(shmfd < 0){
		printf("error creating shared mem C. Error code %d\n", shmfd);
		exit(1);
	}
	ftruncate(shmfd, size * size * sizeof(double));
	c = (double (*))mmap (NULL, size*size*sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED, shmfd,0);
	if(c == NULL){
		printf("shared mem C failed");
		exit(1);
	}
	close(shmfd);
	shm_unlink("/mcdermottsjSharedC");
}



////////////////////////////////////////////////////////////////////////////////
//
// Function     : initSems 
// Description  : intitalizes semaphores 
//
// Inputs       : int size for the dimensions of the matrix
// Outputs      : N/A
void initSems(int size){
	int i;
	for(i = 0; i < size*size; i++){
		if(sem_init(semaphores+i, 1, 1)<0){
			printf("semaphore init error.");
			exit(1);
		}
	}
}



////////////////////////////////////////////////////////////////////////////////
//
// Function     : threadedBlockMult
// Description  : splits matrices into blocks and multiplies the smaller chunks using multiple processes
//
// Inputs       : two double pointers to matrices, a matrix to output to, the block size, and the size of the matrices
// Outputs      : N/A
void threadedBlockMult(double **matrixA, double **matrixB, double *output, int blockSize, int size){
    int i, j, k, l, m, n, imin=0, jmin=0, kmin=0;
	int count =0;

    //Fill output matrices with zeros
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            output[i*size+j]=0;
        }
    }

    //Blocking
    for(i=0; i < size; i += blockSize){

        imin = min(i + blockSize, size);

        for(j=0; j < size; j += blockSize){

            jmin = min(j + blockSize, size);

            for(k=0; k < size; k += blockSize){

                kmin = min( k + blockSize, size);
                
		count++;
		int childPid = fork();
		if(childPid < 0){
			printf("fork faild at %d\n", i);
			printf("Use larger block size\n");
			exit(1);
		}
		else if (childPid > 0){//parent func
			//something for parent
		}
		else { //child block run this
		//for each block
                	for (l=i; l < imin; l++){
                    		for(m=j; m < jmin; m++){
                        		for(n=k; n < kmin; n++){
						sem_wait(&semaphores[l*size+m]);
	                            		output[l*size+m] = output[l*size+m] + matrixA[l][n] * matrixB[n][m];
						sem_post(&semaphores[l*size+m]);
                	        	}
                   		 }
               		}
			exit(0); 
		}
            }
        }
    }
    for(i = 0; i < count; i++){
    	wait(NULL);
    }	
}



int main (void)
{
     int i, j, n, blockSize;
	//srand(time(NULL));TODO uncomment
	

     printf ( "Enter the value of n: ");
     scanf ( "%d", &n);
     printf("Enter the desired block size: ");
     scanf("%d", &blockSize);

	printf("\nn size: %d\nblocksize: %d\n", n, blockSize);

     //Populate arrays....
     a= (double**)malloc(n*sizeof(double));
     b= (double**)malloc(n*sizeof(double));

	setupSharedMem(n);


     for (i=0; i<n; i++)
     {
       a[i]= malloc(sizeof(double)*n);
       b[i]= malloc(sizeof(double)*n);
      }

     for (i=0; i<n; i++)
     {
        for (j=0; j<n; j++)
         a[i][j]=(rand() % (99 + 1 - 0) +0 );
      }

     for (i=0; i<n; i++)
     {
        for (j=0; j<n; j++)
          b[i][j]=(rand() % (99 + 1 - 0) +0 );
      }
      
      //todo remove

     for (i=0; i<n; i++)
     {
        for (j=0; j<n; j++)
         c[i*n+j]=0;
      }


      //standard multiplication
      start = ftime();
      multiply (a,b,c,n);
      stop = ftime();
      used = stop - start;
      mf = (n*n*n *2.0) / used / 1000000.0;
      printf ("\n");
      printf ("Standard Multiplication:\n");
      printf ( "Elapsed time:   %10.2f \n", used);
      printf ( "DP MFLOPS:       %10.2f \n", mf);

      //transposed multiplication
      start = ftime();

      transpose(b, n);
      transposeMult(a,b,c,n);
	transpose(b,n);
      stop = ftime();
      used = stop - start;
      mf = (n*n*n * 2.0) / used / 1000000.0;
      printf ("\n");
      printf ("Transposed Multiplication:\n");
      printf ( "Elapsed time:   %10.2f \n", used);
      printf ( "DP MFLOPS:       %10.2f \n", mf);
	

      //Blocked multiplication
	start = ftime();
      blockMult(a,b,c,blockSize, n);
      stop = ftime();
      used = stop - start;
      mf = (n*n*n *2.0) / used / 1000000.0;//TODO check this
      printf ("\n");
      printf ("Block Multiplication:\n");
      printf ( "Elapsed time:   %10.2f \n", used);
      printf ( "DP MFLOPS:       %10.2f \n", mf);
	

	//threaded mult
	start = ftime();
	initSems(n);
      threadedBlockMult(a,b,c,blockSize, n);
      stop = ftime();
      used = stop - start;
      mf = (n*n*n *2.0) / used / 1000000.0;//TODO check this
      printf ("\n");
      printf ("Threaded Block Multiplication:\n");
      printf ( "Elapsed time:   %10.2f \n", used);
      printf ( "DP MFLOPS:       %10.2f \n", mf);

	printf("\n\n");
      return (0);
}


