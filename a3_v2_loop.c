#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

#define print 0
#define debug 0
#define verify 0

double MedOfThree(double *A, int left, int right);
void seqQsort(double *A, int length);
void printArray(double *A, int length, int rank);
int verify_sort(double *A, int length);


int main(int argc, char **argv)
{
	/* Initialization*/
	int p, i, rank, length, pal, group_size, array_size, mod, ii, loopy, loops;
	double *A_glob, *A, *buff, piv, begin_time, end_time, fastest_time;
	int *send_count, *displs;
	p = atoi(argv[1]);
	if (p != 1)
	{
		if (p % 2 != 0)
		{
			printf("Number of processes must be a power of 2");
			return -1;
		}
	}
	length = atoi(argv[2]);
	loops = atoi(argv[3]);
	MPI_Request send_req;
	MPI_Status status;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &group_size);
	


	/* Arranging individual array sizes*/
	mod = length % group_size;
	array_size = (int)(length/group_size);
	if (mod != 0)
	{
		for (i = 0; i < mod; ++i)
		{
			if (rank == i)
			{
				array_size++;
			}
		}
		#if debug
		printf("My rank is = %d and my array_size is = %d \n", rank, array_size);
		#endif	
	}
	for (loopy = 0; loopy < loops; loopy++)
	{
		if (rank == 0)
		{
			/* Fixing offsets to scatter array to different processes */
			send_count = (int*)malloc(group_size*sizeof(int));
			displs = (int*)malloc(group_size*sizeof(int));
			for (i = 0; i < group_size; ++i)
			{
				{
					displs[i] = i*array_size;
					send_count[i] = array_size;
				}
				if (i >= mod && mod > 0)
				{
					send_count[i] = (array_size-1);
					displs[i] = i*(array_size - 1) + mod;
				}
			}
		
			/*Initializing a random array*/
			A_glob = (double*)malloc(length*sizeof(double));
			srand(time(NULL));
			for (i = 0; i < length; ++i)
			{
				A_glob[i] = (rand() % 1000);//(RAND_MAX*0.001));
			}
		}
		begin_time = MPI_Wtime();
	
		/* Scattering large array to  group_size individual smaller arrays */
		A = (double*)malloc(array_size*sizeof(double));
		MPI_Scatterv(A_glob, send_count, displs, MPI_DOUBLE, A, array_size, \
							MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
		/* Sorting individual arrays */
		seqQsort(A, array_size);
	#if print
			printArray(A, array_size, rank);
	#endif
	
		/*Picking first common median and broadcasting it and freeing initial root array */
		if (rank == 0)
		{
			free(send_count);
			free(displs);
			free(A_glob);
			piv = MedOfThree(A, 0, array_size-1);
		}
		MPI_Bcast(&piv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#if debug
			printf("My rank is = %d pivot is = %lf \n", rank, piv);
		#endif	
	
		/*Initializing variables needed in distribution and merging loop*/
		int index, send_size, recv_size, start, end, new_size, j, k, l, cntr;
		double *tmp_array;
		cntr = 0;
	
		/* Looping to send and merge sorted arrays */
		for ( pal = group_size/2; pal > 0; pal = pal >> 1)
		{

			/*Checking what index values in sorted array are larger than pivot*/
			index = 0;
			for (i = 0; i < array_size; ++i)
			{
				if (A[i] < piv)

					index++;
			}
			#if debug
			cntr++;
			printf("iteration = %d, rank/pal mod 2 =  %d, rank = %d, pal  = %d \n", cntr,((rank/pal) % 2), \
						rank, pal);
			#endif
		
			/*If on the lower (left) of two partners send elements larger than
			* pivot to the higher (right) else send elements smaller than pivot
			* to the left */
			if ((rank/pal) % 2 == 0)
			{
				send_size = 0;
				send_size = array_size - index;


					MPI_Send(&A[index], send_size, MPI_DOUBLE, rank + pal, rank, \
								MPI_COMM_WORLD);

				recv_size = 0;
				MPI_Probe(rank + pal, rank + pal, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_DOUBLE, &recv_size);

					free(buff);
					buff = (double*)malloc(recv_size*sizeof(double));
					MPI_Recv(buff, recv_size, MPI_DOUBLE, rank+pal, rank+pal, \
								MPI_COMM_WORLD, &status);


			
				#if debug
					printf("RECEIVED \n");
					printArray(buff, recv_size, rank);
				#endif
			}
			else
			{
				send_size = 0;
				send_size = index;

				recv_size = 0;
				MPI_Probe(rank - pal, rank - pal, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_DOUBLE, &recv_size);


					free(buff);
					buff = (double*)malloc(recv_size*sizeof(double));
					MPI_Recv(buff, recv_size, MPI_DOUBLE, rank-pal, rank-pal, \
								MPI_COMM_WORLD, &status);



					MPI_Send(&A[0], send_size, MPI_DOUBLE, rank - pal, rank, \
								MPI_COMM_WORLD);							

				#if debug
					printf("RECEIVED \n");
					printArray(buff, recv_size, rank);
				#endif
			}
		
		
			/*if on the left keep smaller elements and merge whats kept
			* with recieved. Else keep large and to the same.*/
			if ((rank / pal) % 2 == 0)
			{
				start = 0;
				end = index;
				new_size = end-start + recv_size;
				j = start;
				#if debug
					printf("My rank is = %d start = %d, end = %d, new size = %d, Array size = %d, " \
								"recv_size = %d, index = %d \n", rank, start, \
					end, new_size, array_size, recv_size, index);
				#endif
				k = 0;
				l = 0;
				tmp_array = (double*)malloc(new_size*sizeof(double));
				while (j < end && l < recv_size)
				{
					if (A[j] <= buff[l])
					{
						tmp_array[k] = A[j];
						j++;
					}
					else
					{
						tmp_array[k] = buff[l];
						l++;
					}
					k++;
				}
				if (j == end)
				{
					while (l < recv_size)
					{
							tmp_array[k] = buff[l];
							k++;
							l++;				
					}
				}
				else 
				{
					while (j < end)
					{
							tmp_array[k] = A[j];
							k++;
							j++;
					}
				}
				free(A);
				A = tmp_array;
				tmp_array = NULL;
				array_size = new_size;
		
			}
			else
			{
				start = index;
				end = array_size;
				new_size = end - start + recv_size;
				j = start;
				k = 0;
				l = 0;
				tmp_array = (double*)malloc(new_size*sizeof(double));
				#if debug
					printf("My rank is = %d start = %d, end = %d, new size = %d, Array size = %d, " \
								"recv_size = %d, index = %d \n", rank, start, \
					end, new_size, array_size, recv_size, index);
				#endif	
				while (j < end && l < recv_size)
				{
					if (A[j] <= buff[l])
					{
						tmp_array[k] = A[j];
						j++;
					}
					else
					{
						tmp_array[k] = buff[l];
						l++;
					}
					k++;
				}
				if (j == end)
				{
					while (l < recv_size)
					{
							tmp_array[k] = buff[l];
							k++;
							l++;
					}
				}
				else
				{
					while (j < end)
					{
							tmp_array[k] = A[j];
							k++;
							j++;
					}
				}

				free(A);
				A = tmp_array;
				tmp_array = NULL;
				array_size = new_size;
			}
			/*If evenly divisible by pal send pivot until rank + partner-1
			* is reached. Else recieve pivot.*/
			if ((rank % pal) == 0)
			{
				
					piv = MedOfThree(A, 0, array_size-1);
				for (ii = 1; ii < pal; ii++)
				{
					MPI_Send(&piv, 1, MPI_DOUBLE, rank + ii, pal +1, \
								MPI_COMM_WORLD);
				}

			
			}
			else
			{
			
					MPI_Recv(&piv, 1, MPI_DOUBLE, (int)((rank/pal)*pal), pal +1, \
								MPI_COMM_WORLD, &status);
			

			}
	
			#if print
				printArray(A, array_size, rank);
			#endif

		}
	
		/*Allocate memory for new full array and displacements and size vectors*/
		int *displ2, *array_size2;
		double *array;
		if (rank == 0)
		{
			displ2 = (int*)malloc(group_size*sizeof(int));
			array_size2 = (int*)malloc(group_size*sizeof(int));
			array = (double*)malloc(length*sizeof(double));
		}
	
		/*Gather final indiviual array sizes*/
	
		/*Fixing displacements*/
		MPI_Gather(&array_size, 1, MPI_INT, array_size2, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			displ2[0] = 0;
			for (i = 1; i < group_size; ++i)
				displ2[i] = displ2[i-1] + array_size2[i-1];
		}
	
		/*Gathering data from individual processes into full sorted array.*/
		MPI_Gatherv(A, array_size, MPI_DOUBLE, array, array_size2, displ2, \
						MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
		free(A);
		free(buff);
		end_time = MPI_Wtime();
		if (rank  == 0)
		{
			fastest_time = 999999999;
			/*Checking if sorted*/
			#if verify
				int sorted = verify_sort(array, length);
			#endif
			if (fastest_time > (end_time - begin_time))
			{
				fastest_time = end_time - begin_time;
			}
			free(displ2);
			free(array_size2);
			free(array);
		}
		
	}
	if (rank == 0)
	{
		printf("Array size = %d, Number of processes = %d, Lowest time (of %d runs) to sort = %lf seconds. \n", length, group_size, loops, fastest_time);
	}
	MPI_Finalize();

	
}

double MedOfThree(double *A, int left, int right)
{
	int center, tmp;
	center = (int)(left+right)/2;
	if (A[left] > A[right])
	{
		tmp = left;
		left = right;
		right = tmp;
	}
	if (A[center] < A[left])
	{
		tmp = center;
		center = left;
		left = tmp;
	}
	if (A[right] < A[center])
	{
		tmp = right;
		right = center;
		center = tmp;
	}
	return A[center];
}

void seqQsort(double *A, int length)
{
	if (length < 2) return;
	
	int i,j;
	double tmp, pivot;
	pivot = (A[(int)((rand()/((float)RAND_MAX+1))*length)] \
				+ A[(int)((rand()/((float)RAND_MAX+1))*length)])*0.5;
	for (i = 0, j = length-1; ; i++, j--)
	{
		while (A[i] < pivot) i++;
		while (A[j] > pivot) j--;
		if (i >= j) break;
		
		tmp = A[i];
		A[i] = A[j];
		A[j] = tmp;
	}

	seqQsort(A, i);
	seqQsort(A+i, length -i);
}

void printArray(double *A, int length, int rank)
{
	int i;
	printf("This array is %d long, my rank is: %d, printing: \n", length, rank);
	printf("{ ");
	for (i = 0; i < length; ++i)
	{
		printf("%lf ", A[i]);
	}
	printf("} \n");
	printf("\n");
	printf("\n");
	printf("\n");
}


int verify_sort(double *A, int length)
{
	if (length == 0) return 0;
	int i;
	for (i = 1; i < length; ++i)
	{
		if (A[i-1] > A[i])
		{
			printf("Array unsorted from index: %d \n", i);
			return -1;
		}
	}
	printf("**** ARRAY SORTED ***** \n");
	return 1;
}
