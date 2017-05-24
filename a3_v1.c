#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

#define print 1
double* merge(double *Arr, double *Brr, int recieve_size, \
					int index, int start,int *blk_size);
double getPivots(double *P, int length);
double MedOfThree(double *A, int left, int right);
void seqQsort(double *A, int length);
void printArray(double *A, int length, int rank);

int main(int argc, char **argv)
{
	/*Declaring variables.*/
	int p, i, j, ndim, rank, glob_rank, reorder, length, iter, group_size;
	double *A_glob, *A, *P, *recvbuffer, pivot;
	
	p = atoi(argv[1]);
	length = atoi(argv[2]);
	MPI_Request send_req[2], recv_req[2];
	MPI_Status status;
	MPI_Init(&argc, &argv);
	ndim = 1;
	reorder = 0;
	int cyclic[] = {0, 0};
	MPI_Comm crt;
	
	
	MPI_Comm_rank(MPI_COMM_WORLD, &glob_rank);
	MPI_Cart_create(MPI_COMM_WORLD, ndim, &p, cyclic, reorder, &crt);
	MPI_Comm_size(crt, &group_size);

	MPI_Comm_rank(crt, &rank);
	int size;
	size = length/p;
	if (rank == 0)
	{
		printf("groupsize = %d \n", group_size);
		A_glob = (double*)malloc(length*sizeof(double));
		srand(time(NULL));
		for (i = 0; i < length; ++i)
		{
			A_glob[i] = (rand()/(RAND_MAX*0.1));
		}
#if print
		printArray(A_glob, length, rank);
#endif
	}
/*#if print	
	printf("before Scatter rank = %d \n", rank);
#endif */
	A = (double*)malloc(size*sizeof(double));
	MPI_Scatter(A_glob, size, MPI_DOUBLE, A, size, MPI_DOUBLE, 0, crt);
/*#if print
	printf("after Scatter, rank = %d \n", rank);
#endif	*/

	seqQsort(A, size);
#if print
		printArray(A, size, rank);
#endif
	int brdcster;
	brdcster = (int)(rand()/((double)RAND_MAX)*rank);
	if (rank == 0)
	{
/*#if print
		printf("Going to broadcast Pivot, my rank is %d \n", rank);
#endif*/
		pivot = MedOfThree(A, (int)(0.25*size), (int)(0.75*size));

	}

	MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, crt);
	//MPI_Wait(&send_req[1], &status);
#if print
	int cntr = 0;
	printf("Broadcasted pivot recieved at rank = %d, pivot = %lf \n", rank, pivot);
#endif
	int index, start, end, sendsize, recvsize;
	index = 0;
	for (iter = group_size/2; iter > 0; iter >>= 1)
	{
#if print
		cntr++;
		
		printf("In loop number = %d rank = %d \n", cntr, rank);
#endif	

			for (i = 0; i < size; ++i)
			{
				if (A[i] < pivot)
				{
					index++;
				}
			}

		
/*#if print
		printf("Before sending \n");
#endif*/	
		if ((rank/iter) % 2 == 0 )
		{
#if print
			printf("I am rank  = %d Sending to %d \n", rank, rank+iter);
#endif	
			sendsize = size - index;
			if (sendsize > 0)
			{
				MPI_Isend(&A[index], sendsize, MPI_DOUBLE, rank+iter, rank, crt, &send_req[0]);
			}
			MPI_Probe(rank+iter, rank+iter, crt, &status);
			MPI_Get_count(&status, MPI_DOUBLE, &recvsize);
			if (recvsize > 0)
			{
				recvbuffer = (double*)malloc(recvsize*sizeof(double));
				MPI_Recv(recvbuffer, recvsize, MPI_DOUBLE, rank+iter, \
					rank+iter, crt, &status);
				//MPI_Wait(&recv_req[0], &status);
			}
#if print
			printf("I am rank  = %d Receiving from %d \n", rank, rank+iter);
#endif	
#if print
			printArray(recvbuffer, recvsize, rank);
#endif

		}
		else
		{
#if print
			printf("I am rank  = %d sending to %d \n", rank, rank-iter);
#endif	
			sendsize = index;
			if (sendsize > 0)
			{
				MPI_Isend(&A[0], sendsize, MPI_DOUBLE, rank-iter, rank, crt, &send_req[1]);
			}
			MPI_Probe(rank-iter, rank-iter, crt, &status);
			MPI_Get_count(&status, MPI_DOUBLE, &recvsize);

			if (recvsize > 0)
			{
				recvbuffer = (double*)malloc(recvsize*sizeof(double));
				MPI_Recv(recvbuffer, recvsize, MPI_DOUBLE, rank-iter, \
					rank-iter, crt, &status);
				//MPI_Wait(&recv_req[1], &status);
			}
#if print
			printf("I am rank  = %d receiving from %d \n", rank, rank-iter);
#endif
#if print
			printArray(recvbuffer, recvsize, rank);
#endif
		}
//		MPI_Barrier(crt);
		
		if ((rank/iter) % 2 == 0)
		{
/*#if print
			printf("Merging in even %d \n", rank);
#endif*/
			start = 0;
			end = index;
		//	A = merge(&A[0], &recvbuffer[0], recvsize, end, start, &size);
	double *tmp;
	int i, j, k, new_size, remain_size, blk_size_tmp;
	//blk_size_tmp = *blk_size;
	i = 0;
	j = start;
	k = 0; 
	remain_size = end - start;
	new_size = (remain_size + recvsize);
	tmp = (double*)malloc((new_size)*sizeof(double));

#if print
				printf("I am rank  = %d In merging loop part 1\n", rank);
#endif
		while (j < (end) && i <  recvsize)
		{
			if (A[j] <= recvbuffer[i])
			{
				tmp[k] = A[j];

#if print
				printf("while... j = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",j, k, j, A[j] \
					, i, recvbuffer[i], k, tmp[k]);
#endif
				j++;
			}
			else
			{
				tmp[k] = recvbuffer[i];

#if print
				printf("while... i = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",i, k, j, A[j] \
					, i, recvbuffer[i], k, tmp[k]);
#endif
				i++;
			}
			 

			k++;
		}
		if (i < recvsize)
		{
			for (int m = i; m < recvsize; m++)
			{
				tmp[k] = recvbuffer[m];

#if print
				printf("while... m = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",m, k, j, A[j] \
					, i, recvbuffer[m], k, tmp[k]);
#endif
				k++;
			}
		}
		else
		{
			for (int n = j; n < (end); n++)
			{
				tmp[k] = A[n];

#if print
				printf("while... n = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",n, k, j, A[n] \
					, i, recvbuffer[i], k, tmp[k]);
#endif
				k++;
			}
		}

	free(A);
	free(recvbuffer);
	
	size = new_size;
	A = tmp;

		}
		else
		{
/*#if print
			printf("Merging in odd %d \n", rank);
#endif*/
			start = index;
			end = size;
		//	A = merge(&A[0], &recvbuffer[0], recvsize, end, start, &size);
	double *tmp;
	int i, j, k, new_size, remain_size, blk_size_tmp;
	//blk_size_tmp = *blk_size;
	i = 0;
	j = start;
	k = 0; 
	remain_size = end - start;
	new_size = (remain_size + recvsize);
	tmp = (double*)malloc((new_size)*sizeof(double));

#if print
				printf("I am rank  = %d  In merging loop part 1\n", rank);
#endif
		while (j < (end) && i <  recvsize)
		{
			if (A[j] <= recvbuffer[i])
			{
				tmp[k] = A[j];

#if print
				printf("while... j = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",j, k, j, A[j] \
					, i, recvbuffer[i], k, tmp[k]);
#endif
				j++;
			}
			else
			{
				tmp[k] = recvbuffer[i];

#if print
				printf("while... i = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",i, k, j, A[j] \
					, i, recvbuffer[i], k, tmp[k]);
#endif
				i++;
			}
			 

			k++;
		}
		if (i < recvsize)
		{
			for (int m = i; m < recvsize; m++)
			{
				tmp[k] = recvbuffer[m];

#if print
				printf("while... m = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",m, k, j, A[j] \
					, i, recvbuffer[m], k, tmp[k]);
#endif
				k++;
			}
		}
		else
		{
			for (int n = j; n < (end); n++)
			{
				tmp[k] = A[n];

#if print
				printf("while... n = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",n, k, j, A[n] \
					, i, recvbuffer[i], k, tmp[k]);
#endif
				k++;
			}
		}

	free(A);
	free(recvbuffer);
	
	size = new_size;
	A = tmp;

		}
/*#if print
				printf("After merge, rank %d \n", rank);
#endif	*/	
#if print
		printArray(A, size, rank);
#endif	
		
		if (iter > 1)
		{
			double piv_cand_send, piv_cand_recv;
			piv_cand_send = MedOfThree(A, (int)(0.25*(size)), (int)(0.75*(size)));
			if ((rank/iter) % 2 == 0)
			{
#if print
				printf("sendig pivot from even %d \n", rank);
#endif
				MPI_Isend(&piv_cand_send, 1, MPI_DOUBLE, rank+(iter>>1), rank, \
					crt, &send_req[0]);
					
				MPI_Recv(&piv_cand_recv, 1, MPI_DOUBLE, rank+(iter>>1), rank+(iter>>1), \
					crt, &status);
				//MPI_Wait(&recv_req[0], &status);
			}
			else
			{
#if print
				printf("sendig pivot from odd %d \n", rank);
#endif
				MPI_Isend(&piv_cand_send, 1, MPI_DOUBLE, rank-(iter>>1), rank, \
					crt, &send_req[1]);
					
				MPI_Recv(&piv_cand_recv, 1, MPI_DOUBLE, rank-(iter>>1), rank-(iter>>1), \
					crt, &status);
				//MPI_Wait(&recv_req[1], &status);
			}
			pivot = (piv_cand_send + piv_cand_recv)*0.5;
		}
	}
	

	MPI_Finalize();
}

double* merge(double *Arr, double *Brr, int recieve_size, \
					int end, int start, int *blk_size)
{
	double *tmp;
	int i, j, k, size, remain_size, blk_size_tmp;
	blk_size_tmp = *blk_size;
	i = 0;
	j = start;
	k = 0; 
	remain_size = end - start;
	size = (remain_size + recieve_size);
	tmp = (double*)malloc((size)*sizeof(double));

#if print
				printf("In merging loop part 1\n");
#endif
		while (j < (end) && i <  recieve_size)
		{
			if (Arr[j] <= Brr[i])
			{
				tmp[k] = Arr[j];

#if print
				printf("while... j = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",j, k, j, Arr[j] \
					, i, Brr[i], k, tmp[k]);
#endif
				j++;
			}
			else
			{
				tmp[k] = Brr[i];

#if print
				printf("while... i = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",i, k, j, Arr[j] \
					, i, Brr[i], k, tmp[k]);
#endif
				i++;
			}
			 

			k++;
		}
		if (i < recieve_size)
		{
			for (int m = i; m < recieve_size; m++)
			{
				tmp[k] = Brr[m];

#if print
				printf("while... m = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",m, k, j, Arr[j] \
					, i, Brr[m], k, tmp[k]);
#endif
				k++;
			}
		}
		else
		{
			for (int n = j; n < (end); n++)
			{
				tmp[k] = Arr[n];

#if print
				printf("while... n = %d, k = %d, A[%d] = %lf, B[%d] = %lf, tmp[%d] = %lf \n",n, k, j, Arr[n] \
					, i, Brr[i], k, tmp[k]);
#endif
				k++;
			}
		}

	free(Arr);
	free(Brr);
	
	blk_size = &size;
	return tmp;
}

double getPivots(double *P, int length)
{
	return P[(int)(length/2)];
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
		printf("%lf ", A[i]);
	printf("} \n");
}