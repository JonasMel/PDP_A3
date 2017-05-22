#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

void qs(double *A, int length);
double pivot(double *A, length);

int main(int *argc, char **argv)
{
	/*Declaring variables.*/
	int p, i, j, ndim, rank, glob_rank, reorder, length, iter;
	double *A_glob, *A; *P, *recvbuffer, pivot;
	
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
	MPI_Cart_create(MPI_COMM_WORLD, ndim, p, cyclic, reorder, &crt);
	MPI_Comm_rank(crt, &rank);
	int size;
	size = length/p
	if (rank == 0)
	{
		A_glob = (double*)malloc(length*sizeof(double));
		srand(time(NULL));
		for (i = 0; i < length; ++i)
		{
			A_glob[i] = (rand()/(RAND_MAX))*1000;
		}
	}
	A = (double*)malloc(size*sizeof(double));
	MPI_Scatter(A_glob, size, MPI_DOUBLE, A, size, MPI_DOUBLE 0, crt);
	seqQsort(A, size);
	if (round(rand()/(RAND_MAX+1)*rank) == rank)
	{
		pivot = getPivot(A, size);
		MPI_Bcast(&pivot, 1, MPI_DOUBLE, rank, crt);
	}
	int start, end, sendsize, recvsize;
	start = 0;
	end = size;
	for (iter = p/2; iter > 0; iter >>= 1)
	{
		if (rank < iter)
		{
			for (i = 0; i < size; ++i)
			{
				if (A[i] >= pivot)
				{
					start = i;
					break;
				}
			}
			sendsize = size - start;
		}
		else
		{
			for (i = size; i >= 0; --i)
			{
				if (A[i] < pivot)
				{
					end = i;
					break;
				}
			}
			sendsize = end;
		}
		if ((rank/iter) % 2 == 0 && sendsize)
		{
			MPI_Isend(A+start, sendsize, MPI_DOUBLE, rank+iter, rank, crt, &send_req);
			MPI_Probe(rank+iter, rank+iter, crt, &status);
			MPI_Get_count(&status, MPI_DOUBLE, &recvsize);
			if (recvsize)
			{
				recvbuffer = (double*)malloc(recvsize*sizeof(double));
				MPI_Irecv(recvbuffer, recvsize, MPI_DOUBLE, rank+iter, \
					rank+iter, crt, &recv_req);
				MPI_Wait(&recv_req, &status);
			}
		}
		else if (sendsize)
		{
			MPI_Isend(A, sendsize, MPI_DOUBLE, rank-iter, rank, crt, &send_req);
			MPI_Probe(rank-iter, rank-iter, crt, &status);
			MPI_Get_count(&status, MPI_DOUBLE, &recvsize);
			if (recvsize)
			{
				recvbuffer = (double*)malloc(recvsize*sizeof(double));
				MPI_Irecv(recvbuffer, recvsize, MPI_DOUBLE, rank-iter, \
					rank-iter, crt, &recv_req);
				MPI_Wait(&recv_req, &status);
			}
		}
		MPI_Barrier(crt);
		if ((rank/iter) % 2 == 0)
		{
			merge(A, recvbuffer, recvsize, sendsize, &size);
		}
		else
		{
			merge(A+sendsize, recvbuffer, recvsize, sendsize, &size);
		}
		if (iter > 1)
		{
			double piv_cand_send, piv_cand_recv;
			piv_cand_send = MedOfThree(A, (int)(0.25*(size)), (int)(0.75*(size)));
			if ((rank/iter) % 2 == 0)
			{
				MPI_Isend(&piv_cand_send, 1, MPI_DOUBLE, rank+iter*0.5, rank+iter*0.5, \
					crt, &send_req);
					
				MPI_Irecv(&piv_cand_recv, 1, MPI_DOUBLE, rank+iter*0.5, rank+iter*0.5, \
					crt, &recv_req);
				MPI_Wait(&recv_req, &status);
			}
			else
			{
				MPI_Isend(&piv_cand_send, 1, MPI_DOUBLE, rank-iter*0.5, rank-iter*0.5, \
					crt, &send_req);
					
				MPI_Irecv(&piv_cand_recv, 1, MPI_DOUBLE, rank-iter*0.5, rank-iter*0.5, \
					crt, &recv_req);
				MPI_Wait(&recv_req, &status);
			}
			pivot = (piv_cand_send + piv_cand_recv)*0.5;
		}
	}
	
	
	
}

double* merge(double *Arr, double *Brr, int recieve_size, int sent_size, int &blk_size)
{
	double *tmp;
	int i, j, k;
	i = 0;
	tmp = (double*)malloc((blk_size - sent_size + recievesize)*sizeof(double))
	if (recieve_size)
	{
		for (k = 0; k < (blk_size - sent_size); ++k)
		{
			if (Arr[j] <= Brr[i])
			{
				tmp[k] = Arr[j]
				j++
			}
			else if (i < recieve_size)
			{
				tmp[k] = Brr[i];
				i++;
			}
		}
	}
	else
	{
		for (i = 0; i < (blk_size - sent_size); ++i)
			tmp[i] = Arr[i]
	}
 	for (j = i; j < recieve_size; ++j)
	{
		tmp[j] = Brr[j];
	}
	free(Arr);
	free(Brr);
	blk_size = (blk_size - sent_size + recievesize);
	return tmp;
}

double getPivots(double *P, length)
{
	return A[floor(length/2)];
}

double MedOfThree(double *A, double left, double right)
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
	pivot = (A[round((rand()/(RAND_MAX+1))*length)] \
				+ A[round((rand()/(RAND_MAX+1))*length)])*0.5;
	for (i = 0, j = length-1; ; i++, j--)
	{
		while (A[i] < pivot) i++;
		while (A[j] > pivot) j--;
		if (i => j) break;
		
		tmp = A[i];
		A[i] = A[j];
		A[j] = tmp;
	}

	seqQsort(A, i);
	seqQsort(A+i, length -i);
}

