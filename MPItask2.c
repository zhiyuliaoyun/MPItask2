#pragma warning(disable:4996)
#include <stdio.h>
#include <mpi.h>
#include <malloc.h>
#include <math.h>
float f(float x)
{
	return (x);
}
float Trap(float a, float b, int j)
{
	float h = (b - a) / j, T = 0;
	for (int i = 1; i < j; i++)  T += f(a + i * h);  
	return h * (f(a) + 2 * T + f(b)) / 2;
}
main(int argc, char** argv) {
	int numProcs, my_rank;
	int master = 0;
	int i=0;
	int j = 10;
	int k=0;
	int count=0;
	float total = 0;
	float a=10, b=100, n=1024;
	float part_trap;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	if (my_rank = master) {
		for (int i = 1; i < numProcs; i++) {
			MPI_Send(&k, 1, MPI_FLOAT,i, MPI_ANY_TAG, MPI_COMM_WORLD);
			k = k + 1;
		}
		for (int i = 1; i < numProcs; i++){
			MPI_Recv(&part_trap, 1, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			total += part_trap;
			if (k < 10) {
				k = k + 1;
				MPI_Send(&k, 1, MPI_FLOAT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD);
			}
		}
		printf(" the integral of f(x)=x from 0 to 10 is %1.3f\t", total);
	}
	else {
		while (1) {
			MPI_Recv(&k, 1, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if (k < 10) {
				float h = (b - a) / n;
				float part_a = a + k * h;
				float part_b = a + (k + 1) * h;
				float part_trap = Trap(part_a, part_b, j);
				MPI_Send(&part_trap, 1, MPI_FLOAT, 0, i, MPI_COMM_WORLD);
			}
			else {
				break;
			}
		}
	}
	MPI_Finalize();
}
