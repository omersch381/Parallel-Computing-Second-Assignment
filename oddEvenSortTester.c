#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#define MASTER 0
#define TAG 0
#define MIN_TO_MAX 0

// Hamami aya po

void oddEvenSort(float *myValue, int n, int pLeft, int pRight, int location, int direction, int my_rank);
void compare(float *myValue, float *otherValue, int direction);

int main(int argc, char *argv[])
{

    int my_rank;       /* rank of process */
    int p;             /* number of processes */
    int source;        /* rank of sender */
    int dest;          /* rank of receiver */
    int tag = 0;       /* tag for messages */
    char message[100]; /* storage for message */
    MPI_Status status; /* return status for receive */

    /* start up MPI */

    MPI_Init(&argc, &argv);

    /* find out process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    float array[] = {4.0, 3.0, 2.0, 1.0};
    float null = 999;
    float myValue;
    switch (my_rank)
    {
    case 0:
        myValue = array[0];
        oddEvenSort(&myValue, p, null, my_rank + 1, 0, MIN_TO_MAX, my_rank);
        printf("I am process num #%d and myValue is %f\n", my_rank, myValue);
        break;
    case 1:
        myValue = array[1];
        oddEvenSort(&myValue, p, my_rank - 1, my_rank + 1, 1, MIN_TO_MAX, my_rank);
        printf("I am process num #%d and myValue is %f\n", my_rank, myValue);
        break;
    case 2:
        myValue = array[2];
        oddEvenSort(&myValue, p, my_rank - 1, my_rank + 1, 2, MIN_TO_MAX, my_rank);
        printf("I am process num #%d and myValue is %f\n", my_rank, myValue);
        break;
    default:
        myValue = array[3];
        oddEvenSort(&myValue, p, my_rank - 1, null, 3, MIN_TO_MAX, my_rank);
        printf("I am process num #%d and myValue is %f\n", my_rank, myValue);
        break;
    }

    /* shut down MPI */
    MPI_Finalize();

    return 0;
}

void oddEvenSort(float *myValue, int n, int pLeft, int pRight, int location, int direction, int my_rank)
{
    MPI_Status status; /* return status for receive */
    float *otherValue = (float *)malloc(sizeof(float));
    ;
    int i;
    for (i = 0; i < n; i++)
    {
        if (i % 2 == 0) // An Even Iteration
        {
            if (location % 2 == 0) // If the location is even
            {
                // As for right now - the even location on the even iteration sorts and sends back.
                MPI_Recv(otherValue, 1, MPI_FLOAT, pRight, TAG, MPI_COMM_WORLD, &status);
                compare(myValue, otherValue, direction);
                MPI_Send(otherValue, 1, MPI_FLOAT, pRight, TAG, MPI_COMM_WORLD);
            }
            else // if the location is odd
            {
                MPI_Send(myValue, 1, MPI_FLOAT, pLeft, TAG, MPI_COMM_WORLD);
                MPI_Recv(myValue, 1, MPI_FLOAT, pLeft, TAG, MPI_COMM_WORLD, &status);
            }
        }
        else // if odd iteration
        {
            if (!(location == 0 || location == n - 1)) // edge cases do not participate in odd iterations
            {
                if (location % 2 == 0) // if location is even
                {
                    MPI_Send(myValue, 1, MPI_FLOAT, pLeft, TAG, MPI_COMM_WORLD);
                    MPI_Recv(myValue, 1, MPI_FLOAT, pLeft, TAG, MPI_COMM_WORLD, &status);
                }
                else // if location is odd
                {
                    // As for right now - the odd location on the odd iteration sorts and sends.
                    MPI_Recv(otherValue, 1, MPI_FLOAT, pRight, TAG, MPI_COMM_WORLD, &status);
                    compare(myValue, otherValue, direction);
                    MPI_Send(otherValue, 1, MPI_FLOAT, pRight, TAG, MPI_COMM_WORLD);
                }
            }
        }
    }
}

void compare(float *myValue, float *otherValue, int direction)
{
    float temp;
    if (direction == MIN_TO_MAX)
    {
        if (*myValue > *otherValue) // which is not what we want
        {
            temp = *otherValue;
            *otherValue = *myValue;
            *myValue = temp;
        }
        //else (if myValue is already smaller than otherValue) do nothing
    }
    else // if direction == MAX_TO_MIN
    {
        if (*myValue < *otherValue) // which is not what we want
        {
            temp = *otherValue;
            *otherValue = *myValue;
            *myValue = temp;
        }
    }
}
