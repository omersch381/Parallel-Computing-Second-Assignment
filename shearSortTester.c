#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #define ROWS 4
// #define COLUMNS 3
#define TAG 0
#define MIN_TO_MAX 0
#define NUM_OF_DIMENSIONS 2
#define MAX_NUM_OF_DIMENSIONS 2
#define HORIZONTAL_DIRECTION 1
#define VERTICAL_DIRECTION 0
#define DIRECT_NEIGHBORS 1

int floorSqrt(int x);
void oddEvenSort(float *myValue, int n, int pLeft, int pRight, int location, int direction);
void compare(float *myValue, float *otherValue, int direction);
void sanityCheck(int numOfProcesses, int n);
MPI_Comm getNewCommunicator(int n);
void getCoordinatesFromGivenRank(MPI_Comm newCommunicator, int rank, int *rowIndex, int *colIndex);
int getRankFromGivenCoordinates(MPI_Comm newCommunicator, int *rowIndex, int *colIndex);
void getNeighborsRanksFromGivenRank(MPI_Comm newCommunicator, int rank, int *pRightRank, int *pLeftRank);

int main(int argc, char *argv[])
{
    int rank, numOfProcesses;
    MPI_Comm newCommunicator;
    int numOfDimensionsArray[2], isPeriodicArray[2], reorder;
    int coordinatesArray[2], newId;
    int rowIndex, colIndex;
    int pLeftRank, pRightRank;
    int n;
    int volume, height;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);

    n = floorSqrt(numOfProcesses);

    sanityCheck(numOfProcesses, n);

    height = rank;
    volume = 2 * rank;
    newCommunicator = getNewCommunicator(n);

    getCoordinatesFromGivenRank(newCommunicator, rank, &rowIndex, &colIndex);
    printf("The coordinates are: %d %d\n", rowIndex, colIndex);

    rank = getRankFromGivenCoordinates(newCommunicator, &rowIndex, &colIndex);

    getNeighborsRanksFromGivenRank(newCommunicator, rank, &pRightRank, &pLeftRank);

    MPI_Finalize();
    return 0;
}

void getNeighborsRanksFromGivenRank(MPI_Comm newCommunicator, int rank, int *pRightRank, int *pLeftRank)
{
    MPI_Cart_shift(newCommunicator, HORIZONTAL_DIRECTION, DIRECT_NEIGHBORS, pLeftRank, pRightRank);
    // printf("Rank = %d pLeftRank = %d  pRightRank %d\n", rank, *pLeftRank, *pRightRank);
    fflush(stdout);
}

int getRankFromGivenCoordinates(MPI_Comm newCommunicator, int *rowIndex, int *colIndex)
{
    int newRank, coordinatesArray[2];
    coordinatesArray[0] = *rowIndex;
    coordinatesArray[1] = *colIndex;
    MPI_Cart_rank(newCommunicator, coordinatesArray, &newRank);
    // printf("The processor at position (%d, %d) has rank %d\n", coordinatesArray[0], coordinatesArray[1], newId);
    fflush(stdout);
    return newRank;
}

void getCoordinatesFromGivenRank(MPI_Comm newCommunicator, int rank, int *rowIndex, int *colIndex)
{
    int coordinatesArray[2];
    MPI_Cart_coords(newCommunicator, rank, MAX_NUM_OF_DIMENSIONS, coordinatesArray);
    // printf("Rank %d coordinates are %d %d\n", rank, coordinatesArray[0], coordinatesArray[1]);
    *(rowIndex) = coordinatesArray[0];
    *(colIndex) = coordinatesArray[1];
    fflush(stdout);
}

MPI_Comm getNewCommunicator(int n)
{
    MPI_Comm newCommunicator;
    int numOfDimensionsArray[2], isPeriodicArray[2], reorder;
    numOfDimensionsArray[0] = n;
    numOfDimensionsArray[1] = n;
    isPeriodicArray[0] = 0;
    isPeriodicArray[1] = 0;
    reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, NUM_OF_DIMENSIONS, numOfDimensionsArray, isPeriodicArray, reorder, &newCommunicator);
    return newCommunicator;
}

void sanityCheck(int numOfProcesses, int n)
{
    if (numOfProcesses != 16)
    {
        printf("Please run with 16 processes.\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int isNPowerOfTwo = (n / 2) % 2;
    if (!isNPowerOfTwo == 0)
    {
        printf("Please run with (n = power of two) processes.\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
void oddEvenSort(float *myValue, int n, int pLeftRank, int pRightRank, int location, int direction)
{
    MPI_Status status;
    float *otherValue = (float *)malloc(sizeof(float));
    int i;
    for (i = 0; i < n; i++)
    {
        if (i % 2 == 0) // An Even Iteration
        {
            if (location % 2 == 0) // If the location is even
            {
                // As for right now - the even location on the even iteration sorts and sends back.
                MPI_Recv(otherValue, 1, MPI_FLOAT, pRightRank, TAG, MPI_COMM_WORLD, &status);
                compare(myValue, otherValue, direction);
                MPI_Send(otherValue, 1, MPI_FLOAT, pRightRank, TAG, MPI_COMM_WORLD);
            }
            else // if the location is odd
            {
                MPI_Send(myValue, 1, MPI_FLOAT, pLeftRank, TAG, MPI_COMM_WORLD);
                MPI_Recv(myValue, 1, MPI_FLOAT, pLeftRank, TAG, MPI_COMM_WORLD, &status);
            }
        }
        else // if odd iteration
        {
            if (!(location == 0 || location == n - 1)) // edge cases do not participate in odd iterations
            {
                if (location % 2 == 0) // if location is even
                {
                    MPI_Send(myValue, 1, MPI_FLOAT, pLeftRank, TAG, MPI_COMM_WORLD);
                    MPI_Recv(myValue, 1, MPI_FLOAT, pLeftRank, TAG, MPI_COMM_WORLD, &status);
                }
                else // if location is odd
                {
                    // As for right now - the odd location on the odd iteration sorts and sends.
                    MPI_Recv(otherValue, 1, MPI_FLOAT, pRightRank, TAG, MPI_COMM_WORLD, &status);
                    compare(myValue, otherValue, direction);
                    MPI_Send(otherValue, 1, MPI_FLOAT, pRightRank, TAG, MPI_COMM_WORLD);
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
int floorSqrt(int x)
{
    // Base cases
    if (x == 0 || x == 1)
        return x;

    // Staring from 1, try all numbers until
    // i*i is greater than or equal to x.
    int i = 1, result = 1;
    while (result <= x)
    {
        i++;
        result = i * i;
    }
    return i - 1;
}