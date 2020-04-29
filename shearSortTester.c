#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define TAG 0
#define MIN_TO_MAX 0
#define MAX_TO_MIN 1
#define NUM_OF_DIMENSIONS 2
#define MAX_NUM_OF_DIMENSIONS 2
#define HORIZONTAL_DIRECTION 1
#define VERTICAL_DIRECTION 0
#define DIRECT_NEIGHBORS 1
#define MASTER 0

int floorSqrt(int x);
void oddEvenSort(float *myValue, int n, int pLeftRank, int pRightRank, int location, int direction);
void compare(float *myValue, float *otherValue, int direction);
void sanityCheck(int numOfProcesses, int n);
MPI_Comm getNewCommunicator(int n);
void getCoordinatesFromGivenRank(MPI_Comm newCommunicator, int rank, int *rowIndex, int *colIndex);
int getRankFromGivenCoordinates(MPI_Comm newCommunicator, int *rowIndex, int *colIndex);
void getNeighborsRanksFromGivenRank(MPI_Comm newCommunicator, int rank, int *pRightRank, int *pLeftRank);
void shearSort(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, float myStruct[], int n);
void sortRows(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, float myStruct[], int n);
void sortRow(MPI_Comm newCommunicator, int rank, float myStruct[], int n, int row, int direction);
void testPrint(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, float myStruct[], int n);
void swapFloats(float *myValue, float *otherValue);

int main(int argc, char *argv[])
{
    int rank, numOfProcesses;
    MPI_Comm newCommunicator;
    int numOfDimensionsArray[2], isPeriodicArray[2], reorder;
    int coordinatesArray[2], newId;
    int rowIndex, colIndex;
    int pLeftRank, pRightRank;
    int n;
    float volume, height;
    float myStruct[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);

    n = floorSqrt(numOfProcesses);

    sanityCheck(numOfProcesses, n);

    height = rank;     //test
    volume = 2 * rank; //test
    myStruct[0] = volume;
    myStruct[1] = height;
    // until here - we received the structs
    newCommunicator = getNewCommunicator(n);

    getCoordinatesFromGivenRank(newCommunicator, rank, &rowIndex, &colIndex);

    rank = getRankFromGivenCoordinates(newCommunicator, &rowIndex, &colIndex);

    shearSort(newCommunicator, rank, rowIndex, colIndex, myStruct, n);

    testPrint(newCommunicator, rank, rowIndex, colIndex, myStruct, n); //DELETE_ME

    MPI_Finalize();
    return 0;
}

void testPrint(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, float myStruct[], int n)
{
    int a = 0, b = 1, c = 2, d = 3;
    int *ap = &a, *bp = &b, *cp = &c, *dp = &d;
    int firstRank = getRankFromGivenCoordinates(newCommunicator, ap, ap);
    int secondRank = getRankFromGivenCoordinates(newCommunicator, ap, bp);
    int thirdRank = getRankFromGivenCoordinates(newCommunicator, ap, cp);
    int fourthRank = getRankFromGivenCoordinates(newCommunicator, ap, dp);
    int fifthRank = getRankFromGivenCoordinates(newCommunicator, bp, ap);
    int sixthRank = getRankFromGivenCoordinates(newCommunicator, bp, bp);
    int seventhRank = getRankFromGivenCoordinates(newCommunicator, bp, cp);
    int eigthRank = getRankFromGivenCoordinates(newCommunicator, bp, dp);

    if (rank == firstRank)
        printf("my rank is %d and I am at the coordinates (0,0), my value is %f\n", rank, myStruct[0]);
    if (rank == secondRank)
        printf("my rank is %d and I am at the coordinates (0,1), my value is %f\n", rank, myStruct[0]);
    if (rank == thirdRank)
        printf("my rank is %d and I am at the coordinates (0,2), my value is %f\n", rank, myStruct[0]);
    if (rank == fourthRank)
        printf("my rank is %d and I am at the coordinates (0,3), my value is %f\n", rank, myStruct[0]);
    if (rank == fifthRank)
        printf("my rank is %d and I am at the coordinates (1,0), my value is %f\n", rank, myStruct[0]);
    if (rank == sixthRank)
        printf("my rank is %d and I am at the coordinates (1,1), my value is %f\n", rank, myStruct[0]);
    if (rank == seventhRank)
        printf("my rank is %d and I am at the coordinates (1,2), my value is %f\n", rank, myStruct[0]);
    if (rank == eigthRank)
        printf("my rank is %d and I am at the coordinates (1,3), my value is %f\n", rank, myStruct[0]);
}

void shearSort(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, float myStruct[], int n)
{
    // manual insertion - so we can make sure it works
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    int a = 0, b = 1, c = 2, d = 3;
    int *ap = &a, *bp = &b, *cp = &c, *dp = &d;
    int firstRank = getRankFromGivenCoordinates(newCommunicator, ap, ap);
    int secondRank = getRankFromGivenCoordinates(newCommunicator, ap, bp);
    int thirdRank = getRankFromGivenCoordinates(newCommunicator, ap, cp);
    int fourthRank = getRankFromGivenCoordinates(newCommunicator, ap, dp);
    int fifthRank = getRankFromGivenCoordinates(newCommunicator, bp, ap);
    int sixthRank = getRankFromGivenCoordinates(newCommunicator, bp, bp);
    int seventhRank = getRankFromGivenCoordinates(newCommunicator, bp, cp);
    int eigthRank = getRankFromGivenCoordinates(newCommunicator, bp, dp);
    if (rank == firstRank)
        myStruct[0] = 8.0;
    if (rank == secondRank)
        myStruct[0] = 6.0;
    if (rank == thirdRank)
        myStruct[0] = 7.0;
    if (rank == fourthRank)
        myStruct[0] = 10.0;
    if (rank == fifthRank)
        myStruct[0] = 1.0;
    if (rank == sixthRank)
        myStruct[0] = 4.0;
    if (rank == seventhRank)
        myStruct[0] = 3.0;
    if (rank == eigthRank)
        myStruct[0] = 2.0;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    int nothing;
    int arrayOfNothing[16]; // we gather only to synchronize them - we don't need anything else
    MPI_Gather(&nothing, sizeof(int), MPI_BYTE, arrayOfNothing, sizeof(int), MPI_BYTE, MASTER, newCommunicator);
    sortRows(newCommunicator, rank, rowIndex, colIndex, myStruct, n);
}

void sortRows(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, float myStruct[], int n)
{
    int nothing, row;
    int arrayOfNothing[16]; // we gather only to synchronize them - we don't need anything else
    int direction;
    MPI_Gather(&nothing, sizeof(int), MPI_BYTE, arrayOfNothing, sizeof(int), MPI_BYTE, MASTER, newCommunicator);

    for (row = 0; row < n; row++)
    {
        direction = row % 2 == 0 ? MIN_TO_MAX : MAX_TO_MIN;
        sortRow(newCommunicator, rank, myStruct, n, row, direction);
        MPI_Gather(&nothing, sizeof(int), MPI_BYTE, arrayOfNothing, sizeof(int), MPI_BYTE, MASTER, newCommunicator);
    }
}

void sortRow(MPI_Comm newCommunicator, int rank, float myStruct[], int n, int row, int direction)
{
    int pRightRank, pLeftRank;
    int rowArray[n];
    int col;
    for (col = 0; col < n; col++)
    {
        int currentRank = getRankFromGivenCoordinates(newCommunicator, &row, &col); //(0,0)
        if (rank == currentRank)
        {
            getNeighborsRanksFromGivenRank(newCommunicator, currentRank, &pRightRank, &pLeftRank);
            int location = row * n + col; // converting to the location paramerter that oddEvenSort accept
            oddEvenSort(&myStruct[0], n, pLeftRank, pRightRank, location, direction);
        }
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
            if (!(location == 0 || location % n == n - 1)) // edge cases do not participate in odd iterations
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
    if (direction == MIN_TO_MAX && *myValue > *otherValue) // which is not what we want
        swapFloats(myValue, otherValue);
    else if (direction == MAX_TO_MIN && *myValue < *otherValue) // which is not what we want
        swapFloats(myValue, otherValue);
}

void swapFloats(float *myValue, float *otherValue)
{
    float temp = *otherValue;
    *otherValue = *myValue;
    *myValue = temp;
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
void getNeighborsRanksFromGivenRank(MPI_Comm newCommunicator, int rank, int *pRightRank, int *pLeftRank)
{
    MPI_Cart_shift(newCommunicator, HORIZONTAL_DIRECTION, DIRECT_NEIGHBORS, pLeftRank, pRightRank);
    fflush(stdout);
}

int getRankFromGivenCoordinates(MPI_Comm newCommunicator, int *rowIndex, int *colIndex)
{
    int newRank, coordinatesArray[2];
    coordinatesArray[0] = *rowIndex;
    coordinatesArray[1] = *colIndex;
    MPI_Cart_rank(newCommunicator, coordinatesArray, &newRank);
    fflush(stdout);
    return newRank;
}

void getCoordinatesFromGivenRank(MPI_Comm newCommunicator, int rank, int *rowIndex, int *colIndex)
{
    int coordinatesArray[2];
    MPI_Cart_coords(newCommunicator, rank, MAX_NUM_OF_DIMENSIONS, coordinatesArray);
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