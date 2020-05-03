#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define INPUT_FILE_NAME "cuboids.dat"
#define TAG 0
#define MIN_TO_MAX 0
#define MAX_TO_MIN 1
#define NUM_OF_DIMENSIONS 2
#define MAX_NUM_OF_DIMENSIONS 2
#define DIRECT_NEIGHBORS 1
#define MASTER 0
#define ROW 1
#define COL 0

struct
{
    int id;
    int length;
    int width;
    int height;
} typedef Box;

int floorSqrt(int x);
void compare(Box *myBox, Box *otherBox, int direction);
void sanityCheck(int numOfProcesses);
MPI_Comm getNewCommunicator(int n);
void getCoordinatesFromGivenRank(MPI_Comm newCommunicator, int rank, int *rowIndex, int *colIndex);
int getRankFromGivenCoordinates(MPI_Comm newCommunicator, int *rowIndex, int *colIndex);
void getNeighborsRanksFromGivenRank(MPI_Comm newCommunicator, int rank, int *pRightRank, int *pLeftRank, int direction);
void shearSort(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, Box *myStruct, int n);
void swapBoxes(Box *myBox, Box *otherBox);
void sortUnit(MPI_Comm newCommunicator, int rank, Box *myStruct, int n, int row, int direction, int unit);
void sortUnits(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, Box *myStruct, int n, int unit);
void oddEvenSortGeneral(Box *myStruct, int n, int pSourceRank, int pTargetRank, int location, int direction, int rank, int unit);
void readFromFile(int rank, Box boxesArray[]);
void writeToFile(int rank, Box boxesArray[], int numOfProcesses);

int main(int argc, char *argv[])
{
    int rank, numOfProcesses;
    int rowIndex, colIndex;
    int n;
    Box myStruct, boxesArray[numOfProcesses];
    MPI_Comm newCommunicator;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);

    sanityCheck(numOfProcesses);

    n = floorSqrt(numOfProcesses);

    readFromFile(rank, boxesArray);

    MPI_Scatter(boxesArray, sizeof(myStruct), MPI_BYTE, &myStruct, sizeof(myStruct), MPI_BYTE, MASTER, MPI_COMM_WORLD);

    newCommunicator = getNewCommunicator(n);

    getCoordinatesFromGivenRank(newCommunicator, rank, &rowIndex, &colIndex);

    rank = getRankFromGivenCoordinates(newCommunicator, &rowIndex, &colIndex);

    shearSort(newCommunicator, rank, rowIndex, colIndex, &myStruct, n);

    MPI_Gather(&myStruct, sizeof(myStruct), MPI_BYTE, boxesArray, sizeof(myStruct), MPI_BYTE, MASTER, MPI_COMM_WORLD);

    writeToFile(rank, boxesArray, numOfProcesses);

    MPI_Finalize();

    return 0;
}

void writeToFile(int rank, Box boxesArray[], int numOfProcesses)
{
    if (rank == 0)
    {
        int i;
        for (i = 0; i < numOfProcesses; i++)
            printf("%d\n", boxesArray[i].id); //TODO
    }
}

void readFromFile(int rank, Box boxesArray[])
{
    Box input;
    int i;
    FILE *input_file;
    if (rank == 0) //master reads the file
    {
        input_file = fopen(INPUT_FILE_NAME, "r");
        if (input_file == NULL)
        {
            fprintf(stderr, "\nError opening file\n");
            exit(1);
        }
        while (fread(&input, sizeof(Box), 1, input_file))
        {
            boxesArray[i] = input;
            i++;
        }
        fclose(input_file);
    }
}

void shearSort(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, Box *myStruct, int n)
{
    int i;
    for (i = 0; i < n / 2; i++)
    {
        sortUnits(newCommunicator, rank, rowIndex, colIndex, myStruct, n, ROW);
        sortUnits(newCommunicator, rank, rowIndex, colIndex, myStruct, n, COL);
    }
    sortUnits(newCommunicator, rank, rowIndex, colIndex, myStruct, n, ROW);
}

void sortUnits(MPI_Comm newCommunicator, int rank, int rowIndex, int colIndex, Box *myStruct, int n, int unit)
{
    int row;
    int direction;
    for (row = 0; row < n; row++)
    {
        if (unit == ROW)
            direction = row % 2 == 0 ? MIN_TO_MAX : MAX_TO_MIN;
        else // unit == COL
            direction = MIN_TO_MAX;
        sortUnit(newCommunicator, rank, myStruct, n, row, direction, unit);
    }
}
void sortUnit(MPI_Comm newCommunicator, int rank, Box *myStruct, int n, int row, int direction, int unit)
{
    int pTargetRank, pSourceRank;
    int rowArray[n];
    int col;
    for (col = 0; col < n; col++)
    {
        int currentRank = getRankFromGivenCoordinates(newCommunicator, &row, &col); //(0,0)
        if (rank == currentRank)
        {
            getNeighborsRanksFromGivenRank(newCommunicator, currentRank, &pTargetRank, &pSourceRank, unit);
            int location = row * n + col; // converting to the location paramerter that oddEvenSort accept
            oddEvenSortGeneral(myStruct, n, pSourceRank, pTargetRank, location, direction, rank, unit);
        }
    }
}

void oddEvenSortGeneral(Box *myStruct, int n, int pSourceRank, int pTargetRank, int location, int direction, int rank, int unit)
{

    MPI_Status status;
    Box *otherStruct = (Box *)malloc(sizeof(Box));
    int i;
    for (i = 0; i < n; i++)
    {
        if (i % 2 == 0) // An Even Iteration
        {
            int sortRowsAndEvenLocation = ((location % 2 == 0) && unit == ROW);
            int sortColsAndNotOnEdges = (((location < n) || (location / n) % 2 == 0) && unit == COL);
            if (sortRowsAndEvenLocation || sortColsAndNotOnEdges)
            {
                MPI_Recv(otherStruct, sizeof(Box), MPI_BYTE, pTargetRank, TAG, MPI_COMM_WORLD, &status);
                compare(myStruct, otherStruct, direction);
                MPI_Send(otherStruct, sizeof(Box), MPI_BYTE, pTargetRank, TAG, MPI_COMM_WORLD);
            }
            else // if the location is odd
            {
                MPI_Send(myStruct, sizeof(Box), MPI_BYTE, pSourceRank, TAG, MPI_COMM_WORLD);
                MPI_Recv(myStruct, sizeof(Box), MPI_BYTE, pSourceRank, TAG, MPI_COMM_WORLD, &status);
            }
        }
        else // if odd iteration
        {
            int edgeCaseCol = (location < n || location >= n * (n - 1));
            int edgeCaseRow = (location == 0 || location % n == n - 1);
            if (((!edgeCaseCol) && unit == COL) || ((!edgeCaseRow) && unit == ROW)) // edge cases do not participate in odd iterations
            {
                int sortRowsAndEvenLocation = ((location % 2 == 0) && unit == ROW);
                int sortColsAndNotOnEdges = (((location < n) || (location / n) % 2 == 0) && unit == COL);
                if (sortColsAndNotOnEdges || sortRowsAndEvenLocation)
                {
                    MPI_Send(myStruct, sizeof(Box), MPI_BYTE, pSourceRank, TAG, MPI_COMM_WORLD);
                    MPI_Recv(myStruct, sizeof(Box), MPI_BYTE, pSourceRank, TAG, MPI_COMM_WORLD, &status);
                }
                else // if location is odd
                {
                    // The odd location on the odd iteration sorts and sends.
                    MPI_Recv(otherStruct, sizeof(Box), MPI_BYTE, pTargetRank, TAG, MPI_COMM_WORLD, &status);
                    compare(myStruct, otherStruct, direction);
                    MPI_Send(otherStruct, sizeof(Box), MPI_BYTE, pTargetRank, TAG, MPI_COMM_WORLD);
                }
            }
        }
    }
    free(otherStruct);
}

void compare(Box *myBox, Box *otherBox, int direction)
{
    float myVolume = myBox->height * myBox->length * myBox->width;
    float otherVolume = otherBox->height * otherBox->length * otherBox->width;

    if (direction == MIN_TO_MAX)
    {
        if (myVolume > otherVolume)
            swapBoxes(myBox, otherBox);
        else if (myVolume == otherVolume && myBox->height > otherBox->height)
            swapBoxes(myBox, otherBox);
    }
    else // direction == MAX_TO_MIN
    {
        if (myVolume < otherVolume) // which is not what we want
            swapBoxes(myBox, otherBox);
        else if (myVolume == otherVolume && myBox->height < otherBox->height)
            swapBoxes(myBox, otherBox);
    }
}

void swapBoxes(Box *myBox, Box *otherBox)
{
    Box temp = *otherBox;
    *otherBox = *myBox;
    *myBox = temp;
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
void sanityCheck(int numOfProcesses)
{
    if (numOfProcesses < 4)
    {
        printf("Please Run With 4 Or More Processes! Aborting....\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int areTheNumberOfProcessesPowerOfTwo = ceil(log2(numOfProcesses)) == floor(log2(numOfProcesses));
    if (!areTheNumberOfProcessesPowerOfTwo) // which means that numOfProcesses is not a power of 2
    {
        printf("Please Run With a Number of Processes That is a Power of 2! Aborting....\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
void getNeighborsRanksFromGivenRank(MPI_Comm newCommunicator, int rank, int *pRightRank, int *pLeftRank, int direction)
{
    MPI_Cart_shift(newCommunicator, direction, DIRECT_NEIGHBORS, pLeftRank, pRightRank);
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
