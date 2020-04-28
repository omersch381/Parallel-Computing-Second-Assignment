#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define FILE_NAME "input.dat"
#define MASTER 0
#define TAG 0
#define MIN_TO_MAX 0
#define BUFFER_SIZE 100

typedef struct box
{
    int id;
    int length;
    int width;
    int height;
};

int floorSqrt(int x);
void initFloatsRead(int n, int my_rank /*, TODO*/);
void scatterValuesToOtherProcesses(struct box boxesArray[], int my_rank);
void printBox(struct box givenBox);
struct box handleFileRead(int n, int my_rank /*TODO*/);
void oddEvenSort(float *myValue, int n, int pLeft, int pRight, int location, int direction);
void compare(float *myValue, float *otherValue, int direction);
// int pack(char *buffer, int *id, int *length, int *width, int *height);
// void unpack(char *buffer, int *id, int *length, int *width, int *height);

int main(int argc, char *argv[])
{
    // Variables
    int my_rank;       /* rank of process */
    int size;          /* number of processes */
    int source;        /* rank of sender */
    int dest;          /* rank of receiver */
    char message[100]; /* storage for message */
    MPI_Status status; /* return status for receive */
    MPI_Comm newCommunicator;
    int numOfDimensionsArray[2], isPeriodicArray[2], reorder;
    int coordinatesArray[2], newId;
    int row, column;
    int pLeftRank, pRightRank;
    int n;

    /* start up MPI */
    MPI_Init(&argc, &argv);

    /* find out process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* create a type for struct box */
    const int nitems = 4;
    int blocklengths[4] = {1, 1, 1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_box_type;
    MPI_Aint offsets[4];

    offsets[0] = offsetof(struct box, id);
    offsets[1] = offsetof(struct box, length);
    offsets[2] = offsetof(struct box, width);
    offsets[3] = offsetof(struct box, height);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_box_type);
    MPI_Type_commit(&mpi_box_type);

    n = floorSqrt(size);
    initFloatsRead(my_rank, n);

    MPI_Type_free(&mpi_box_type);
    /* shut down MPI */
    MPI_Finalize();

    return 0;
}

void initFloatsRead(int n, int my_rank /*, TODO*/)
{
    struct box boxesArray[n];
    boxesArray = handleFileRead(n, my_rank);
    scatterValuesToOtherProcesses(boxesArray, my_rank /*TODO*/);
    // handleNewCommunicatorVariables(/*TODO*/);
}

void scatterValuesToOtherProcesses(struct box boxesArray[], int my_rank)
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    const int nitems = 4;
    int blocklengths[4] = {1, 1, 1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_box_type;
    MPI_Aint offsets[4];

    offsets[0] = offsetof(struct box, id);
    offsets[1] = offsetof(struct box, length);
    offsets[2] = offsetof(struct box, width);
    offsets[3] = offsetof(struct box, height);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_box_type);
    MPI_Type_commit(&mpi_box_type);
    //////////////////////////////////////////////////////////////////////////////////////////////////

    struct box input;
    MPI_Scatter(boxesArray, 1, mpi_box_type, &input, 1, mpi_box_type, MASTER, MPI_COMM_WORLD);
    // data = data + 1;
    // data = 5;
    MPI_Gather(&input, 1, mpi_box_type, boxesArray, 1, mpi_box_type, MASTER, MPI_COMM_WORLD);
    if (my_rank == 0)
        for (int i = 0; i < 16; i++)
            printBox(boxesArray[i]);
}

void printBox(struct box givenBox)
{
    printf("My ID is %d\n", givenBox.id);
    printf("My length is %d\n", givenBox.length);
    printf("My width is %d\n", givenBox.width);
    printf("My heigth is %d\n", givenBox.height);
}

struct box handleFileRead(int n, int my_rank /*TODO*/)
{
    FILE *input_file;
    struct box input;
    struct box boxesArray[n];
    int data = 0, i;

    if (my_rank == 0) //master reads the file
    {
        input_file = fopen(FILE_NAME, "r");
        if (input_file == NULL)
        {
            fprintf(stderr, "\nError opening file\n");
            exit(1);
        }
        while (fread(&input, sizeof(struct body), 1, input_file))
        {
            // printf("id = %d length = %f width = %f height = %f \n", input.id, input.length, input.width, input.height);
            boxesArray[i] = input;
            i++;
        }
        fclose(input_file);
    }
    return boxesArray;
}

// void handleNewCommunicatorVariables(/*TODO*/)
// {
//     numOfDimensionsArray[0] = COLUMNS;
//     numOfDimensionsArray[1] = ROWS;
//     isPeriodicArray[0] = 0;
//     isPeriodicArray[1] = 0;
//     reorder = 1; // probably lets it to reorder
// }

// int pack(char *buffer, int *id, float *length, float *width, float *height)
// {
//     int position = 0;
//     MPI_Pack(id, 1, MPI_FLOAT, buffer, BUFFER_SIZE, &position, MPI_COMM_WORLD);
//     MPI_Pack(length, 1, MPI_FLOAT, buffer, BUFFER_SIZE, &position, MPI_COMM_WORLD);
//     MPI_Pack(width, 1, MPI_FLOAT, buffer, BUFFER_SIZE, &position, MPI_COMM_WORLD);
//     MPI_Pack(height, 1, MPI_FLOAT, buffer, BUFFER_SIZE, &position, MPI_COMM_WORLD);
//     return position;
// }

// void unpack(char *buffer, int *id, float *length, float *width, float *height)
// {
//     int position = 0;
//     MPI_Unpack(buffer, BUFFER_SIZE, &position, id, 1, MPI_FLOAT, MPI_COMM_WORLD);
//     MPI_Unpack(buffer, BUFFER_SIZE, &position, length, 1, MPI_FLOAT, MPI_COMM_WORLD);
//     MPI_Unpack(buffer, BUFFER_SIZE, &position, width, 1, MPI_FLOAT, MPI_COMM_WORLD);
//     MPI_Unpack(buffer, BUFFER_SIZE, &position, height, 1, MPI_FLOAT, MPI_COMM_WORLD);
// }

// int slaveJob(/*TODO*/)
// {
//     MPI_Status status;
//     int id, length, width, height, volume;
//     char buffer[BUFFER_SIZE];
//     MPI_Recv(buffer, BUFFER_SIZE, MPI_PACKED, MASTER, TAG, MPI_COMM_WORLD, &status); // receive
//     unpack(buffer, &id, &length, &width, &height);
//     volume = length * width * height;
// }

// from here and above IDK

void oddEvenSort(float *myValue, int n, int pLeft, int pRight, int location, int direction)
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