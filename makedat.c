#include <stdio.h> 
#include <stdlib.h> 

struct body
{
    int id;
    int length;
    int width;
    int height;
};
  
int main () 
{ 
    FILE *outfile; 
      
    // open file for writing 
    outfile = fopen ("input2.dat", "w"); 
    if (outfile == NULL) 
    { 
        fprintf(stderr, "\nError opend file\n"); 
        exit (1); 
    } 
  
    struct body input1  = {1,1, 1, 5}; 
    struct body input2  = {2,1, 1, 2};
    struct body input3  = {3,1, 1, 9}; 
    struct body input4  = {4,1, 1, 7}; 
    struct body input5  = {5,1, 1, 8}; 
    struct body input6  = {6,1, 1, 11}; 
    struct body input7  = {7,1, 1, 14}; 
    struct body input8  = {8,1, 1, 13}; 
    struct body input9  = {9,1, 1, 1}; 
    struct body input10 = {10,1, 1, 10}; 
    struct body input11 = {11,1, 1, 4}; 
    struct body input12 = {12,1, 1, 3}; 
    struct body input13 = {13,1, 1, 15}; 
    struct body input14 = {14,1, 1, 16}; 
    struct body input15 = {15,1, 1, 12}; 
    struct body input16 = {16,1, 1, 6}; 

      
    // write struct to file 
    fwrite (&input1, sizeof(struct body), 1, outfile); 
    fwrite (&input2, sizeof(struct body), 1, outfile);
    fwrite (&input3, sizeof(struct body), 1, outfile); 
    fwrite (&input4, sizeof(struct body), 1, outfile); 
    fwrite (&input5, sizeof(struct body), 1, outfile); 
    fwrite (&input6, sizeof(struct body), 1, outfile); 
    fwrite (&input7, sizeof(struct body), 1, outfile); 
    fwrite (&input8, sizeof(struct body), 1, outfile); 
    fwrite (&input9, sizeof(struct body), 1, outfile); 
    fwrite (&input10, sizeof(struct body), 1, outfile); 
    fwrite (&input11, sizeof(struct body), 1, outfile); 
    fwrite (&input12, sizeof(struct body), 1, outfile); 
    fwrite (&input13, sizeof(struct body), 1, outfile); 
    fwrite (&input14, sizeof(struct body), 1, outfile); 
    fwrite (&input15, sizeof(struct body), 1, outfile); 
    fwrite (&input16, sizeof(struct body), 1, outfile);  
      
    if(fwrite != 0)  
        printf("contents to file written successfully !\n"); 
    else 
        printf("error writing file !\n"); 
  
    // close file 
    fclose (outfile); 
  
    return 0; 
} 