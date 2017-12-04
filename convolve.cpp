/**
 * CPSC 501 Assignment 4 - Optimizing Program Performance
 * Baseline program:
 * Time Domain Convolution of .wav audio files via Input Side Algorithm.
 * Program takes command-line for names of input recording and impulse response files, 
 * and produces a convolved signal
 * 
 * Created by Johnny Chung - based on demo code from Leonard Manzara
 */

//compile all files: g++ *.cpp -o convolve
//run: ./convolve <inputfile.wav> <IRfile.wav> <outputfile.wav>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "WavFile.h"

using namespace std;

#define DEBUG_MODE


/*Function declarations */
int checkExtension(char* fileName);
int setExtensionFlag(int* fileExtensions);
void createOutputWAV();
void convolve(short x[], int N, short h[], int M, short y[], int P);

/*Instance variables */
WavFile *inputFile;
WavFile *impulseFile;

int main(int argc, char* argv[]){

    char* inputFileName;
    char* impulseFileName;
    char* outputFileName;

    //check command-line arguments
    if(argc <= 3){
        printf("Usage: ./convolve <inputfile.wav> <IRfile.wav> <outputfile.wav>\n");
        return EXIT_FAILURE;
    }
    else{
        inputFileName = argv[1];
        impulseFileName = argv[2];
        outputFileName = argv[3];

        //check file extensions of command-line arguments
        int inputFileExt = checkExtension(inputFileName);
        int impulseFileExt = checkExtension(impulseFileName);
        int outputFileExt = checkExtension(outputFileName);

        int fileExtensions[3] = {inputFileExt, impulseFileExt, outputFileExt};
        int fileExtFlag = setExtensionFlag(fileExtensions);

        if(fileExtFlag == 1){
            printf("Usage: ./convolve <inputfile.wav> <IRfile.wav> <outputfile.wav>\n");       
            printf("returning %d\n", EXIT_FAILURE);
            return EXIT_FAILURE;         
        }
        else{
            //read input file
            inputFile = new WavFile();
            inputFile->readWAV(inputFileName);

            //read IR file
            impulseFile = new WavFile();
            impulseFile->readWAV(impulseFileName);

            printf("\nInput Size: %d, Impulse Size: %d\n", inputFile->signalSize, impulseFile->signalSize);

            createOutputWAV();
        }
    }
}


void convolve(short x[], int N, short h[], int M, short y[], int P){
    int n, m;

    //check size of output buffer for size P = N + M - 1
    if(P != (N + M - 1)){
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;        
    }

    //clear output buffer
    for(n = 0; n < P; n++){
        y[n] = 0;
    }

    printf("Starting convolution loops...\n");
    clock_t startTime = clock();

    //Outer loop: process each x[n]
    for(n = 0; n < N; n++){
        //inner loop: process each h[m] for current x[n]
        for(m = 0; m < M; n++){
            y[n+m] += x[n] * h[m];

            //provide periodic printout for long convolution
            clock_t currentTime = clock();
            clock_t elapsedTime = currentTime - startTime;
            
            if((m == 100000) && ((n%200)==0) && ((elapsedTime%30) == 0)){
                printf("Convolving %d...\n", (n+m));
            }

        } 
    }


}

int checkExtension(char* fileName){
    
    printf("Checking file extension for %s...", fileName);
    const char* fileExtension = strrchr(fileName, '.');
    
    if(!fileExtension){
        //no file extension
        printf("No file extension found!\n");
        return EXIT_FAILURE;
    }
    else if((strcmp(fileExtension, ".wav")) != 0){
        //wrong file extension
        printf("Wrong file extension found!\n");    
        return EXIT_FAILURE;
    }
    else {
        printf("correct file extension\n");
        return EXIT_SUCCESS;
    }
       
}

int setExtensionFlag(int* fileExtensions){
    int flag = 0;  
    for(int i = 0; i < sizeof(fileExtensions); i++){
        if(fileExtensions[i] == EXIT_FAILURE)
            flag = 1;
    }

    return flag;
}

void createOutputWAV(){
    //P = N + M -1
    int outputSize = (inputFile->signalSize) + (impulseFile->signalSize) - 1;
    short* outputSignal = new short[outputSize];

    printf("made it here...\n");
    //TODO: normalize impulse before convolving
    
    //TODO: convolve the signals

    //TODO: open file stream 

    //TODO: write header for output WAV file
        //use fputs for big endian?
        //write separate functions for little endians (shorts and ints)?

    //TODO: write convolved output signal into output WAV file

    //TODO: close file stream
}


