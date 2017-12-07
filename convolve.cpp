/**
 * CPSC 501 Assignment 4 - Optimizing Program Performance
 * Baseline program:
 * Time Domain Convolution of .wav audio files via Input Side Algorithm.
 * Program takes command-line for names of input recording and impulse response files, 
 * and produces a convolved signal
 * 
 * Created by Johnny Chung - based on demo code from Leonard Manzara
 */

/*
compile all files: g++ *.cpp -o convolve
run: ./convolve <inputfile.wav> <IRfile.wav> <outputfile.wav>
profiler: g++ -p *.cpp -o convolve
          ./convolve <inputfile.wav> <IRfile.wav> <outputfile.wav>
          gprof convolve
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <fstream>
#include <iostream>
#include "WavFile.h"

using namespace std;

/*CONSTANTS */
#define DEBUG_MODE


/*Function declarations */
int checkExtension(char* fileName);
int setExtensionFlag(int* fileExtensions);
void convolve(double x[], int N, double h[], int M, double y[], int P);
void createOutputWAV(char* fileName);
void writeWAVHeader(int numChannels, int numSamples, int bitsPerSample, int sampleRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *fileStream);
size_t fwriteShortLSB(short data, FILE* fileStream);
void signalToDouble(WavFile* wav,double signalDouble[]);
void scaleOutputSignal(double* outputSignal, WavFile* inputFile, int outputSize);

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

            //convolve and generate output WAV file
            createOutputWAV(outputFileName);

            return 0;
        }
    }
}


/** 
 * Function: checkExtension
 * Description: check file name contains the .wav extension
 * Parameters:
 *          char* fileName - file name from command-line arg
 */
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


/**
 * Function: setExtensionFlag
 * Description: set flag in main for valid file extensions for command-line args
 * Parameters:
 *          int* fileExtensions - pointer to list of file extension flags
 */
int setExtensionFlag(int* fileExtensions){
    int flag = 0;  
    for(int i = 0; i < sizeof(fileExtensions); i++){
        if(fileExtensions[i] == EXIT_FAILURE)
            flag = 1;
    }

    return flag;
}


/**
 * Function: createOutputWAV
 * Description: write header and convolved sound data to output WAV file
 * Parameters:
 *          char* fileName - name of output file
 */
void createOutputWAV(char* fileName){
    //P = N + M -1
    int outputSize = (inputFile->signalSize) + (impulseFile->signalSize) - 1;
    short* outputSignal = new short[outputSize];
    double* yDouble = new double[outputSize];
   
    //convert input and impulse signals to double before convolving
    double* xDouble = new double[inputFile->signalSize];
    signalToDouble(inputFile, xDouble);

    double* hDouble = new double[impulseFile->signalSize];
    signalToDouble(impulseFile, hDouble);

    //convolve the signals
    time_t startTime, endTime;
    time(&startTime);
    
    convolve(xDouble, inputFile->signalSize, hDouble, impulseFile->signalSize, yDouble, outputSize);

    time(&endTime);
    double elapsed = difftime(endTime, startTime);
    printf("DONE convolution in %.2f seconds!\n\n", elapsed);


    scaleOutputSignal(yDouble, inputFile, outputSize);

    for(int i = 0; i < outputSize; i++){
        
        //double ySample = yDouble[i] * 32767;
        //outputSignal[i] = (short)ySample ;
        outputSignal[i] = (short) yDouble[i];
    }
   
    //open file stream 
    FILE* outputFile = fopen(fileName, "wb");


    //write header for output WAV file
    printf("Writing WAV header for %s...\n", fileName);
    writeWAVHeader(inputFile->numChannels, outputSize, inputFile->bitsPerSample, inputFile->sampleRate, outputFile);


    
    //write convolved output signal into output WAV file
    //use LSB method since signal is short array
    printf("Writing convolved signal to %s\n", fileName);
    for(int i = 0; i < outputSize; i++){
        fwriteShortLSB(outputSignal[i], outputFile);
    }

    printf("Done writing %s!\n", fileName);

    //close file stream
    fclose(outputFile);
}


/** 
 * Function: writeWAVHeader
 * Description: writes the header for WAV file format to output file
 * Parameters:
 *          int numChannels - number of audio channels
 *          int numSamples - number of sound data samples
 *          int bitsPerSample - number of bits per sound data sample
 *          int sampleRate - sampling frequency
 *          FILE* outputFile - file stream for output     
 */
void writeWAVHeader(int numChannels, int numSamples, int bitsPerSample, int sampleRate, FILE *outputFile){

    //calculations for header fields

    //subchunk sizes
    int subChunk2Size = numChannels * numSamples * (bitsPerSample/8); 
    int chunkSize = 36 + subChunk2Size; //36 is sum of all field sizes before data
    
    //number of bytes for one sample after accounting all channels (frame size)
    short blockAlign = numChannels * (bitsPerSample/8);
    
    //bytes per second  = sample rate * total bytes per sample
    int byteRate = (int) sampleRate * blockAlign;

    //write header to file
    //use fputs for big endian fields, use respective LSB methods for little endian

    //RIFF chunk descriptor
    fputs("RIFF", outputFile);
    fwriteIntLSB(chunkSize, outputFile);
    fputs("WAVE", outputFile);
    
    //fmt subchunk
    fputs("fmt ", outputFile);
    fwriteIntLSB(16, outputFile); //subchunk1size should be fixed 16 bytes
    fwriteShortLSB(1, outputFile); // AudioFormat = 1 for PCM
    fwriteShortLSB(numChannels, outputFile);
    fwriteIntLSB(sampleRate, outputFile);
    fwriteIntLSB(byteRate, outputFile);
    fwriteShortLSB(blockAlign, outputFile);
    fwriteShortLSB(bitsPerSample, outputFile);

    //data subchunk
    fputs("data", outputFile);
    fwriteIntLSB(subChunk2Size, outputFile);

}


/** 
 * Function: fwriteIntLSB
 * Description: Writes a 4-byte integer to the file stream, starting
 *              with the least significant byte (i.e. writes the int
 *              in little-endian form).  This routine will work on both
 *              big-endian and little-endian architecture
 * Parameters:
 *          int data - int to write to file stream
 *          FILE* fileStream - file stream
 */
size_t fwriteIntLSB(int data, FILE *fileStream){

    unsigned char charArray[4];

    //write int (4 bytes) into fileStream in little-endian
    //little endian writes from least significant byte (LSB) first
    charArray[3] = (unsigned char)((data >> 24) & 0xFF);
    charArray[2] = (unsigned char)((data >> 16) & 0xFF);
    charArray[1] = (unsigned char)((data >> 8) & 0xFF);
    charArray[0] = (unsigned char)(data & 0xFF);

    //use charArray to write values as characters to file
    return fwrite(charArray, sizeof(unsigned char), 4, fileStream);
}


/** 
 * Function: fwriteShortLSB
 * Description: Writes a 2-byte integer to the file stream, starting
 *              with the least significant byte (i.e. writes the int
 *              in little-endian form).  This routine will work on both
 *              big-endian and little-endian architectures.
 * Parameters:
 *          short data - short to write to file stream
 *          FILE* fileStream - file stream
 */
size_t fwriteShortLSB(short data, FILE* fileStream){

    unsigned char charArray[2];

    //write short (2 bytes) into fileStream in little-endian
    //little endian writes from least significant byte (LSB) first
    charArray[1] = (unsigned char)((data >> 8) & 0xFF);
    charArray[0] = (unsigned char)(data & 0xFF);

    //use charArray to write values as characters to file
    return fwrite(charArray, sizeof(unsigned char), 2, fileStream);
}


/**
 * Function: convolve
 * Description: convolve input and impulse signals via input side algorithm
 * Parameters:
 *          double x[] - input signal array x[n]
 *          int N -  size of x[n]
 *          double h[]  - impulse signal array h[n]
 *          int M - size of h[n]
 *          double y[] - output signal y[n]
 *          int P - size of y[n] (must be P = N + M -1)    
 */
void convolve(double x[], int N, double h[], int M,double y[], int P){

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
        y[n] = 0.0;
    }
       
    printf("Starting convolution loops...\n");
    time_t startTime = time(NULL);
    
    //Outer loop: process each x[n]
    for(n = 0; n < N; n++){
        //inner loop: process each h[m] for current x[n]
        for(m = 0; m < M; m++){
                
            y[n+m] += x[n] * h[m];
            
            //provide periodic printout for long convolution
            time_t currentTime = time(NULL);
            time_t elapsedTime = currentTime - startTime;
                
            if((m == 100000) && ((n%200)==0) && ((elapsedTime%30) == 0)){
                printf("Convolving %d...\n", (n+m));
            }
    
        } 
    }
    
    
}


/** 
 * Function: signalToDouble
 * Description: Convert wav file sound data from short to double array
 * Parameters:
 *          WavFile* wav - object instance of WAV file
 *          double signalDouble[] - output double array for wav file sound data 
 */
void signalToDouble(WavFile* wav,double signalDouble[]){
    for(int i  = 0; i < (wav->signalSize); i++){
        signalDouble[i] = ((double) wav->signal[i])/32678.0;
    }
}


/** 
 * Function: scaleOutputSignal
 * Description: scale output signal to remove overflows/clipping in resulting audio
 * Parameters: 
 *          double* outputSignal - convolved signal for output sound data
 *          WavFile* inputFile - object for input signal WAV file
 *          int outputSize - number of sound data samples in output signal
 */
void scaleOutputSignal(double* outputSignal, WavFile* inputFile, int outputSize){
    double inputMaxValue = 0.0;
    double outputMaxValue = 0.0;

    //check for max value in both original and output signals
    for(int i = 0; i < outputSize; i++){
        if(inputFile->signal[i] > inputMaxValue)
            inputMaxValue = inputFile->signal[i];

        if(outputSignal[i] > outputMaxValue)
            outputMaxValue = outputSignal[i];
    }
    
    for(int i  = 0; i < outputSize; i++){
        outputSignal[i] = outputSignal[i] / outputMaxValue * inputMaxValue;
    }

}


