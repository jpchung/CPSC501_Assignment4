/**
 * CPSC 501 Assignment 4 - Optimizing Program Performance
 * Algorithm-Based Optimization of Baseline program:
 * Frequency Domain Convolution of .wav audio files via Fast Fourier Transform (FFT)
 * Program takes command-line for names of input recording and impulse response files, 
 * and produces a convolved signal, but faster than Baseline Program
 * 
 * Created by Johnny Chung - based on demo code from Leonard Manzara
 */

//compile all files: g++ *.cpp -o convolveFFT
//run: ./convolveFFT <inputfile.wav> <IRfile.wav> <outputfile.wav>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "WavFile.h"

using namespace std;

/*CONSTANTS*/
#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr 

/*Function declarations*/
int checkExtension(char* fileName);
int setExtensionFlag(int* fileExtensions);
void writeWAVHeader(int numChannels, int numSamples, int bitsPerSample, int sampleRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *fileStream);
size_t fwriteShortLSB(short data, FILE* fileStream);
void signalToDouble(WavFile* wav,double signalDouble[]);
int isPowerOfTwo(int num);

int main(int argc, char* argv[]){

    double complexData[SIZE*2];

    char* inputFileName;
    char* impulseFileName;
    char* outputFileName;

    //check command-line arguments
    if(argc <= 3){
        printf("Usage: ./convolveFFT <inputfile.wav> <IRfile.wav> <outputfile.wav>\n");
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
        else {
            //read input and IR wav files
            WavFile* inputFile = new WavFile();
            inputFile->readWAV(inputFileName);

            WavFile* impulseFile = new WavFile();
            impulseFile->readWAV(impulseFileName);

            printf("\nInput Size: %d, Impulse Size: %d\n", inputFile->signalSize, impulseFile->signalSize);
            
            /*
            CONVOLUTION USING FFT
            1. turn Time Domain signals x[n], h[n] into Freq Domain arrays X[k], H[k]
                - ZERO-PAD so that:
                    - x[n] and h[n] have same length
                    - length is a power of 2 (needed for FFT)
                    - length is long enough to avoid circular convolution
            2. multiply X[k] and H[k] point by point (complex multiplication)
            3. convert result Y[k] back to time domain with IFFT
             */

            //gonna use four1 algorithm for FFT/IFFT

            //TODO: Time to Freq Domain transformation (four1 FFT)
            //TODO: Freq Domain convolution (complex multiplication)
            //TODO: inverse FFT (four1)
            //TODO: normalize/scale result
            //TODO: write to file
            
            /* 1.turn Time Domain signals x[n], h[n] into Freq Domain arrays X[k], H[k] */
             
            //turn x[n], h[n] signals into double arrays (need for FFT)
            double x[inputFile->signalSize];
            double h[impulseFile->signalSize];

            signalToDouble(inputFile, x);
            signalToDouble(impulseFile, h);

            printf("converted signals to double\n");

            //zero-pad arrays so that:
            //x[n] and h[n] have same length
            //length is a power of 2 (need for FFT)
            //length is long enough to avoid circular convolution
            int sizeFreqX = inputFile->signalSize;
            int sizeFreqH = impulseFile->signalSize;
            
            int maxLength = 0;
            if(sizeFreqX <= sizeFreqH)
                maxLength = sizeFreqH;
            else
                maxLength = sizeFreqX;

            printf("Max length of signals: %d\n", maxLength);
            //TODO: pad to equal length, check for power of 2

            //check if maxLength is power of 2
            int maxLengthPow2 = 0;
            int pow2 = 0; //minimum power of 2 to produce length at least as large as maxLength
            if(isPowerOfTwo(maxLength) != 1){
                printf("%d NOT a power of 2\n", maxLength);
                //log2 truncated when cast from double to int, so add 1 for next closest power
                pow2 = (int) log2(maxLength) + 1;
                maxLengthPow2 = (int) pow(2,pow2);
                printf("minimum power: 2^(%d) = %d\n", pow2, maxLengthPow2);
            }
            else{
                printf("%d is power of 2\n", maxLength);
                pow2 = (int) log2(maxLength);                
                maxLengthPow2 = maxLength;
            }

            //X[k],H[k] have real & imaginary parts, so length  = maxLength * 2
            //real part is signal, imaginary part should be set to 0's
            double* freqX = new double[maxLengthPow2 *2];
            double* freqH = new double[maxLengthPow2 *2];
            
            //TODO: initialize X[k],H[k] to all zeroes for imaginary part

            //TODO: rewrite every other index with respective signal for real part



        }
    }
}


//check if filename has .wav extension
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

//set flag for valid file extensions
int setExtensionFlag(int* fileExtensions){
    int flag = 0;  
    for(int i = 0; i < sizeof(fileExtensions); i++){
        if(fileExtensions[i] == EXIT_FAILURE)
            flag = 1;
    }

    return flag;
}



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


size_t fwriteShortLSB(short data, FILE* fileStream){

    unsigned char charArray[2];

    //write short (2 bytes) into fileStream in little-endian
    //little endian writes from least significant byte (LSB) first
    charArray[1] = (unsigned char)((data >> 8) & 0xFF);
    charArray[0] = (unsigned char)(data & 0xFF);

    //use charArray to write values as characters to file
    return fwrite(charArray, sizeof(unsigned char), 2, fileStream);
}


void signalToDouble(WavFile* wav,double signalDouble[]){
    for(int i = 0; i < wav->signalSize; i++){
        signalDouble[i] = ((double) wav->signal[i])/32678.0;
    }
}


//check if a number is a power of 2 (divisible by 2 n times)
int isPowerOfTwo(int num){
    if(num == 0)
        return 0;

    //keep dividing num by 2 until num is 1 or num not divisible by 2
    while(num != 1){
        if((num%2) != 0)
            return 0;
        else
            num = num/2;
    }

    return 1;
}