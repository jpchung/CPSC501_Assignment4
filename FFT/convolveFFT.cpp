/**
 * CPSC 501 Assignment 4 - Optimizing Program Performance
 * Algorithm-Based Optimization of Baseline program:
 * Frequency Domain Convolution of .wav audio files via Fast Fourier Transform (FFT)
 * Program takes command-line for names of input recording and impulse response files, 
 * and produces a convolved signal, but faster than Baseline Program
 * 
 * Created by Johnny Chung - based on demo code from Leonard Manzara
 */

/*
compile all files: g++ *.cpp -o convolveFFT
run: ./convolveFFT <inputfile.wav> <IRfile.wav> <outputfile.wav>
profiler: g++ -p *.cpp -o convolveFFT
          ./convolveFFT <inputfile.wav> <IRfile.wav> <outputfile.wav>
          gprof convolveFFT
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
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
void four1(double data[], int nn, int isign);
void convolveFreqs(double* freqX, double *freqH, double* freqY, int arrayLength);
void scaleOutputFreq(double* outFreq, WavFile* inputFile, int numSamples);
void createOutputWAV(char* fileName, double* outputSignal, int numSamples, WavFile* inputFile);


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
            
            /* 1.Time Domain to Freq Domain Transformation (FFT)*/

            //Turn time domain signals x[n], h[n] into freq domain X[k], H[k]
            //zero-pad arrays so that:
            //x[n] and h[n] have same length
            //length is a power of 2 (need for FFT)
            //length is long enough to avoid circular convolution

            //turn x[n], h[n] signals into double arrays (need for FFT)           
            double* x = new double[inputFile->signalSize];
            double* h = new double[impulseFile->signalSize];

            signalToDouble(inputFile, x);
            signalToDouble(impulseFile, h);

            printf("converted signals to double\n");

            int sizeFreqX = inputFile->signalSize;
            int sizeFreqH = impulseFile->signalSize;
            
            int maxLength = 0;
            if(sizeFreqX <= sizeFreqH)
                maxLength = sizeFreqH;
            else
                maxLength = sizeFreqX;

            printf("Max length of signals: %d\n", maxLength);


            //find next largest power of 2 that's at least as large as maxLength
            int pow2 = 1;
            while(pow2 < maxLength){
                pow2 *= 2;
            }

            //double that power of 2 for length of freq arrays X[k],H[k]
            //i.e. make sure freq arrays long enough for real and imaginary parts
            int maxLengthPow2 = pow2*2;
            printf("pow2: %d, maxLengthPow2: %d\n", pow2, maxLengthPow2);

            double* freqX = new double[maxLengthPow2];
            double* freqH = new double[maxLengthPow2];
            
            //zero pad X[k],H[k] to all zeroes for imaginary part
            //for(int i = 0; i < maxLengthPow2*2; i++){
            for(int i  = 0; i < maxLengthPow2; i++){
               freqX[i] = 0.0;
            }

            for(int i = 0; i < maxLengthPow2; i++){
                freqH[i] = 0.0;
            }
            //rewrite every other index with respective signal for real part
            //index i = real, index i+1 = imaginary
            for(int i = 0; i < sizeFreqX; i++){
                freqX[i*2] = x[i];
            }

            for(int i = 0; i < sizeFreqH; i++){
                freqH[i*2] = h[i];
            }

            //complete freq domain transformation using four1 algorithm
            printf("Performing FFT for %s signal...\n", inputFileName);
            four1(freqX-1, pow2, 1);

            printf("Performing FFT for %s signal...\n", impulseFileName);
            four1(freqH-1, pow2, 1);

            printf("FFT for input and impulse signals complete!\n");

            /*2. Frequency Domain Convolution */

            //Create Y[k] by multiplying X[k] and H[k] point by point (complex multiplication)

            clock_t startTime;
            double* freqY = new double[maxLengthPow2];
            startTime = clock();

            convolveFreqs(freqX, freqH, freqY, maxLengthPow2);
            
            double elapsed = clock() - startTime;
            printf("DONE FFT convolution in %.3f seconds!\n\n", elapsed/CLOCKS_PER_SEC);
            
            /*3. Freq Domain to Time Domain Transformation (IFFT)*/
            
            //use inverse four1 to turn Y[k] to y[n]
            printf("Performing IFFT on output frequency...\n");
            four1(freqY-1, pow2, -1);
            printf("IFFT on output frequency complete!\n");
            

            
            //separate convolved signal (i.e. real indices) from freqY
            int outputSize = (inputFile->signalSize) + (impulseFile->signalSize) - 1;   
            double outputMaxValue = 0.0;

            double* freqYSignal = new double[outputSize];
            for(int i = 0; i < outputSize; i++){
                double signalSample = freqY[i*2];
                freqYSignal[i] = signalSample;

                if(outputMaxValue < abs(signalSample)){
                    outputMaxValue = abs(signalSample);
                }
            }

            //scale convolved output signal
            scaleOutputFreq(freqYSignal, inputFile, outputSize);
            
            //write to file
            createOutputWAV(outputFileName, freqYSignal, outputSize, inputFile);
            
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
 * Function: signalToDouble
 * Description: Convert wav file sound data from short to double array
 * Parameters:
 *          WavFile* wav - object instance of WAV file
 *          double signalDouble[] - output double array for wav file sound data 
 */
void signalToDouble(WavFile* wav,double signalDouble[]){   
    for(int i = 0; i < (wav->signalSize); i++){
        signalDouble[i] = ((double) wav->signal[i])/32678.0;

    }
}


//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).
void four1(double data[], int nn, int isign){

    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}


/** 
 * Function: convolveFreqs
 * Description: convolve freq domain signals X[k], H[k] using complex multiplication to produce Y[k]
 * Parameters:
 *          double* freqX - freq domain signal X[k] for input sound data
 *          double* freqH - freq domain signal H[k] for impulse sound data
 *          double* freqY - freq domain signal Y[k] for output sound data
 *          int arrayLength - padded length of  input/impulse freq domain arrays 
 */
void convolveFreqs(double* freqX, double *freqH, double* freqY, int arrayLength){
    printf("Starting convolution loops...\n");
    //complex multiplication:
    //realY[k] = realX[k]*realH[k] - imX[k]*imH[k]
    //imY[k] = imX[k]*realH[k] + realX[k]*imH[k]
    for(int i = 0; i < arrayLength; i+= 2){
        //index i = real, index i+1 = imaginary
        freqY[i] = freqX[i] * freqH[i] - freqX[i+1] * freqH[i+1]; //real
        freqY[i+1] = freqX[i+1] * freqH[i] + freqX[i] * freqH[i+1]; //imaginary
    
        if((i%100000) == 0)
            printf("Convolving %d...\n", i);
    }
}


/** 
 * Function: scaleOutputFreq
 * Description: scale output freq signal to remove overflows/clipping in resulting audio
 * Parameters: 
 *          double* outFreq - convolved signal for output sound data
 *          WavFile* inputFile - object for input signal WAV file
 *          int numSamples - number of sound data samples in output signal
 */
void scaleOutputFreq(double* outFreq, WavFile* inputFile, int numSamples){
    double inputMaxValue = 0.0;
    double outputMaxValue = 0.0;

    //check for max value in both original and output signals
    for(int i = 0; i < numSamples; i++){
        if(inputFile->signal[i] > inputMaxValue)
            inputMaxValue = inputFile->signal[i];

        if(outFreq[i] > outputMaxValue)
            outputMaxValue = outFreq[i];
    }

    //scale output frequency to maximums (divide by outputMax, multiply by inputMax)
    //should hold that outputMax > inputMax, so will scale down
    for(int i  = 0; i < numSamples; i++){
        outFreq[i] = outFreq[i] / outputMaxValue * inputMaxValue;
    }

}


/**
 * Function: createOutputWAV
 * Description: write header and convolved sound data to output WAV file
 * Parameters:
 *          char* fileName - name of output file
 *          double* outputSignal - convolved and scaled output sound data
 *          int numSamples - number of sound data samples in output signal
 *          WavFile* inputFile - object for input signal WAV file
 */
void createOutputWAV(char* fileName, double* outputSignal, int numSamples, WavFile* inputFile){
    
    //open file stream
    FILE* outputFile = fopen(fileName, "wb");

    //write header for WAV output file
    printf("Writing WAV Header for %s...\n", fileName);
    writeWAVHeader(inputFile->numChannels, numSamples, inputFile->bitsPerSample, inputFile->sampleRate, outputFile);

    //write convolved output signal into output WAV file
    //use LSB method since signal needs to be short
    printf("Writing convolved signal to %s\n", fileName);
    for(int i = 0; i < numSamples; i++){
        fwriteShortLSB((short) outputSignal[i], outputFile);
    }
    
    //close file stream
    printf("Done writing %s!\n", fileName);
    fclose(outputFile);
}
