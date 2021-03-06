/**
 * CPSC 501 Assignment 4 - Optimizing Program Performance
 * Regression Testing program:
 * This program takes command-line for names of output WAV files 
 * generated by Baseline program and FFT algorithm optimized program,
 * and compares the header fields to see if the two files match
 * 
 * Created by Johnny Chung
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

/*Function Declarations */
int checkExtension(char* fileName);
int compareRIFF(WavFile* baseFile, WavFile* fftFile);
int compareSubchunk1(WavFile* baseFile, WavFile* fftFile);
int compareSubchunk2(WavFile* baseFile, WavFile* fftFile);
int compareSignalData(WavFile* baseFile, WavFile* fftFile);


int main(int argc, char* argv[]){
    char* baseFileName;
    char* fftFileName;

    //check number of command-line args
    if(argc <= 2){
        printf("Usage: ./regressionTest <outputFileBase.wav> <outputFileFFT.wav>\n");
        return EXIT_FAILURE;
    }
    else{

        baseFileName = argv[1];
        fftFileName = argv[2];

        //check file extensions for file names
        int extBase = checkExtension(baseFileName);
        int extFFT = checkExtension(fftFileName);

        if((extBase == EXIT_FAILURE) || (extFFT == EXIT_FAILURE)){
            printf("Usage: ./regressionTest <outputFileBase.wav> <outputFileFFT.wav>\n");
            return EXIT_FAILURE;
        }
        else{
            //read the output files
            WavFile* baseFile = new WavFile();
            baseFile->readWAV(baseFileName);

            WavFile* fftFile = new WavFile();
            fftFile->readWAV(fftFileName);


            printf("\nStarting WAV header comparisons...\n");
            int equivalentFiles = 0;
            
            //compare RIFF fields for both files
            if(compareRIFF(baseFile, fftFile) != 0){
                printf("RIFF chunk descriptors NOT equal!\n");
                equivalentFiles = 1;
            }
            else
                printf("RIFF chunk descriptors equal...\n");
            
            //compare subchunk1 fields for both files
            if(compareSubchunk1(baseFile, fftFile) != 0){
                printf("fmt subchunks NOT equal!\n");
                equivalentFiles = 1;
            }
            else
                printf("fmt subchunks equal...\n");

            //compare subchunk2 fields (except data) for both files
            if(compareSubchunk2(baseFile, fftFile) != 0){
                printf("data subchunks NOT equal!\n");
                equivalentFiles = 1;
            }
            else
                printf("data subchunks equal...\n");


            //compare the signal data for both files
            printf("\nStarting signal data comparisons...\n");
            if(compareSignalData(baseFile, fftFile) != 0){
                printf("output signals NOT equal!\n");
                equivalentFiles = 1;
            }
            else
                printf("output signals equal!\n");

            if(equivalentFiles != 0)
                printf("\n%s and %s NOT equivalent files!\n", baseFileName, fftFileName);
            else{
                printf("\n%s equivalent to %s!\n", baseFileName, fftFileName);
            }

            return EXIT_SUCCESS;
        
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
 * Function: compareRIFF
 * Description: compare the RIFF chunk descriptor fields of the base and fft WAV files
 * Parameters:
 *          WavFile* baseFile - object of baseline output WAV file
 *          WavFile* fftFile - object of FFT algorithm optimized output WAV file
 */
int compareRIFF(WavFile* baseFile, WavFile* fftFile){
    int flagRIFF = 0;

    //compare fields
    if(strcmp(baseFile->chunkID, fftFile->chunkID) != 0)
        flagRIFF = 1;
     
    if(baseFile->chunkSize != fftFile->chunkSize)
        flagRIFF = 1;

    if(strcmp(baseFile->format, fftFile->format) != 0)
        flagRIFF = 1;

    return flagRIFF;
}


/**
 * Function: compareSubChunk1
 * Description: compare the fmt subchunk fields of the base and fft WAV files
 * Parameters:
 *          WavFile* baseFile - object of baseline output WAV file
 *          WavFile* fftFile - object of FFT algorithm optimized output WAV file
 */
int compareSubchunk1(WavFile* baseFile, WavFile* fftFile){
    int flagSubChunk1 = 0;

    //compare fields
    if(strcmp(baseFile->subChunk1ID, fftFile->subChunk1ID) != 0)
        flagSubChunk1 = 1;

    if(baseFile->subChunk1Size != fftFile->subChunk1Size)
        flagSubChunk1 = 1;

    if(baseFile->audioFormat != fftFile->audioFormat)
        flagSubChunk1 = 1;

    if(baseFile->numChannels != fftFile->numChannels)
        flagSubChunk1 = 1;
    
    if(baseFile->sampleRate != fftFile->sampleRate)
        flagSubChunk1 = 1;

    if(baseFile->byteRate != fftFile->byteRate)
        flagSubChunk1 = 1;

    if(baseFile->blockAlign != fftFile->blockAlign)
        flagSubChunk1 = 1;

    if(baseFile->bitsPerSample != fftFile->bitsPerSample)
        flagSubChunk1 = 1;

    return flagSubChunk1;
}


/**
 * Function: compareSubChunk2
 * Description: compare the data subchunk fields of the base and fft WAV files
 * Parameters:
 *          WavFile* baseFile - object of baseline output WAV file
 *          WavFile* fftFile - object of FFT algorithm optimized output WAV file
 */
int compareSubchunk2(WavFile* baseFile, WavFile* fftFile){
    int flagSubChunk2 = 0;

    //compare fields
    if(strcmp(baseFile->subChunk2ID, fftFile->subChunk2ID) != 0)
        flagSubChunk2 = 1;

    if(baseFile->subChunk2Size != fftFile->subChunk2Size)
        flagSubChunk2 = 1;

    return flagSubChunk2;
}


/**
 * Function: compareSignalData
 * Description: compare the sound signal data of the base and fft WAV files
 * Parameters:
 *          WavFile* baseFile - object of baseline output WAV file
 *          WavFile* fftFile - object of FFT algorithm optimized output WAV file
 */
int compareSignalData(WavFile* baseFile, WavFile* fftFile){
    int flagSignalData = 0;

    int baseSignalSize = baseFile->signalSize;
    int fftSignalSize = fftFile->signalSize;

    if(baseSignalSize != fftSignalSize){
        printf("baseFile signal size (%d) NOT equal to fftFile signal size (%d)\n", baseSignalSize, fftSignalSize);
        flagSignalData = 1;
    }
    else {
        for(int i = 0; i < baseSignalSize; i++){
            if(baseFile->signal[i] != fftFile->signal[i])
                flagSignalData = 1;
        }

        if(flagSignalData == 1)
            printf("base signalData NOT identical to fft signalData!\n");
    }

    return flagSignalData;
}
