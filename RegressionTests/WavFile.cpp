#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "WavFile.h"

using namespace std;

void WavFile::readWAV(char* fileName){

    printf("\nReading file %s...\n", fileName);

    //open file stream in binary mode for input operations
    ifstream wav(fileName, ios::in | ios:: binary);

    //read(char* storageArray, streamsize n)
    //RIFF chunk descriptor
    wav.seekg(0, ios::beg);
    wav.read(chunkID, 4);
    wav.read((char*) &chunkSize, 4);
    wav.read(format, 4);
    //fmt subchunk
    wav.read(subChunk1ID, 4);
    wav.read((char*) &subChunk1Size, 4);
    wav.read((char*) &audioFormat, sizeof(short));
    wav.read((char*) &numChannels, sizeof(short));
    wav.read((char*) &sampleRate, sizeof(int));
    wav.read((char*) &byteRate, sizeof(int));
    wav.read((char*) &blockAlign, sizeof(short));
    wav.read((char*) &bitsPerSample, sizeof(short));

    
    printf("ChunkID: %s\n", chunkID);
    printf("ChunkSize: %d\n", chunkSize);
    printf("Format: %s\n", format);

    printf("SubChunk1ID: %s\n", subChunk1ID);
    printf("SubChunk1Size: %d\n", subChunk1Size);
    printf("AudioFormat: %d\n", audioFormat);
    printf("NumChannels: %d\n", numChannels);
    printf("SampleRate: %d\n", sampleRate);
    printf("ByteRate: %d\n", byteRate);
    printf("BitsPerSample: %d\n", bitsPerSample);

    //read extra bytes if necessary
    if(subChunk1Size == 18){
        short emptyBytes;
        wav.read((char*) &emptyBytes, sizeof(short));
    }

    //data subchunk
    wav.read(subChunk2ID, 4);
    wav.read((char*) &subChunk2Size, sizeof(int));

    printf("SubChunk2ID: %s\n", subChunk2ID);
    printf("SubChunk2Size: %d\n", subChunk2Size);

    int dataSize = subChunk2Size; //subchunk2size = size of data following

    dataArray = new char[dataSize]; 
    wav.read(dataArray, dataSize);

    //close file stream
    wav.close();

    signal = NULL;

    //convert signal based on bits per sample
    if(bitsPerSample == 8){
        signalSize = dataSize;
        signal = new short[signalSize];
        for(int i = 0; i < dataSize; i++){
            signal[i] = (short) ((unsigned char) dataArray[i]);
        }
    }
    else {
        signalSize = dataSize/2;
        signal = new short[signalSize];
        //every 2 chars (i.e. a short) is one signal sample
        short shortData;
        for(int i = 0; i < dataSize; i+= 2){
            shortData = (short) ((unsigned char) dataArray[i]);
            shortData = (short) ((unsigned char) dataArray[i+1]) *256; //shift to next 8 bits
            signal[i/2] = shortData;
        }
    }

    printf("Done reading %s!\n", fileName);

    
}