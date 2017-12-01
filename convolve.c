/**
 * CPSC 501 Assignment 4 - Optimizing Program Performance
 * Baseline program:
 * Time Domain Convolution of .wav audio files via Input Side Algorithm.
 * Program takes command-line for names of input recording and impulse response files, 
 * and produces a convolved signal
 * 
 * Created by: Johnny Chung -  Modified from demo code by Leonard Manzara
 * /


/*  Include files  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/*Data structure definition for WAVE file format*/
typedef struct wavFile {
    char chunkID[4];
    int chunkSize;
    char format[4];
    char subChunk1_ID[4];
    int subChunk1_Size;
    short audioFormat;
    short numChannels;
    int sampleRate;
    int byteRate;
    short blockAlign;
    short bitsPerSample;
    char subChunk2_ID[4];
    int subChunk2_Size;
} wavFile;


/*Instance variable declarations */

//array for holding data from input WAV file
short* inputWAVdata;

//RIFF Chunk Descriptor (IRfile.wav)
char impulseChunkID[4];
int impulseChunkSize;
char impulseFormat[4];
//fmt subchunk (IRfile.wav)
char impulseSubChunk1ID[4];
int impulseSubChunk1Size;
short impulseAudioFormat;
short impulseNumChannels;
int impulseSampleRate;
int impulseByteRate;
short impulseBlockAlign;
short impulseBitsPerSample;
//data subchunk (IRfile.wav)
char impulseSubChunk2ID[4];
int impulseSubChunk2Size;

//RIFF Chunk Descriptor (outputfile.wav)
char outputChunkID[4];
int outputChunkSize;
char outputFormat[4];
//fmt subchunk (outputfile.wav)
char outputSubChunk1ID[4];
int outputSubChunk1Size;
short outputAudioFormat;
short outputNumChannels;
int outputSampleRate;
int outputByteRate;
short outputBlockAlign;
short outputBitsPerSample;
//data subchunk (outputfile.wav)
char outputSubChunk2ID[4];
int outputSubChunk2Size;

/*  Function declarations  */
void convolve(float x[], int N, float h[], int M, float y[], int P);
void print_vector(char *title, float x[], int N);
int check_fileExtension(char* fileName);
int readInputWAV(char* fileName);


/*****************************************************************************
*
*    Function:     main
*
*    Description:  Tests the convolve function with various input signals
*
*****************************************************************************/

int main(int argc, char *argv[])
{

    char* inputFileName;
    char* impulseFileName;
    char* outputFileName;

    //check number of command-line arguments
    if(argc <= 3){
        printf("Usage: ./convolve <inputfile.wav> <IRfile.wav> <outputfile.wav>\n");
        return EXIT_FAILURE;
    }
    else {
        printf("Correct number of arguments...\n");


        inputFileName = argv[1];
        impulseFileName = argv[2];
        outputFileName = argv[3];


        int inputExtension = check_fileExtension(inputFileName);
        int impulseExtension = check_fileExtension(impulseFileName);
        int outputExtension = check_fileExtension(outputFileName);
        
        int fileExtensionFlag = 0;        
        int fileExtensions[3] = {inputExtension, impulseExtension, outputExtension};
        for(int i = 0; i < sizeof(fileExtensions); i++){
            if(fileExtensions[i] == EXIT_FAILURE){
                fileExtensionFlag = 1;
            }
        }

        if(fileExtensionFlag == 1){
            printf("Usage: ./convolve <inputfile.wav> <IRfile.wav> <outputfile.wav>\n");       
            printf("returning %d\n", EXIT_FAILURE);
            return EXIT_FAILURE;
        }
        else {
            //read input wav file
            readInputWAV(inputFileName);
        }
        
        
    }


    float input_signal[100], impulse_response[20], output_signal[120];
    int input_size, impulse_size, output_size;

    /*  Create an example input signal  */
    input_signal[0] = 1.0;
    input_signal[1] = 0.5;
    input_signal[2] = 0.25;
    input_signal[3] = 0.125;
    input_size = 4;
        
    /*  Print out the input signal to the screen  */
    //print_vector("Original input signal", input_signal, input_size);

    /*  Create an "identity" impulse response.  The output should be
        the same as the input when convolved with this  */
    impulse_response[0] = 1.0;
    impulse_size = 1;

    /*  Set the expected size of the output signal  */
    output_size = input_size + impulse_size - 1;

    /*  Do the convolution, and print the output signal  */
    convolve(input_signal, input_size, impulse_response, impulse_size,
        output_signal, output_size);
    //print_vector("Output signal using identity IR", output_signal, output_size);

  
    /*  End of program  */
    return 0;
}



/*****************************************************************************
*
*    Function:     convolve
*
*    Description:  Convolves two signals, producing an output signal.
*                  The convolution is done in the time domain using the
*                  "Input Side Algorithm" (see Smith, p. 112-115).
*
*    Parameters:   x[] is the signal to be convolved
*                  N is the number of samples in the vector x[]
*                  h[] is the impulse response, which is convolved with x[]
*                  M is the number of samples in the vector h[]
*                  y[] is the output signal, the result of the convolution
*                  P is the number of samples in the vector y[].  P must
*                       equal N + M - 1
*
*****************************************************************************/

void convolve(float x[], int N, float h[], int M, float y[], int P)
{
    int n, m;

    /*  Make sure the output buffer is the right size: P = N + M - 1  */
    if (P != (N + M - 1)) {
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;
    }

    /*  Clear the output buffer y[] to all zero values  */  
    for (n = 0; n < P; n++)
        y[n] = 0.0;

    /*  Do the convolution  */
    /*  Outer loop:  process each input value x[n] in turn  */
    for (n = 0; n < N; n++) {
        /*  Inner loop:  process x[n] with each sample of h[]  */
        for (m = 0; m < M; m++)
        y[n+m] += x[n] * h[m];
    }
}



/*****************************************************************************
*
*    Function:     print_vector
*
*    Description:  Prints the vector out to the screen
*
*    Parameters:   title is a string naming the vector
*                  x[] is the vector to be printed out
*                  N is the number of samples in the vector x[]
*
*****************************************************************************/

void print_vector(char *title, float x[], int N)
{
  int i;

  printf("\n%s\n", title);
  printf("Vector size:  %-d\n", N);
  printf("Sample Number \tSample Value\n");
  for (i = 0; i < N; i++)
    printf("%-d\t\t%f\n", i, x[i]);
}

/**
 * Function: check_fileExtension
 * Description: checks for .wav file extension in file names
 * Parameter(s): fileName is string for name of file
*/
int check_fileExtension(char* fileName){

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


int readInputWAV(char* fileName){
    int fileFlag = 1;

    //open file
    FILE* fp = fopen(fileName, "rb");
    if(fp != NULL){

        struct wavFile input;

        //read WAV file by fields
        //fread(storage pointer, size of bytes, n times to read, filepointer)
        printf("Reading file: %s...\n", fileName); 
        
        //RIFF chunk descriptor
        fread(input.chunkID,4,1,fp);
        fread(&(input.chunkSize),4,1,fp);
        printf("input chunk size:  %d\n", input.chunkSize);
        fread(input.format,4,1,fp);
        //fmt subchunk
        fread(input.subChunk1_ID,4,1,fp);
        fread(&(input.subChunk1_Size),4,1,fp); 
        fread(&(input.audioFormat),2,1,fp); 
        fread(&(input.numChannels),2,1,fp);
        fread(&(input.sampleRate),4,1,fp); 
        fread(&(input.byteRate),4,1,fp); 
        fread(&(input.blockAlign),2,1,fp); 
        fread(&(input.bitsPerSample),2,1,fp); 

        printf("subchunk 1 size: %d\n", input.subChunk1_Size);
        printf("audio format: %d\n", input.audioFormat);        
        printf("num channels: %d\n", input.numChannels);
        printf("sample rate: %d\n", input.sampleRate);

        //check size of data subchunk, read extra bytes if necessary
        if(input.subChunk1_Size == 18){
            short emptyBytes;
            fread(&emptyBytes,1,2, fp);
        }

        //data subchunk
        fread(input.subChunk2_ID,4,1,fp);
        fread(&(input.subChunk2_Size),4,1,fp);
        printf("input subchunk2 size: %d\n", input.subChunk2_Size);

        //calculate bytes per sample, number of samples         
        int inputBytesPerSample = (input.bitsPerSample)/8;
        printf("bytes per sample: %d\n", inputBytesPerSample);
        int inputNumberOfSamples = (input.subChunk2_Size)/inputBytesPerSample; //since subchunk2size = sample data size in bytes
        printf("number of samples: %d\n", inputNumberOfSamples);

        //allocate memory for input sample data
        inputWAVdata = (short*) malloc(sizeof(short) * inputNumberOfSamples);
        
        //read the sound data a sample at a time, checking for expected bytes per sample
        short sampleData = 0;
        int i = 0;
        while(fread(&sampleData,1, inputBytesPerSample,fp) == inputBytesPerSample){
            inputWAVdata[i++] = sampleData;
            sampleData = 0; //clear for next sample
            if((i % 100000) == 0)
                printf("Reading data sample: %d...\n", i);
        }
        printf("last read sample: %d\n", i-1);        
        
        //close file
        fclose(fp);
        printf("Finished reading %s!\n", fileName);

        fileFlag = 0;
    }
    else{
        printf("Error opening file %s\n", fileName);
    }

    return fileFlag;
}