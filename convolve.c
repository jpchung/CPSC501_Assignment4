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
#include <time.h>



/*Data structure definition for WAVE file format*/
typedef struct wavHeader {
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
} wavHeader;


/*Instance variable declarations */

//arrays for holding data from input, IR and output WAV files
short* inputWAVdata;
double* xDouble;
int input_size;
short* impulseWAVdata;
double* hDouble;
int impulse_size;
short* outputWAVdata;
double* yDouble;
int output_size;




/*  Function declarations  */
int check_fileExtension(char* fileName);
int readWAV(char* fileName, char* signalType);
int writeWAV(char* fileName);
void convolve(short x[], int N, short h[], int M, short y[], int P);
void print_vector(char *title, float x[], int N);




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
            readWAV(inputFileName, "input");
            //printf("\n");
            readWAV(impulseFileName, "impulse");

            //print_vector("Output signal using identity IR", output_signal, output_size);
            
            
            //int input_size = sizeof(inputWAVdata);
            //int impulse_size = sizeof(impulseWAVdata);
            printf("input size: %d, impulse size: %d\n", input_size, impulse_size);

          
            output_size = input_size + impulse_size - 1;
            printf("output size: %d\n", output_size);

            outputWAVdata = (short*) malloc(sizeof(short) * output_size);
            
            time_t timeStart,timeEnd;
            time(&timeStart);
            convolve(inputWAVdata, input_size, impulseWAVdata, impulse_size, outputWAVdata, output_size);
            time(&timeEnd);
            double elapsed = difftime(timeEnd, timeStart);
            printf("DONE convolution in %.2f seconds!\n\n", elapsed);


            writeWAV(outputFileName);
            


        }
        
        
    }




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

void convolve(short x[], int N, short h[], int M, short y[], int P)
{
    int n, m;

    
    /*  Make sure the output buffer is the right size: P = N + M - 1  */
    if (P != (N + M - 1)) {
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;
    }


    xDouble = malloc(sizeof(double) * N);
    for(int i = 0; i < N; i++){
        xDouble[i] = (double) x[i]/32768.0;
    }

    hDouble =  malloc(sizeof(double) * M);
    for (int i = 0; i < M; i++){
        hDouble[i] = (double) h[i]/32768.0;
    }

    yDouble = malloc(sizeof(double) * P);


    printf("converted to float?\n");

    //printf("Made it here\n");    

    /*  Clear the output buffer y[] to all zero values  */  
    for (n = 0; n < P; n++)
        y[n] = 0;
        
    for(n = 0; n < P; n++)
        yDouble[n] = 0.0;

    printf("about to start convolution loops...\n");
    time_t start = time(NULL);
    /*  Do the convolution  */
    /*  Outer loop:  process each input value x[n] in turn  */
    for (n = 0; n < N; n++) {
        /*  Inner loop:  process x[n] with each sample of h[]  */
        for (m = 0; m < M; m++){
            yDouble[n+m] = xDouble[n] * hDouble[m];
            //y[n+m] += x[n] * h[m];
            y[n+m] = (short) (yDouble[n+m] * 32768); 

            time_t now = time(NULL);
            time_t diff =  now - start;
    
            if((m == 100000) && ((n%200)==0) && ((diff%30) == 0))
                printf("Convolving %d...\n", (n+m));
        }
        
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


int readWAV(char* fileName, char* signalType){
    int fileFlag = 1;

    //open file
    FILE* fp = fopen(fileName, "rb");
    if(fp != NULL){

        struct wavHeader wav;

        //read WAV file by fields
        //fread(storage pointer, size of bytes, n times to read, filepointer)
        printf("Reading file: %s...\n", fileName); 
        
        //RIFF chunk descriptor
        fread(wav.chunkID,4,1,fp);
        fread(&(wav.chunkSize),4,1,fp);
        fread(wav.format,4,1,fp);

        printf("chunk size:  %d\n", wav.chunkSize);
        
        //fmt subchunk
        fread(wav.subChunk1_ID,4,1,fp);
        fread(&(wav.subChunk1_Size),4,1,fp); 
        fread(&(wav.audioFormat),2,1,fp); 
        fread(&(wav.numChannels),2,1,fp);
        fread(&(wav.sampleRate),4,1,fp); 
        fread(&(wav.byteRate),4,1,fp); 
        fread(&(wav.blockAlign),2,1,fp); 
        fread(&(wav.bitsPerSample),2,1,fp); 

        printf("subchunk1 size: %d\n", wav.subChunk1_Size);
        printf("audio format: %d\n", wav.audioFormat);        
        printf("num channels: %d\n", wav.numChannels);
        printf("sample rate: %d\n", wav.sampleRate);

        //check size of data subchunk, read extra bytes if necessary
        if(wav.subChunk1_Size == 18){
            short emptyBytes;
            fread(&emptyBytes,1,2, fp);
        }

        //data subchunk
        fread(wav.subChunk2_ID,4,1,fp);
        fread(&(wav.subChunk2_Size),4,1,fp);
        printf("subchunk2 size: %d\n", wav.subChunk2_Size);

                
        //calculate bytes per sample, number of samples         
        int wavBytesPerSample = (wav.bitsPerSample)/8;
        printf("bytes per sample: %d\n", wavBytesPerSample);
        int wavNumberOfSamples = (wav.subChunk2_Size)/wavBytesPerSample; //since subchunk2size = sample data size in bytes
        printf("number of samples: %d\n", wavNumberOfSamples);

        //allocate memory for input sample data
        if(signalType == "input"){
            //input_size = wav.subChunk2_Size;
            input_size = wavNumberOfSamples;
            
            printf("%s\n", signalType);
            inputWAVdata = (short*) malloc(sizeof(short) * wavNumberOfSamples);

            //read the sound data a sample at a time, checking for expected bytes per sample
            short sampleData = 0;
            int i = 0;
            while(fread(&sampleData,1, wavBytesPerSample,fp) == wavBytesPerSample){
                inputWAVdata[i++] = sampleData;
                sampleData = 0; //clear for next sample
                //if((i % 100000) == 0)
                //    printf("Reading data sample: %d...\n", i);
            }
            printf("last read sample: %d\n", i-1);     
            
        }
        else if(signalType == "impulse"){
            //impulse_size = wav.subChunk2_Size;
            impulse_size = wavNumberOfSamples;
            
            printf("%s\n", signalType);
            impulseWAVdata = (short*) malloc(sizeof(short) * wavNumberOfSamples);
            
            //read the sound data a sample at a time, checking for expected bytes per sample
            short sampleData = 0;
            int i = 0;
            while(fread(&sampleData,1, wavBytesPerSample,fp) == wavBytesPerSample){
                impulseWAVdata[i++] = sampleData;
                sampleData = 0; //clear for next sample
                //if((i % 100000) == 0)
                //    printf("Reading data sample: %d...\n", i);
            }
            printf("last read sample: %d\n", i-1); 

        }

         
        //close file
        fclose(fp);
        printf("Finished reading %s!\n\n", fileName);

        fileFlag = 0;
    }
    else{
        printf("Error opening file %s\n", fileName);
    }

    return fileFlag;
}

int writeWAV(char* fileName){
    int fileFlag = 1;

    FILE* fp = fopen(fileName, "wb");
    if(fp != NULL){

        //struct wavHeader wavOut;

        //write WAV file by fields
        //fwrite(storage pointer, size of elements in bytes, number of elements, filepointer)
        printf("Writing file: %s...\n", fileName); 

        //RIFF chunk descriptor
        fwrite("RIFF",4,1,fp);
        fwrite(&(wavOut.chunkSize),4,1,fp);
        fwrite(wavOut.format,4,1,fp);

        printf("file write check 1\n");
        //fmt subchunk
        fwrite(wavOut.subChunk1_ID,4,1,fp);
        fwrite(&(wavOut.subChunk1_Size),4,1,fp);
        fwrite(&(wavOut.audioFormat),2,1,fp);
        fwrite(&(wavOut.numChannels),2,1,fp);
        fwrite(&(wavOut.sampleRate),4,1,fp);
        fwrite(&(wavOut.byteRate),4,1,fp);
        fwrite(&(wavOut.blockAlign),2,1,fp);
        fwrite(&(wavOut.bitsPerSample),2,1,fp);

        printf("file write check 2\n");
        if(wavOut.subChunk1_Size == 18){
            short emptyBytes;
            fwrite(&emptyBytes,1,2,fp);
        }

        //data subchunk
        fwrite(wavOut.subChunk2_ID,4,1,fp);
        fwrite(&(wavOut.subChunk2_Size),4,1,fp);

        printf("file write check 3\n");
        //calculate bytes per sample, number of samples
        int wavOutBytesPerSample = (wavOut.bitsPerSample)/8;
        int wavOutNumberOfSamples = output_size/wavOutBytesPerSample; //since output size is P = N + M - 1
        


        //write the output sound data a sample at a time
        short outSampleData;
        for(int i = 0; i < output_size; i++){
        //for(int i = 0; i < wavOutNumberOfSamples; i++){
            printf("pls work...\n");
            //outSampleData = (short) (outputWAVdata[i] * 32767);
            outSampleData = (short) outputWAVdata[i];
            
            fwrite((char*) &outSampleData, 1, 2, fp);
        }

        fclose(fp);
        printf("Finished writing %s!\n", fileName);


        fileFlag = 0;
    }
    else{
        printf("Error writing file %s\n", fileName);
    }
    
    return fileFlag;

}