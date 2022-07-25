/*
    mfe.c is part of the NUPACK software suite
    Copyright (c) 2007 Caltech. All rights reserved.
    Coded by: Robert Dirks, 6/2006 and Justin Bois 1/2007


    This function will calculate and print all mfe structures (if the
    -degenerate flag is selected) or one mfe structure (if not),
    taking into account symmetry corrections.  Consequently, if
    -degenerate is chosen, this could scale as poorly as exponential
    with regard to space and time.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include <thermo/core.h>
#include <shared.h>

typedef struct ThreadInfo {
	pthread_t thread;
	char **seqs;
	int seqNum;
	DBL_TYPE mfeSum;
} ThreadInfo;

int threadNum = 2;
int complexity = 3;

void *mfeproc(void *a);
 
/* ************************************************ */

int main(int argc, char *argv[])
{

  char **seqs;
//  int seqNum[MAXLINE];
  int nStrands;
  int strandLen;
  DBL_TYPE mfeSum = 0;

//  int tmpLength;
//  DBL_TYPE mfe;
  int vs;
  char inputFile[MAXLINE];
  char outFile[MAXLINE];
  int inputFileSpecified;
  FILE *fp;
  int seqPerThread;
  ThreadInfo *tinfo;

//  dnaStructures mfeStructs = {NULL, 0, 0, 0, NAD_INFINITY};

  // Get the command line arguments
  strcpy(inputFile, "");

  inputFileSpecified = ReadCommandLineNPK(argc, argv, inputFile);
  if (NupackShowHelp)
  {
    printf("Usage: mfes [OPTIONS] PREFIX\n");
    printf("Compute and store the minimum free energy and the MFE\n");
    printf("secondary structure(s) of the input sequence(s).\n");
    printf("Example: mfes -T 25 -material dna -multi example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
    exit(1);
  }

  if (!inputFileSpecified)
  {
    printf("No input file specified, aborting.");
    abort();
  }


  // Create output file path
  strncpy(outFile, inputFile, strlen(inputFile) - 3);
  outFile[strlen(inputFile) - 3] = '\0';
  strcat(outFile, ".mfes");

  // Populate output file with basic information
  header(argc, argv, "mfes", outFile);
  printInputs(argc, argv, NULL, 1, NULL, NULL, outFile);

  if (!DO_PSEUDOKNOTS)
  {
    complexity = 3;
  }
  else
  {
    complexity = 5;
  }

  // Read file
  if (ReadInputFileSizes(inputFile, &nStrands, &strandLen) == 0)
  {
    printf("Unable to read input file sizes, aborting.");
    abort();
  }

  // Allocate memory
  seqs = (char **)malloc(nStrands * sizeof(char *));
  for (int i = 0; i < nStrands; i++)
  {
    seqs[i] = (char *)malloc((strandLen + 1) * sizeof(char));
  }

//  printf("Allocated memory for sequences.");

  // Read in sequences
  if (ReadInputFileADSCustom(inputFile, seqs, nStrands, strandLen, &vs, NULL, NULL, NULL) == 0)
  {
    printf("Unable to read input file sequences, aborting.");
    abort();
  }

  seqPerThread = nStrands / threadNum;
  if (nStrands%threadNum != 0)
	seqPerThread++;

  // Get mfe for each sequence
  tinfo = calloc(threadNum, sizeof(ThreadInfo));
  for(int i = 0; i < threadNum; i++) {
	ThreadInfo *ti = &tinfo[i];
	int b, e;

	b = i*seqPerThread;
	e = b + seqPerThread;
	if (b >= nStrands) {
		// this can happen if we have too few strands, less than number of threads
		threadNum = i;
		break;
	}

	if (e > nStrands)
		e = nStrands;

	ti->seqs = &seqs[b];
	ti->seqNum = e - b;
	pthread_create(&ti->thread, NULL, mfeproc, ti);
  }

  for(int i = 0; i < threadNum; i++) {
	ThreadInfo *ti = &tinfo[i];
	pthread_join(ti->thread, NULL);
	mfeSum += ti->mfeSum;
  }

/*
  for (int i = 0; i < nStrands; i++)
  {
    
    // get MFE of sequence
    tmpLength = strlen(seqs[i]);
    convertSeq(seqs[i], seqNum, tmpLength);

    mfe = mfeFullWithSym(seqNum, tmpLength, &mfeStructs, complexity, DNARNACOUNT,
                         DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, vs,
                         1, SODIUM_CONC, MAGNESIUM_CONC,
                         USE_LONG_HELIX_FOR_SALT_CORRECTION);

    printf("%Lf\n", mfe);

    clearDnaStructures(&mfeStructs);

    mfeSum += mfe;
  }
*/

  // Write to output file
  fp = fopen(outFile, "a");

  fprintf(fp, "%.3Lf\n", mfeSum);
  printf("final mfe sum is %Lf\n", mfeSum);

  fclose(fp);

  return 0;
}

/* ****** */

void *mfeproc(void *a) {
	int slen;
	ThreadInfo *ti = a;
	int seqNum[MAXLINE];
	DBL_TYPE mfe;
	dnaStructures mfeStructs = {NULL, 0, 0, 0, NAD_INFINITY};

	for(int i = 0; i < ti->seqNum; i++) {
		slen = strlen(ti->seqs[i]);
		convertSeq(ti->seqs[i], seqNum, slen);

		printf("%s ", ti->seqs[i]);
		mfe = mfeFullWithSym(seqNum, slen, &mfeStructs, complexity, DNARNACOUNT,
                         DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, 0,
                         1, SODIUM_CONC, MAGNESIUM_CONC,
                         USE_LONG_HELIX_FOR_SALT_CORRECTION);

		printf("%Lf\n", mfe);
		clearDnaStructures(&mfeStructs);
		ti->mfeSum += mfe;
	}

	return NULL;
}
