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

#include <thermo/core.h>
#include <shared.h>

/* ************************************************ */

int main(int argc, char *argv[])
{

  // TODO find a better way to reading in a large file  to the sequence string
  // char seq[MAXSEQLENGTH_ADS];
  // int seqNum[MAXSEQLENGTH_ADS + 1]; // last index contains -1 to indicate end of sequences

  // DBL_TYPE strandMfes[MAXSTRANDS_ADS + 1]; // contains mfe values for each strand -- output array

  printf("in mfes init vars");
  char **seqs;
  int seqNum[MAXLINE];
  int nStrands;
  int strandLen;
  DBL_TYPE mfeSum = 0;

  int complexity = 3;
  int tmpLength;
  DBL_TYPE mfe;
  int vs;
  char inputFile[MAXLINE];
  char outFile[MAXLINE];
  int inputFileSpecified;
  FILE *fp;

  dnaStructures mfeStructs = {NULL, 0, 0, 0, NAD_INFINITY};

  printf("in mfes");

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

  /*
    if (!inputFileSpecified)
    {
      printf("Enter output file prefix: ");
      scanf("%s", inputFile);
      strcat(inputFile, ".in"); // Here, .in is just a placeholder
    }

    // printf("line 63 mfes.c\n");

    // Read the input file
    if (!inputFileSpecified ||
        !ReadInputFileADSCustom(inputFile, seq, &vs, NULL, NULL, NULL))
    {
      if (inputFileSpecified == 0)
        getUserInput(seq, &vs, NULL, NULL);
      else
        abort();
    }
  */

  // Create output file path
  strncpy(outFile, inputFile, strlen(inputFile) - 3);
  outFile[strlen(inputFile) - 3] = '\0';
  strcat(outFile, ".mfes");

  // Populate output file with basic information
  header(argc, argv, "mfes", outFile);
  printInputs(argc, argv, NULL, vs, NULL, NULL, outFile);

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

  printf("Allocated memory for sequences.");

  // Read in sequences
  if (ReadInputFileADSCustom(inputFile, seqs, nStrands, strandLen, &vs, NULL, NULL, NULL) == 0)
  {
    printf("Unable to read input file sequences, aborting.");
    abort();
  }

  // Get mfe for each sequence
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

  // Write to output file
  fp = fopen(outFile, "a");

  fprintf(fp, "%.3Lf\n", mfeSum);

  fclose(fp);

  return 0;
}

/* ****** */
