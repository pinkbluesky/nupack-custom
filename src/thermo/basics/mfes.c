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

/* ************************************************ */

int main(int argc, char *argv[])
{

  char seq[MAXSEQLENGTH];
  int seqNum[MAXSEQLENGTH + 1];

  int complexity = 3;
  int tmpLength;
  DBL_TYPE mfe;
  int vs;
  char inputFile[MAXLINE];
  char outFile[MAXLINE];
  int inputFileSpecified;
  FILE *fp;

  dnaStructures mfeStructs = {NULL, 0, 0, 0, NAD_INFINITY};

  // Get the command line arguments
  strcpy(inputFile, "");

  inputFileSpecified = ReadCommandLineNPK(argc, argv, inputFile);
  if (NupackShowHelp)
  {
    printf("Usage: mfes [OPTIONS] PREFIX\n");
    printf("Compute and store the minimum free energy and the MFE\n");
    printf("secondary structure(s) of the input sequence(s).\n");
    printf("Example: mfes -T 25 -material dna example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
    exit(1);
  }

  if (!inputFileSpecified)
  {
    printf("Enter output file prefix: ");
    scanf("%s", inputFile);
    strcat(inputFile, ".in"); // Here, .in is just a placeholder
  }

  printf("line 63 mfes.c\n");

  // Read the input file
  if (!inputFileSpecified ||
      !ReadInputFileADSCustom(inputFile, seq, &vs, NULL, NULL, NULL))
  {
    if (inputFileSpecified == 0)
      getUserInput(seq, &vs, NULL, NULL);
    else
      abort();
  }
  strncpy(outFile, inputFile, strlen(inputFile) - 3);
  outFile[strlen(inputFile) - 3] = '\0';
  strcat(outFile, ".mfes");

  // Populate output file with basic information
  header(argc, argv, "mfes", outFile);
  printInputs(argc, argv, seq, vs, NULL, NULL, outFile);

  if (!DO_PSEUDOKNOTS)
  {
    complexity = 3;
  }
  else
  {
    complexity = 5;
  }

  fp = fopen(outFile, "a");

  // Iterate through each strand
  printf("\n\nfull seq: %s\n\n", seq);
  char *token = strtok(seq, "+");

  while (token)
  {
    printf("token: %s\n", token);
    // printf("token is alpha: %d\n\n", isalpha(token));
    if (isalpha(token[0]) == 0) // if strand is not alpha, exit
    {
      printf("Error in input file, sequence is not alpha.");
      return 1;
    }

    // get MFE of sequence
    tmpLength = strlen(token);
    convertSeq(token, seqNum, tmpLength);

    //int seqNumCopy[MAXSEQLENGTH + 1];
    //strncpy(seqNumCopy, seqNum, tmpLength); 


    mfe = mfeFullWithSym(seqNum, tmpLength, &mfeStructs, complexity, DNARNACOUNT,
                         DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, vs,
                         1, SODIUM_CONC, MAGNESIUM_CONC,
                         USE_LONG_HELIX_FOR_SALT_CORRECTION);


    // Get next strand
    token = strtok(NULL, "+");

    //mfe = 0;
    printf("[test] mfe of token: %Lf\n", mfe);

    // write to output file
    fprintf(fp, "%.3Lf\n", mfe);

    // Reset vars
    clearDnaStructures(&mfeStructs);
  }

  fclose(fp);

  return 0;
}
/* ****** */
