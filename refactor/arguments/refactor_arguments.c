/*
 * Definition file for refactored version of ONeSAMP that handles command line
 * arguments.
 */

// HEADERS

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "refactor_arguments.h"

// INCLUDE GUARD

#ifndef REFACTOR_ARGUMENTS_C
#define REFACTOR_ARGUMENTS_C

// VARIABLES

// Stored variables from the command line
int randomFlag;
int num_loci;
int input_individuals_count;
int *bottleneck_individuals_count = NULL;
int *bottleneck_length = NULL;
int final_individuals_count;
int isMicrosats;
double minAlleleFrequency;
double omitThreshold;
double *theta = NULL;
double *mutation_rate = NULL;
int repetitions;
char *programName = NULL;
int syntax_check;
int example;
int raw_stats;
int single_generation;
int example_pop;
int absentDataExtrapolate;
int input_individuals_count_allocation = -1;
int num_loci_allocation = -1;
int *motif_lengths = NULL;

// Stored arrays from the command line
int *bottleneck_individuals_count_random_choices = NULL;
int *bottleneck_length_random_choices = NULL;
double *theta_random_choices = NULL;
double *mutation_rate_random_choices = NULL;

// Proportion of missing data
double proportionMissingData;

// FUNCTIONS

/*! \brief Reports an error generated by the parser from OneSamp.
 */
void reportParseError(char *message, long row, long column){
  fprintf(stderr, "ONESAMP PARSE ERROR, line %ld, column %ld \n", row, column);
  fprintf(stderr, message, parseProgramName());
  fprintf(stderr, "\nExiting...\n");
  fflush(stderr);
  exit(1);
}

/*! \brief Reports an error (in general) in OneSamp.
 */
void reportError(char *message){
  fprintf(stderr, "ONESAMP ERROR\n");
  fprintf(stderr, message, parseProgramName());
  fprintf(stderr, "\nExiting...\n");
  fflush(stderr);
  exit(1);
}

/*! \brief Exits the program when an error with the arguments is encountered.
 */
void reportArgumentError(char *message){
  fprintf(stderr, "ONESAMP ERROR\n");
  fprintf(stderr, message, parseProgramName());
  fprintf(stderr, "\nError reading inputs. Please see the README for details.\n");
  fprintf(stderr, "\nExiting...\n");
  fflush(stderr);
  exit(1);
}

/*! \brief Exits the program when an error with theta value is encountered.
 */
void reportThetaArgumentError(char *message, double thetaMin, double thetaMax){
  fprintf(stderr, "ONESAMP ERROR\n");
  fprintf(stderr, message, parseProgramName(), thetaMin, thetaMax);
  fprintf(stderr, "\nError reading inputs. Please see the README for details.\n");
  fprintf(stderr, "\nExiting...\n");
  fflush(stderr);
  exit(1);
}

/*! \brief Resets memory that stores arguments.
 */
void resetArguments(){
  randomFlag = -1; // True for C random, False for GFSR
  num_loci = -1;
  input_individuals_count = -1;
  final_individuals_count = -1;
  isMicrosats = -1;
  repetitions = -1;
  minAlleleFrequency = -1;
  omitThreshold = -1;
  syntax_check = FALSE;
  example = FALSE;
  raw_stats = FALSE;
  single_generation = FALSE;
  absentDataExtrapolate = FALSE;
  // Program Name
  if(programName != NULL) free(programName);
  programName = NULL;
  // Population bottleneck length
  if(bottleneck_length != NULL) free(bottleneck_length);
  bottleneck_length = NULL;
  if(bottleneck_length_random_choices) free(bottleneck_length_random_choices);
  bottleneck_length_random_choices = NULL;
  // Population bottleneck size
  if(bottleneck_individuals_count != NULL) free(bottleneck_individuals_count);
  bottleneck_individuals_count = NULL;
  if(bottleneck_individuals_count_random_choices != NULL) free(bottleneck_individuals_count_random_choices);
  bottleneck_individuals_count_random_choices = NULL;
  // Population mutation rate
  if(mutation_rate != NULL) free(mutation_rate);
  mutation_rate = NULL;
  if(mutation_rate_random_choices != NULL) free(mutation_rate_random_choices);
  mutation_rate_random_choices = NULL;
  // Evolutionary parameter theta
  if(theta != NULL) free(theta);
  theta = NULL;
  if(theta_random_choices != NULL) free(theta_random_choices);
  theta_random_choices = NULL;
  if(motif_lengths != NULL) free(motif_lengths);
}

/*! \brief Allocates memory for simulated bottlenecks.
 */
void allocateStruct7(int **bottleneck_individuals_count_random_choices_ptr){
  int j;
  *bottleneck_individuals_count_random_choices_ptr = (int *)malloc(parseIterations() * sizeof(int *));
  for(j = 0; j < parseIterations(); j++) (*bottleneck_individuals_count_random_choices_ptr)[j] = randomQuantizedIntervalSelection(parseBottleneckMin(), parseBottleneckMax(), 1);
}

/*! \brief Allocates memory for simulated bottleneck lengths.
 */
void allocateStruct8(int **bottleneck_length_random_choices_ptr){
  int j;
  *bottleneck_length_random_choices_ptr = (int *)malloc(parseIterations() * sizeof(int *));
  for(j = 0; j < parseIterations(); j++) (*bottleneck_length_random_choices_ptr)[j] = randomQuantizedIntervalSelection(parseBottleneckLengthMin(), parseBottleneckLengthMax(), 1);
}

/*! \brief Allocates memory for simulated theta values.
 */
void allocateStruct9(double **theta_random_choices_ptr){
  int j;
  *theta_random_choices_ptr = (double *)malloc(parseIterations() * sizeof(double *));
  for(j = 0; j < parseIterations(); j++) (*theta_random_choices_ptr)[j] = randomQuantizedIntervalSelection(parseThetaMin(), parseThetaMax(), 0.00000001);
}

/*! \brief Allocates memory for simulated mutation rates.
 */
void allocateStruct10(double **mutation_rate_random_choices_ptr){
  int j;
  *mutation_rate_random_choices_ptr = (double *)malloc(parseIterations() * sizeof(double *));
  for(j = 0; j < parseIterations(); j++) (*mutation_rate_random_choices_ptr)[j] = randomQuantizedIntervalSelection(parseMRateMin(), parseMRateMax(), 0.00000001);
}

/*! \brief Returns the name of the executable.
 */
char *parseProgramName(){
  return programName;
}

/* \brief Returns whether or not to interpolate missing data.
 */
int parseFillInAbsentData(){
  return absentDataExtrapolate;
}

/* \brief Returns the minimum allele frequency of mutated alleles in populations
 */
double parseMinAlleleFrequency(){
  if(minAlleleFrequency < 0 || minAlleleFrequency > 0.5) reportArgumentError((char *) "%s: argument -f, minimum proportion of mutated alleles, is either missing or not a floating point number between 0 and 0.5");
  return minAlleleFrequency;
}

/* \brief Returns the threshold by which loci should be omitted.
 */
double parseOmitLocusThreshold(){
  if(omitThreshold == -1) reportArgumentError((char *) "%s: argument -o, minimum proportion of individuals with completely specified genotypes for loci to be included in computation, is not specified");
  return omitThreshold;
}

/*! \brief Returns whether the input data is composed of microsatellites or SNPs.
 */
int parseFormFlag(){
  if(isMicrosats == -1) reportArgumentError((char *) "%s: argument -s or -m, SNPs or microsatellites loci, not specified.");
  return isMicrosats;
}

/*! \brief Returns the number of trials of bottleneck generations to simulate.
 */
int parseIterations(){
  if(repetitions <= 0) reportArgumentError((char *) "%s: argument -t, number of repetitions, must be a positive integer.");
  return repetitions;
}

/*! \brief Returns the number of polymorphic loci
 */
int parseNLoci(){
  if(num_loci <= 0) reportArgumentError((char *) "%s: argument -l, num of unlinked polymorphic loci, must be a positive integer.");
  return num_loci;
}

/*! \brief Sets the number of loci after some may have been thrown out.
 */
void setNLoci(int size){
  num_loci = size;
}

/*! \brief Returns the size of the input samples.
 */
int parseNLociAllocation(){
  return num_loci_allocation;
}

/*! \brief Returns the size of the input samples.
 */
int parseInputSamples(){
  if(input_individuals_count <= 0) reportArgumentError((char *) "%s: argument -i, num of input samples, must be a positive integer.");
  return input_individuals_count;
}

/*! \brief Sets the size of the allocation.
 */
void setInputSamples(int size){
  input_individuals_count = size;
  final_individuals_count = size;
}

/*! \brief Returns the size of the input samples.
 */
int parseInputSamplesAllocation(){
  return input_individuals_count_allocation;
}

/*! \brief Returns the size of the bottleneck of the current sample.
 */
int parseBottleneck(int samp){
  if(bottleneck_individuals_count == NULL) reportArgumentError((char *) "%s: argument -b, num of individuals in bottleneck generation, must be a positive even integer");
  return bottleneck_individuals_count_random_choices[samp];
}

/*! \brief Returns the minimum possible size of the bottleneck specified by the data from the command line.
 */
int parseBottleneckMin(){
  if(bottleneck_individuals_count == NULL) reportArgumentError((char *) "%s: argument -b, num of individuals in bottleneck generation, must be a positive even integer");
  return bottleneck_individuals_count[0];
}

/*! \brief Returns the maximum possible size of the bottleneck specified by the data from the command line.
 */
int parseBottleneckMax(){
  if(bottleneck_individuals_count == NULL) reportArgumentError((char *) "%s: argument -b, num of individuals in bottleneck generation, must be a positive even integer");
  return bottleneck_individuals_count[1];
}

/*! \brief Returns the value of the duration of the bottleneck of the specified simulation over the sample.
 */
int parseBottleneckLength(int samp){
  if(bottleneck_length == NULL || bottleneck_length[0] < 0 || bottleneck_length[1] < 0) reportArgumentError((char *) "%s: argument -d, duration of bottleneck generations, must be a nonnegative integer");
  return bottleneck_length_random_choices[samp];
}

/*! \brief Returns the minimum duration of the bottleneck specified by the data from the command line.
 */
int parseBottleneckLengthMin(){
  if(bottleneck_length == NULL || bottleneck_length[0] < 0 || bottleneck_length[1] < 0) reportArgumentError((char *) "%s: argument -d, duration of bottleneck generations, must be a nonnegative integer");
  return bottleneck_length[0];
}

/*! \brief Returns the maximum duration of the bottleneck specified by the data from the command line.
 */
int parseBottleneckLengthMax(){
  if(bottleneck_length == NULL || bottleneck_length[0] < 0 || bottleneck_length[1] < 0) reportArgumentError((char *) "%s: argument -d, duration of bottleneck generations, must be a nonnegative integer");
  return bottleneck_length[1];
}

/*! \brief Returns the value of theta of the specified sample.
 */
double parseTheta(int samp){
  if(theta == NULL || theta[0] <= 0 || theta[1] <= 0) reportArgumentError((char *) "%s: argument -v, theta value, must be a positive real number");
  return theta_random_choices[samp];
}

/*! \brief Returns the minimum value of theta specified by the data from the command line.
 */
double parseThetaMin(){
  parseMRateMin(); // Helps this function issue a recommendation for theta
  parseMRateMax(); // Helps this function issue a recommendation for theta
  parseBottleneckMin(); // Helps this function issue a recommendation for theta
  parseBottleneckMax(); // Helps this function issue a recommendation for theta
  if(theta == NULL || theta[0] <= 0 || theta[1] <= 0) reportThetaArgumentError((char *) "%s: argument -v, theta value, must be a positive real number. Recommended input based on choices of mutation rate and bottleneck min and max: -v%e,%e", 4000 * parseMRateMin(), 400000 * parseMRateMax());
  return theta[0];
}

/*! \brief Returns the maximum value of theta specified by the data from the command line.
 */
double parseThetaMax(){
  parseMRateMin(); // Helps this function issue a recommendation for theta
  parseMRateMax(); // Helps this function issue a recommendation for theta
  parseBottleneckMin(); // Helps this function issue a recommendation for theta
  parseBottleneckMax(); // Helps this function issue a recommendation for theta
  if(theta == NULL || theta[0] <= 0 || theta[1] <= 0) reportThetaArgumentError((char *) "%s: argument -v, theta value, must be a positive real number. Recommended input based on choices of mutation rate and bottleneck min and max: -v%e,%e", 4000 * parseMRateMin(), 400000 * parseMRateMax());
  return theta[1];
}

/*! \brief Returns the mutation rate of the specified sample.
 */
double parseMRate(int samp){
  if(mutation_rate == NULL || mutation_rate[0] < 0 || mutation_rate[1] < 0) reportArgumentError((char *) "%s: argument -u, mutation rate during simulation, must be a nonnegative real number");
  return mutation_rate_random_choices[samp];
}

/*! \brief Returns the minimum value of the mutation rate specified by the data from the command line.
 */
double parseMRateMin(){
  if(mutation_rate == NULL || mutation_rate[0] < 0 || mutation_rate[1] < 0) reportArgumentError((char *) "%s: argument -u, mutation rate during simulation, must be a nonnegative real number");
  return mutation_rate[0];
}

/*! \brief Returns the maximum value of the mutation rate specified by the data from the command line.
 */
double parseMRateMax(){
  if(mutation_rate == NULL || mutation_rate[0] < 0 || mutation_rate[1] < 0) reportArgumentError((char *) "%s: argument -u, mutation rate during simulation, must be a nonnegative real number");
  return mutation_rate[1];
}

/*! \brief Specifies whether the program is to use random numbers from the GFSR algorithm or C's library.
 */
int parseRFlag(){
  if(randomFlag == -1) reportArgumentError((char *) "%s: argument -r, flag to determine source of random numbers, must be -rGFSR (for GFSR values), -rRESET (for GFSR values with a reset of the GFSR register), or -rC (for values from C's random number generator)");
  return randomFlag;
}

/*! \brief Returns true if only verifying the syntax of the input data.
 */
int parseSyntaxCheck(){
  return syntax_check;
}

/*! \brief Returns true if generating data initialized from a coalescent.
 */
int parseExample(){
  return example;
}

/*! \brief Returns true if generating an example population of known effective population size.
 */
int parseExamplePop(){
  return example_pop;
}

/*! \brief Returns true if computing statistics from input data without simulating intermediate generations.
 */
int parseRawSample(){
  return raw_stats;
}

/*! \brief Returns true if computing statistics from input data without simulating intermediate generations.
 */
int parseSingleGeneration(){
  return single_generation;
}

/*! \brief Get microsatellite motif lengths
 */
int *getMotifLengths(){
  return motif_lengths;
}

void setMotifLengths(int *new_motif_lengths){
  motif_lengths = new_motif_lengths;
}

void freeMotifLengths(){
  free(motif_lengths);
}

/*! \brief Reads microsatellite motif lengths.
 */
int readMicrosatelliteMotifLengths(int argc, char **argv){
  int i, j, inputs;
  int *inputsList;
  for(j = 1; j < argc; j++){
    if(argv[j][1] != 'm') continue;
    inputs = (strlen(argv[j]) - 1) / 2;
    inputsList = (int *) malloc(inputs * sizeof(int));
    for(i = 0; i < inputs; i++){
      sscanf(argv[j] + 2 * i + 2, "%d", inputsList + i);
      if(!(inputsList[i] == 2 || inputsList[i] == 3 || inputsList[i] == 4 || inputsList[i] == 6)) return -1;
    }
  }
  motif_lengths = inputsList;
  return inputs;
}

/*! \brief Parses command line arguments passed in to the program.
 */
void parseArguments(int argc, char **argv){
  resetArguments();
  int i;

  // Parse the program name
  int lengthName = 1 + strlen(argv[0]);
  programName = (char *) malloc(sizeof(char) * lengthName);
  strcpy(programName, argv[0]);

  // Read in the remaining arguments
  for(i = 1; i < argc; i++){
    char *currentArg = argv[i];
    if(currentArg[0] != '-') reportError("Arguments to OneSamp must start with hypens.");
    if(currentArg[1] == 'r'){
      if(randomFlag == 1 || randomFlag == 0){
        reportError("Duplicate flag: -r");
      }
      if(strlen(currentArg) < 3){
        reportError("Mangled command line argument under -r: check documentation.");
      }
      // Randomness
      if(strcmp(currentArg + 2, (char *) "C\0") == 0){
        randomFlag = 1; // Use C random numbers
        opengfsr();
      }
      else if(strcmp(currentArg + 2, (char *) "GFSR\0") == 0){
        randomFlag = 0; // Use GFSR random numbers
        opengfsr();
      }
      else if(strcmp(currentArg + 2, (char *) "RESET\0") == 0){
        randomFlag = 0; // Use GFSR random numbers
        resetgfsr();
      }
      else{
        reportError("Mangled command line argument under -r: check documentation.");
      }
    }
    else if(currentArg[1] == 'l') {
      // Loci
      if(num_loci != -1) reportError("Duplicate flag: -l");
      num_loci = parsePositiveInt(i, argv);
      num_loci_allocation = num_loci;
    }
    else if(currentArg[1] == 'i') {
      // Initial individuals
      if(input_individuals_count != -1) reportError("Duplicate flag: -i");
      input_individuals_count = parsePositiveInt(i, argv);
      input_individuals_count_allocation = input_individuals_count;

      final_individuals_count = input_individuals_count;
    }
    else if(currentArg[1] == 'b') {
      // Bottleneck
      if(bottleneck_individuals_count != NULL) reportError("Duplicate flag: -b");
      bottleneck_individuals_count = parsePositiveIntPair(i, argv);
      if(bottleneck_individuals_count[0] <= 0 || bottleneck_individuals_count[0] & 1 != 0 || bottleneck_individuals_count[1] <= 0 || bottleneck_individuals_count[1] & 1 != 0) reportArgumentError((char *) "%s: argument -b, num of individuals in bottleneck generation, must be a positive even integer at least 2");
      bottleneck_individuals_count[0] /= 2;
      bottleneck_individuals_count[1] /= 2;
    }
    else if(currentArg[1] == 'd') {
      // Bottleneck
      if(bottleneck_length != NULL) reportError("Duplicate flag: -d");
      bottleneck_length = parsePositiveIntPair(i, argv);
    }
    else if(currentArg[1] == 'm' || currentArg[1] == 's') {
      // SNPs versus Microsatellites
      if(isMicrosats != -1){
        reportError("Duplicate flag: -m and/or -s");
      } else if (currentArg[1] == 'm') {isMicrosats = TRUE;}
      else isMicrosats = FALSE;
    }
    else if(currentArg[1] == 't') {
      // Number of repetitions
      if(repetitions != -1) reportError("Duplicate flag: -t");
      repetitions = parsePositiveInt(i, argv);
    }
    else if(currentArg[1] == 'u') {
      // Mutation rate
      if(mutation_rate != NULL) reportError("Duplicate flag: -u");
      mutation_rate = parsePositiveDoublePair(i, argv);
    }
    else if(currentArg[1] == 'v') {
      // Theta
      if(theta != NULL) reportError("Duplicate flag: -v");
      theta = parsePositiveDoublePair(i, argv);
    }
    else if(currentArg[1] == 'x' || currentArg[1] == 'e' || currentArg[1] == 'w' || currentArg[1] == 'g' || currentArg[1] == 'p') {
      // Check syntax
      if(syntax_check != FALSE || example != FALSE || raw_stats != FALSE || single_generation != FALSE)
        reportError("Duplicate flag: -x and/or -e and/or -w and/or -g and/or -p");
      if(currentArg[1] == 'e') example = TRUE;
      if(currentArg[1] == 'x') syntax_check = TRUE;
      if(currentArg[1] == 'w') raw_stats = TRUE;
      if(currentArg[1] == 'g') single_generation = TRUE;
      if(currentArg[1] == 'p') example_pop = TRUE;
    }
    else if(currentArg[1] == 'f') {
      if(minAlleleFrequency != -1) reportError("Duplicate flag: -f");
      minAlleleFrequency = parsePositiveDouble(i, argv);
    }
    else if(currentArg[1] == 'a') {
      // Interpolate absent data
      if(absentDataExtrapolate != FALSE) reportError("Duplicate flag: -a");
      absentDataExtrapolate = TRUE;
    }
    else if(currentArg[1] == 'o') {
      // Threshold to omit loci
      if(omitThreshold != -1) reportError("Duplicate flag: -o");
      omitThreshold = parsePositiveDouble(i, argv);
    }
    else {
      reportError("Unknown flag passed in to OneSamp.");
    }
  }
  // Make sure all of the arguments are valid by fetching them.
  // If we're computing stats directly, avoid allocating the wrong amount of
  // memory. No need to simulate intermediate generations here.
  if(!parseSyntaxCheck()){
    if(parseSingleGeneration()) repetitions = 1;
    allocateStruct7(&bottleneck_individuals_count_random_choices);
    if(parseRawSample()){
      final_individuals_count = input_individuals_count;
      bottleneck_length = 0;
      bottleneck_individuals_count[0] = 2;
      bottleneck_individuals_count[1] = 2;
      repetitions = 1;
    } else if(!parseSyntaxCheck()){
      if(!parseSingleGeneration()) {
        allocateStruct8(&bottleneck_length_random_choices);
        allocateStruct9(&theta_random_choices);
      }
      allocateStruct10(&mutation_rate_random_choices);
    }
    if(parseSingleGeneration()){
      parseNLoci();
      parseInputSamples();
      parseMRate(0);
      parseBottleneck(0);
      parseFormFlag();
      parseRFlag();
    }
    else
    {
      parseFormFlag();
      if(!parseRawSample()) parseIterations();
      parseNLoci();
      parseInputSamples();
      parseFormFlag();
      parseOmitLocusThreshold();
      if(!parseRawSample()) parseBottleneck(0);
      if(!parseRawSample()) parseBottleneckLength(0);
      if(!parseRawSample()) parseMRate(0);
      if(!parseRawSample()) parseRFlag();
      if(!parseRawSample()) parseTheta(0);
    }
  }
  if(!parseSyntaxCheck() && !parseRawSample() && !parseExample() && !parseSingleGeneration() && !parseExamplePop()){
    reportArgumentError((char *) "%s: missing an operation to perform on the input file: -x (syntax check operation), -w (compute statistics of input sample), -g (simulate a single generation from an input population and display to standard out), -p (dump out an example population with known effective population size), or -e (compute stats of coalescent sample after a few generations have passed)");
  }
}

/*! \brief Deallocates memory for data from command line arguments.
 */
void flushArguments(){
  free(programName); programName = NULL;
  free(bottleneck_individuals_count); bottleneck_individuals_count = NULL;
  free(bottleneck_individuals_count_random_choices); bottleneck_individuals_count_random_choices = NULL;
  if(parseRawSample()){
    free(bottleneck_length); bottleneck_length = NULL;
    free(theta); theta = NULL;
    free(mutation_rate); mutation_rate = NULL;
    free(bottleneck_length_random_choices); bottleneck_length_random_choices = NULL;
    free(theta_random_choices); theta_random_choices = NULL;
    free(mutation_rate_random_choices); mutation_rate_random_choices = NULL;
  }

}

/*! \brief Retrieves positive double from command line argument.
 */
double parsePositiveDouble(int argNum, char **argv){
  char *arg = argv[argNum];
  double result = -1.0;
  sscanf(arg + 2, "%lf", &result);
  if(result < 0) result = -1;
  return result;
}

/*! \brief Retrieves positive double pair from command line argument.
 */
double *parsePositiveDoublePair(int argNum, char **argv){
  double temp;
  char *arg = argv[argNum];
  double *result =(double *) malloc(2 * sizeof(double));
  result[0] = -1;
  result[1] = -1;
  sscanf(arg + 2, "%lf,%lf", result, result + 1);
  if(result[0] < 0) result[0] = -1;
  if(result[1] < 0) result[1] = -1;
  if(result[1] == -1) result[1] = result[0];
  if(result[0] > result[1])
  {
    temp = result[0]; result[0] = result[1]; result[1] = temp;
  }
  return result;
}

/*! \brief Retrieves positive int from command line argument.
 */
int parsePositiveInt(int argNum, char **argv){
  char *arg = argv[argNum];
  int result = -1;
  sscanf(arg + 2, "%d", &result);
  if(result < 0) result = -1;
  return result;
}

/*! \brief Retrieves positive double pair from command line argument.
 */
int *parsePositiveIntPair(int argNum, char **argv){
  int temp;
  char *arg = argv[argNum];
  int *result = (int *) malloc(2 * sizeof(int));
  result[0] = -1;
  result[1] = -1;
  sscanf(arg + 2, "%d,%d", result, result + 1);
  if(result[0] < 0) result[0] = -1;
  if(result[1] < 0) result[1] = -1;
  if(result[1] == -1) result[1] = result[0];
  if(result[0] > result[1])
  {
    temp = result[0]; result[0] = result[1]; result[1] = temp;
  }
  return result;
}

/*! \brief Returns the proportion of missing data in the population.
 */
double getProportionMissingData(){
  return proportionMissingData;
}

/*! \brief Sets the proportion of missing data in the population.
 */
void setProportionMissingData(double val){
  proportionMissingData = val;
}

/* \brief Initialize a genotype at a locus
 */
int initializeMicrosat1(int locus){
  int i;
  int val = locus % parseNLoci();
  if(parseExamplePop()) return 192; // A sample microsatellite length
  int totalMicrosatLength = 0;
  int gene1;
  int gene2;
  for(i = 0; i < parseInputSamples(); i++){
    loadInitialGenotype(i, val, &gene1, &gene2);
    totalMicrosatLength += gene1 + gene2;
  }
  // Return microsatellite length rounded to nearest multiple of motif
  return (totalMicrosatLength * getMotifLengths()[val] + (getMotifLengths()[val] / 2)) / (2 * getMotifLengths()[val] * parseInputSamples());
}

/* \brief Initialize a distinct genotype at a locus
 */
int initializeMicrosat2(int locus){
  return initializeMicrosat1(locus) + getMotifLengths()[locus % parseNLoci()];
}

// INCLUDE GUARD

#endif
