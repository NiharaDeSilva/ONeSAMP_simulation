#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>



/* \brief Initialize a genotype at a locus
 */
int initializeMicrosat1(int locus){
  return 200;
}

/* \brief Initialize a distinct genotype at a locus
 */
int initializeMicrosat2(int locus){
  return 202;
}


#ifndef TYPEDEFS
#define TYPEDEFS
#define GFSR_STYPE_UNSIGNED_MAX UINT_MAX
#define STRUCT_GTYPE_SIZE (2 * sizeof(int *))
#define ALLELE_TYPE short
struct gtype_type {
    ALLELE_TYPE *pgtype;
    ALLELE_TYPE *mgtype;
};
typedef struct gtype_type gtype_type;
#define gfsr4() (((float) rand()) / INT_MAX ) / (GFSR_STYPE_UNSIGNED_MAX + 1.0)
#define disrand(l, t) ((int) ((((unsigned int) intrand())%((t) - (l) + 1)) + (l)))
#define intrand() rand()
#define randBit() disrand(0, 1)
#endif

int lociDirectInput = 0;
double parseNLoci(){
  return lociDirectInput;
}

double thetaDirectInput = 0;
double parseTheta(int i){
  return thetaDirectInput;
}

int individualsDirectInput = 0;
int parseBottleneck(){
  return individualsDirectInput / 2;
}

int formFlag = 0;
int parseFormFlag(){
  return formFlag;
}

double minAlleleFrequencyDirectInput = 0;
double parseMinAlleleFrequency(){
  return minAlleleFrequencyDirectInput;
}

/*! \def fallingQuotient(s, t1, t2, c)
 *  \brief Returns s * ((t1)/(t2)) * ((t1-1)/(t2-1)) * ... ((t1-c+1)/(t2-c+1))
 *  The value c must be nonnegative.
 */
double fallingQuotient(double s, double t1, double t2, int c){
  return (c == 0) ? s : fallingQuotient(s * t1 / t2, t1 - 1, t2 - 1, c - 1);
}

/*! \def allelePr()
 *  \brief Returns likelihood of finding val1 alleles in one locus and val2 mutated alleles
 */
double allelePr(int val1, int val2, double theta){
  int n = val1 + val2;
  double result = fallingQuotient(1, n, theta + n - 1, n);
  int j;
  if(val1 != 0) result *= theta / (double) val1;
  if(val2 != 0) result *= theta / (double) val2;
  if(val1 == val2) result /= 2;
  return result;
}

int main(int argc, char **argv) {
  printf("hope  not");
  if(argc != 6){
    printf("usage: ./csim SNP <theta> <loci> <diploid individuals> <min allele threshold frequency>\n");
    printf("usage: ./csim MICROSAT <theta> <loci> <diploid individuals> <min allele threshold frequency>\n");
    exit(1);
  }

  // Seed random number generator
  srand(time(NULL));

  // Load command line arguments
  int offset = 1;

  // Formflag is 1 if we are working with microsatellite data, 0 otherwise
  formFlag = argv[offset][0] == 'M' ? 1 : 0;

  // Theta is the specified evolutionary parameter
  thetaDirectInput = (double) atof(argv[offset + 1]);

  // Loci is the numver of loci we have in the simulation
  lociDirectInput = (int) atoi(argv[offset + 2]);

  // Individuals is the number of diploid individuals present in the sample
  individualsDirectInput = (int) atoi(argv[offset + 3]);

  // Min allele frequency is the minimum frequency of the mutated allele
  minAlleleFrequencyDirectInput = (double) atof(argv[offset + 4]);

  // Load dummy variables
  int i = 0;
  int k = 0;
  int extraProportionOfBufferLoci = 1;

  // These data structures are targets adapted to attach to the main code base
  struct gtype_type *females[1], *males[1];
  females[0] = malloc(parseBottleneck() * STRUCT_GTYPE_SIZE);
  males[0] = malloc(parseBottleneck() * STRUCT_GTYPE_SIZE);
  for(k = 0; k < parseBottleneck(); k++){
    females[0][k].pgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci()*sizeof(ALLELE_TYPE));
    males[0][k].pgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci()*sizeof(ALLELE_TYPE));
    females[0][k].mgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci()*sizeof(ALLELE_TYPE));
    males[0][k].mgtype = (ALLELE_TYPE *)malloc(extraProportionOfBufferLoci * parseNLoci()*sizeof(ALLELE_TYPE));
  }

  #include "../engine/refactor_coalescent_engine.txt"
  
  int num_loci = extraProportionOfBufferLoci * parseNLoci();
  int num_samples = 1;
  int final_indivs_count = parseBottleneck();
  int l;
  printf("Auto-generated genotype output.\n");
  for (l = 0; l < num_loci;++l) {
    printf("%1d\n",l+1);
  }
  printf("Pop\n");
  for(l = 0; l < final_indivs_count; l++){
    printf("%d , ", l+1);
    for(i = 0; i < num_loci;++i){
      printf("%02d%02d ", females[0][l].mgtype[i], females[0][l].pgtype[i]);
    }
    printf("\n");
  }
  for(l = 0; l < final_indivs_count; l++){
    printf("%d , ", final_indivs_count + l+1);
    for(i = 0; i < num_loci;++i){
      printf("%02d%02d ", males[0][l].mgtype[i], males[0][l].pgtype[i]);
    }
    printf("\n");
  }

  // Deallocate intermediate genotype arrays
  for(j = 0; j < parseBottleneck(); j++){
    free(females[0][j].pgtype);
    free(males[0][j].pgtype);
    free(females[0][j].mgtype);
    free(males[0][j].mgtype);
  }
  free(females[0]);
  free(males[0]);
  return 0;
}
