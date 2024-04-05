#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FRANDOM drand48()  //put here your favorite random number generator

#define Q 5

struct vector {
  int n[Q];
} zero, sum, sum_flag, group[Q], *perm;

int N, NoverQ, M, *graph, **neigh, **flag, *deg, *color, numUnsat, fact[Q+1];
double B;



void error(char *string) {
  fprintf(stderr, "ERROR: %s\n", string);
  exit(EXIT_FAILURE);
}

void allocateMem(void) {
  int i;

  fact[0] = 1;
  for ( i = 0; i < Q; i++) {
    zero.n[i] = 0;
    fact[i+1] = (i+1) * fact[i];
  }
  graph = (int*)calloc(2*M, sizeof(int));
  deg = (int*)calloc(N, sizeof(int));
  neigh = (int**)calloc(N, sizeof(int*));
  flag = (int**)calloc(N, sizeof(int*));
  color = (int*)calloc(N, sizeof(int));
  perm = (struct vector*)calloc(fact[Q], sizeof(struct vector));
}


struct vector rotate(struct vector input, int modulo) {
  struct vector output;
  int i;

  for (i = 0; i < modulo; i++)
    output.n[i] = (input.n[i] + 1) % modulo;
  return output;
}

void initPerm(int max) { //permutations needed to compute the overlap
  int i;
  
  if (max == 1)
    perm[0].n[0] = 0;
  else {
    initPerm(max-1);
    for (i = 0; i < fact[max-1]; i++)
      perm[i].n[max-1] = max-1;
    for (; i < fact[max]; i++)
      perm[i] = rotate(perm[i-fact[max-1]], max);
  }
}

void makePlantedGraph(void) {
  int i, var1, var2;

  for (i = 0; i < N; i++)
    deg[i] = 0;
  for (i = 0; i < M; i++) {
    var1 = (int)(FRANDOM * N);
    do {
      var2 = (int)(FRANDOM * N);
    } while ((int)(var1/NoverQ) == (int)(var2/NoverQ)); 
    graph[2*i] = var1;
    graph[2*i+1] = var2;
    deg[var1]++;
    deg[var2]++;
  }        //we do not check double links, but in the large N limit the probability of these events is very low
  for (i = 0; i < N; i++) {
    neigh[i] = (int*)calloc(deg[i], sizeof(int));
    flag[i] = (int*)calloc(deg[i], sizeof(int));
    deg[i] = 0;
  }
  for (i = 0; i < M; i++) {
    var1 = graph[2*i];
    var2 = graph[2*i+1];
    neigh[var1][deg[var1]++] = var2;
    neigh[var2][deg[var2]++] = var1;
  }
}

void initColors(void) {
  int i;
  
  for (i = 0; i < N; i++)
    color[i] = (int)(FRANDOM * Q);
}

void initFlag(){
  int i,j;

  for (i = 0; i < N; i++){
    for (j = 0; j < deg[i]; j++){
      if(FRANDOM>B){
	flag[i][j]=0;
      }else{
	flag[i][j]=1;
      }
    }
  }
  
}

void quenchT0() {
  int i, j, newColor, deltaUnsat, deltaUnsat_flag;

  for (i = 0; i < N; i++) {
    do {
      newColor = (int)(FRANDOM * Q);
    } while (newColor == color[i]);
    sum = zero;
    sum_flag = zero;
    for (j = 0; j < deg[i]; j++){
      if(flag[i][j]==1){
	sum_flag.n[color[neigh[i][j]]]++;
      }
      sum.n[color[neigh[i][j]]]++;
    }
    
    deltaUnsat = sum.n[newColor] - sum.n[color[i]];
    deltaUnsat_flag = sum_flag.n[newColor] - sum_flag.n[color[i]];
    if (deltaUnsat_flag <= 0) {
      color[i] = newColor;
      numUnsat += deltaUnsat;
    }
  }
  
}

int computeUnsat(void) {
  int i, res = 0;

  for (i = 0; i < M; i++)
    res += (color[graph[2*i]] == color[graph[2*i+1]]);
  return res;
}

double overlapPlanted(void) {
  int i, j, overlap, maxOver=0;

  for (i = 0; i < Q; i++)
    group[i] = zero;
  for (i = 0; i < N; i++)
    group[(int)(i/NoverQ)].n[color[i]]++;
  for (i = 0; i < fact[Q]; i++) {
    overlap = 0;
    for (j = 0; j < Q; j++)
      overlap += group[j].n[perm[i].n[j]];
    if (overlap > maxOver) maxOver = overlap;
  }
  return (double)(Q*maxOver-N)/(Q-1)/N;
}

int main(int argc, char *argv[]) {
  int i, nIter;
  double c;
  long int seed;

  FILE *devran = fopen("/dev/urandom","r");
  fread(&seed, 4, 1, devran);
  fclose(devran);

  if (argc != 5) {
    fprintf(stderr, "usage: %s <N> <c> <max number of iterations> <B: Fraction of edges in minibatch> \n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  N = atoi(argv[1]);
  c = atof(argv[2]);
  nIter = atoi(argv[3]);
  B = atof(argv[4]);
  
  if (Q * (int)(N/Q) != N) error("Q must divide N");
  NoverQ = (int)(N/Q);
  M = (int)(0.5 * c * N + 0.5);
  printf("# Q = %i   N = %i   M = %i   c = %f   nIter = %i   B: %g   seed = %ld\n", Q, N, M, c, nIter, B, seed);
  printf("# 1:t 2:intensive_energy  3:overlapPlanted\n");

  srand48(seed);

  allocateMem();
  initPerm(Q);

  makePlantedGraph();
  initColors();
  numUnsat = computeUnsat();

  for (i = nIter; i > 0 && numUnsat ; i--) {
    initFlag();
    quenchT0();
    printf("%d %g %g\n",nIter-i, numUnsat/((double)N), overlapPlanted());
  }

  return EXIT_SUCCESS;
}
