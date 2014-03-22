/**************************************************************************************************
 * Jean-Baptiste Pettit
 * European Bioinformatics Institute
 * 
 * NOTICE OF LICENSE
 *
 * This source file is subject to the Academic Free License (AFL 3.0)
 * that is bundled with this package in the file LICENSE_AFL.txt.
 * It is also available through the world-wide-web at this URL:
 * http://opensource.org/licenses/afl-3.0.php.
****************************************************************************************************/


#define LG_LIG_MAX 1000000
#define LG_GENES_MAX 1000
#define LG_NEI_MAX 200
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#ifdef linux
#include <sys/signal.h>
#endif

/****************************START Defining structures***********************************/
	/*Model parameters we want to estimate or that are set*/	
	typedef struct {
		double ** theta; /* k*p float table Be parameters*/
	} params;
	
	/*Data point expression values and properties*/
	typedef struct {
		int * expVect; /* expression vector*/
		int * nei; /*Indexes of neighbours in dataSet*/
		int numNei; /*Number of neighbours*/
	} dataPoint;

	typedef struct {
		dataPoint * obs; /*observation x*/
		int num; /*number of Points*/
		int length; /*length of expression vectors*/
	} dataSet;


	
	/*Classification Z*/
	typedef struct {
		int * clust;
		params parameters;
		dataSet set;
		int numClust;
		double ** tihm;
		double ** cellDensities;
		double likelihood;
		double fullLikelihood;
		double * beta; /*smoothness param */
	} classif;

	

/****************************END Defining structures***********************************/

/****************************START function prototypes***********************************/
/*Parses one line, places numbers un int vector*/
int parse_vect(char * line, int * res);


/*Reading and storing data
returns : void
*/
void load_data(FILE * data /*I*/,dataSet * myData/*I\O*/);

/*Parses neighbouring files, first number is number of neighbours then indexes*/
int parse_nei(char * line, int * res);


/*Reading and storing spatial information
returns : void
*/
void load_nei(FILE * data /*I*/,dataSet * myData/*I\O*/);

/*Initialize classification randomly
returns void*/
void initClassifRand(dataSet set/*I*/,int numClust/*I*/, classif* myClassif/*I\O*/,double beta);

/*Initialize classification from file
returns void*/
void initClassifFile(dataSet set/*I*/,FILE * data, classif* myClassif/*I\O*/,double beta);


/*function to check if all elements of the count vector > 0
* returns 1 if vector contains 0s 0 otherwise
*/
int zerosInVector(int * vector,int length);
/*Checks if current classif has at least one point in each cluster
*If not it displays a warning.
*Empty classes usually lead to NaN final likelihood, if you get that error a lot try initializing the the clutering differently
*
*/
void noEmptyClass(classif * myClassif);
/*Computes logLikelihood*/
double computeFullLogLikelihood(classif * myClassif/*I*/);
/*Computes the beta part of the expected likelihood*/
double computeBetaExpectation(classif * myClassif/*I*/);
/*Computes expectation*/
double computeFullExpectation(classif * myClassif/*I*/);

double computePseudoLogLikelihood(classif * myClassif/*I*/);
/*Set model pseudo-logLikelihood*/
void computeCellDensities(classif * myClassif/*I/O*/);
double cellDensity(classif * myClassif/*i*/,int clust/*I*/,int cell/*I*/);
double logCellDensity(classif * myClassif/*i*/,int clust/*I*/,int cell/*I*/);


double computeNeiCoef(classif * myClassif/*i\O*/,double ** currentT/*I*/,int clust/*I*/,int cell/*I*/);


/*Computes current thims*/
void computeThims(classif * myClassif /*I/O*/, int numIterFixed /*I*/);

/*Gradient descent algorithm to maximize an unique beta*/
void gradientAscent(classif * myClassif /*I\O*/);

/*E Step of the EM algorithm*/
int eStep(classif * myClassif /*I\O*/);

int numCellsAlike(int exp,int indexGene,int clust,classif * myClassif);

/*Maximize thetas*/
void maxThetas(classif * myClassif /*I/O*/);

/*M Step of the algorithm*/
void mStep(classif * myClassif,int type_beta);



/****************************END function prototypes***********************************/
