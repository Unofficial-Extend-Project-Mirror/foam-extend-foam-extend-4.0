/*
 * IMlib.h
 *
 * Irene's library of most frequently used routines
 *
 */

#ifndef _IMLIB_H_
#define _IMLIB_H_


/* Undefine the following #define in order to use short int as the idxtype */
#define IDXTYPE_INT
/* Undefine the following #define in order to use float as the realtype */
/*#define TYPE_REAL*/

/* Indexes are as long as integers for now */
#ifdef IDXTYPE_INT
typedef int idxtype;
#else
typedef short idxtype;
#endif

/* floats for now */
#ifdef TYPE_REAL
typedef float realtype;
#else
typedef double realtype;
#endif


/*************************************************************************
* Header file inclusion section
**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>
#include "rand48.h"

#ifdef DMALLOC
#include <dmalloc.h>
#else
#include <malloc.h>
#endif

/*************************************************************************
* Data structure definition section
**************************************************************************/
/*-------------------------------------------------------------
 * The following data structure stores int key-value pairs
 *-------------------------------------------------------------*/
struct IKeyValueType {
  int key;
  int val;
};

typedef struct IKeyValueType IKeyValueType;


/*-------------------------------------------------------------
 * The following data structure stores int key-value pairs
 *-------------------------------------------------------------*/
struct idxKeyValueType {
  idxtype key;
  idxtype val;
};

typedef struct idxKeyValueType idxKeyValueType;


/*-------------------------------------------------------------
 * The following data structure stores int-key - double-value pairs
 *-------------------------------------------------------------*/
struct FKeyValueType {
  double key;
  int val, val1, val2;
};

typedef struct FKeyValueType FKeyValueType;


/*-------------------------------------------------------------
 * The following data structure stores int-key - double-value pairs
 *-------------------------------------------------------------*/
struct realKeyValueType {
  realtype key;
  int val, val1, val2;
};

typedef struct realKeyValueType realKeyValueType;

/*************************************************************************
* Definition Section
**************************************************************************/
#define LTERM                   (void **) 0     /* List terminator for IMfree() */



/*************************************************************************
* Macros Section
**************************************************************************/
/*-------------------------------------------------------------
 * Usefull commands
 *-------------------------------------------------------------*/
#define sign(a, b) ((b) >= 0 ? ((a) >= 0.0 ? a : -a) : ((a) >= 0.0 ? -a : a))
#define amax(a, b) ((a) >= (b) ? (a) : (b))
#define amin(a, b) ((a) >= (b) ? (b) : (a))
#define RandomInRange(u) ((int)(drand48()*((double)(u))))
#define RandomInRangeFast(u) ((rand()>>3)%(u))
#define SWAP(a, b, tmp) do {(tmp) = (a); (a) = (b); (b) = (tmp);} while(0)
#define INC_DEC(a, b, val) do {(a) += (val); (b) -= (val);} while(0)
#define icopy(n, a, b) (int *)memcpy((void *)(b), (void *)(a), sizeof(int)*(n))
#define idxcopy(n, a, b) (idxtype *)memcpy((void *)(b), (void *)(a), sizeof(idxtype)*(n))
#define scopy(n, a, b) (double *)memcpy((void *)(b), (void *)(a), sizeof(double)*(n))
#define fcopy(n, a, b) (double *)memcpy((void *)(b), (void *)(a), sizeof(double)*(n))
#define realcopy(n, a, b) (realtype *)memcpy((void *)(b), (void *)(a), sizeof(realtype)*(n))


/*-------------------------------------------------------------
 * Timing macros
 *-------------------------------------------------------------*/
#define cleartimer(tmr) (tmr = 0.0)
#define starttimer(tmr) (tmr -= seconds())
#define stoptimer(tmr) (tmr += seconds())
#define gettimer(tmr) (tmr)


/*-------------------------------------------------------------
 * Debuging memory leaks
 *-------------------------------------------------------------*/
#ifdef DMALLOC
#define imalloc(n, msg) (malloc(sizeof(int)*(n)))
#define fmalloc(n, msg) (malloc(sizeof(double)*(n)))
#define idxmalloc(n, msg) ((idxtype *)malloc(sizeof(idxtype)*(n)))
#define realmalloc(n, msg) ((realtype *)malloc(sizeof(realtype)*(n)))
#define ismalloc(n, val, msg) (iset((n), (val), malloc(sizeof(int)*(n))))
#define idxsmalloc(n, val, msg) (idxset((n), (val), (idxtype *)malloc(sizeof(idxtype)*(n))))
#define fsmalloc(n, val, msg) (fset((n), (val), malloc(sizeof(double)*(n))))
#define realsmalloc(n, val, msg) (realset((n), (val), (realtype *)malloc(sizeof(realtype)*(n))))
#define IMmalloc(a, b) (malloc((a)))
#endif

#ifdef DMALLOC
#   define MALLOC_CHECK(ptr)                                          \
    if (malloc_verify((ptr)) == DMALLOC_VERIFY_ERROR) {  \
        printf("***MALLOC_CHECK failed on line %d of file %s: " #ptr "\n", \
              __LINE__, __FILE__);                               \
        abort();                                                \
    }
#else
#   define MALLOC_CHECK(ptr) ;
#endif


/*-------------------------------------------------------------
 * CSR conversion macros
 *-------------------------------------------------------------*/
#define MAKECSR(i, n, a) \
   do { \
     for (i=1; i<n; i++) a[i] += a[i-1]; \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0)


/*-------------------------------------------------------------
 * Program Assertions
 *-------------------------------------------------------------*/
#ifdef DEBUG
#   define ASSERT(expr)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        abort();                                                \
    }
#else
#   define ASSERT(expr) ;
#endif

#ifdef DEBUG
#   define ASSERTP(expr,msg)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        printf msg ; \
        printf("\n"); \
        abort();                                                \
    }
#else
#   define ASSERTP(expr,msg) ;
#endif


/*************************************************************************
* Function prototypes
**************************************************************************/
/*-------------------------------------------------------------
 * blas.c
 *-------------------------------------------------------------*/
int *iset(int, int, int *);
idxtype *idxset(int, idxtype, idxtype *);
double *fset(int, double, double *);
realtype *realset(int, realtype, realtype *);
int iamax(int, int *);
int idxamax(int, idxtype *);
int famax(int, double *);
int iamin(int, int *);
int idxamin(int, idxtype *);
int famin(int, double *);
int charsum(int, char *);
int isum(int, int *);
int idxsum(int, idxtype *);
double ssum(int, double *);
double ssum_strd(int, double *, int);
void sscale(int, double, double *);
double snorm2(int, double *);
double sdot(int, double *, double *);
void saxpy(int, double, double *, int, double *, int);


/*-------------------------------------------------------------
 * file.c
 *-------------------------------------------------------------*/
FILE *IMfopen(char *, char *, char *);
void IMfclose(FILE *);

/*-------------------------------------------------------------
 * memory.c
 *-------------------------------------------------------------*/
#ifndef DMALLOC
int *imalloc(int, char *);
idxtype *idxmalloc(int, char *);
double *fmalloc(int, char *);
realtype *realmalloc(int, char *);
int *ismalloc(int, int, char *);
idxtype *idxsmalloc(int, idxtype, char *);
double *fsmalloc(int, double, char *);
realtype *realsmalloc(int, realtype, char *);
void *IMmalloc(int, char *);
#endif
/* void IMfree(void **, ...); */


/*-------------------------------------------------------------
 * util.c
 *-------------------------------------------------------------*/
void *errexit(char *,...);
int IMlog2(int);
double flog2(double);
int ispow2(int);
double seconds(void);


/*-------------------------------------------------------------
 * Sorting routines
 *-------------------------------------------------------------*/
void dfkeysort(int, FKeyValueType *);
void dkeysort(int, IKeyValueType *);
void ifkeysort(int, FKeyValueType *);
void ifkeysort2(int, FKeyValueType *);
void ifloatsort(int, double *);
void iintsort(int, int *);
void ikeysort(int, IKeyValueType *);
void idxkeysort(int, idxKeyValueType *);
void ikeysort2(int, IKeyValueType *);
void idxkeysort2(int, idxKeyValueType *);

/*-------------------------------------------------------------
 * sort.c
 *-------------------------------------------------------------*/
void ikeyvalsort_org(int, IKeyValueType *);
int IncKeyValueCmp(const void *, const void *);
void dkeyvalsort(int, IKeyValueType *);
void BucketSortKeysInc(int, idxtype, idxtype *, int *, int *);
int DecKeyValueCmp(const void *, const void *);
int BSearch(int, idxtype *, int);
void RandomPermute(int, idxtype *, int);
void RandomPermuteFine(int, int *, int);
void FastRandomPermute(int, idxtype *, int);

#endif
