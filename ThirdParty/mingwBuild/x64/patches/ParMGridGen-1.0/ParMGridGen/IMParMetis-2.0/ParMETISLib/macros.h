/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * macros.h
 *
 * This file contains macros used in multilevel
 *
 * George
 */

#include "rand48.h"

/*************************************************************************
* The following macro returns a random number in the specified range
**************************************************************************/
#define RandomInRange(u) ((int)(drand48()*((double)(u))))
#define RandomInRangeFast(u) ((rand()>>2)%(u))

#define amax(a, b) ((a) >= (b) ? (a) : (b))
#define amin(a, b) ((a) >= (b) ? (b) : (a))

#define AND(a, b) ((a) < 0 ? ((-(a))&(b)) : ((a)&(b)))
#define OR(a, b) ((a) < 0 ? -((-(a))|(b)) : ((a)|(b)))
#define XOR(a, b) ((a) < 0 ? -((-(a))^(b)) : ((a)^(b)))

#define SWAP(a, b, tmp)  \
                 do {(tmp) = (a); (a) = (b); (b) = (tmp);} while(0)

#define INC_DEC(a, b, val) \
                 do {(a) += (val); (b) -= (val);} while(0)


#define scopy(n, a, b) memcpy((b), (a), sizeof(float)*(n))
#define idxcopy(n, a, b) memcpy((b), (a), sizeof(idxtype)*(n))
#define realcopy(n, a, b) memcpy((b), (a), sizeof(realtype)*(n))


#define HASHFCT(key, size) ((key)%(size))


/*************************************************************************
* Timer macros
**************************************************************************/
#define cleartimer(tmr) (tmr = 0.0)
#define starttimer(tmr) (tmr -= MPI_Wtime())
#define stoptimer(tmr) (tmr += MPI_Wtime())
#define gettimer(tmr) (tmr)


/*************************************************************************
* This macro is used to handle dbglvl
**************************************************************************/
#define IFSET(a, flag, cmd) if ((a)&(flag)) (cmd);

/*************************************************************************
* These macros are used for debuging memory leaks
**************************************************************************/
#ifdef DMALLOC
#define imalloc(n, msg) (malloc(sizeof(int)*(n)))
#define fmalloc(n, msg) (malloc(sizeof(float)*(n)))
#define idxmalloc(n, msg) (malloc(sizeof(idxtype)*(n)))
#define realmalloc(n, msg) ((realtype *)malloc(sizeof(realtype)*(n)))
#define ismalloc(n, val, msg) (iset((n), (val), malloc(sizeof(int)*(n))))
#define idxsmalloc(n, val, msg) (idxset((n), (val), malloc(sizeof(idxtype)*(n))))
#define realsmalloc(n, val, msg) (realset((n), (val), (realtype *)malloc(sizeof(realtype)*(n))))
#define IMmalloc(a, b) (malloc(a))
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

/*************************************************************************
* This macro converts a length array in a CSR one
**************************************************************************/
#define MAKECSR(i, n, a) \
   do { \
     for (i=1; i<n; i++) a[i] += a[i-1]; \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0)




#ifdef DEBUG
#   define ASSERT(ctrl, expr)                                          \
    if (!(expr)) {                                               \
        myprintf(ctrl, "***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        abort();                                                \
    }
#else
#   define ASSERT(ctrl, expr) ;
#endif

#ifdef DEBUG
#   define ASSERTP(ctrl, expr,msg)                                          \
    if (!(expr)) {                                               \
        myprintf(ctrl, "***ASSERTION failed on line %d of file %s:" #expr "\n", \
              __LINE__, __FILE__);                               \
        myprintf msg ; \
        abort();                                                \
    }
#else
#   define ASSERTP(ctrl, expr,msg) ;
#endif

#ifdef DEBUGS
#   define ASSERTS(expr)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        abort();                                                \
    }
#else
#   define ASSERTS(expr) ;
#endif

#ifdef DEBUGS
#   define ASSERTSP(expr, msg)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        printf msg ; \
        abort();                                                \
    }
#else
#   define ASSERTSP(expr, msg) ;
#endif
