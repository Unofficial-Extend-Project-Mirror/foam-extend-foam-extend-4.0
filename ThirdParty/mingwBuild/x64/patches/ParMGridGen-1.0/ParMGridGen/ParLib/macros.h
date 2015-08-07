/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * macros.h
 *
 * This file contains macros used in multilevel
 *
 * George Irene
 */

#include "rand48.h"

/*************************************************************************
* This macro is used to normalize the weights of two nodes
**************************************************************************/
#define ARATIO1(dim, surf, vol) ((dim == 2) ? (pow((surf), 2)/(vol)) : (pow((surf), 1.5)/(vol)))
#define ARATIO(dim, surf, vol) ((dim == 2) ? ((surf)*(surf)/(vol)) : (sqrt((surf)*(surf)*(surf))/(vol)))
#define ARATIO2(dim, surf, vol) ((dim == 2) ? ((surf)*(surf)*(surf)*(surf)/(vol)*(vol)) : ((surf)*(surf)*(surf)/((vol)*(vol))))


/*************************************************************************
* The following macro returns a random number in the specified range
**************************************************************************/
#define RandomInRange(u) ((int)(drand48()*((double)(u))))

#define amax(a, b) ((a) >= (b) ? (a) : (b))
#define amin(a, b) ((a) >= (b) ? (b) : (a))

#define AND(a, b) ((a) < 0 ? ((-(a))&(b)) : ((a)&(b)))
#define OR(a, b) ((a) < 0 ? -((-(a))|(b)) : ((a)|(b)))
#define XOR(a, b) ((a) < 0 ? -((-(a))^(b)) : ((a)^(b)))

#define SWAP(a, b, tmp)  \
                 do {(tmp) = (a); (a) = (b); (b) = (tmp);} while(0)

#define INC_DEC(a, b, val) \
                 do {(a) += (val); (b) -= (val);} while(0)


#define HASHFCT(key, size) ((key)%(size))


/*************************************************************************
* Timer macros
**************************************************************************/
#undef cleartimer
#undef starttimer
#undef stoptimer
#undef gettimer
#define cleartimer(tmr) (tmr = 0.0)
#define starttimer(tmr) (tmr -= MPI_Wtime())
#define stoptimer(tmr) (tmr += MPI_Wtime())
#define gettimer(tmr) (tmr)


/*************************************************************************
* This macro is used to handle dbglvl
**************************************************************************/
#define IFSET(a, flag, cmd) if ((a)&(flag)) (cmd);

#undef ASSERT
#undef ASSERTP

#ifdef DEBUG
#   define ASSERT(ctrl, expr)                                          \
    if (!(expr)) {                                               \
        MGridmyprintf(ctrl, "***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        abort();                                                \
    }
#else
#   define ASSERT(ctrl, expr) ;
#endif

#ifdef DEBUG
#   define ASSERTP(ctrl, expr,msg)                                          \
    if (!(expr)) {                                               \
        MGridmyprintf(ctrl, "***ASSERTION failed on line %d of file %s:" #expr "\n", \
              __LINE__, __FILE__);                               \
        MGridmyprintf msg ; \
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
