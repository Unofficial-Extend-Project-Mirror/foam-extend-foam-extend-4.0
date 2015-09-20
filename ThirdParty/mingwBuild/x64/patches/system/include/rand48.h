#ifndef RAND48_H
#define RAND48_H

#define drand48()       (rand()*(1./RAND_MAX))
static long _rand = 1;

static __inline__ void srand48(long seed)
{
        _rand = seed;
}

static __inline__ long lrand48(void)
{
        long val = (int)(abs(10000.0*sin(_rand)));
        _rand++;
        return val;
}

#endif //RAND48_H
