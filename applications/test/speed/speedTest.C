#include "primitiveFields.H"
#include "cpuTime.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    const label nIter = 10;
    const label size = 10000000;

    double* f1 = new double[size];
    double* f2 = new double[size];
    double* f3 = new double[size];
    double* f4 = new double[size];

    for (register label i=0; i<size; i++)
    {
        f1[i] = 1.0;
        f2[i] = 1.0;
        f3[i] = 1.0;
    }

    cpuTime executionTime1;

    for (int j=0; j<nIter; j++)
    {
        for (register label i=0; i<size; i++)
        {
            f4[i] = f1[i] + f2[i] - f3[i];
        }
    }

    Info<< "ExecutionTime = "
        << executionTime1.elapsedCpuTime()
        << " s\n" << endl;

    Info << f4[1] << endl << endl;


    scalarField sf1(size, 1.0), sf2(size, 1.0), sf3(size, 1.0), sf4(size);

    cpuTime executionTime2;

    for (register int j=0; j<nIter; j++)
    {
        sf4 = sf1 + sf2 - sf3;
        //sf4 = sf1;
        //sf4 += sf2;
        //sf4 -= sf3;
    }

    Info<< "ExecutionTime = "
        << executionTime2.elapsedCpuTime()
        << " s\n" << endl;

    Info << sf4[1] << endl << endl;


    vectorField 
        vf1(size, vector::one),
        vf2(size, vector::one),
        vf3(size, vector::one),
        vf4(size);

    cpuTime executionTime3;

    for (register int j=0; j<nIter; j++)
    {
        vf4 = vf1 + vf2 - vf3;
    }

    Info<< "ExecutionTime = "
        << executionTime3.elapsedCpuTime()
        << " s\n" << endl;

    Info << vf4[1] << endl << endl;

    cpuTime executionTime4;

    scalarField sf11(size, 1.0), sf12(size, 1.0), sf13(size, 1.0), sf14(size);
    scalarField sf21(size, 1.0), sf22(size, 1.0), sf23(size, 1.0), sf24(size);
    scalarField sf31(size, 1.0), sf32(size, 1.0), sf33(size, 1.0), sf34(size);

    for (register int j=0; j<nIter; j++)
    {
        sf14 = sf11 + sf12 - sf13;
        sf24 = sf21 + sf22 - sf23;
        sf34 = sf31 + sf32 - sf33;
    }

    Info<< "ExecutionTime = "
        << executionTime4.elapsedCpuTime()
        << " s\n" << endl;

    Info << sf14[1] << sf24[1] << sf34[1] << endl << endl;


}
