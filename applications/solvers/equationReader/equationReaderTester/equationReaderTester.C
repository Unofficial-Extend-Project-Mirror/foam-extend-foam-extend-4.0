/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    equationReaderTester

Description
    Sample application testing the equationReader in a finite volume solver
    environment.
    
Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "IFstream.H"
#include "IOEquationReader.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fileName path(args.rootPath()/args.caseName());

    // Create the equationReader object
    Info << "Creating the equationReader object" << token::NL << endl;
    IOEquationReader eqns
    (
        IOobject
        (
            "eqns",
            runTime.timeName(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        true
    );

    // Create dictionaries
    Info << "Reading testDict dictionaries" << token::NL << endl;
    IFstream tfIF1(path/"testDict1");
    const dictionary testDict1(tfIF1);
    IFstream tfIF2(path/"testDict2");
    const dictionary testDict2(tfIF2);
    IFstream tfIF3(path/"testDict3");
    const dictionary testDict3(tfIF3);
    IFstream tfIF4(path/"testDict4");
    const dictionary testDict4(tfIF4);
    IFstream tfIF5(path/"testDict5");
    const dictionary testDict5(tfIF5);
    IFstream tfIF6(path/"testDict6");
    const dictionary testDict6(tfIF6);
    IFstream tfIF7(path/"testDict7");
    const dictionary testDict7(tfIF7);
    IFstream tfIF8(path/"testDict8");
    const dictionary testDict8(tfIF8);
    IFstream tfIF9(path/"testDict9");
    const dictionary testDict9(tfIF9);
    IFstream tfIF10(path/"testDict10");
    const dictionary testDict10(tfIF10);
    scalar Sa(0.1);
    scalar Sb(0.2);
    scalar Sc(0.3);
    scalar Sd(0.4);
    scalar Se(0.5);
    scalar Sf(0.6);
    scalar Sg(0.7);
    scalar Sh(0.8);
    scalar Si(0.9);
    scalar Sj(0.10);
    scalar Sk(0.11);
    scalar Sl(0.12);
    dimensionedScalar DSa("DSa", dimless, 1);
    dimensionedScalar DSb("DSb", dimless, 2);
    dimensionedScalar DSc("DSc", dimless, 3);
    dimensionedScalar DSd("DSd", dimless, 4);
    dimensionedScalar DSe("DSe", dimless, 5);
    dimensionedScalar DSf("DSf", dimless, 6);
    dimensionedScalar DSg("DSg", dimless, 7);
    dimensionedScalar DSh("DSh", dimless, 8);
    dimensionedScalar DSi("DSi", dimless, 9);
    dimensionedScalar DSj("DSj", dimless, 10);
    dimensionedScalar DSk("DSk", dimless, 11);
    dimensionedScalar DSl("DSl", dimless, 12);
    dimensionedScalar DSm("DSm", dimless, 13);
    dimensionedScalar DSn("DSn", dimless, 14);
    dimensionedScalar DSo("DSo", dimless, 15);
    dimensionedScalar DStime("DStime", dimless, 0);
    dimensionedScalar Aa("Aa", dimless, 0);
    dimensionedScalar Ab("Ab", dimless, 0);
    dimensionedScalar Ac("Ac", dimless, 0);
    dimensionedScalar Ad("Ad", dimless, 0);
    dimensionedScalar Ae("Ae", dimless, 0);
    dimensionedScalar Af("Af", dimless, 0);
    dimensionedScalar Pa("Pa", dimless, 0);
    dimensionedScalar Pb("Pb", dimless, 0);
    dimensionedScalar Pc("Pc", dimless, 0);
    dimensionedScalar Pd("Pd", dimless, 0);
    dimensionedScalar Pe("Pe", dimless, 0);
    dimensionedScalar Pf("Pf", dimless, 0);

    eqns.addDataSource(testDict1);
    eqns.addDataSource(testDict2);
    eqns.addDataSource(testDict3);
    eqns.addDataSource(testDict4);
    eqns.addDataSource(testDict5);
    eqns.addDataSource(testDict6);
    eqns.addDataSource(testDict7);
    eqns.addDataSource(testDict8);
    eqns.addDataSource(testDict9);
    eqns.addDataSource(testDict10);

    eqns.addDataSource(DSa);
    eqns.addDataSource(DSb);
    eqns.addDataSource(DSc);
    eqns.addDataSource(DSd);
    eqns.addDataSource(DSe);
    eqns.addDataSource(DSf);
    eqns.addDataSource(DSg);
    eqns.addDataSource(DSh);
    eqns.addDataSource(DSi);
    eqns.addDataSource(DSj);
    eqns.addDataSource(DSk);
    eqns.addDataSource(DSl);
    eqns.addDataSource(DSm);
    eqns.addDataSource(DSn);
    eqns.addDataSource(DSo);
    eqns.addDataSource(DStime);

    eqns.addDataSource(Sa, "Sa");
    eqns.addDataSource(Sb, "Sb");
    eqns.addDataSource(Sc, "Sc");
    eqns.addDataSource(Sd, "Se");
    eqns.addDataSource(Se, "Sd");
    eqns.addDataSource(Sf, "Sf");
    eqns.addDataSource(Sg, "Sg");
    eqns.addDataSource(Sh, "Sh");
    eqns.addDataSource(Si, "Si");
    eqns.addDataSource(Sj, "Sj");
    eqns.addDataSource(Sk, "Sk");
    eqns.addDataSource(Sl, "Sl");
    
    label listIndex(0);
    eqns.addDataSource(p);
    eqns.addDataSource(dummy);
    eqns.setListIndex(listIndex);

    eqns.readEquation(testDict1, "Pa");
    eqns.readEquation(testDict1, "Pb");
    eqns.readEquation(testDict1, "Pc");
    eqns.readEquation(testDict1, "Pd");
    eqns.readEquation(testDict1, "Pe");
    eqns.readEquation(testDict1, "Pf");

    eqns.readEquation(testDict1, "Aa", Aa);
    eqns.readEquation(testDict1, "Ab", Ab);
    eqns.readEquation(testDict1, "Ac", Ac);
    eqns.readEquation(testDict1, "Ad", Ad);
    eqns.readEquation(testDict1, "Ae", Ae);
    eqns.readEquation(testDict1, "Af", Af);
    eqns.readEquation(testDict1, "nu", nu);

    scalar saA(readScalar(testDict1.lookup("saA")));
    scalar saB(readScalar(testDict1.lookup("saB")));
    scalar saC(readScalar(testDict1.lookup("saC")));
    scalar saD(readScalar(testDict1.lookup("saD")));
    scalar saE(readScalar(testDict1.lookup("saE")));
    scalar saF(readScalar(testDict1.lookup("saF")));
    
    dimensionedScalar dsaA(testDict1.lookup("dsaA"));
    dimensionedScalar dsaB(testDict1.lookup("dsaB"));
    dimensionedScalar dsaC(testDict1.lookup("dsaC"));
    dimensionedScalar dsaD(testDict1.lookup("dsaD"));
    dimensionedScalar dsaE(testDict1.lookup("dsaE"));
    dimensionedScalar dsaF(testDict1.lookup("dsaF"));

    Info << "Stand-alone test:" << endl;
    Info << "saA = " << saA << endl;
    Info << "saB = " << saB << endl;
    Info << "saC = " << saC << endl;
    Info << "saD = " << saD << endl;
    Info << "saE = " << saE << endl;
    Info << "saF = " << saF << endl;

    Info << "dsaA = " << dsaA << endl;
    Info << "dsaB = " << dsaB << endl;
    Info << "dsaC = " << dsaC << endl;
    Info << "dsaD = " << dsaD << endl;
    Info << "dsaE = " << dsaE << endl;
    Info << "dsaF = " << dsaF << endl;

    Info << "Pa is at index " << eqns.lookup("Pa") << endl;
    Info << "Pb is at index " << eqns.lookup("Pb") << endl;
    Info << "Pc is at index " << eqns.lookup("Pc") << endl;
    Info << "Pd is at index " << eqns.lookup("Pd") << endl;
    Info << "Pe is at index " << eqns.lookup("Pe") << endl;
    Info << "Pf is at index " << eqns.lookup("Pf") << endl;
    Info << "Aa is at index " << eqns.lookup("Aa") << endl;
    Info << "Ab is at index " << eqns.lookup("Ab") << endl;
    Info << "Ac is at index " << eqns.lookup("Ac") << endl;
    Info << "Ad is at index " << eqns.lookup("Ad") << endl;
    Info << "Ae is at index " << eqns.lookup("Ae") << endl;
    Info << "Af is at index " << eqns.lookup("Af") << endl;
    Info << "nu is at index " << eqns.lookup("nu") << endl;

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        DStime.value() = runTime.value();
        
        Info << "Moving p index to ";
        listIndex++;
        if (listIndex == p.size())
        {
            listIndex = 0;
        }
        Info << listIndex << "..." << endl;
        eqns.setListIndex(listIndex);

        Info << "Passive reread..." << endl;
        eqns.readEquation(testDict1, "Pa");
        eqns.readEquation(testDict1, "Pb");
        eqns.readEquation(testDict1, "Pc");

        Info << "Active reread..." << endl;
        eqns.readEquation(testDict1, "Aa", Aa);
        eqns.readEquation(testDict1, "Ab", Ab);
        eqns.readEquation(testDict1, "Ac", Ac);

        Info << "Updating active equations..." << endl;
        eqns.update();
        Info << "Aa = " << Aa << endl;
        Info << "Ab = " << Ab << endl;
        Info << "Ac = " << Ac << endl;
        Info << "Ad = " << Ad << endl;
        Info << "Ae = " << Ae << endl;
        Info << "Af = " << Af << endl;
        Info << "nu = " << nu << endl;

        Info << "Evaluating passive equations: Pa, ";
        Pa = eqns.evaluate("Pa");
        Info << "Pb, ";
        Pb = eqns.evaluate("Pb");
        Info << "Pc, ";
        Pc = eqns.evaluate("Pc");
        Info << "Pd, ";
        Pe = eqns.evaluate("Pd");
        Info << "Pe, ";
        Pd = eqns.evaluate("Pe");
        Info << "Pf." << endl;
        Pf = eqns.evaluate("Pf");

        Info << "Pa = " << Pa << endl;
        Info << "Pb = " << Pb << endl;
        Info << "Pc = " << Pc << endl;
        Info << "Pd = " << Pd << endl;
        Info << "Pe = " << Pe << endl;
        Info << "Pf = " << Pf << endl;

#       include "readPISOControls.H"
#       include "CourantNo.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        solve(UEqn == -fvc::grad(p));

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField rUA = 1.0/UEqn.A();

            U = rUA*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf()) 
                + fvc::ddtPhiCorr(rUA, U, phi);

            adjustPhi(phi, U, p);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

