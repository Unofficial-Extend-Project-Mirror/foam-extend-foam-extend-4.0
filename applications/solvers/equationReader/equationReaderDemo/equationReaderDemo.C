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
    equationReaderDemo

Description
    Sample application testing the equationReader extension, and demonstrating
    its use.
    
Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "IFstream.H"
#include "OFstream.H"
#include "equationReader.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
#   include "setRootCase.H"

    fileName path(args.rootPath()/args.caseName());

    // Create dictionary
    Info << "Reading testDict dictionary" << token::NL << endl;
    IFstream tfIF(path/"testDict");
    const dictionary testDict(tfIF);

    // Demonstrate stand-alone operation
    Info << "Begin stand-alone operation..." << endl;

    Info << "Reading a scalar... ";
    scalar readSa(readScalar(testDict.lookup("standAloneScalar")));
    Info << "done.  Result = " << readSa << endl;
    Info << "Reading a dimensionedScalar ... ";
    dimensionedScalar readDSa(testDict.lookup("standAloneDScalar"));
    Info << "done.  Result = " << readDSa << token::NL << endl;

    // Create the equationReader object
    Info << "Creating the equationReader object" << token::NL << endl;
    equationReader eqns;

    // Demonstrate giving data sources to equationReader
    // -First create the data sources
    Info << "Creating data sources: dictionary ptrs... ";
    IFstream tfIF2(path/"testDict2");
    const dictionary testDict2(tfIF2);
    IFstream tfIF3(path/"testDict3");
    const dictionary testDict3(tfIF3);
    Info << "scalars... ";
    scalar Sa(0.1);
    scalar Sb(0.2);
    scalar Sc(0.3);
    Info << "dimensionedScalar ptrs... ";
    dimensionedScalar DSa("DSa", dimless, 1);
    dimensionedScalar DSb("DSb", dimless, 2);
    dimensionedScalar DSc("DSc", dimless, 3);
    Info << "output dimensionedScalar ptrs... ";
    dimensionedScalar passiveOutA("passiveOutA", dimless, 0);
    dimensionedScalar activeOutB("activeOutB", dimless, 0);
    dimensionedScalar passiveOutC("passiveOutC", dimless, 0);
    dimensionedScalar passiveOutD("passiveOutD", dimless, 0);
    dimensionedScalar passiveOutF("passiveOutF", dimless, 0);
    Info << "done." << endl;

    Info << "Linking in the data sources: dictionary ptrs... ";
    eqns.addDataSource(testDict);
    eqns.addDataSource(testDict2);
    eqns.addDataSource(testDict3);

    Info << "dimensionedScalar ptrs... ";
    eqns.addDataSource(DSa);
    eqns.addDataSource(DSb);
    eqns.addDataSource(DSc);

    Info << "scalar ptrs... ";
    eqns.addDataSource(Sa, "Sa");
    eqns.addDataSource(Sb, "Sb");
    eqns.addDataSource(Sc, "Sc");
    Info << "done." << token::NL << endl;

    // Demonstrate passive output
    Info << "Reading equation a from testDict with no output variable" << endl;
    eqns.readEquation(testDict, "a");

    Info << "Evaluating equation a ... ";
    passiveOutA = eqns.evaluate("a");
    Info << "done.  Result = " << passiveOutA << token::NL << endl;

    // Demonstrate active output
    Info << "Reading equation b from testDict, linking an output variable"
        << endl;
    eqns.readEquation(testDict, "b", activeOutB);
    
    Info << "Output variable before update() = " << activeOutB << endl;
    Info << "Begining .update() - this evaluates all equations with active "
        << "output..." << endl;
    eqns.update();
    Info << "Done.  Output variable after update() = " << activeOutB
        << token::NL << endl;

    // Demonstrating variable dependence
    Info << "Equation c depends on equation a.  Reading it will link them."
        << endl;
    Info << "Reading equation c from testDict... ";
    eqns.readEquation(testDict, "c");
    Info << "done." << endl;
    Info << "Evaluating c will force an evaluate of a." << endl;
    Info << "Evaluating ... ";
    passiveOutC = eqns.evaluate("c");
    Info << "done.  Result = " << passiveOutC << token::NL << endl;

    // Demonstrate on-the-fly equation creation
    Info << "Equation d depends on equation e, but equation e is never "
        << "explicitly " << endl << "read by equationReaderDemo.  Reading equation"
        << " d will automatically " << endl << "create equation e on-the-fly. "
        << endl;
    Info << "Reading equation d from testDict ... ";
    eqns.readEquation(testDict, "d");
    Info << "done." << endl << "Again, evaluating d will force an evaluate of "
        << "e." << endl;
    Info << "Evaluating d ... ";
    passiveOutD = eqns.evaluate("d");
    Info << "done.  The result is = " << passiveOutD << token::NL << endl;

    // Demonstrate dependence
    Info << "Equations can draw from any sources added to equationReader." << endl;
    Info << "Equation f is very complex, drawing from numerous sources." << endl;
    Info << "Reading equation f ... ";
    eqns.readEquation(testDict, "f");
    Info << "done.  Evaluating equation f ... ";
    passiveOutF = eqns.evaluate("f");
    Info << "done." << token::NL << "The result is: " << passiveOutF << endl;
    
    Info << token::NL << "Creating output..." << endl;
    OFstream os(path/"outputDict");
    os << eqns;
    eqns.dataSourceStatus(os);
    
    return(0);
}
