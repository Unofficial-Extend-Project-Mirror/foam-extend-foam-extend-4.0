/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Converts CHEMKINIII thermodynamics and reaction data files into FOAM format

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "chemkinReader.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.clear();
    argList::validArgs.append("CHEMKINFile");
    argList::validArgs.append("CHEMKINThermodynamicsFile");
    argList::validArgs.append("FOAMChemistryFile");
    argList::validArgs.append("FOAMThermodynamicsFile");
    argList args(argc, argv);

    fileName CHEMKINFileName(args.additionalArgs()[0]);
    fileName thermoFileName(args.additionalArgs()[1]);
    fileName FOAMChemistryFileName(args.additionalArgs()[2]);
    fileName FOAMThermodynamicsFileName(args.additionalArgs()[3]);

    chemkinReader cr(CHEMKINFileName, thermoFileName);

    OFstream reactionsFile(FOAMChemistryFileName);
    reactionsFile
        << "species" << cr.species() << ';' << endl << endl
        << "reactions" << cr.reactions() << ';' << endl;

    OFstream thermoFile(FOAMThermodynamicsFileName);
    thermoFile<< cr.sThermo() << endl;

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
