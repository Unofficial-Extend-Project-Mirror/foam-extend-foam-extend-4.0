/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Application
    spatialMeanValue

Description
    Calculates the spatial mean of the specified volScalarField over the
    whole domain (for regular grid).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"

    timeSelector::addOptions();
    argList::validArgs.append("fieldName");

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedMesh.H"

    word fieldName(args.additionalArgs()[0]);

    Info << "Calculating avarage pressure heads:" << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check field exists
        if (fieldHeader.headerOk())
        {
            mesh.readUpdate();

            // Read field and calc mean, for regular grid
            if (fieldHeader.headerClassName() == volScalarField::typeName)
            {
                volScalarField field(fieldHeader, mesh);

                int nbMesh;
                nbMesh = 0;

                forAll(field, cellI)
                {
                    nbMesh++;
                }

                Info<< runTime.timeName()<< " "
                    << sum(field).value()/nbMesh<< " "
                    << endl;
            }
            else
            {
                FatalError
                    << "Only possible for volScalarField "<< "s "
                    << nl << exit(FatalError);
            }
        }
        else
        {
            Info<< "    No field " << fieldName << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
