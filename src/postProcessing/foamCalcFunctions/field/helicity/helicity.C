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

\*---------------------------------------------------------------------------*/

#include "helicity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace calcTypes
    {
        defineTypeNameAndDebug(helicity, 0);
        addToRunTimeSelectionTable(calcType, helicity, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::helicity::helicity()
:
    calcType()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::helicity::~helicity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::helicity::init()
{
    argList::validArgs.append("helicity");
    argList::validArgs.append("fieldName");
}


void Foam::calcTypes::helicity::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{}


void Foam::calcTypes::helicity::calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    const word& fieldName = args.additionalArgs()[1];

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
        bool processed = false;

        writeHelicityField(fieldHeader, mesh, processed);

        if (!processed)
        {
            FatalError
                << "Unable to process " << fieldName << nl
                << "No call to helicity for fields of type "
                << fieldHeader.headerClassName() << nl << nl
                << exit(FatalError);
        }
    }
    else
    {
        Info<< "    No " << fieldName << endl;
    }
}


void Foam::calcTypes::helicity::writeHelicityField
(
    const IOobject& header,
    const fvMesh& mesh,
    bool& processed
)
{
    if (header.headerClassName() == volVectorField::typeName)
    {
        Info<< "    Reading " << header.name() << endl;
        volVectorField field(header, mesh);

        Info<< "    Calculating helicity" << header.name() << endl;
        volScalarField helicityField
        (
            IOobject
            (
                "helicity" + header.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            field & (fvc::curl(field))
        );

        Info << "helicity(" << header.name() << "): max: "
            << gMax(helicityField.internalField())
            << " min: " << gMin(helicityField.internalField()) << endl;

        helicityField.write();

        processed = true;
    }
}


// ************************************************************************* //

