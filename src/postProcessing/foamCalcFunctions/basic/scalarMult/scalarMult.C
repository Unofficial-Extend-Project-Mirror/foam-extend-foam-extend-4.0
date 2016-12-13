/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "scalarMult.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace calcTypes
    {
        defineTypeNameAndDebug(scalarMult, 0);
        addToRunTimeSelectionTable(calcType, scalarMult, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::calcTypes::scalarMult::writeScalarMultValues
(
    const Time& runTime,
    const fvMesh& mesh,
    const IOobject& baseFieldHeader
)
{
    bool processed = false;

    writeScalarMultValue<scalar>
    (
        baseFieldHeader,
        scalarMultValueStr_,
        mesh,
        processed
    );
    writeScalarMultValue<vector>
    (
        baseFieldHeader,
        scalarMultValueStr_,
        mesh,
        processed
    );
    writeScalarMultValue<sphericalTensor>
    (
        baseFieldHeader,
        scalarMultValueStr_,
        mesh,
        processed
    );
    writeScalarMultValue<symmTensor>
    (
        baseFieldHeader,
        scalarMultValueStr_,
        mesh,
        processed
    );
    writeScalarMultValue<tensor>
    (
        baseFieldHeader,
        scalarMultValueStr_,
        mesh,
        processed
    );

    if (!processed)
    {
        FatalErrorIn("calcTypes::scalarMult::writeScalarMultValue()")
            << "Unable to process " << baseFieldName_
            << " + " << scalarMultValueStr_ << nl
            << "No call to scalarMult for fields of type "
            << baseFieldHeader.headerClassName() << nl << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::scalarMult::scalarMult()
:
    calcType(),
    baseFieldName_(""),
    scalarMultValueStr_(""),
    resultName_("")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::scalarMult::~scalarMult()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::scalarMult::init()
{
    argList::validArgs.append("scalarMult");
    argList::validArgs.append("baseField");
    argList::validOptions.insert("value", "valueString");
    argList::validOptions.insert("resultName", "fieldName");
}


void Foam::calcTypes::scalarMult::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    baseFieldName_ = args.additionalArgs()[1];

    if (args.optionFound("value"))
    {
        scalarMultValueStr_ = args.option("value");
    }
    else
    {
        FatalErrorIn("calcTypes::scalarMult::preCalc")
            << "scalarMult requires -value option"
            << nl << exit(FatalError);
    }

    if (args.optionFound("resultName"))
    {
        resultName_ = args.option("resultName");
    }
}


void Foam::calcTypes::scalarMult::calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    IOobject baseFieldHeader
    (
        baseFieldName_,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (baseFieldHeader.headerOk())
    {
        writeScalarMultValues(runTime, mesh, baseFieldHeader);
    }
    else
    {
        FatalErrorIn("calcTypes::scalarMult::calc")
            << "Unable to read base field: " << baseFieldName_
            << nl << exit(FatalError);
    }
}


// ************************************************************************* //
