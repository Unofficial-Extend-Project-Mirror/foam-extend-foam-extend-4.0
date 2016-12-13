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

#include "calcType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(calcType, 0);

defineRunTimeSelectionTable(calcType, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcType::calcType()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcType::~calcType()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::calcType::init()
{
    // Do nothing
}


void Foam::calcType::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    // Do nothing
}


void Foam::calcType::calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    // Do nothing
}


void Foam::calcType::postCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    // Do nothing
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcType::tryInit()
{
    FatalIOError.throwExceptions();

    try
    {
        init();
    }
    catch(IOerror& err)
    {
        Warning<< err << endl;
    }
}


void Foam::calcType::tryPreCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    FatalIOError.throwExceptions();

    try
    {
        preCalc(args, runTime, mesh);
    }
    catch(IOerror& err)
    {
        Warning<< err << endl;
    }
}


void Foam::calcType::tryCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    FatalIOError.throwExceptions();

    try
    {
        calc(args, runTime, mesh);
    }
    catch(IOerror& err)
    {
        Warning<< err << endl;
    }
}


void Foam::calcType::tryPostCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    FatalIOError.throwExceptions();

    try
    {
        postCalc(args, runTime, mesh);
    }
    catch(IOerror& err)
    {
        Warning<< err << endl;
    }
}


// ************************************************************************* //
