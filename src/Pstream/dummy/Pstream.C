/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "Pstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Pstream::addValidParOptions(HashTable<string>& validParOptions)
{}


bool Foam::Pstream::init(int& argc, char**& argv)
{
    FatalErrorIn("Pstream::init(int& argc, char**& argv)")
        << "Trying to use the dummy Pstream library." << nl
        << "This dummy library cannot be used in parallel mode"
        << Foam::exit(FatalError);

    return false;
}


void Foam::Pstream::exit(int errnum)
{
    notImplemented("Pstream::exit(int errnum)");
}


void Foam::Pstream::abort()
{
    notImplemented("Pstream::abort()");
}


void Foam::reduce(scalar&, const sumOp<scalar>&)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
