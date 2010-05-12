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

\*---------------------------------------------------------------------------*/

#include "solution.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int Foam::solution::debug(Foam::debug::debugSwitch("solution", false));


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solution::solution(const objectRegistry& obr, const fileName& dictName)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            IOobject::READ_IF_PRESENT,  // Allow default dictionary creation
            IOobject::NO_WRITE
        )
    ),
    relaxationFactors_(ITstream("relaxationFactors", tokenList())()),
    solvers_(ITstream("solvers", tokenList())())
{
    if (!headerOk())
    {
        if (debug)
        {
            InfoIn
            (
                "Foam::solution::solution(const objectRegistry& obr, "
                "const fileName& dictName)"
            )   << "Solution dictionary not found.  Adding default entries"
                << endl;
        }
    }

    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solution::read()
{
    bool readOk = false;

    if (headerOk())
    {
        readOk = regIOobject::read();
    }

    if (readOk)
    {
        const dictionary& dict = solutionDict();

        if (dict.found("relaxationFactors"))
        {
            relaxationFactors_ = dict.subDict("relaxationFactors");
        }

        if (dict.found("solvers"))
        {
            solvers_ = dict.subDict("solvers");
        }

        return true;
    }

    return readOk;
}


const Foam::dictionary& Foam::solution::solutionDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }
    else
    {
        return *this;
    }
}


bool Foam::solution::relax(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup relax for " << name << endl;
    }

    return relaxationFactors_.found(name);
}


Foam::scalar Foam::solution::relaxationFactor(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup relaxationFactor for " << name << endl;
    }

    return readScalar(relaxationFactors_.lookup(name));
}


const Foam::dictionary& Foam::solution::solverDict(const word& name) const
{
    if (debug)
    {
        InfoIn("solution::solverDict(const word& name)")
            << "Lookup solver for " << name << endl;
    }

    return solvers_.subDict(name);
}


Foam::ITstream& Foam::solution::solver(const word& name) const
{
    if (debug)
    {
        InfoIn("solution::solver(const word& name)")
            << "Lookup solver for " << name << endl;
    }

    return solvers_.lookup(name);
}


bool Foam::solution::writeData(Ostream& os) const
{
    // Write dictionaries
    os << nl << "solvers";
    solvers_.write(os, true);

    os << nl << "relaxationFactors";
    relaxationFactors_.write(os, true);

    // Write direct entries of the solution dictionary
    // HJ, 16/Feb/2010
    dictionary::write(os, false);

    return true;
}


// ************************************************************************* //
