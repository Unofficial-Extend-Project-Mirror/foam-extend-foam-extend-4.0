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

\*---------------------------------------------------------------------------*/

#include "speciesTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::speciesTable::setIndices()
{
    forAll (*this, i)
    {
        specieIndices_.insert(wordList::operator[](i), i);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from list of specie names
Foam::speciesTable::speciesTable(const wordList& specieNames)
:
    wordList(specieNames)
{
    setIndices();
}


// Construct from number of species and list of specie names
Foam::speciesTable::speciesTable(const label nSpecies, const char** specieNames)
:
    wordList(nSpecies)
{
    forAll (*this, i)
    {
        wordList::operator[](i) = specieNames[i];
    }

    setIndices();
}


// Construct from Istream
Foam::speciesTable::speciesTable(Istream& is)
:
    wordList(is)
{
    setIndices();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, speciesTable& st)
{
    is >> static_cast<wordList&>(st);
    st.setIndices();

    return is;
}

// ************************************************************************* //
