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

Description
    Constructor of cellModeller: just sets the cellModeller's params.

\*---------------------------------------------------------------------------*/

#include "cellModeller.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

cellModeller::~cellModeller()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Returns a pointer to a model which matches the string symbol
// supplied. A null pointer is returned if there is no suitable match.

const cellModel* cellModeller::lookup(const word& symbol)
{
    HashTable<const cellModel*>::iterator iter = modelDictionary_.find(symbol);

    if (iter != modelDictionary_.end())
    {
        return iter();
    }
    else
    {
        return NULL;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
