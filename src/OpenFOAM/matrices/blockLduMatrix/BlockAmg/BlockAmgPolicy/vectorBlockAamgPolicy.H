/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    BlockAamgPolicy

Description
    Specialisation of the BlockAamgPolicy for vectors. 

Author
    Klas Jareteg, 2013-01-31

\*---------------------------------------------------------------------------*/

#ifndef vectorBlockAamgPolicy_H
#define vectorBlockAamgPolicy_H

#include "blockAmgPolicies.H"
#include "blockAamgPolicies.H"
#include "BlockAamgPolicy.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class BlockAamgPolicy Declaration
\*---------------------------------------------------------------------------*/

//- Restrict matrix

template<>
inline autoPtr<BlockLduMatrix<vector> > BlockAamgPolicy<vector>::restrictMatrix
() const
{
    FatalErrorIn
    (
        "autoPtr<BlockLduMatrix<vector> > BlockAamgPolicy<Type>::"
        "restrictMatrix() const"
    )   << "Function not implemented for Type=vector. " << endl
        << abort(FatalError);

    // Dummy return to keep compiler happy
    return autoPtr<BlockLduMatrix<vector> >(NULL);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
