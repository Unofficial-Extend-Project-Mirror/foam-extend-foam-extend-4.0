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
    BlockCoeffComponentNorm

\*---------------------------------------------------------------------------*/

#include "BlockCoeffComponentNorm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockCoeffComponentNorm<Type>::BlockCoeffComponentNorm
(
    const dictionary& dict
)
:
    BlockCoeffNorm<Type>(dict),    
    dict_(dict),
    cmpt_(readInt(this->dict().lookup("normComponent")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::BlockCoeffComponentNorm<Type>::normalize
(
    const Foam::BlockCoeff<Type>& a
)
{ 
    return mag(a.component(cmpt_));
}


template<class Type>
void Foam::BlockCoeffComponentNorm<Type>::coeffMag
(
    const Foam::CoeffField<Type>& a,
    Foam::Field<scalar>& b
)
{
    b = mag(a.component(cmpt_));
}


// ************************************************************************* //
