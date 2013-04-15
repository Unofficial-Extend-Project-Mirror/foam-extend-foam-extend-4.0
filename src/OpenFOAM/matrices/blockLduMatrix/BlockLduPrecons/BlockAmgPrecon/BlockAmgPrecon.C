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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    BlockAmgPrecon

Description
    AMG preconditioner for BlockLduMatrix

Author
    Klas Jareteg, 2013-04-15

\*---------------------------------------------------------------------------*/

#include "BlockAmgPrecon.H"
#include "fineBlockAmgLevel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
Foam::BlockAmgPrecon<Type>::BlockAmgPrecon
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduPrecon<Type>
    (
        matrix
    ),
    cycle_(BlockAmgCycle<Type>::cycleNames_.read(dict.lookup("cycle"))),
    nPreSweeps_(readLabel(dict.lookup("nPreSweeps"))),
    nPostSweeps_(readLabel(dict.lookup("nPostSweeps"))),
    nMaxLevels_(readLabel(dict.lookup("nMaxLevels"))),
    scale_(dict.lookup("scale")),
    amgPtr_
    (
        new BlockAmgCycle<Type>
        (
            autoPtr<BlockAmgLevel<Type> >
            (
                new fineBlockAmgLevel<Type>
                (
                    matrix,
                    dict,
                    dict.lookup("policy"),
                    readLabel(dict.lookup("groupSize")),
                    readLabel(dict.lookup("minCoarseEqns")),
                    dict.lookup("smoother")
                )
            )
        )
    ),
    xBuffer_(matrix.lduAddr().size())
{
    // Make coarse levels
    amgPtr_->makeCoarseLevels(nMaxLevels_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockAmgPrecon<Type>::~BlockAmgPrecon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::BlockAmgPrecon<Type>::nLevels() const
{
    return amgPtr_->nLevels();
}


template<class Type>
const Foam::Field<Type>& Foam::BlockAmgPrecon<Type>::residual
(
    const Field<Type>& x,
    const Field<Type>& b
) const
{
    // Calculate residual
    amgPtr_->residual(x, b, xBuffer_);

    return xBuffer_;
}


template<class Type>
void Foam::BlockAmgPrecon<Type>::cycle
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    amgPtr_->fixedCycle
    (
        x,
        b,
        xBuffer_,
        cycle_,
        nPreSweeps_,
        nPostSweeps_,
        scale_
    );
}

template<class Type>
void Foam::BlockAmgPrecon<Type>::precondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    // Execute preconditioning
    residual(x, b);
    cycle(x, b);
}


// ************************************************************************* //
