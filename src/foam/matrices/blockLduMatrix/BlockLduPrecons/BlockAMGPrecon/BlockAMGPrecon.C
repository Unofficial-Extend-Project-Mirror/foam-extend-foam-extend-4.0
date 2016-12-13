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

Class
    BlockAMGPrecon

Description
    AMG preconditioner for BlockLduMatrix

Author
    Klas Jareteg, 2013-04-15

\*---------------------------------------------------------------------------*/

#include "BlockAMGPrecon.H"
#include "fineBlockAMGLevel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockAMGPrecon<Type>::BlockAMGPrecon
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduPrecon<Type>
    (
        matrix
    ),
    cycle_(BlockAMGCycle<Type>::cycleNames_.read(dict.lookup("cycle"))),
    nPreSweeps_(readLabel(dict.lookup("nPreSweeps"))),
    nPostSweeps_(readLabel(dict.lookup("nPostSweeps"))),
    nMaxLevels_(readLabel(dict.lookup("nMaxLevels"))),
    scale_(dict.lookup("scale")),
    amgPtr_
    (
        new BlockAMGCycle<Type>
        (
            autoPtr<BlockAMGLevel<Type> >
            (
                new fineBlockAMGLevel<Type>
                (
                    matrix,
                    dict,
                    dict.lookup("coarseningType"),
                    readLabel(dict.lookup("groupSize")),
                    readLabel(dict.lookup("minCoarseEqns"))
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
Foam::BlockAMGPrecon<Type>::~BlockAMGPrecon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::BlockAMGPrecon<Type>::nLevels() const
{
    return amgPtr_->nLevels();
}


template<class Type>
const Foam::Field<Type>& Foam::BlockAMGPrecon<Type>::residual
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
void Foam::BlockAMGPrecon<Type>::cycle
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
void Foam::BlockAMGPrecon<Type>::precondition
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
