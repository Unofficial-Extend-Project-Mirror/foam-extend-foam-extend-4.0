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
    BlockCoeffNorm

Description
    Abstract class for calculating norm of a BlockCoeff

Author
    Klas Jareteg, 2013-01-28

SourceFiles
    BlockCoeffNorm.C

\*---------------------------------------------------------------------------*/

#ifndef BlockCoeffNorm_H
#define BlockCoeffNorm_H

#include "blockCoeffs.H"
#include "CoeffField.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class BlockCoeffNorm Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockCoeffNorm
{
    // Private data

        //- Control data dictionary
        dictionary dict_;


protected:

    // Protected Member Functions

        //- Return dictionary
        const dictionary& dict() const
        {
            return dict_;
        }


public:

    //- Runtime type information
    TypeName("BlockCoeffNorm");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockCoeffNorm,
            dictionary,
            (
                const dictionary& dict
            ),
            (
                dict
            )
        );


    // Constructors

        //- Construct from dictionary
        BlockCoeffNorm
        (
            const dictionary& dict
        );


    // Selectors

        //- Return a new norm calculator
        static autoPtr<BlockCoeffNorm<Type> > New
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~BlockCoeffNorm()
        {}


    // Member functions

        //- Solve
        virtual scalar normalize
        (
            const BlockCoeff<Type>& a
        ) = 0;

        virtual void coeffMag
        (
            const CoeffField<Type>& a,
            Field<scalar>& b
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockCoeffNorm.C"
#endif

#define makeBlockCoeffNorm(BlockCoeffNormType, typeBlockCoeffNormType)        \
                                                                              \
defineNamedTemplateTypeNameAndDebug(typeBlockCoeffNormType, 0);               \
                                                                              \
addToRunTimeSelectionTable(BlockCoeffNormType, typeBlockCoeffNormType, dictionary);

#define makeBlockCoeffNorms(blockCoeffNormType)                               \
                                                                              \
makeBlockCoeffNorm(blockScalarCoeffNorm, blockCoeffNormType##Scalar);         \
makeBlockCoeffNorm(blockVectorCoeffNorm, blockCoeffNormType##Vector);         \
makeBlockCoeffNorm(blockTensorCoeffNorm, blockCoeffNormType##Tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
