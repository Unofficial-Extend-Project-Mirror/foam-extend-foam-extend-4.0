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

Description
    Class for BlockCoeff norm based on a specified component.

Author
    Klas Jareteg, 2013-03-07

SourceFiles
    BlockCoeffComponentNorm.C

\*---------------------------------------------------------------------------*/

#ifndef BlockCoeffComponentNorm_H
#define BlockCoeffComponentNorm_H

#include "blockCoeffs.H"
#include "blockCoeffNorms.H"
#include "BlockCoeffNorm.H"
#include "CoeffField.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class BlockCoeffComponentNorm Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockCoeffComponentNorm
:
    public BlockCoeffNorm<Type>
{
    // Private data

        //- Control data dictionary
        dictionary dict_;
        
        //- Component from which the norm is determined
        const direction cmpt_;


protected:

    // Protected Member Functions

        //- Return dictionary
        const dictionary& dict() const
        {
            return dict_;
        }


public:

    //- Runtime type information
    TypeName("componentNorm");


    // Constructors

        //- Construct from dictionary
        BlockCoeffComponentNorm
        (
            const dictionary& dict
        );
        
    // Destructor

        virtual ~BlockCoeffComponentNorm()
        {}


    // Member functions

        //- Solve
        virtual scalar normalize
        (
            const BlockCoeff<Type>& a
        );
        
        virtual void coeffMag
        (
            const CoeffField<Type>& a,
            Field<scalar>& b
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockCoeffComponentNorm.C"
#endif
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
