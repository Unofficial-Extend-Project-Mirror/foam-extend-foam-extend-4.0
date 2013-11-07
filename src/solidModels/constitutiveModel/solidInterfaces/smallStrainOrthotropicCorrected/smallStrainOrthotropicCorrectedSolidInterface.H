/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    smallStrainOrthotropicCorrectedSolidInterface

Description
    solidInterface procedure for small strain.
    no non-orthogonal correction.
    Interface tangential gradients are extrapolted from adjacent
    cell centres

SourceFiles
    smallStrainOrthotropicCorrectedSolidInterface.C

Author
    Philip Cardiff UCD
    Original procedure by Tukovic et al 2012
\*---------------------------------------------------------------------------*/

#ifndef smallStrainOrthotropicCorrectedSolidInterface_H
#define smallStrainOrthotropicCorrectedSolidInterface_H

#include "solidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class smallStrainOrthotropicCorrectedSolidInterface Declaration
\*---------------------------------------------------------------------------*/

class smallStrainOrthotropicCorrectedSolidInterface
:
    public solidInterface
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        smallStrainOrthotropicCorrectedSolidInterface
        (const smallStrainOrthotropicCorrectedSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const smallStrainOrthotropicCorrectedSolidInterface&);

public:

    //- Runtime type information
    TypeName("smallStrainOrthotropicCorrected");

    // Static data members


    // Constructors

        //- Construct from dictionary
        smallStrainOrthotropicCorrectedSolidInterface
        (
            const word& name,
            const fvMesh& mesh,
            const constitutiveModel& rheology
        );


    // Destructor

        virtual ~smallStrainOrthotropicCorrectedSolidInterface();


    // Member Functions

        // correct UEqn
        virtual void correct(fvVectorMatrix& UEqn);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
