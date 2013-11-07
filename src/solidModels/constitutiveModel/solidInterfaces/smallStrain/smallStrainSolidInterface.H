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
    smallStrainSolidInterface

Description
    solidInterface procedure for small strain.
    no non-orthogonal correction.
    Interface tangential gradients are extrapolted from adjacent
    cell centres

SourceFiles
    smallStrainSolidInterface.C

Author
    Philip Cardiff UCD
    Original procedure by Tukovic et al 2012
\*---------------------------------------------------------------------------*/

#ifndef smallStrainSolidInterface_H
#define smallStrainSolidInterface_H

#include "solidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class smallStrainSolidInterface Declaration
\*---------------------------------------------------------------------------*/

class smallStrainSolidInterface
:
    public solidInterface
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        smallStrainSolidInterface(const smallStrainSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const smallStrainSolidInterface&);

public:

    //- Runtime type information
    TypeName("smallStrain");

    // Static data members


    // Constructors

        //- Construct from dictionary
        smallStrainSolidInterface
        (
            const word& name,
            const fvMesh& mesh,
            const constitutiveModel& rheology
        );


    // Destructor

        virtual ~smallStrainSolidInterface();


    // Member Functions

        // correct UEqn
        virtual void correct(fvVectorMatrix& UEqn);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
