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
    largeStrainULCorrectedSolidInterface

Description
    solidInterface procedure for large strain updated Lagrangian.
    with non-orthogonal correction.
    Interface tangential gradients are extrapolted from adjacent
    cell centres

SourceFiles
    largeStrainULCorrectedSolidInterface.C

Author
    Philip Cardiff UCD
    Original procedure by Tukovic et al 2012
\*---------------------------------------------------------------------------*/

#ifndef largeStrainULCorrectedSolidInterface_H
#define largeStrainULCorrectedSolidInterface_H

#include "solidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class largeStrainULCorrectedSolidInterface Declaration
\*---------------------------------------------------------------------------*/

class largeStrainULCorrectedSolidInterface
:
    public solidInterface
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        largeStrainULCorrectedSolidInterface(const largeStrainULCorrectedSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const largeStrainULCorrectedSolidInterface&);

public:

    //- Runtime type information
    TypeName("largeStrainULCorrected");

    // Static data members


    // Constructors

        //- Construct from dictionary
        largeStrainULCorrectedSolidInterface
        (
	 const word& name,
	 const fvMesh& mesh,
	 const constitutiveModel& rheology
        );


    // Destructor

        virtual ~largeStrainULCorrectedSolidInterface();


    // Member Functions

        // correct UEqn
        virtual void correct(fvVectorMatrix& UEqn);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
