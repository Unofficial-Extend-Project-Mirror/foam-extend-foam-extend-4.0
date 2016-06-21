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
    largeStrainULSolidInterface

Description
    solidInterface procedure for large strain updated Lagrangian.
    no non-orthogonal correction.
    Interface tangential gradients are extrapolted from adjacent
    cell centres

SourceFiles
    largeStrainULSolidInterface.C

Author
    Philip Cardiff UCD
    Original procedure by Tukovic et al 2012
\*---------------------------------------------------------------------------*/

#ifndef largeStrainULSolidInterface_H
#define largeStrainULSolidInterface_H

#include "solidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class largeStrainULSolidInterface Declaration
\*---------------------------------------------------------------------------*/

class largeStrainULSolidInterface
:
    public solidInterface
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        largeStrainULSolidInterface(const largeStrainULSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const largeStrainULSolidInterface&);

public:

    //- Runtime type information
    TypeName("largeStrainUL");

    // Static data members


    // Constructors

        //- Construct from dictionary
        largeStrainULSolidInterface
        (
            const word& name,
            const fvMesh& mesh,
            const constitutiveModel& rheology
        );


    // Destructor

        virtual ~largeStrainULSolidInterface();


    // Member Functions

        // correct UEqn
        virtual void correct(fvVectorMatrix& UEqn);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
