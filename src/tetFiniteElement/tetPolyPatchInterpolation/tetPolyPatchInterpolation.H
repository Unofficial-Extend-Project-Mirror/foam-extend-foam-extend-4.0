/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    tetPolyPatchInterpolation

Description
    Interpolates fields defined on faces into points on a tetPolyPatch.

SourceFiles
    tetPolyPatchInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyPatchInterpolation_H
#define tetPolyPatchInterpolation_H

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class faceTetPolyPatch;

/*---------------------------------------------------------------------------*\
              Class tetPolyPatchInterpolation Declaration
\*---------------------------------------------------------------------------*/

class tetPolyPatchInterpolation
{
    // Private data

        //- Reference to patch
        const faceTetPolyPatch& patch_;

        //- Primitive interpolation
        primitivePatchInterpolation interpolator_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyPatchInterpolation
        (
            const tetPolyPatchInterpolation&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetPolyPatchInterpolation&);


public:

    // Constructors

        //- Construct from patch
        tetPolyPatchInterpolation
        (
            const faceTetPolyPatch& p
        );


    // Destructor - default


    // Member Functions

        //- Interpolate from faces to points
        template<class Type>
        tmp<Field<Type> > faceToPointInterpolate
        (
            const Field<Type>& ff
        ) const;

        template<class Type>
        tmp<Field<Type> > faceToPointInterpolate
        (
            const tmp<Field<Type> >& tff
        ) const;


        //- Interpolate from mesh points to tet FEM points
        template<class Type>
        tmp<Field<Type> > pointToPointInterpolate
        (
            const Field<Type>& ff
        ) const;

        template<class Type>
        tmp<Field<Type> > pointToPointInterpolate
        (
            const tmp<Field<Type> >& tff
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetPolyPatchInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
