/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    tetPolyPatchInterpolationFaceDecomp

Description
    Interpolates fields defined on faces into points on a tetPolyPatch.

SourceFiles
    tetPolyPatchInterpolationFaceDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyPatchInterpolationFaceDecomp_H
#define tetPolyPatchInterpolationFaceDecomp_H

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class faceTetPolyPatchFaceDecomp;

/*---------------------------------------------------------------------------*\
              Class tetPolyPatchInterpolationFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyPatchInterpolationFaceDecomp
{
    // Private data

        //- Reference to patch
        const faceTetPolyPatchFaceDecomp& patch_;

        //- Primitive interpolation
        primitivePatchInterpolation interpolator_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyPatchInterpolationFaceDecomp
        (
            const tetPolyPatchInterpolationFaceDecomp&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetPolyPatchInterpolationFaceDecomp&);


public:

    // Constructors

        //- Construct from patch
        tetPolyPatchInterpolationFaceDecomp
        (
            const faceTetPolyPatchFaceDecomp& p
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
#   include "tetPolyPatchInterpolateFaceDecomp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
