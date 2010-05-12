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
    tetPolyPatchInterpolationCellDecomp

Description
    Interpolates fields defined on faces into points on a tetPolyPatch.

SourceFiles
    tetPolyPatchInterpolationCellDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyPatchInterpolationCellDecomp_H
#define tetPolyPatchInterpolationCellDecomp_H

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class faceTetPolyPatchCellDecomp;

/*---------------------------------------------------------------------------*\
              Class tetPolyPatchInterpolationCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyPatchInterpolationCellDecomp
{
    // Private data

        //- Reference to patch
        const faceTetPolyPatchCellDecomp& patch_;

        //- Primitive interpolation
        primitivePatchInterpolation interpolator_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyPatchInterpolationCellDecomp
        (
            const tetPolyPatchInterpolationCellDecomp&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetPolyPatchInterpolationCellDecomp&);


public:

    // Constructors

        //- Construct from patch
        tetPolyPatchInterpolationCellDecomp
        (
            const faceTetPolyPatchCellDecomp& p
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
#   include "tetPolyPatchInterpolateCellDecomp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
