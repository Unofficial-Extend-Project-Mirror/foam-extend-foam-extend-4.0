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
    Foam::PrimitivePatchInterpolation

Description
    Interpolation class within a primitive patch. Allows interpolation from
    points to faces and vice versa

SourceFiles
    PrimitivePatchInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef PrimitivePatchInterpolation_H
#define PrimitivePatchInterpolation_H

#include "scalarList.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class PrimitivePatchInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Patch>
class PrimitivePatchInterpolation
{
    // Private data

        //- Reference to patch
        const Patch& patch_;


    // Demand-driven data

        //- Face-to-point weights
        mutable scalarListList* faceToPointWeightsPtr_;

        //- Face-to-edge weights
        mutable scalarList* faceToEdgeWeightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PrimitivePatchInterpolation(const PrimitivePatchInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const PrimitivePatchInterpolation&);


        //- Face-to-point weights
        const scalarListList& faceToPointWeights() const;

        //- Make face-to-point weights
        void makeFaceToPointWeights() const;

        //- Face-to-edge weights
        const scalarList& faceToEdgeWeights() const;

        //- Make face-to-edge weights
        void makeFaceToEdgeWeights() const;

        //- Clear weights
        void clearWeights();


public:

    // Constructors

        //- Construct from PrimitivePatch
        PrimitivePatchInterpolation(const Patch& p);


    // Destructor

        ~PrimitivePatchInterpolation();


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

        //- Interpolate from points to faces
        template<class Type>
        tmp<Field<Type> > pointToFaceInterpolate
        (
            const Field<Type>& pf
        ) const;

        template<class Type>
        tmp<Field<Type> > pointToFaceInterpolate
        (
            const tmp<Field<Type> >& tpf
        ) const;

        //- Interpolate from faces to edges
        template<class Type>
        tmp<Field<Type> > faceToEdgeInterpolate
        (
            const Field<Type>& ff
        ) const;

        template<class Type>
        tmp<Field<Type> > faceToEdgeInterpolate
        (
            const tmp<Field<Type> >& tff
        ) const;

        //- Do what is neccessary if the mesh has moved
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PrimitivePatchInterpolationTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
