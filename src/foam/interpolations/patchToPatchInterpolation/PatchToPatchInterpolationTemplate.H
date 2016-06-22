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
    Foam::PatchToPatchInterpolation

Description
    Interpolation class dealing with transfer of data between two
    primitivePatches

SourceFiles
    PatchToPatchInterpolation.C
    PatchToPatchInterpolate.C
    CalcPatchToPatchWeights.C

\*---------------------------------------------------------------------------*/

#ifndef PatchToPatchInterpolation_H
#define PatchToPatchInterpolation_H

#include "className.H"
#include "labelList.H"
#include "scalarField.H"
#include "pointField.H"
#include "FieldFields.H"
#include "faceList.H"
#include "intersection.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class PatchToPatchInterpolationName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(PatchToPatchInterpolation);


/*---------------------------------------------------------------------------*\
                  Class PatchToPatchInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class FromPatch, class ToPatch>
class PatchToPatchInterpolation
:
    public PatchToPatchInterpolationName
{
    // Private data

        //- Reference to the source patch
        const FromPatch& fromPatch_;

        //- Reference to the target patch
        const ToPatch& toPatch_;

        //- Type of intersection algorithm to use in projection
        intersection::algorithm alg_;

        //- Direction projection to use in projection
        intersection::direction dir_;


    // Static data

        //- Relative merge tolerance for projected points missing the target
        //  Expressed as the fraction of min involved edge size
        static const debug::tolerancesSwitch projectionTol_;

        //- Direct hit tolerance
        static const debug::tolerancesSwitch directHitTol_;


        // Point addressing

            //- Face into which each point of target patch is projected
            mutable labelList* pointAddressingPtr_;

            //- Weighting factors
            mutable FieldField<Field, scalar>* pointWeightsPtr_;

            //- Distance to intersection for patch points
            mutable scalarField* pointDistancePtr_;


        // Face addressing

            //- Face into which each face centre of target patch is projected
            mutable labelList* faceAddressingPtr_;

            //- Weighting factors
            mutable FieldField<Field, scalar>* faceWeightsPtr_;

            //- Distance to intersection for patch face centres
            mutable scalarField* faceDistancePtr_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const PatchToPatchInterpolation&);

        //- Calculate point weights
        void calcPointAddressing() const;

        //- Calculate face weights
        void calcFaceAddressing() const;

        //- Clear all geometry and addressing
        void clearOut();


protected:

    // Protected member functions

        //- Set weights when read from file
        void setWeights
        (
            labelList* paPtr,
            FieldField<Field, scalar>* pwPtr,
            scalarField* pdPtr,
            labelList* faPtr,
            FieldField<Field, scalar>* fwPtr,
            scalarField* fdPtr
        );


public:

    // Constructors

        //- Construct from components
        PatchToPatchInterpolation
        (
            const FromPatch& fromPatch,
            const ToPatch& toPatch,
            const intersection::algorithm alg = intersection::FULL_RAY,
            const intersection::direction dir = intersection::VECTOR
        );

        //- Construct as copy
        PatchToPatchInterpolation
        (
            const PatchToPatchInterpolation<FromPatch, ToPatch>&
        );


    // Destructor

        ~PatchToPatchInterpolation();


    // Member Functions

        // Access

            //- Return ype of intersection algorithm to use in projection
            intersection::algorithm projectionAlgo() const
            {
                return alg_;
            }

            //- Return direction projection to use in projection
            intersection::direction projectionDir() const
            {
                return dir_;
            }

            //- Return reference to point addressing
            const labelList& pointAddr() const;

            //- Return reference to point weights
            const FieldField<Field, scalar>& pointWeights() const;

            //- Return reference to face addressing
            const labelList& faceAddr() const;

            //- Return reference to face weights
            const FieldField<Field, scalar>& faceWeights() const;

            //- Return distance to intersection for patch points
            const scalarField& pointDistanceToIntersection() const;

            //- Return distance to intersection for patch face centres
            const scalarField& faceDistanceToIntersection() const;


        // Interpolation functions

            //- Interpolate point field
            template<class Type>
            tmp<Field<Type> > pointInterpolate(const Field<Type>& pf) const;

            template<class Type>
            tmp<Field<Type> > pointInterpolate
            (
                const tmp<Field<Type> >& tpf
            ) const;

            //- Interpolate face field
            template<class Type>
            tmp<Field<Type> > faceInterpolate(const Field<Type>& pf) const;

            template<class Type>
            tmp<Field<Type> > faceInterpolate
            (
                const tmp<Field<Type> >& tpf
            ) const;


        // Edit

            //- Correct weighting factors for moving mesh.
            bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "PatchToPatchInterpolationTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
