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
    edgeCorrectedVolPointInterpolation

Description
    Volume-to-point interpolation class which introduces the
    correction for the points shared between the patches.  This is
    done by gradient-based extrapolation into the vertices from the
    boundary faces around them using inverse-distance weighting.  This
    kind of interpolation is useful when trying to use the volume
    field as a basis of the mesh motion.

SourceFiles
    edgeCorrectedVolPointInterpolation.C
    edgeCorrectedVolPointInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef edgeCorrectedVolPointInterpolation_H
#define edgeCorrectedVolPointInterpolation_H

#include "volPointInterpolation.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class edgeCorrectedVolPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

class edgeCorrectedVolPointInterpolation
:
    public volPointInterpolation
{
    // Private data

        //- Extrapolation vectors
        mutable vectorListList extrapolationVectors_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        edgeCorrectedVolPointInterpolation
        (
            const edgeCorrectedVolPointInterpolation&
        );

        //- Disallow default bitwise assignment
        void operator=(const edgeCorrectedVolPointInterpolation&);


        //- Construct extrapolation vectors
        void makeExtrapolationVectors() const;


public:

    // Constructors

        //- Constructor given fvMesh and pointMesh.
        edgeCorrectedVolPointInterpolation(const fvMesh&, const pointMesh&);


    // Destructor

        virtual ~edgeCorrectedVolPointInterpolation();


    // Member Functions

        //- Update mesh topology using the morph engine
        virtual void updateMesh();

        //- Correct weighting factors for moving mesh.
        virtual bool movePoints();


    // Interpolation functions

        //- Interpolate from volField to pointField
        //  using inverse distance weighting with boundary correction
        //  given the field and its gradient
        template<class Type>
        void interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const GeometricField
            <
                typename outerProduct<vector, Type>::type,
                fvPatchField,
                volMesh
            >&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate volField returning pointField
        //  using inverse distance weighting with boundary correction
        //  given the field and its gradient
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const GeometricField
            <
                typename outerProduct<vector, Type>::type,
                fvPatchField,
                volMesh
            >&
        ) const;

        //- Interpolate tmp<volField> returning pointField
        //  using inverse distance weighting with boundary correction
        //  given the field and its gradient
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
            const GeometricField
            <
                typename outerProduct<vector, Type>::type,
                fvPatchField,
                volMesh
            >&
        ) const;

        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const tmp<GeometricField
            <
                typename outerProduct<vector, Type>::type,
                fvPatchField,
                volMesh
            > >&
        ) const;

        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
            const tmp<GeometricField
            <
                typename outerProduct<vector, Type>::type,
                fvPatchField,
                volMesh
            > >&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeCorrectedVolPointInterpolate.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
