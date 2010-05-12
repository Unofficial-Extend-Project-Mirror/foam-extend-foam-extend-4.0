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
    Foam::surfaceSlipDisplacementPointPatchVectorField

Description
    Displacement follows a triSurface. Use in a displacement fvMotionSolver.
    Following is either
    - NEAREST : nearest
    - POINTNORMAL : intersection with point normal
    - FIXEDNORMAL : intersection with fixed vector

    Optionally (intersection only) removes a component ("wedgePlane") to
    stay in 2D.

    Needs:
    - projectSurfaces : names of triSurfaceMeshes (in constant/triSurface)
    - followMode : see above
    - projectDirection : if followMode = fixedNormal
    - wedgePlane : -1 or component to knock out of intersection normal
    - frozenPointsZone : empty or name of pointZone containing points
                         that do not move

SourceFiles
    surfaceSlipDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceSlipDisplacementPointPatchVectorField_H
#define surfaceSlipDisplacementPointPatchVectorField_H

#include "pointPatchFields.H"
#include "searchableSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class surfaceSlipDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class surfaceSlipDisplacementPointPatchVectorField
:
    public pointPatchVectorField
{

public:

    // Public data types

        enum followMode
        {
            NEAREST,
            POINTNORMAL,
            FIXEDNORMAL
        };

private:

    // Private data

        //- follow mode names
        static const NamedEnum<followMode, 3> followModeNames_;

        //- names of surfaces
        const dictionary surfacesDict_;

        //- How to follow/project onto surface
        const followMode projectMode_;

        //- direction to project
        const vector projectDir_;

        //- plane for 2D wedge case or -1.
        const label wedgePlane_;

        //- pointZone with frozen points
        const word frozenPointsZone_;

        //- Demand driven: surface to follow
        mutable autoPtr<searchableSurfaces> surfacesPtr_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const surfaceSlipDisplacementPointPatchVectorField&);


public:

    //- Runtime type information
    TypeName("surfaceSlipDisplacement");


    // Constructors

        //- Construct from patch and internal field
        surfaceSlipDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceSlipDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        surfaceSlipDisplacementPointPatchVectorField
        (
            const surfaceSlipDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy
        surfaceSlipDisplacementPointPatchVectorField
        (
            const surfaceSlipDisplacementPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchVectorField> clone() const
        {
            return autoPtr<pointPatchVectorField>
            (
                new surfaceSlipDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        surfaceSlipDisplacementPointPatchVectorField
        (
            const surfaceSlipDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchVectorField> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchVectorField>
            (
                new surfaceSlipDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }

    // Member Functions

        //- Surfaces to follow
        const dictionary& surfacesDict() const
        {
            return surfacesDict_;
        }

        //- Surface to follow. Demand loads surfaceNames.
        const searchableSurfaces& surfaces() const;

        //- Mode of projection/following
        followMode projectMode() const
        {
            return projectMode_;
        }

        //- Direction to project back onto surface
        const vector& projectDir() const
        {
            return projectDir_;
        }

        //- Normal of wedgeplane (0, 1, 2) or -1. Note: should be obtained
        //  from twoDPointCorrector.
        label wedgePlane() const
        {
            return wedgePlane_;
        }

        //- Zone containing frozen points
        const word& frozenPointsZone() const
        {
            return frozenPointsZone_;
        }

        //- Update the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
