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
    Foam::sampledPlane

Description
    A sampledSurface defined by a cuttingPlane. Triangulated by default.

SourceFiles
    sampledPlane.C

\*---------------------------------------------------------------------------*/

#ifndef sampledPlane_H
#define sampledPlane_H

#include "sampledSurface.H"
#include "cuttingPlane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sampledPlane Declaration
\*---------------------------------------------------------------------------*/

class sampledPlane
:
    public sampledSurface,
    public cuttingPlane
{
    // Private data

        //- zone name (if restricted to zones)
        word zoneName_;

        //- Triangulated faces in terms of intersection points.
        //  Non-triangulated faces are obtained directly from cuttingPlane
        faceList faces_;

        //- For every triangulated face, the original cell in mesh
        labelList meshCells_;

    // Private Member Functions

        //- Do all to create geometry. Triangulate as required
        void createGeometry();

        //- sample field on faces
        template <class Type>
        tmp<Field<Type> > sampleField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;


        template <class Type>
        tmp<Field<Type> >
        interpolateField(const interpolation<Type>&) const;


public:

    //- Runtime type information
    TypeName("sampledPlane");


    // Constructors

        //- Construct from components
        sampledPlane
        (
            const word& name,
            const polyMesh& mesh,
            const plane& planeDesc,
            const word& zoneName = word::null,
            const bool triangulate = true
        );

        //- Construct from dictionary
        sampledPlane
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~sampledPlane();


    // Member Functions

        //- Points of surface
        const pointField& points() const
        {
            return cuttingPlane::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            if (triangulate())
            {
                return faces_;
            }
            else
            {
                return cuttingPlane::faces();
            }
        }

        //- For every face original cell in mesh
        const labelList& meshCells() const
        {
            if (triangulate())
            {
                return meshCells_;
            }
            else
            {
                return cuttingPlane::cells();
            }
        }


        //- Correct for mesh movement and/or field changes
        virtual void correct(const bool meshChanged);


        //- sample field on surface
        virtual tmp<scalarField> sample
        (
            const volScalarField&
        ) const;


        //- sample field on surface
        virtual tmp<vectorField> sample
        (
            const volVectorField&
        ) const;

        //- sample field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const volSphericalTensorField&
        ) const;

        //- sample field on surface
        virtual tmp<symmTensorField> sample
        (
            const volSymmTensorField&
        ) const;

        //- sample field on surface
        virtual tmp<tensorField> sample
        (
            const volTensorField&
        ) const;


        //- interpolate field on surface
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>&
        ) const;


        //- interpolate field on surface
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>&
        ) const;

        //- interpolate field on surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>&
        ) const;

        //- interpolate field on surface
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>&
        ) const;

        //- interpolate field on surface
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>&
        ) const;

        //- Write
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sampledPlaneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
