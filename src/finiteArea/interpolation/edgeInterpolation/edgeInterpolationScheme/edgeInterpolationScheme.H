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
    edgeInterpolationScheme

Description
    Abstract base class for edge interpolation schemes.

SourceFiles
    edgeInterpolationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef edgeInterpolationScheme_H
#define edgeInterpolationScheme_H

#include "tmp.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faMesh;

/*---------------------------------------------------------------------------*\
                    Class edgeInterpolationScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class edgeInterpolationScheme
:
    public refCount
{
    // Private data

        //- Hold reference to mesh
        const faMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const edgeInterpolationScheme&);


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            edgeInterpolationScheme,
            Mesh,
            (
                const faMesh& mesh,
                Istream& schemeData
            ),
            (mesh, schemeData)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            edgeInterpolationScheme,
            MeshFlux,
            (
                const faMesh& mesh,
                const edgeScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );


    // Constructors

        //- Construct from mesh
        edgeInterpolationScheme(const faMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<edgeInterpolationScheme<Type> > New
        (
            const faMesh& mesh,
            Istream& schemeData
        );

        //- Return new tmp interpolation scheme
        static tmp<edgeInterpolationScheme<Type> > New
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream& schemeData
        );


    // Destructor

        virtual ~edgeInterpolationScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return mesh_;
        }


        //- Return the face-interpolate of the given cell field
        //  with the given owner and neighbour weigting factors
        static tmp<GeometricField<Type, faePatchField, edgeMesh> >
        interpolate
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const tmp<edgeScalarField>&,
            const tmp<edgeScalarField>&
        );

        //- Return the face-interpolate of the given cell field
        //  with the given weigting factors
        static tmp<GeometricField<Type, faePatchField, edgeMesh> >
        interpolate
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const tmp<edgeScalarField>&
        );


        //- Return the euclidian edge-interpolate of the given area field
        //  with the given weigting factors
        static tmp<GeometricField<Type, faePatchField, edgeMesh> >
        euclidianInterpolate
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const tmp<edgeScalarField>&
        );


        //- Return the interpolation weighting factors for the given field
        virtual tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const = 0;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the face-interpolate
        //  for the given field
        virtual tmp<GeometricField<Type, faePatchField, edgeMesh> >
        correction(const GeometricField<Type, faPatchField, areaMesh>&) const
        {
            return tmp<GeometricField<Type, faePatchField, edgeMesh> >(NULL);
        }

        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        virtual tmp<GeometricField<Type, faePatchField, edgeMesh> >
        interpolate(const GeometricField<Type, faPatchField, areaMesh>&) const;

        //- Return the euclidian edge-interpolate of the given area field
        //  without explicit correction
        virtual tmp<GeometricField<Type, faePatchField, edgeMesh> >
        euclidianInterpolate
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        //- Return the face-interpolate of the given tmp cell field
        //  with explicit correction
        tmp<GeometricField<Type, faePatchField, edgeMesh> >
        interpolate
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh> >&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeEdgeInterpolationTypeScheme(SS, Type)                              \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
edgeInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type> >            \
    add##SS##Type##MeshConstructorToTable_;                                    \
                                                                               \
edgeInterpolationScheme<Type>::addMeshFluxConstructorToTable<SS<Type> >        \
    add##SS##Type##MeshFluxConstructorToTable_;

#define makeEdgeInterpolationScheme(SS)                                        \
                                                                               \
makeEdgeInterpolationTypeScheme(SS, scalar)                                    \
makeEdgeInterpolationTypeScheme(SS, vector)                                    \
makeEdgeInterpolationTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "edgeInterpolationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
