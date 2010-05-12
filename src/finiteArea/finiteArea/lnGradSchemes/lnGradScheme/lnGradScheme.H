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
    lnGradScheme

Description
    Abstract base class for lnGrad schemes.

SourceFiles
    lnGradScheme.C

\*---------------------------------------------------------------------------*/

#ifndef lnGradScheme_H
#define lnGradScheme_H

#include "tmp.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                 Class lnGradScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class lnGradScheme
:
    public refCount
{
    // Private data

        //- Hold reference to mesh
        const faMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const lnGradScheme&);


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            lnGradScheme,
            Mesh,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        lnGradScheme(const faMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<lnGradScheme<Type> > New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~lnGradScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the lnGrad of the given cell field
        //  with the given weigting factors
        static tmp<GeometricField<Type, faePatchField, edgeMesh> >
        lnGrad
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const tmp<edgeScalarField>&,
            const word& snGradName = "lnGrad"
        );

        //- Return the interpolation weighting factors for the given field
        virtual tmp<edgeScalarField> deltaCoeffs
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const = 0;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the lnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, faePatchField, edgeMesh> >
        correction(const GeometricField<Type, faPatchField, areaMesh>&) const
        {
            return tmp<GeometricField<Type, faePatchField, edgeMesh> >(NULL);
        }

        //- Return the lnGrad of the given cell field
        //  with explicit correction
        virtual tmp<GeometricField<Type, faePatchField, edgeMesh> >
        lnGrad(const GeometricField<Type, faPatchField, areaMesh>&) const;

        //- Return the lnGrad of the given tmp cell field
        //  with explicit correction
        tmp<GeometricField<Type, faePatchField, edgeMesh> >
        lnGrad
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh> >&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeLnGradTypeScheme(SS, Type)                                         \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
lnGradScheme<Type>::addMeshConstructorToTable<SS<Type> >                       \
    add##SS##Type##MeshConstructorToTable_;

#define makeLnGradScheme(SS)                                                   \
                                                                               \
makeLnGradTypeScheme(SS, scalar)                                               \
makeLnGradTypeScheme(SS, vector)                                               \
makeLnGradTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "lnGradScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
