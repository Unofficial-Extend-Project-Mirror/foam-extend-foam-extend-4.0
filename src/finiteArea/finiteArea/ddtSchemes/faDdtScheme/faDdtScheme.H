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
    faDdtScheme

Description
    Abstract base class for ddt schemes.

SourceFiles
    faDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef faDdtScheme_H
#define faDdtScheme_H

#include "tmp.H"
#include "dimensionedType.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class faMatrix;

class faMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                           Class faDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faDdtScheme
:
    public refCount
{

protected:

    // Protected data

        const faMesh& mesh_;


    // Private Member Functions

        //- Disallow copy construct
        faDdtScheme(const faDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const faDdtScheme&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            faDdtScheme,
            Istream,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        faDdtScheme(const faMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Construct from mesh and Istream
        faDdtScheme(const faMesh& mesh, Istream&)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new faDdtScheme created on freestore
        static tmp<faDdtScheme<Type> > New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~faDdtScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const dimensioned<Type>
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const dimensioned<Type>
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faePatchField, edgeMesh> > facDdt0
        (
            const GeometricField<Type, faePatchField, edgeMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type> > famDdt
        (
            GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type> > famDdt
        (
            const dimensionedScalar&,
            GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type> > famDdt
        (
            const areaScalarField&,
            GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFaDdtTypeScheme(SS, Type)                                          \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
faDdtScheme<Type>::addIstreamConstructorToTable<SS<Type> >                     \
    add##SS##Type##IstreamConstructorToTable_;


#define makeFaDdtScheme(SS)                                                    \
                                                                               \
makeFaDdtTypeScheme(SS, scalar)                                                \
makeFaDdtTypeScheme(SS, vector)                                                \
makeFaDdtTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
