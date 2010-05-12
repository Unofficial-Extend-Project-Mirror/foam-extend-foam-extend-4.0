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
    convectionScheme

Description
    Abstract base class for convection schemes.

SourceFiles
    faConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef faConvectionScheme_H
#define faConvectionScheme_H

#include "tmp.H"
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
                           Class convectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class convectionScheme
:
    public refCount
{
    // Private data

        const faMesh& mesh_;


    // Private Member Functions

        //- Disallow copy construct
        convectionScheme(const convectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const convectionScheme&);


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            convectionScheme,
            Istream,
            (
                const faMesh& mesh,
                const edgeScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );


    // Constructors

        //- Construct from mesh, flux and Istream
        convectionScheme
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux
        )
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new convectionScheme created on freestore
        static tmp<convectionScheme<Type> > New
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream& schemeData
        );


    // Destructor

        virtual ~convectionScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<GeometricField<Type, faePatchField, edgeMesh> > flux
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const = 0;

        virtual tmp<faMatrix<Type> > famDiv
        (
            const edgeScalarField&,
            GeometricField<Type, faPatchField, areaMesh>&
        ) const = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDiv
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFaConvectionTypeScheme(SS, Type)                                   \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
convectionScheme<Type>::addIstreamConstructorToTable<SS<Type> >                \
    add##SS##Type##IstreamConstructorToTable_;


#define makeFaConvectionScheme(SS)                                             \
                                                                               \
makeFaConvectionTypeScheme(SS, scalar)                                         \
makeFaConvectionTypeScheme(SS, vector)                                         \
makeFaConvectionTypeScheme(SS, tensor)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
