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
    Foam::fv::adjConvectionScheme

Description
    Abstract base class for adjConvection schemes.

SourceFiles
    adjConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef adjConvectionScheme_H
#define adjConvectionScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "multivariateSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class fvMatrix;

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class adjConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class adjConvectionScheme
:
    public refCount
{
    // Private data

        const fvMesh& mesh_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            adjConvectionScheme,
            Istream,
            (
                const fvMesh& mesh,
                const volVectorField& Up,
                Istream& schemeData
            ),
            (mesh, Up, schemeData)
        );


    // Constructors

        //- Copy construct
        adjConvectionScheme(const adjConvectionScheme&);

        //- Construct from mesh, flux and Istream
        adjConvectionScheme
        (
            const fvMesh& mesh,
            const volVectorField&
        )
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new adjConvectionScheme created on freestore
        static tmp<adjConvectionScheme<Type> > New
        (
            const fvMesh& mesh,
            const volVectorField& Up,
            Istream& schemeData
        );


    // Destructor

        virtual ~adjConvectionScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<fvMatrix<Type> > fvmAdjDiv
        (
            const volVectorField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcAdjDiv
        (
            const volVectorField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;


    // Member operators

        void operator=(const adjConvectionScheme<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvAdjConvectionTypeScheme(SS, Type)                               \
                                                                              \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                             \
                                                                              \
adjConvectionScheme<Type>::addIstreamConstructorToTable<SS<Type> >            \
    add##SS##Type##IstreamConstructorToTable_;

#define makeFvAdjConvectionScheme(SS)                                         \
                                                                              \
makeFvAdjConvectionTypeScheme(SS, vector)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "adjConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
