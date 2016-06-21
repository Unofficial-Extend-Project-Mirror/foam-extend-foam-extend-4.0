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
    Foam::fv::gradScheme

Description
    Abstract base class for gradient schemes.

SourceFiles
    gradScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gradScheme_H
#define gradScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "BlockLduSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class gradScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gradScheme
:
    public refCount
{
public:

    // Public typedefs

    typedef Field<Type> FieldType;
    typedef GeometricField<Type, fvPatchField, volMesh> GeoFieldType;
    typedef typename GeoFieldType::GeometricBoundaryField GeoBoundaryFieldType;

    typedef Field<typename outerProduct<vector, Type>::type> GradFieldType;
    typedef GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    > GeoGradFieldType;
    typedef BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
    GradMatrixType;


    // Private data

        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow copy construct
        gradScheme(const gradScheme&);

        //- Disallow default bitwise assignment
        void operator=(const gradScheme&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            gradScheme,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        gradScheme(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new gradScheme created on freestore
        static tmp<gradScheme<Type> > New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~gradScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //  Moved from gaussGrad into base class.  HJ, 14/Jun/2013
        //- Correct the boundary values of the gradient using the patchField
        //  snGrad functions
        static void correctBoundaryConditions
        (
            const GeoFieldType&,
            GeoGradFieldType&
        );


        // Gradient functions

            //- Calculate and return the grad of the given field
            //  which may have been cached
            tmp<GeoGradFieldType> grad
            (
                const GeoFieldType&,
                const word& name
            ) const;

            //- Calculate and return the grad of the given field
            //  with the default name which may have been cached
            tmp<GeoGradFieldType> grad
            (
                const GeoFieldType&
            ) const;

            //- Calculate and return the grad of the given field
            //  with the default name which may have been cached
            tmp<GeoGradFieldType> grad
            (
                const tmp<GeoFieldType>&
            ) const;


        // Virtual function interface

            //- Calculate and return the grad of the given field.
            //  Used by grad either to recalculate the cached gradient when
            //  out of date with respect to the field or when it is not cached.
            virtual tmp<GeoGradFieldType> calcGrad
            (
                const GeoFieldType&,
                const word& name
            ) const = 0;

            //- Return the BlockLduSystem corresponding to the implicit grad
            // discretization. For block coupled systems.
            virtual tmp<GradMatrixType> fvmGrad
            (
                const GeoFieldType&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvGradTypeScheme(SS, Type)                                        \
                                                                              \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                             \
                                                                              \
gradScheme<Type>::addIstreamConstructorToTable<SS<Type> >                     \
    add##SS##Type##IstreamConstructorToTable_;


#define makeFvGradScheme(SS)                                                  \
                                                                              \
makeFvGradTypeScheme(SS, scalar)                                              \
makeFvGradTypeScheme(SS, vector)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gradScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
