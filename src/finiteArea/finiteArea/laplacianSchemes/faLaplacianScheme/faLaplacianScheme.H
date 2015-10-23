/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    faLaplacianScheme

Description
    Abstract base class for laplacian schemes.

SourceFiles
    faLaplacianScheme.C
    faLaplacianSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef faLaplacianScheme_H
#define faLaplacianScheme_H

#include "tmp.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "linearEdgeInterpolation.H"
#include "correctedLnGrad.H"
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
                           Class laplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class laplacianScheme
:
    public refCount
{

protected:

    // Protected data

        const faMesh& mesh_;
        tmp<edgeInterpolationScheme<scalar> > tinterpGammaScheme_;
        tmp<lnGradScheme<Type> > tlnGradScheme_;


    // Private Member Functions

        //- Disallow copy construct
        laplacianScheme(const laplacianScheme&);

        //- Disallow default bitwise assignment
        void operator=(const laplacianScheme&);


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            laplacianScheme,
            Istream,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        laplacianScheme(const faMesh& mesh)
        :
            mesh_(mesh),
            tinterpGammaScheme_(new linearEdgeInterpolation<scalar>(mesh)),
            tlnGradScheme_(new correctedLnGrad<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        laplacianScheme(const faMesh& mesh, Istream& is)
        :
            mesh_(mesh),
            tinterpGammaScheme_(NULL),
            tlnGradScheme_(NULL)
        {
            if (is.eof())
            {
                tinterpGammaScheme_ = tmp<edgeInterpolationScheme<scalar> >
                (
                    new linearEdgeInterpolation<scalar>(mesh)
                );

                tlnGradScheme_ = tmp<lnGradScheme<Type> >
                (
                    new correctedLnGrad<Type>(mesh)
                );
            }
            else
            {
                tinterpGammaScheme_ = tmp<edgeInterpolationScheme<scalar> >
                (
                    edgeInterpolationScheme<scalar>::New(mesh, is)
                );

                tlnGradScheme_ = tmp<lnGradScheme<Type> >
                (
                    lnGradScheme<Type>::New(mesh, is)
                );
            }
        }


    // Selectors

        //- Return a pointer to a new laplacianScheme created on freestore
        static tmp<laplacianScheme<Type> > New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~laplacianScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<faMatrix<Type> > famLaplacian
        (
            const edgeScalarField&,
            GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type> > famLaplacian
        (
            const areaScalarField&,
            GeometricField<Type, faPatchField, areaMesh>&
        );

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facLaplacian
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facLaplacian
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facLaplacian
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFaLaplacianTypeScheme(SS, Type)                                    \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
laplacianScheme<Type>::addIstreamConstructorToTable<SS<Type> >                 \
    add##SS##Type##IstreamConstructorToTable_;


#define makeFaLaplacianScheme(SS)                                              \
                                                                               \
makeFaLaplacianTypeScheme(SS, scalar)                                          \
makeFaLaplacianTypeScheme(SS, vector)                                          \
makeFaLaplacianTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
