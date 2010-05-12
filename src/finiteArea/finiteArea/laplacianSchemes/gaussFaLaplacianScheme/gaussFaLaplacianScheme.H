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
    gaussFaLaplacianScheme

Description
    Basic second-order laplacian using face-gradients and Gauss' theorem.

SourceFiles
    gaussFaLaplacianScheme.C
    gaussFaLaplacianSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef gaussFaLaplacianScheme_H
#define gaussFaLaplacianScheme_H

#include "faLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class gaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussLaplacianScheme
:
    public fa::laplacianScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        gaussLaplacianScheme(const gaussLaplacianScheme&);

        //- Disallow default bitwise assignment
        void operator=(const gaussLaplacianScheme&);


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct null
        gaussLaplacianScheme(const faMesh& mesh)
        :
            laplacianScheme<Type>(mesh)
        {}

        //- Construct from Istream
        gaussLaplacianScheme(const faMesh& mesh, Istream& is)
        :
            laplacianScheme<Type>(mesh, is)
        {}


    // Member Functions

        tmp<faMatrix<Type> > famLaplacian
        (
            const edgeScalarField&,
            GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh> > facLaplacian
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh> > facLaplacian
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gaussFaLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
