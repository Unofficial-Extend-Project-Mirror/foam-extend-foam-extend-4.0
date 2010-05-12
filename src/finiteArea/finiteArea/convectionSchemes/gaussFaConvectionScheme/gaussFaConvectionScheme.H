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
    gaussConvectionScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    gaussFaConvectionScheme.C
    gaussFaConvectionSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef gaussFaConvectionScheme_H
#define gaussFaConvectionScheme_H

#include "faConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class gaussConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussConvectionScheme
:
    public fa::convectionScheme<Type>
{
    // Private data

        tmp<edgeInterpolationScheme<Type> > tinterpScheme_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        gaussConvectionScheme(const gaussConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const gaussConvectionScheme&);


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct from flux and interpolation scheme
        gaussConvectionScheme
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            const tmp<edgeInterpolationScheme<Type> >& scheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_(scheme)
        {}

        //- Construct from flux and Istream
        gaussConvectionScheme
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_
            (
                edgeInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        tmp<GeometricField<Type, faePatchField, edgeMesh> > flux
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        tmp<faMatrix<Type> > famDiv
        (
            const edgeScalarField&,
            GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        tmp<GeometricField<Type, faPatchField, areaMesh> > facDiv
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gaussFaConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
