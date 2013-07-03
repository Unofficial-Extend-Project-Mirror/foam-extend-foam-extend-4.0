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
    Foam::fv::noAdjConvectionScheme

Description
    Dummy adjConvection scheme, used to switch off adjConvection without
    recompilation

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    noAdjConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef noAdjConvectionScheme_H
#define noAdjConvectionScheme_H

#include "adjConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class noAdjConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class noAdjConvectionScheme
:
    public fv::adjConvectionScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        noAdjConvectionScheme(const noAdjConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const noAdjConvectionScheme&);


public:

    //- Runtime type information
    TypeName("off");


    // Constructors

        //- Construct from flux and interpolation scheme
        noAdjConvectionScheme
        (
            const fvMesh& mesh,
            const volVectorField& Up
        )
        :
            adjConvectionScheme<Type>(mesh, Up)
        {}

        //- Construct from flux and Istream
        noAdjConvectionScheme
        (
            const fvMesh& mesh,
            const volVectorField& Up,
            Istream& is
        )
        :
            adjConvectionScheme<Type>(mesh, Up)
        {}


    // Member Functions

        tmp<fvMatrix<Type> > fvmAdjDiv
        (
            const volVectorField&,
            GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcAdjDiv
        (
            const volVectorField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "noAdjConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
