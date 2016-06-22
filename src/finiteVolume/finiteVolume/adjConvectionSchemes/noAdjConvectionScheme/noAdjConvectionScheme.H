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
            const GeometricField<Type, fvPatchField, volMesh>&
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
