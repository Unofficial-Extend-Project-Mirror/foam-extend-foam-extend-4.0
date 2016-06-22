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
    Foam::fv::explicitAdjConvectionScheme

Description
    AdjConvection scheme used to make adjConvection explicit without
    recompilation

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    explicitAdjConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef explicitAdjConvectionScheme_H
#define explicitAdjConvectionScheme_H

#include "adjConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class explicitAdjConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class explicitAdjConvectionScheme
:
    public fv::adjConvectionScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        explicitAdjConvectionScheme(const explicitAdjConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const explicitAdjConvectionScheme&);


public:

    //- Runtime type information
    TypeName("explicit");


    // Constructors

        //- Construct from flux and interpolation scheme
        explicitAdjConvectionScheme
        (
            const fvMesh& mesh,
            const volVectorField& Up
        )
        :
            adjConvectionScheme<Type>(mesh, Up)
        {}

        //- Construct from flux and Istream
        explicitAdjConvectionScheme
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
#   include "explicitAdjConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
