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
    Foam::fa::leastSquaresFaGrad

Description
    Second-order gradient scheme using least-squares.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    leastSquaresFaGrad.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresFaGrad_H
#define leastSquaresFaGrad_H

#include "faGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                      Class leastSquaresFaGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class leastSquaresFaGrad
:
    public fa::gradScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        leastSquaresFaGrad(const leastSquaresFaGrad&);

        //- Disallow default bitwise assignment
        void operator=(const leastSquaresFaGrad&);


public:

    //- Runtime type information
    TypeName("leastSquares");


    // Constructors

        //- Construct from mesh
        leastSquaresFaGrad(const faMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from Istream
        leastSquaresFaGrad(const faMesh& mesh, Istream&)
        :
            gradScheme<Type>(mesh)
        {}


    // Member Functions

        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > grad
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "leastSquaresFaGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
