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
    Foam::fv::leastSquaresSolidInterfaceGrad

Description
    Second-order gradient scheme using least-squares.
    philipc-
    Modified to use face centre values at multi-material interfaces.

SourceFiles
    leastSquaresSolidInterfaceGrad.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresSolidInterfaceGrad_H
#define leastSquaresSolidInterfaceGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class leastSquaresSolidInterfaceGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class leastSquaresSolidInterfaceGrad
:
    public fv::gradScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        leastSquaresSolidInterfaceGrad(const leastSquaresSolidInterfaceGrad&);

        //- Disallow default bitwise assignment
        void operator=(const leastSquaresSolidInterfaceGrad&);


public:

    //- Runtime type information
    TypeName("leastSquaresSolidInterface");


    // Constructors

        //- Construct from mesh
        leastSquaresSolidInterfaceGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from Istream
        leastSquaresSolidInterfaceGrad(const fvMesh& mesh, Istream&)
        :
            gradScheme<Type>(mesh)
        {}


    // Member Functions

        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > grad
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "leastSquaresSolidInterfaceGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
