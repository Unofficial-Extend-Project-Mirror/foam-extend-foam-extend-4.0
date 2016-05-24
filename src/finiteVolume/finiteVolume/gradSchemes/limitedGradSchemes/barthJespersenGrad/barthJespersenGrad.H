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
    Foam::fv::barthJespersenGrad

Description
    Barth-Jespersen gradient limiter applied to a runTime selected
    base gradient scheme.

SourceFiles
    barthJespersenGrad.C

\*---------------------------------------------------------------------------*/

#ifndef barthJespersenGrad_H
#define barthJespersenGrad_H

#include "gradScheme.H"
#include "LimitedGrad.H"
#include "BarthJespersenLimiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class barthJespersenGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class barthJespersenGrad
:
    public fv::gradScheme<Type>,
    public LimitedGrad<Type, BarthJespersenLimiter>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        barthJespersenGrad(const barthJespersenGrad&);

        //- Disallow default bitwise assignment
        void operator=(const barthJespersenGrad&);


public:

    //- RunTime type information
    TypeName("barthJespersen");


    // Constructors

        //- Construct from mesh and schemeData
        barthJespersenGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            LimitedGrad<Type, BarthJespersenLimiter>(mesh, schemeData)
        {}


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const word& name
        ) const
        {
            return LimitedGrad<Type, BarthJespersenLimiter>::gradientField
            (
                vf,
                name
            );
        }

        //- Return the BlockLduSystem corresponding to the implicit cell
        //  limited grad discretization.  For block coupled systems.
        virtual tmp
        <
            BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
        > fvmGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return LimitedGrad<Type, BarthJespersenLimiter>::gradientMatrix
            (
                vf
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
