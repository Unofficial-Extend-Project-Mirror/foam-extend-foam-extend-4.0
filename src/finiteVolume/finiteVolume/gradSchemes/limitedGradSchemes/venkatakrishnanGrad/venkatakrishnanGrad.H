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
    Foam::fv::venkatakrishnanGrad

Description
    Venkatakrishnan gradient limiter applied to a runTime selected
     base gradient scheme.

SourceFiles
    venkatakrishnanGrad.C

\*---------------------------------------------------------------------------*/

#ifndef venkatakrishnanGrad_H
#define venkatakrishnanGrad_H

#include "gradScheme.H"
#include "LimitedGrad.H"
#include "VenkatakrishnanLimiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class venkatakrishnanGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class venkatakrishnanGrad
:
    public fv::gradScheme<Type>,
    public LimitedGrad<Type, VenkatakrishnanLimiter>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        venkatakrishnanGrad(const venkatakrishnanGrad&);

        //- Disallow default bitwise assignment
        void operator=(const venkatakrishnanGrad&);


public:

    //- RunTime type information
    TypeName("venkatakrishnan");


    // Constructors

        //- Construct from mesh and schemeData
        venkatakrishnanGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            LimitedGrad<Type, VenkatakrishnanLimiter>(mesh, schemeData)
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
            return LimitedGrad<Type, VenkatakrishnanLimiter>::gradientField
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
            return LimitedGrad<Type, VenkatakrishnanLimiter>::gradientMatrix
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
