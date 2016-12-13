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
    Foam::fv::LimitedGrad

Description
    Second-order gradient scheme with a templated limiter.

SourceFiles
    LimitedGrad.C

\*---------------------------------------------------------------------------*/

#ifndef LimitedGrad_H
#define LimitedGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class LimitedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GradientLimiter>
class LimitedGrad
{
public:

    // Public typedefs
    typedef typename outerProduct<vector, Type>::type OuterProductType;

    typedef Field<Type> FieldType;
    typedef GeometricField<Type, fvPatchField, volMesh> GeoFieldType;
    typedef typename GeoFieldType::GeometricBoundaryField GeoBoundaryFieldType;

    typedef Field<typename outerProduct<vector, Type>::type> GradFieldType;
    typedef GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    > GeoGradFieldType;
    typedef BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
    GradMatrixType;


private:

    // Private Data

        //- Base gradient scheme
        tmp<fv::gradScheme<Type> > basicGradScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        LimitedGrad(const LimitedGrad&);

        //- Disallow default bitwise assignment
        void operator=(const LimitedGrad&);


        //- Calculate limiter
        tmp<GeoFieldType> limiter
        (
            const GeoFieldType& vf,
            const GeoGradFieldType& gradVf
        ) const;


public:

    // Constructors

        //- Construct from mesh and schemeData
        LimitedGrad(const fvMesh& mesh, Istream& schemeData)
        :
            basicGradScheme_(fv::gradScheme<Type>::New(mesh, schemeData))
        {}


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        tmp<GeoGradFieldType> gradientField
        (
            const GeoFieldType& vf,
            const word& name
        ) const;

        //- Return the BlockLduSystem corresponding to the implicit gradient
        //  discretization. For block coupled systems.
        tmp<GradMatrixType> gradientMatrix
        (
            const GeoFieldType& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LimitedGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
