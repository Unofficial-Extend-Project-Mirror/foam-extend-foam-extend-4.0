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
    gaussGrad

Description
    Basic second-order gradient scheme using face-interpolation
    and Gauss' theorem.

SourceFiles
    gaussFaGrad.C
    gaussFaGrads.C

\*---------------------------------------------------------------------------*/

#ifndef gaussFaGrad_H
#define gaussFaGrad_H

#include "faGradScheme.H"
#include "edgeInterpolationScheme.H"
#include "linearEdgeInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                           Class gaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussGrad
:
    public fa::gradScheme<Type>
{
    // Private data

        tmp<edgeInterpolationScheme<Type> > tinterpScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        gaussGrad(const gaussGrad&);

        //- Disallow default bitwise assignment
        void operator=(const gaussGrad&);


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct from mesh
        gaussGrad(const faMesh& mesh)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(new linearEdgeInterpolation<Type>(mesh))
        {}

        //- Construct from Istream
        gaussGrad(const faMesh& mesh, Istream& is)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(NULL)
        {
            if (is.eof())
            {
                tinterpScheme_ =
                    tmp<edgeInterpolationScheme<Type> >
                    (
                        new linearEdgeInterpolation<Type>(mesh)
                    );
            }
            else
            {
                tinterpScheme_ =
                    tmp<edgeInterpolationScheme<Type> >
                    (
                        edgeInterpolationScheme<Type>::New(mesh, is)
                    );
            }
        }


    // Member Functions

        //- Return the gradient of the given field calculated
        //  using Gauss' theorem on the interpolated field
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > grad
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        //- Correct the boundary values of the gradient using the patchField
        // snGrad functions
        static void correctBoundaryConditions
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gaussFaGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
