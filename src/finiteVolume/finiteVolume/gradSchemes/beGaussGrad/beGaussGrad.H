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
    beGaussGrad

Description
    Second-order gradient scheme using face-interpolation
    and Gauss' theorem.  Boundary faces without specified gradient are treated
    through implicit extrapolation.

SourceFiles
    beGaussGrad.C

\*---------------------------------------------------------------------------*/

#ifndef beGaussGrad_H
#define beGaussGrad_H

#include "gradScheme.H"
#include "surfaceInterpolationScheme.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class beGaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class beGaussGrad
:
    public fv::gradScheme<Type>
{
    // Private data

        tmp<surfaceInterpolationScheme<Type> > tinterpScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        beGaussGrad(const beGaussGrad&);

        //- Disallow default bitwise assignment
        void operator=(const beGaussGrad&);


public:

    //- Runtime type information
    TypeName("beGauss");


    // Constructors

        //- Construct from mesh
        beGaussGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(new linear<Type>(mesh))
        {}

        //- Construct from Istream
        beGaussGrad(const fvMesh& mesh, Istream& is)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(NULL)
        {
            if (is.eof())
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type> >
                    (
                        new linear<Type>(mesh)
                    );
            }
            else
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type> >
                    (
                        surfaceInterpolationScheme<Type>::New(mesh, is)
                    );
            }
        }


    // Member Functions

        //- Return the gradient of the given field calculated
        //  using Gauss' theorem on the interpolated field
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
#   include "beGaussGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
