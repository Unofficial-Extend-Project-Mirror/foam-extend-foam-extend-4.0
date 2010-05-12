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
    Foam::fa::limitedLnGrad

Description
    Central-difference lnGrad scheme with limited non-orthogonal correction.

    The limiter is controlled by a coefficient with a value between 0 and 1
    which when 0 switches the correction off and the scheme behaves as
    uncorrectedSnGrad, when set to 1 the full correction is applied and the
    scheme behaves as correctedSnGrad and when set to 0.5 the limiter is
    calculated such that the non-orthogonal contribution does not exceed the
    orthogonal part.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    limitedLnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef limitedLnGrad_H
#define limitedLnGrad_H

#include "lnGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                           Class limitedLnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class limitedLnGrad
:
    public lnGradScheme<Type>
{
    // Private data

        //- Limiter.  0 = no limiting, 1 = full limiting
        scalar limitCoeff_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const limitedLnGrad&);


public:

    //- Runtime type information
    TypeName("limited");


    // Constructors

        //- Construct from mesh
        limitedLnGrad(const faMesh& mesh)
        :
            lnGradScheme<Type>(mesh)
        {}


        //- Construct from mesh and data stream
        limitedLnGrad(const faMesh& mesh, Istream& is)
        :
            lnGradScheme<Type>(mesh),
            limitCoeff_(readScalar(is))
        {
            if (limitCoeff_ < 0 || limitCoeff_ > 1)
            {
                FatalIOErrorIn
                (
                    "limitedLnGrad(const faMesh& mesh, Istream& is) : ",
                    is
                )   << "limitCoeff is specified as " << limitCoeff_
                    << " but should be >= 0 && <= 1"
                    << exit(FatalIOError);
            }
        }


    // Destructor

        virtual ~limitedLnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<edgeScalarField> deltaCoeffs
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return this->mesh().deltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return !this->mesh().orthogonal();
        }

        //- Return the explicit correction to the limitedLnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, faePatchField, edgeMesh> >
        correction(const GeometricField<Type, faPatchField, areaMesh>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "limitedLnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
