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
    Foam::fv::skewCorrectedSnGrad

Description
    Simple central-difference snGrad scheme with skewness and
    non-orthogonal correction.

SourceFiles
    skewCorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef skewCorrectedSnGrad_H
#define skewCorrectedSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class skewCorrectedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewCorrectedSnGrad
:
    public snGradScheme<Type>
{
    // Private data

        scalar limitCoeff_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const skewCorrectedSnGrad&);


public:

    //- Runtime type information
    TypeName("skewCorrected");


    // Constructors

        //- Construct from mesh
        skewCorrectedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}


        //- Construct from mesh and data stream
        skewCorrectedSnGrad(const fvMesh& mesh, Istream& is)
        :
            snGradScheme<Type>(mesh),
            limitCoeff_(readScalar(is))
        {
            if (limitCoeff_ < 0 || limitCoeff_ > 1)
            {
                FatalIOErrorIn
                (
                    "skewCorrectedSnGrad(const fvMesh& mesh, Istream& is) : ",
                    is
                )   << "limitCoeff is specified as " << limitCoeff_
                    << " but should be >= 0 && <= 1"
                    << exit(FatalIOError);
            }
        }


    // Destructor

        virtual ~skewCorrectedSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().deltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
//             return !this->mesh().orthogonal();
        }

        //- Return the explicit correction to the skewCorrectedSnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >
skewCorrectedSnGrad<scalar>::correction
(
    const GeometricField<scalar, fvPatchField, volMesh>&
) const;

template<>
tmp<GeometricField<vector, fvsPatchField, surfaceMesh> >
skewCorrectedSnGrad<vector>::correction
(
    const GeometricField<vector, fvPatchField, volMesh>&
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "skewCorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
