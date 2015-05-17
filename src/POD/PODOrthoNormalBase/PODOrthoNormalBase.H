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
    PODOrthoNormalBase

Description
    Establish POD ortho-normal base and interpolation coefficients give a list
    of fields. Size of ortho-normal base is calculated from the desired
    accuracy, e.g. 0.99-0.99999 (in energy terms)

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    PODOrthoNormalBase.C

\*---------------------------------------------------------------------------*/

#ifndef PODOrthoNormalBase_H
#define PODOrthoNormalBase_H

#include "volFields.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PODOrthoNormalBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PODOrthoNormalBase
{
public:

    // Public typedefs

        typedef GeometricField<Type, fvPatchField, volMesh> GeoTypeField;

        typedef Field<Type> TypeField;


private:

    // Private data

        //- List of ortho-normal fields
        PtrList<GeoTypeField> orthoFields_;

        //- Interpolation coefficients
        RectangularMatrix<Type>* interpolationCoeffsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PODOrthoNormalBase(const PODOrthoNormalBase&);

        //- Disallow default bitwise assignment
        void operator=(const PODOrthoNormalBase&);

        //- Calculate ortho-normal base
        void calcOrthoBase
        (
            const PtrList<GeoTypeField>& snapshots,
            const scalar accuracy
        );


public:

    // Constructors

        //- Construct given list of snapshots and accuracy
        PODOrthoNormalBase
        (
            const PtrList<GeoTypeField>& snapshots,
            const scalar accuracy
        );


    // Destructor

        ~PODOrthoNormalBase();


    // Member Functions

        //- Return base size
        label baseSize() const
        {
            return orthoFields_.size();
        }

        //- Return n-th ortho-normal base field
        const GeoTypeField& orthoField(const label n) const
        {
            return orthoFields_[n];
        }

        //- Return interpolation coefficients (snapshots x orthoVectors)
        const RectangularMatrix<Type>& interpolationCoeffs() const
        {
            return *interpolationCoeffsPtr_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PODOrthoNormalBase.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
