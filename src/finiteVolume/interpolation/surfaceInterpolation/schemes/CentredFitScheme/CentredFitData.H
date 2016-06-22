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
    Foam::CentredFitData

Description
    Data for the quadratic fit correction interpolation scheme

SourceFiles
    CentredFitData.C

\*---------------------------------------------------------------------------*/

#ifndef CentredFitData_H
#define CentredFitData_H

#include "FitData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class extendedCentredCellToFaceStencil;

/*---------------------------------------------------------------------------*\
                        Class CentredFitData Declaration
\*---------------------------------------------------------------------------*/

template<class Polynomial>
class CentredFitData
:
    public FitData
    <
        CentredFitData<Polynomial>,
        extendedCentredCellToFaceStencil,
        Polynomial
    >
{
    // Private data

        //- For each cell in the mesh store the values which multiply the
        //  values of the stencil to obtain the gradient for each direction
        List<scalarList> coeffs_;


    // Private member functions

        //- Calculate the fit for the all the mesh faces
        //  and set the coefficients
        void calcFit();


public:

    TypeName("CentredFitData");


    // Constructors

        //- Construct from components
        CentredFitData
        (
            const fvMesh& mesh,
            const extendedCentredCellToFaceStencil& stencil,
            const scalar linearLimitFactor,
            const scalar centralWeight
        );


    //- Destructor

        virtual ~CentredFitData()
        {}


    // Member functions

        //- Return reference to fit coefficients
        const List<scalarList>& coeffs() const
        {
            return coeffs_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CentredFitData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
