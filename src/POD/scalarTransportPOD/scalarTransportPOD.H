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
    scalarTransportPOD

Description
    POD solver for a transport equation for a passive scalar

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    scalarTransportPOD.C

\*---------------------------------------------------------------------------*/

#ifndef scalarTransportPOD_H
#define scalarTransportPOD_H

#include "PODODE.H"
#include "PODOrthoNormalBases.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class scalarTransportPOD Declaration
\*---------------------------------------------------------------------------*/

class scalarTransportPOD
:
    public PODODE
{
    // Private data

        //- Field name
        word phiName_;

        //- POD coefficients
        scalarField coeffs_;

        //- Derivative coefficient matrix
        mutable scalarSquareMatrix* derivativeMatrixPtr_;

        //- Ortho-normal base pointer
        mutable scalarPODOrthoNormalBase* orthoBasePtr_;

        //- Reconstructed field pointer
        mutable volScalarField* fieldPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        scalarTransportPOD(const scalarTransportPOD&);

        //- Disallow default bitwise assignment
        void operator=(const scalarTransportPOD&);

        //- Calculate ortho-normal base
        void calcOrthoBase() const;

        //- Calculate coefficients
        void calcDerivativeCoeffs() const;


public:

    //- Runtime type information
    TypeName("scalarTransport");


    // Constructors

        //- Construct from components
        scalarTransportPOD
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~scalarTransportPOD();


    // Member Functions

        // ODE parameters

            //- Return number of equations
            virtual label nEqns() const;

            //- Return reference to interpolation coefficients
            virtual scalarField& coeffs();

            //- Return reference to interpolation coefficients
            virtual const scalarField& coeffs() const;

            //- Return derivatives
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const;

            //- Return Jacobian
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const;


    // Orthogonalisation and fields

        //- Return ortho-normal base
        const scalarPODOrthoNormalBase& orthoBase() const;

        //- Return reconstructed field
        const volScalarField& field() const;


        //- Clear ortho-normal base
        virtual void clearBase() const;

        //- Update reconstructed fields
        virtual void updateFields() const;

        //- Clear reconstructed field
        virtual void clearFields() const;


        //- Write reconstructed fields
        virtual void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
