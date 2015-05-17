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
    Foam::compressible::LESModels::dynOneEqEddy

Description
    One Equation Eddy Viscosity Model for compressible flows.

    Eddy viscosity SGS model using a modeled balance equation to simulate
    the behaviour of k.
    Thus
    @verbatim
        d/dt(k) + div(U*k) - div(nuSgs*grad(k))
        =
        -rho*B*L - ce*rho*k^3/2/delta

    and

        B = 2/3*k*I - 2*nuSgs*dev(D)

    where

        D = symm(grad(U));
        nuSgs = ck*sqrt(k)*delta
    @endverbatim

SourceFiles
    dynOneEqEddy.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleDynOneEqEddy_H
#define compressibleDynOneEqEddy_H

#include "GenEddyVisc.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class dynOneEqEddy Declaration
\*---------------------------------------------------------------------------*/

class dynOneEqEddy
:
    public GenEddyVisc
{
    // Private data

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volSymmTensorField& D);

        //- Calculate ck, ce by filtering the velocity field U.
        dimensionedScalar ck_(const volSymmTensorField& D) const;
        dimensionedScalar ce_(const volSymmTensorField& D) const;

        // Disallow default bitwise copy construct and assignment
        dynOneEqEddy(const dynOneEqEddy&);
        dynOneEqEddy& operator=(const dynOneEqEddy&);


public:

    //- Runtime type information
    TypeName("dynOneEqEddy");


    // Constructors

        //- Constructor from components
        dynOneEqEddy
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel
        );


    // Destructor
    virtual ~dynOneEqEddy()
    {}


    // Member Functions

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", muSgs_ + mu())
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
