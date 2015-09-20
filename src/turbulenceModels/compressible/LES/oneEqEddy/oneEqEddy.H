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
    Foam::compressible::LESModels::oneEqEddy

Description
    One Equation Eddy Viscosity Model for incompressible flows

    Eddy viscosity SGS model using a modeled balance equation to simulate the
    behaviour of k, hence,
    @verbatim
        d/dt(rho*k) + div(rho*U*k) - div(muEff*grad(k))
        =
        -rho*D:B - ce*rho*k^3/2/delta

    and

        B = 2/3*k*I - 2*nuSgs*dev(D)

    where

        D = symm(grad(U));
        muSgs = ck*rho*sqrt(k)*delta
    @endverbatim


SourceFiles
    oneEqEddy.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleOneEqEddy_H
#define compressibleOneEqEddy_H

#include "GenEddyVisc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class oneEqEddy Declaration
\*---------------------------------------------------------------------------*/

class oneEqEddy
:
    public GenEddyVisc
{
    // Private data

        dimensionedScalar ck_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields();

        // Disallow default bitwise copy construct and assignment
        oneEqEddy(const oneEqEddy&);
        oneEqEddy& operator=(const oneEqEddy&);


public:

    //- Runtime type information
    TypeName("oneEqEddy");


    // Constructors

        //- Constructor from components
        oneEqEddy
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel
        );


    // Destructor
    virtual ~oneEqEddy()
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
