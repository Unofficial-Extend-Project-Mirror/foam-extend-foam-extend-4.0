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
    devOneEqEddy

Description
<pre>
    One Equation Eddy Viscosity Model
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Eddy viscosity SGS model using a modeled balance equation to simulate the
    behaviour of k, hence,

        d/dt(k) + div(U*k) - div(nuEff*grad(k))
        =
        -B*L - ce*k^3/2/delta

    and

        B = 2/3*k*I - 2*nuEff*dev(D)

    where

        D = symm(grad(U));
        nuSgs = ck*sqrt(k)*delta
        nuEff = nuSgs + nu
</pre>

SourceFiles
    devOneEqEddy.C

\*---------------------------------------------------------------------------*/

#ifndef devOneEqEddy_H
#define devOneEqEddy_H

#include "GenEddyVisc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class devOneEqEddy Declaration
\*---------------------------------------------------------------------------*/

class devOneEqEddy
:
    public GenEddyVisc
{
    // Private data

        volScalarField k_;

        dimensionedScalar ck_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        devOneEqEddy(const devOneEqEddy&);
        devOneEqEddy& operator=(const devOneEqEddy&);


public:

    //- Runtime type information
    TypeName("devOneEqEddy");

    // Constructors

        //- Constructor from components
        devOneEqEddy
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~devOneEqEddy()
        {}


    // Member Functions

        //- Return SGS kinetic energy
        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nuSgs_ + nu())
            );
        }

        //- Correct Eddy-Viscosity and related properties
        void correct(const tmp<volTensorField>& gradU);

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModelsModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
