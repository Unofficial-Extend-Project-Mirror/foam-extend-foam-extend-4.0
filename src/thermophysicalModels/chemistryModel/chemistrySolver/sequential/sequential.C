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

\*---------------------------------------------------------------------------*/

#include "sequential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sequential, 0);
    addToRunTimeSelectionTable
    (
        chemistrySolver,
        sequential,
        dictionary
    );
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::sequential::sequential
(
    const Foam::dictionary& dict,
    Foam::chemistryModel& chemistry
)
:
    chemistrySolver(dict, chemistry),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    cTauChem_(readScalar(coeffsDict_.lookup("cTauChem"))),
    equil_(coeffsDict_.lookup("equilibriumRateLimiter"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sequential::~sequential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::sequential::solve
(
    scalarField &c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{
    scalar tChemInv = SMALL;

    scalar pf, cf, pb, cb;
    label lRef, rRef;

    for(label i=0; i<chemistry_.reactions().size(); i++)
    {
        const chemistryModel::reaction& R = chemistry_.reactions()[i];

        scalar om0 = chemistry_.omega
        (
            R, c, T, p, pf, cf, lRef, pb, cb, rRef
        );

        scalar omeg = 0.0;
        if (!equil_)
        {
            omeg = om0;
        }
        else
        {
            if (om0<0.0)
            {
                omeg = om0/(1.0 + pb*dt);
            }
            else
            {
                omeg = om0/(1.0 + pf*dt);
            }
        }
        tChemInv = max(tChemInv, mag(omeg));


        // update species
        for(label s=0; s<R.lhs().size(); s++)
        {
            label si = R.lhs()[s].index;
            scalar sl = R.lhs()[s].stoichCoeff;
            c[si] -= dt*sl*omeg;
            c[si] = max(0.0, c[si]);
        }        

        for(label s=0; s<R.rhs().size(); s++)
        {
            label si = R.rhs()[s].index;
            scalar sr = R.rhs()[s].stoichCoeff;
            c[si] += dt*sr*omeg;
            c[si] = max(0.0, c[si]);
        }        

    } // end for(label i...

    return cTauChem_/tChemInv;
}


// ************************************************************************* //
