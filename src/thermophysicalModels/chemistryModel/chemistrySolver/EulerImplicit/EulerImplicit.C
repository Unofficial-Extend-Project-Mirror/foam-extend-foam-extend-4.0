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

#include "EulerImplicit.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EulerImplicit, 0);
    addToRunTimeSelectionTable
    (
        chemistrySolver,
        EulerImplicit,
        dictionary
    );
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EulerImplicit::EulerImplicit
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

Foam::EulerImplicit::~EulerImplicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::EulerImplicit::solve
(
    scalarField &c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    label Ns = chemistry_.Ns();
    simpleMatrix<scalar> RR(Ns);
    
    for(label i=0; i<Ns; i++)
    {
        c[i] = max(0.0, c[i]);
    }

    for(label i=0; i<Ns; i++)
    {
        RR.source()[i] = c[i]/dt;
    }

    for(label i=0; i<chemistry_.reactions().size(); i++)
    {
        const chemistryModel::reaction& R = chemistry_.reactions()[i];

        scalar omegai = chemistry_.omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        scalar corr = 1.0;
        if (equil_)
        {
            if (omegai<0.0)
            {
                corr = 1.0/(1.0 + pr*dt);
            }
            else
            {
                corr = 1.0/(1.0 + pf*dt);
            }
        }

        for(label s=0; s<R.lhs().size(); s++)
        {
            label si = R.lhs()[s].index;
            scalar sl = R.lhs()[s].stoichCoeff;
            RR[si][rRef] -= sl*pr*corr;
            RR[si][lRef] += sl*pf*corr;
        }
            
        for(label s=0; s<R.rhs().size(); s++)
        {
            label si = R.rhs()[s].index;
            scalar sr = R.rhs()[s].stoichCoeff;
            RR[si][lRef] -= sr*pf*corr;
            RR[si][rRef] += sr*pr*corr;
        }
        
    } // end for(label i...

    
    for(label i=0; i<Ns; i++)
    {
        RR[i][i] += 1.0/dt;
    }

    c = RR.LUsolve();
    for(label i=0; i<Ns; i++)
    {
        c[i] = max(0.0, c[i]);
    }

    // estimate the next time step
    scalar tMin = GREAT;
    label n = chemistry_.nEqns();
    scalarField c1(n, 0.0);

    for(label i=0; i<Ns; i++)
    {
        c1[i] = c[i];
    }
    c1[Ns] = T;
    c1[Ns+1] = p;

    scalarField dcdt(n, 0.0);
    chemistry_.derivatives(0.0, c1, dcdt);
    
    scalar sumC = sum(c);

    for(label i=0; i<Ns; i++)
    {
        scalar d = dcdt[i];
        if (d < -SMALL)
        {
            tMin = min(tMin, -(c[i]+SMALL)/d);
        }
        else
        {
            d = max(d, SMALL);
            scalar cm = max(sumC - c[i], 1.0e-5);
            tMin = min(tMin, cm/d);
        }
    }    

    return cTauChem_*tMin;

}


// ************************************************************************* //
