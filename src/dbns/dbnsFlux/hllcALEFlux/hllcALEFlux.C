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

\*---------------------------------------------------------------------------*/

#include "hllcALEFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllcALEFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const scalar& RLeft,
    const scalar& RRight,
    const scalar& CvLeft,
    const scalar& CvRight,
    const vector& Sf,
    const scalar& magSf,
    const vector& dotX
) const
{
    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // Ratio of specific heat capacities
    const scalar kappaLeft = (RLeft + CvLeft)/CvLeft;
    const scalar kappaRight = (RRight + CvRight)/CvRight;

    // Compute conservative variables assuming perfect gas law

    // Density
    const scalar rhoLeft = pLeft/(RLeft*TLeft);
    const scalar rhoRight = pRight/(RRight*TRight);

    // DensityVelocity
    const vector rhoULeft = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // DensityTotalEnergy
    const scalar rhoELeft = rhoLeft*(CvLeft*TLeft + 0.5*magSqr(ULeft));
    const scalar rhoERight = rhoRight*(CvRight*TRight + 0.5*magSqr(URight));

    // Compute left and right total enthalpies:
    const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
    const scalar HRight = (rhoERight + pRight)/rhoRight;

    // Compute velocity relative to mesh
    const vector URelLeft = ULeft - dotX;
    const vector URelRight = URight - dotX;

    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = (URelLeft & normalVector);
    const scalar qRight = (URelRight & normalVector);

    // Speed of sound, for left and right side, assuming perfect gas
    const scalar aLeft =
        Foam::sqrt(max(0 , kappaLeft*pLeft/rhoLeft));

    const scalar aRight =
        Foam::sqrt(max(0, kappaRight*pRight/rhoRight));


    // Step 2:
    // needs rho_{l,r}, U_{l,r}, H_{l,r}, kappa_{l,r}, Gamma_{l,r}, q_{l,r}
    // compute Roe weights
    const scalar rhoLeftSqrt = Foam::sqrt(max(0.0,rhoLeft));
    const scalar rhoRightSqrt = Foam::sqrt(max(0.0,rhoRight));

    const scalar wLeft = rhoLeftSqrt
        /stabilise((rhoLeftSqrt + rhoRightSqrt),VSMALL);
    const scalar wRight = 1 - wLeft;

    // Roe averaged velocity
    const vector UTilde = wLeft*ULeft + wRight*URight;

    // Roe averaged relative face velocity
    const scalar contrURelTilde = ((UTilde - dotX) & normalVector);

    // Roe averaged contravariant velocity
    const scalar contrUTilde = (UTilde & normalVector);

    // Roe averaged total enthalpy
    const scalar HTilde = wLeft*HLeft + wRight*HRight;

    // Roe averaged kappa
    // TODO: needs to be verified!
    const scalar kappaTilde = wLeft*kappaLeft + wRight*kappaRight;

    // Speed of sound with Roe reconstruction values
    // TODO: not sure if the correct (flow speed) and kappa is used here
    const scalar aTilde =
        Foam::sqrt(max(0.0,(kappaTilde-1.0)*(HTilde-0.5*sqr(contrUTilde))));

    // Step 3: compute signal speeds for face:
    const scalar SLeft  = min(qLeft-aLeft,   contrURelTilde-aTilde);
    const scalar SRight = max(qRight+aRight, contrURelTilde+aTilde);

    const scalar SStar = (rhoRight*qRight*(SRight-qRight)
    - rhoLeft*qLeft*(SLeft - qLeft) + pLeft - pRight )/
        stabilise((rhoRight*(SRight-qRight)-rhoLeft*(SLeft-qLeft)),VSMALL);

    // Compute pressure in star region from the right side
    const scalar pStarRight =
        rhoRight*(qRight - SRight)*(qRight - SStar) + pRight;

    // Should be equal to the left side
    const scalar pStarLeft  =
        rhoLeft*(qLeft -  SLeft)*(qLeft - SStar) + pLeft;

    // Give a warning if this is not the case
    if ( mag(pStarRight-pStarLeft) > 1e-6 )
    {
        Info << "mag(pStarRight-pStarLeft) > VSMALL " << endl;
    }

    // Use pStarRight for pStar, as in theory, pStarRight == pStarLeft
    const scalar pStar = pStarRight;

    // Step 4: upwinding - compute states:
    scalar convectionSpeed = 0.0;
    scalar rhoState = 0.0;
    vector rhoUState = vector::zero;
    scalar rhoEState = 0.0;
    scalar pState = 0.0;

    // TODO: Maybe one can use pos/neg implementation, but then one has to
    // evaluate all 4 states at each iteration!
    // label A = pos(SLeft);
    // label B = pos(SStar);
    // label C = pos(SRight);
    // please double check the bool operators again, if one want's to
    // implement this!!!
    // scalar convectionSpeed = A*B*C*qLeft+(1-A)*B*C*SStar
    //     +(1-A)*(1-B)*C*SStar+(1-A)*(1-B)*(1-C)*qRight:

    if (pos(SLeft))
    {
        // Compute F_l
        convectionSpeed = qLeft;
        rhoState  = rhoLeft;
        rhoUState = rhoULeft;
        rhoEState = rhoELeft;
        pState = pLeft;
    }
    else if (pos(SStar))
    {
        scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar), VSMALL);

        // Compute left star region
        convectionSpeed = SStar;
        rhoState  = omegaLeft*(SLeft - qLeft)*rhoLeft;
        rhoUState = omegaLeft*((SLeft - qLeft)*rhoULeft
        + (pStar - pLeft)*normalVector);
        rhoEState = omegaLeft*((SLeft - qLeft)*rhoELeft
        - pLeft*qLeft + pStar*SStar);
        pState = pStar;
    }
    else if (pos(SRight))
    {
        scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar), VSMALL);

        // Compute right star region
        convectionSpeed = SStar;
        rhoState  = omegaRight*(SRight - qRight)*rhoRight;
        rhoUState = omegaRight*((SRight - qRight)*rhoURight
        + (pStar - pRight)*normalVector);
        rhoEState = omegaRight*((SRight - qRight)*rhoERight
        - pRight*qRight + pStar*SStar);
        pState = pStar;
    }
    else if (neg(SRight))
    {
        // Compute F_r
        convectionSpeed = qRight;
        rhoState  = rhoRight;
        rhoUState = rhoURight;
        rhoEState = rhoERight;
        pState = pRight;
    }
    else
    {
        Info << "Error in HLLC Riemann solver" << endl;
    }

    rhoFlux  = (convectionSpeed*rhoState)*magSf;
    rhoUFlux = (convectionSpeed*rhoUState+pState*normalVector)*magSf;
    rhoEFlux = (convectionSpeed*(rhoEState+pState)
    + pState*(dotX & normalVector))*magSf;
}


// ************************************************************************* //
