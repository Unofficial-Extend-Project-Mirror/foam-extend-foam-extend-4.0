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

#include "betaFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::betaFlux::evaluateFlux
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
    const scalar& magSf
) const
{
    // Step 1: decode rho left and right:
    scalar rhoLeft = pLeft/(RLeft*TLeft);
    scalar rhoRight = pRight/(RRight*TRight);

    // Decode left and right total energy:
    scalar eLeft = CvLeft*TLeft+0.5*magSqr(ULeft);
    scalar eRight = CvRight*TRight+0.5*magSqr(URight);

    // Adiabatic exponent is constant for ideal gas but if Cp=Cp(T)
    // it must be computed for each cell and evaluated at each face
    // through reconstruction
    const scalar kappaLeft = (CvLeft+RLeft)/CvLeft;
    const scalar kappaRight = (CvRight+RRight)/CvRight;

    // normal vector
    vector normalVector = Sf/magSf;

    // Compute left and right contravariant velocities:
    const scalar contrVLeft  = (ULeft & normalVector);
    const scalar contrVRight = (URight & normalVector);

    // Compute left and right total enthalpies:
    const scalar hLeft = eLeft + pLeft/rhoLeft;
    const scalar hRight = eRight + pRight/rhoRight;

    // Compute left and right velocity square
    const scalar qLeftSquare  = magSqr(ULeft);
    const scalar qRightSquare = magSqr(URight);

    // compute left and right speed of sound
    const scalar cLeft =
        sqrt(max((kappaLeft - 1)*(hLeft - 0.5*qLeftSquare), SMALL));
    const scalar cRight =
        sqrt(max((kappaRight - 1)*(hRight - 0.5*qRightSquare), SMALL));

    const scalar magULeft  = mag(ULeft);
    const scalar magURight = mag(URight);

    const scalar MLeft = mag(magULeft/cLeft);
    const scalar MRight = mag(magURight/cRight);

    // Compute beta parameter - this should be done in multidimensional way
    // similarly to multidimensional limiters
    const scalar Mmax = max(MLeft, MRight);
    const scalar Mmin = min(MLeft, MRight);

    scalar beta = 0;

    if (Mmin > 0)
    {
        scalar alpha = 0.1;
        scalar r = sqrt(sqr(Mmax) - sqr(Mmin))/Mmin;

        if (r >= alpha)
        {
            beta = min(r,1);
        }
    }

    // Step 2: compute Roe averaged quantities for face:
    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // Some temporary variables:
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
    const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1 - wLeft;

    const vector UTilde = ULeft*wLeft + URight*wRight;
    const scalar hTilde = hLeft*wLeft + hRight*wRight;
    const scalar qTildeSquare = magSqr(UTilde);
    const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;

    // Speed of sound
    const scalar cTilde =
        sqrt(max((kappaTilde - 1)*(hTilde - 0.5*qTildeSquare), SMALL));

    // Roe averaged contravariant velocity
    const scalar contrVTilde = (UTilde & normalVector);

    // Step 3: compute primitive differences:
    const scalar deltaP = pRight - pLeft;
    const scalar deltaRho = rhoRight - rhoLeft;
    const vector deltaU = URight - ULeft;
    const scalar deltaContrV = (deltaU & normalVector);

    // Step 4: compute wave strengths:

    // Roe and Pike - formulation
    const scalar r1 =
        (deltaP - rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));
    const scalar r2 = deltaRho - deltaP/sqr(cTilde);
    const scalar r3 =
        (deltaP + rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));

    // Step 5: compute l vectors

    // rho row:
    const scalar l1rho = 1;
    const scalar l2rho = 1;
    const scalar l3rho = 0;
    const scalar l4rho = 1;

    // first U column
    const vector l1U = UTilde - cTilde*normalVector;

    // second U column
    const vector l2U = UTilde;

    // third U column
    const vector l3U = deltaU - deltaContrV*normalVector;

    // fourth U column
    const vector l4U = UTilde + cTilde*normalVector;

    // E row
    const scalar l1e = hTilde - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU) - contrVTilde*deltaContrV;
    const scalar l4e = hTilde + cTilde*contrVTilde;

    // Step 6: compute eigenvalues

    scalar lambda1 = mag(contrVTilde - cTilde);
    scalar lambda2 = mag(contrVTilde);
    scalar lambda3 = mag(contrVTilde + cTilde);

    scalar lambdaMax = max(max(lambda1, lambda2), lambda3);

    scalar lambda1Max = beta*lambda1 + (1 - beta)*lambdaMax;
    scalar lambda2Max = beta*lambda2 + (1 - beta)*lambdaMax;
    scalar lambda3Max = beta*lambda3 + (1 - beta)*lambdaMax;

    // Step 7: check for Harten entropy correction


    // Components of deltaF1
    const scalar diffF11  = lambda1Max*r1*l1rho;
    const vector diffF124 = lambda1Max*r1*l1U;
    const scalar diffF15  = lambda1Max*r1*l1e;

    // Components of deltaF2
    const scalar diffF21  = lambda2Max*(r2*l2rho + rhoTilde*l3rho);
    const vector diffF224 = lambda2Max*(r2*l2U + rhoTilde*l3U);
    const scalar diffF25  = lambda2Max*(r2*l2e + rhoTilde*l3e);

    // Components of deltaF3
    const scalar diffF31  = lambda3Max*r3*l4rho;
    const vector diffF324 = lambda3Max*r3*l4U;
    const scalar diffF35  = lambda3Max*r3*l4e;

    // Step 8: compute left and right fluxes

    // Left flux 5-vector
    const scalar fluxLeft11 = rhoLeft*contrVLeft;
    const vector fluxLeft124 = ULeft*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15 = hLeft*fluxLeft11;

    // Right flux 5-vector
    const scalar fluxRight11 = rhoRight*contrVRight;
    const vector fluxRight124 = URight*fluxRight11 + normalVector*pRight;
    const scalar fluxRight15 = hRight*fluxRight11;

    // Step 9: compute face flux 5-vector
    const scalar flux1 =
        0.5*(fluxLeft11 + fluxRight11 - (diffF11 + diffF21 + diffF31));
    const vector flux24 =
        0.5*(fluxLeft124 + fluxRight124 - (diffF124 + diffF224 + diffF324));
    const scalar flux5 =
        0.5*(fluxLeft15 + fluxRight15 - (diffF15 + diffF25 + diffF35));

    // Compute private data
    rhoFlux  = flux1*magSf;
    rhoUFlux = flux24*magSf;
    rhoEFlux = flux5*magSf;
}

// ************************************************************************* //
