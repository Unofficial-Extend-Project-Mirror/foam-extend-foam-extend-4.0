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

#include "ODESolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::ODESolver, 0);
namespace Foam
{
    defineRunTimeSelectionTable(ODESolver, ODE);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ODESolver::ODESolver(ODE& ode)
:
    ode_(ode),
    yScale_(ode.nEqns()),
    dydx_(ode.nEqns())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ODESolver::solve
(
    const scalar xStart,
    const scalar xEnd,
    const scalar eps,
    scalar& hEst
) const
{
    const label MAXSTP = 10000;

    scalar x = xStart;
    scalar h = hEst;

    const label nEqns = ode_.nEqns();
    scalarField& y = ode_.coeffs();

    for (label nStep = 0; nStep < MAXSTP; nStep++)
    {
        ode_.derivatives(x, y, dydx_);

        for (label i=0; i<nEqns; i++)
        {
            yScale_[i] = mag(y[i]) + mag(dydx_[i]*h) + SMALL;
        }

        if ((x + h - xEnd)*(x + h - xStart) > 0.0)
        {
            h = xEnd - x;
        }

        scalar hNext, hDid;
        solve(x, y, dydx_, eps, yScale_, h, hDid, hNext);

        if ((x - xEnd)*(xEnd - xStart) >= 0.0)
        {
            hEst = hNext;

            // Solution completed.  Update ODE
            ode_.update(xEnd - xStart);
            return;
        }

        h = hNext;
    }

    FatalErrorIn
    (
        "ODESolver::solve"
        "(const scalar xStart, const scalar xEnd,"
        "scalarField& yStart, const scalar eps, scalar& hEst) const"
    )   << "Too many integration steps"
        << exit(FatalError);
}


// ************************************************************************* //
