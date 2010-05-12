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

#include "ode.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ode, 0);
    addToRunTimeSelectionTable
    (
        chemistrySolver,
        ode,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::ode::ode
(
    const Foam::dictionary& dict,
    Foam::chemistryModel& chemistry
)
:
    chemistrySolver(dict, chemistry),
    chemistry_(chemistry),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    solverName_(coeffsDict_.lookup("ODESolver")),
    odeSolver_(ODESolver::New(solverName_, chemistry)),
    eps_(readScalar(coeffsDict_.lookup("eps"))),
    scale_(readScalar(coeffsDict_.lookup("scale")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ode::~ode()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar Foam::ode::solve
(
    scalarField& c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{
    label Ns = chemistry_.Ns();
    scalarField& c1 = chemistry_.coeffs();

    // copy the concentration, T and P to the total solve-vector
    for(label i=0; i<Ns; i++)
    {
        c1[i] = c[i];
    }
    c1[Ns] = T;
    c1[Ns+1] = p;

    scalar dtEst = dt;

    odeSolver_->solve
    (
        t0,
        t0 + dt,
        eps_,
        dtEst
    );

    for(label i=0; i<c.size(); i++)
    {
        c[i] = max(0.0, c1[i]);
    }

    return dtEst;
}


// ************************************************************************* //
