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

Application
    findRoot

Description
    Test root-finding functions

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "BisectionRoot.H"
#include "RiddersRoot.H"
#include "NewtonSecantRoot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class testFunction
{
public:

    testFunction()
    {}

    scalar operator()(const scalar& x) const
    {
        return sqr(x) - 1 - Foam::sin(x);
    }
};


class testFunctionDerivative
{
public:

    testFunctionDerivative()
    {}

    scalar operator()(const scalar& x) const
    {
        return 0.5*Foam::pow(x,-1.5) - Foam::cos(x);
    }
};


int main(int argc, char *argv[])
{
    testFunction           tf;
    testFunctionDerivative df;

    Info<< setprecision(10)
        << "Bisection root "
        << BisectionRoot<testFunction>(tf, 1e-5).root(0, 10) << nl
        << "Ridders root "
        << RiddersRoot<testFunction>(tf, 1e-5).root(0, 10) << nl
        << "NewtonSecant root "
        << NewtonSecantRoot<testFunction,testFunctionDerivative>
           (
               tf,
               df,
               1e-5
           ).root(1.5) << endl;


    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
