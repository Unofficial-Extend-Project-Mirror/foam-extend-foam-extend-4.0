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
        return 0.5*Foam::pow(x, -1.5) - Foam::cos(x);
    }
};


int main(int argc, char *argv[])
{
    testFunction tf;
    testFunctionDerivative df;

    Info<< setprecision(10)
        << "Bisection root "
        << BisectionRoot<testFunction>(tf, 1e-5).root(0, 10) << nl
        << "Ridders root "
        << RiddersRoot<testFunction>(tf, 1e-5).root(0, 10) << nl
        << "NewtonSecant root "
        << NewtonSecantRoot<testFunction, testFunctionDerivative>
           (
               tf,
               df,
               1e-5
           ).root(1.5) << endl;


    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
