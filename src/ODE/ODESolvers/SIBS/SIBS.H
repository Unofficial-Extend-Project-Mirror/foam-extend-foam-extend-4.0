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

Class
    SIBS

Description
    Semi-Implicit Burlisch-Stoer scheme for stiff equation sets
    Numerical Recipes in C, Secn 16.6 page 742-747
    Alternative reference in Numerical Recipes in C++

SourceFiles
    SIBS.C
    SIMPR.C
    polyExtrapolate.C

\*---------------------------------------------------------------------------*/

#ifndef SIBS_H
#define SIBS_H

#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SIBS Declaration
\*---------------------------------------------------------------------------*/

class SIBS
:
    public ODESolver
{
    // Private data

        static const label kMaxX_ = 7, iMaxX_ = kMaxX_ + 1;
        static const label nSeq_[iMaxX_];

        static const scalar safe1, safe2, redMax, redMin, scaleMX;

        mutable scalarField a_;
        mutable scalarSquareMatrix alpha_;
        mutable scalarRectangularMatrix d_p_;
        mutable scalarField x_p_;
        mutable scalarField err_;

        mutable scalarField yTemp_;
        mutable scalarField ySeq_;
        mutable scalarField yErr_;
        mutable scalarField dfdx_;
        mutable scalarSquareMatrix dfdy_;

        mutable label first_, kMax_, kOpt_;
        mutable scalar epsOld_, xNew_;


    // Private Member Functions

        void SIMPR
        (
            const scalar xStart,
            const scalarField& y,
            const scalarField& dydx,
            const scalarField& dfdx,
            const scalarSquareMatrix& dfdy,
            const scalar deltaX,
            const label nSteps,
            scalarField& yEnd
        ) const;

        void polyExtrapolate
        (
            const label iest,
            const scalar xest,
            const scalarField& yest,
            scalarField& yz,
            scalarField& dy,
            scalarField& x_p,
            scalarRectangularMatrix& d_p
        ) const;


public:

    //- Runtime type information
    TypeName("SIBS");


    // Constructors

        //- Construct from ODE
        SIBS(ODE& ode);


    // Destructor

        virtual ~SIBS()
        {}


    // Member Functions

        virtual void solve
        (
            scalar& x,
            scalarField& y,
            scalarField& dydx,
            const scalar eps,
            const scalarField& yScale,
            const scalar hTry,
            scalar& hDid,
            scalar& hNext
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
