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
    Foam::ODESolver

Description
    Selection for ODE solver

SourceFiles
    ODESolver.C

\*---------------------------------------------------------------------------*/

#ifndef ODESolver_H
#define ODESolver_H

#include "ODE.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ODESolver Declaration
\*---------------------------------------------------------------------------*/

class ODESolver
{

protected:

    // Private data

        //- Reference to ODE
        ODE& ode_;

        mutable scalarField yScale_;
        mutable scalarField dydx_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ODESolver(const ODESolver&);

        //- Disallow default bitwise assignment
        void operator=(const ODESolver&);


public:

    //- Runtime type information
    TypeName("ODESolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ODESolver,
            ODE,
            (ODE& ode),
            (ode)
        );


    // Constructors

        //- Construct for given ODE
        ODESolver(ODE& ode);


    // Selectors

        //- Select null constructed
        static autoPtr<ODESolver> New
        (
            const word& ODESolverTypeName,
            ODE& ode
        );


    // Destructor

        virtual ~ODESolver()
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
        ) const = 0;


        virtual void solve
        (
            const scalar xStart,
            const scalar xEnd,
            const scalar eps,
            scalar& hEst
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
