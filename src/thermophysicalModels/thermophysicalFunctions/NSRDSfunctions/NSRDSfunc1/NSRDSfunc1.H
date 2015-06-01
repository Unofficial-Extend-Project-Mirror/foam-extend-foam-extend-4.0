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

Class
    Foam::NSRDSfunc1

Description
    NSRDS function number 101

    Source:
    @verbatim
                      NSRDS - AICHE
                 Data Compilation Tables
                    of Properties of
                     Pure Compounds

        Design Institute for Physical Property Data
          American Institute of Chemical Engineers
                  345 East 47th Street
                New York, New York 10017

         National Standard Reference Data System
         American Institute of Chemical Engineers

          T.E. Daubert       -       R.P. Danner

            Department of Chemical Engineering
            The Pennsylvania State University
                University Park, PA 16802
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef NSRDSfunc1_H
#define NSRDSfunc1_H

#include "thermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class NSRDSfunc1 Declaration
\*---------------------------------------------------------------------------*/

class NSRDSfunc1
:
    public thermophysicalFunction
{
    // Private data

        // NSRDS function 101 coefficients
        scalar a_, b_, c_, d_, e_;


public:

    //- Runtime type information
    TypeName("NSRDSfunc1");


    // Constructors

        //- Construct from components
        NSRDSfunc1(scalar a, scalar b, scalar c, scalar d, scalar e)
        :
            a_(a),
            b_(b),
            c_(c),
            d_(d),
            e_(e)
        {}

        //- Construct from Istream
        NSRDSfunc1(Istream& is)
        :
            a_(readScalar(is)),
            b_(readScalar(is)),
            c_(readScalar(is)),
            d_(readScalar(is)),
            e_(readScalar(is))
        {}


    // Member Functions

        //- Evaluate the function and return the result
        scalar f(scalar, scalar T) const
        {
            return exp(a_ + b_/T + c_*log(T) + d_*pow(T, e_));
        }


        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            os  << a_ << token::SPACE
                << b_ << token::SPACE
                << c_ << token::SPACE
                << d_ << token::SPACE
                << e_;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const NSRDSfunc1& f)
        {
            f.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
