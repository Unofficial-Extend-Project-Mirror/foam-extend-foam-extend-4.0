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
    Foam::mixture

Description

SourceFiles
    mixtureI.H
    mixture.C
    mixtureIO.C

\*---------------------------------------------------------------------------*/

#ifndef mixture_H
#define mixture_H

#include "error.H"

#include "List.H"
#include "substance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixture Declaration
\*---------------------------------------------------------------------------*/

class mixture
:
    public List<substance>
{
    // Private data

        word name_;


public:

    // Constructors

        //- Construct from Istream
        mixture(Istream& is)
        :
            List<substance>(is),
            name_(is)
        {
            scalar volTot = 0;

            for (label i = 0; i < size(); i++)
            {
                volTot += operator[](i).volFrac();
            }

            if (volTot > 1.001 || volTot < 0.999)
            {
                FatalErrorIn("mixture::mixture(istream& is)")
                    << "Sum of volume fractions for Mixture " << name_
                    << " = " << volTot << endl
                    << "Should equal one."
                    << abort(FatalError);
            }
        }


    // Member Functions

        // Access

            const word& name() const
            {
                return name_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
