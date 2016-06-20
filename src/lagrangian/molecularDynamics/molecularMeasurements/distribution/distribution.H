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
    Foam::distribution

Description

SourceFiles
    distributionI.H
    distribution.C

\*---------------------------------------------------------------------------*/

#ifndef distribution_H
#define distribution_H

#include "Map.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class distribution Declaration
\*---------------------------------------------------------------------------*/

class distribution
:
    public Map<label>
{
    // Private data

        scalar binWidth_;


public:

    // Constructors

        //- Construct null
        distribution();

        //- Construct from binWidth
        distribution(const scalar binWidth);

        //- Construct as copy
        distribution(const distribution&);


    // Destructor

        ~distribution();


    // Member Functions

        label totalEntries() const;

        scalar approxTotalEntries() const;

        scalar mean() const;

        scalar median();

        //- Add a value to the appropriate bin of the distribution.
        void add(const scalar valueToAdd);

        void add(const label valueToAdd);

        void insertMissingKeys();

        List<Pair<scalar> > normalised();

        List<Pair<scalar> > normalisedMinusMean();

        List<Pair<scalar> > normalisedShifted(const scalar shiftValue);

        List<Pair<scalar> > raw();


        // Access

            inline scalar binWidth() const;


    // Member Operators

        void operator=(const distribution&);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const distribution&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "distributionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
