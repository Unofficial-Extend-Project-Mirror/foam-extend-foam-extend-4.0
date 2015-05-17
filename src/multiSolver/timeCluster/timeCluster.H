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
    Foam::timeCluster

Description
    A cluster of time instances, used with multiSolver.  Catalogues all the
    time directories within a superLoop, within a solverDomain.  A single
    time cluster describes what would be analogous to all the time directories
    for a solver that does not use multiSolver.

SourceFiles
    timeCluster.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef timeCluster_H
#define timeCluster_H

#include "word.H"
#include "scalar.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class timeCluster;

// Friend IOstream Operators

Istream& operator>>(Istream&, timeCluster&);
Ostream& operator<<(Ostream&, const timeCluster&);


/*---------------------------------------------------------------------------*\
                        Class timeCluster Declaration
\*---------------------------------------------------------------------------*/

class timeCluster
:
    public instantList
{
    // Private data

        scalar globalOffset_;
        label globalIndex_;
        label superLoop_;
        word solverDomainName_;

        // PreConditioned file name - used only for pre and post conditioning
        // a data set for parallel runs
        word preConName_;

public:

    // Public classes

        //- Less function class used in sorting timeClusters
        class less
        {
        public:

            bool operator()(const timeCluster& one, const timeCluster& two) const
            {
                return one.globalMinValue() < two.globalMinValue();
            }
        };


    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct null
        timeCluster();

        //- Construct from components
        timeCluster
        (
            const instantList& times,
            const scalar globalOffset,
            const label globalIndex,
            const label superLoop,
            const word& solverDomainName,
            const word& preConName = word::null
        );

        //- Construct given a timeCluster and an index
        //  This creates a 'timeCluster' that holds a single instant
        //  in time, whose other values match those of tc.
        timeCluster
        (
            const timeCluster& tc,
            const label index
        );

        //- Construct a time cluster given a scalar value.  This constructs
        //  a timeCluster with a single instant of time at value t, and whose
        //  other values are zero or empty.
        explicit timeCluster(const scalar t);

    // Member Functions

        // Access

            //- Times
            inline const instantList& times() const;

            inline instantList& times();

            //- Global offset
            inline scalar globalOffset() const;

            inline scalar& globalOffset();

            //- Global index
            inline label globalIndex() const;

            inline label& globalIndex();

            //- SuperLoop
            inline label superLoop() const;

            inline label& superLoop();

            //- Solver domain name
            inline const word& solverDomainName() const;

            inline word& solverDomainName();

            //- Solver domain name
            inline const word& preConName() const;

            inline word& preConName();

        // Derived values

            //- Global value at index
            scalar globalValue(const label& index) const;

            //- Search for and return global min value.  If empty,
            //  returns VGREAT.
            scalar globalMinValue() const;

            //- Search for and return global max value.  If empty,
            //  returns 0.
            scalar globalMaxValue() const;

            //- Search for and return index of global min value
            label globalMinIndex() const;

            //- Search for and return index of global max value
            label globalMaxIndex() const;

            //- Global closest time
            scalar globalFindClosestTimeValue(const scalar) const;

            //- Global closest time index
            label globalFindClosestTimeIndex(const scalar) const;

            //- Local value at index
            scalar localValue(const label& index) const;

            //- Search for and return local min value
            scalar localMinValue() const;

            //- Search for and return local max value
            scalar localMaxValue() const;

            //- Search for and return index of global min value.  If empty,
            //  returns VGREAT
            label localMinIndex() const;

            //- Search for and return index of global max value.  If empty,
            //  returns 0.
            label localMaxIndex() const;

            //-Local closest time
            scalar localFindClosestTimeValue(const scalar) const;

            //-Local closest time index
            label localFindClosestTimeIndex(const scalar) const;

        // Operators

            //- Chip off a single instant given an index, return as timeCluster
                // This retains superLoop, solverDomain, etc.. with a timevalue
            timeCluster operator()(const label index) const;


    // Friend IOstream Operators

        friend Istream& operator>>(Istream&, timeCluster&);
        friend Ostream& operator<<(Ostream&, const timeCluster&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeClusterI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
