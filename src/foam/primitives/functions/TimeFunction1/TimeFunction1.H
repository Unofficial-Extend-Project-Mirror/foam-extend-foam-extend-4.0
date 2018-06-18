/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::TimeFunction1

Description
    Light wrapper around Function1 to provide a mechanism to update time-based
    entries.

SourceFiles
    TimeFunction1.C

\*---------------------------------------------------------------------------*/

#ifndef TimeFunction1_H
#define TimeFunction1_H

#include "Function1.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<class Type>
class TimeFunction1;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const TimeFunction1<Type>&
);

/*---------------------------------------------------------------------------*\
                       Class TimeFunction1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TimeFunction1
{

protected:

    // Protected data

        //- Reference to the time database
        const Time& time_;

        //- Name of the data entry
        const word name_;

        //- The underlying Function1
        autoPtr<Function1<Type> > entry_;


public:

    // Constructor

        //- Construct from entry name
        TimeFunction1
        (
            const Time& t,
            const word& name,
            const dictionary& dict
        );

        //- Construct null from entry name
        TimeFunction1
        (
            const Time& t,
            const word& entryName
        );

        //- Copy constructor
        TimeFunction1(const TimeFunction1<Type>& tde);


    //- Destructor
    virtual ~TimeFunction1();


    // Member Functions

        // Access

            //- Reset entry by re-reading from dictionary
            void reset(const dictionary& dict);

            //- Return the name of the entry
            const word& name() const;


        // Evaluation

            //- Return value as a function of (scalar) independent variable
            virtual Type value(const scalar x) const;

            //- Integrate between two (scalar) values
            virtual Type integrate(const scalar x1, const scalar x2) const;


    // I/O

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const TimeFunction1<Type>& de
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TimeFunction1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
