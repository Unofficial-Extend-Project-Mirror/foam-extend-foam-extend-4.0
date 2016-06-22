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
    tetMetric

Description
    Base class for tetrahedral mesh-quality metrics.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef tetMetric_H
#define tetMetric_H

#include "point.H"
#include "scalar.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class word;
class dictionary;

class tetMetric
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetMetric(const tetMetric&);

        //- Disallow default bitwise assignment
        void operator=(const tetMetric&);


public:

    // Typedef for tetrahedral metrics

        typedef scalar (*tetMetricReturnType)
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );


    // Selector

        static tetMetricReturnType New
        (
            const dictionary& dict,
            const word& metricName
        );


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            scalar,
            tetMetric,
            metric,
            Point,
            (
                const point& p0,
                const point& p1,
                const point& p2,
                const point& p3
            ),
            (metricName)
        );

        static scalar metric
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );

    // Destructor

        virtual ~tetMetric()
        {}
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
