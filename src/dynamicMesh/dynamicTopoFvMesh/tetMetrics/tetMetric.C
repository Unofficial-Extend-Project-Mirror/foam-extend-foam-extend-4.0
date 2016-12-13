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
    tetMetrics

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "tetMetric.H"
#include "word.H"
#include "dictionary.H"
#include "dlLibraryTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineMemberFunctionSelectionTable(tetMetric, metric, Point);

// * * * * * * * * * * * * * Static Members Functions * * * * * * * * * * *  //

tetMetric::tetMetricReturnType
tetMetric::New
(
    const dictionary& dict,
    const word& metricName
)
{
    Info << "Selecting metric " << metricName << endl;

    dlLibraryTable::open
    (
        dict,
        "tetMetricLibs",
        metricPointMemberFunctionTablePtr_
    );

    if (!metricPointMemberFunctionTablePtr_)
    {
        FatalErrorIn
        (
            "tetMetric::New(const dictionary&, const word&)"
        )   << "tetMetric table is empty"
            << exit(FatalError);
    }

    metricPointMemberFunctionTable::iterator mfIter =
        metricPointMemberFunctionTablePtr_->find(metricName);

    if (mfIter == metricPointMemberFunctionTablePtr_->end())
    {
        FatalErrorIn
        (
            "tetMetric::New(const dictionary&, const word&)"
        )   << "Unknown metric " << metricName
            << endl << endl
            << "Valid metrics are :" << endl
            << metricPointMemberFunctionTablePtr_->toc()
            << exit(FatalError);
    }

    return mfIter();
}

} // End namespace Foam

// ************************************************************************* //
