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

\*---------------------------------------------------------------------------*/

#include "gnuplotSetWriter.H"
#include "clock.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::gnuplotSetWriter<Type>::gnuplotSetWriter()
:
    writer<Type>()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::gnuplotSetWriter<Type>::~gnuplotSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::gnuplotSetWriter<Type>::getFileName
(
    const coordSet& points,
    const wordList& valueSetNames
) const
{
    return this->getBaseName(points, valueSetNames) + ".gplt";
}


template<class Type>
void Foam::gnuplotSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    os  << "set term postscript color" << nl
        << "set output \"" << points.name() << ".ps\"" << nl
        << "plot";

    forAll(valueSets, i)
    {
        if (i != 0)
        {
            os << ',';
        }

        os  << " \"-\" title \"" << valueSetNames[i] << "\" with lines";
    }
    os  << nl;


    forAll(valueSets, i)
    {
        this->writeTable(points, *valueSets[i], os);
        os  << "e" << nl;
    }
}


template<class Type>
void Foam::gnuplotSetWriter<Type>::write
(
    const bool writeTracks,
    const PtrList<coordSet>& trackPoints,
    const wordList& valueSetNames,
    const List<List<Field<Type> > >& valueSets,
    Ostream& os
) const
{
    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorIn("gnuplotSetWriter<Type>::write(..)")
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }
    if (trackPoints.size() > 0)
    {
        os  << "set term postscript color" << nl
            << "set output \"" << trackPoints[0].name() << ".ps\"" << nl;

        forAll(trackPoints, trackI)
        {
            os  << "plot";

            forAll(valueSets, i)
            {
                if (i != 0)
                {
                    os << ',';
                }

                os  << " \"-\" title \"" << valueSetNames[i] << "\" with lines";
            }
            os << nl;

            forAll(valueSets, i)
            {
                this->writeTable(trackPoints[trackI], valueSets[i][trackI], os);
                os  << "e" << nl;
            }
        }
    }
}


// ************************************************************************* //
