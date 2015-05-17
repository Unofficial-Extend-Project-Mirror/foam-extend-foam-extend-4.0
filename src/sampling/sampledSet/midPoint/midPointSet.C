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

#include "midPointSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(midPointSet, 0);
    addToRunTimeSelectionTable(sampledSet, midPointSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Rework faceOnlySet samples.
// Take two consecutive samples
void Foam::midPointSet::genSamples()
{
    // Generate midpoints.
    if(size() == 0)
    {
        // Nothing to do

        return;
    }

    List<point> midPoints(2*size());
    labelList midCells(2*size());
    labelList midSegments(2*size());
    scalarList midCurveDist(2*size());

    label midI = 0;

    label sampleI = 0;

    while(true)
    {
        // calculate midpoint between sampleI and sampleI+1
        // (if in same segment)
        while
        (
            (sampleI < size() - 1)
         && (segments_[sampleI] == segments_[sampleI+1])
        )
        {
            midPoints[midI] =
                0.5*(operator[](sampleI) + operator[](sampleI+1));

            label cell1 = getCell(faces_[sampleI], midPoints[midI]);
            label cell2 = getCell(faces_[sampleI+1], midPoints[midI]);

            if (cell1 != cell2)
            {
                FatalErrorIn("midPointSet::genSamples()")
                    << "  sampleI:" << sampleI
                    << "  midI:" << midI
                    << "  sampleI:" << sampleI
                    << "  pts[sampleI]:" << operator[](sampleI)
                    << "  face[sampleI]:" << faces_[sampleI]
                    << "  pts[sampleI+1]:" << operator[](sampleI+1)
                    << "  face[sampleI+1]:" << faces_[sampleI+1]
                    << "  cell1:" << cell1
                    << "  cell2:" << cell2
                    << abort(FatalError);
            }

            midCells[midI] = cell1;
            midSegments[midI] = segments_[sampleI];
            midCurveDist[midI] = mag(midPoints[midI] - start());

            midI++;
            sampleI++;
        }

        if (sampleI == size() - 1)
        {
            break;
        }
        sampleI++;
    }

    midPoints.setSize(midI);
    midCells.setSize(midI);
    midSegments.setSize(midI);
    midCurveDist.setSize(midI);
    setSamples
    (
        midPoints,
        midCells,
        labelList(midCells.size(), -1),
        midSegments,
        midCurveDist
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::midPointSet::midPointSet
(
    const word& name,
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const word& axis,
    const point& start,
    const point& end
)
:
    faceOnlySet(name, mesh, searchEngine, axis, start, end)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


Foam::midPointSet::midPointSet
(
    const word& name,
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict
)
:
    faceOnlySet(name, mesh, searchEngine, dict)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::midPointSet::~midPointSet()
{}


// ************************************************************************* //
