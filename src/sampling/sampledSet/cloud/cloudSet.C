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

#include "cloudSet.H"
#include "sampledSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cloudSet, 0);
    addToRunTimeSelectionTable(sampledSet, cloudSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cloudSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    forAll(sampleCoords_, sampleI)
    {
        label cellI = searchEngine().findCell(sampleCoords_[sampleI]);

        if (cellI != -1)
        {
            samplingPts.append(sampleCoords_[sampleI]);
            samplingCells.append(cellI);
            samplingFaces.append(-1);
            samplingSegments.append(0);
            samplingCurveDist.append(1.0 * sampleI);
        }
    }
}


void Foam::cloudSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloudSet::cloudSet
(
    const word& name,
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const word& axis,
    const List<point>& sampleCoords
)
:
    sampledSet(name, mesh, searchEngine, axis),
    sampleCoords_(sampleCoords)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


Foam::cloudSet::cloudSet
(
    const word& name,
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    sampleCoords_(dict.lookup("points"))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cloudSet::~cloudSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::cloudSet::getRefPoint(const List<point>& pts) const
{
    if (pts.size())
    {
        // Use first samplePt as starting point
        return pts[0];
    }
    else
    {
        return vector::zero;
    }
}


// ************************************************************************* //
