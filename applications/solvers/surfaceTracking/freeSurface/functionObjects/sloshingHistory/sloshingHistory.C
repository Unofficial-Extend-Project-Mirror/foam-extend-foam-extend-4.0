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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "sloshingHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sloshingHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sloshingHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sloshingHistory::sloshingHistory
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_(NULL),
    freeSurfacePatchID_(-1),
    leftPointID_(-1),
    rightPointID_(-1)
{
    if (Pstream::parRun())
    {
        FatalErrorIn("sloshingHistory::sloshingHistory(...)")
            << "sloshingHistory objec function "
                << "is not implemented for parallel run"
                << abort(FatalError);
    }

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    if (Pstream::master())
    {
        if (!time_.processorCase())
        {
            mkDir
            (
                time_.path()
               /"history"
               /time_.timeName()
            );

            historyFilePtr_ =
                new OFstream
                (
                    time_.path()
                   /"history"
                   /time_.timeName()
                   /"history.dat"
                );
        }
        else
        {
            mkDir
            (
                time_.path()/".."/"history"
               /time_.timeName()
            );

            historyFilePtr_ =
                new OFstream
                (
                    time_.path()/".."
                   /"history"
                   /time_.timeName()
                   /"history.dat"
                );
        }

        (*historyFilePtr_)
            << "Time" << tab
                << "Yleft" << tab
                << "Yright" << endl;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    freeSurfacePatchID_ =
        mesh.boundaryMesh().findPatchID("freeSurface");

    if (freeSurfacePatchID_ == -1)
    {
        FatalErrorIn("dropHistory::dropHistory(...)")
            << "Problem with finding freeSurface patch"
                << abort(FatalError);
    }

    const vectorField& fsPoints =
        mesh.boundaryMesh()[freeSurfacePatchID_].localPoints();

    scalar maxX = 0;
    scalar minX = GREAT;

    forAll(fsPoints, pointI)
    {
        if (fsPoints[pointI].x() < minX)
        {
            minX = fsPoints[pointI].x();
            leftPointID_ = pointI;
        }

        if (fsPoints[pointI].x() > maxX)
        {
            maxX = fsPoints[pointI].x();
            rightPointID_ = pointI;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sloshingHistory::start()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const vectorField& fsPoints =
        mesh.boundaryMesh()[freeSurfacePatchID_].localPoints();

    if (Pstream::master())
    {
        historyFilePtr_->precision(12);

        (*historyFilePtr_) << time_.value() << tab
            << fsPoints[leftPointID_].y() << tab
            << fsPoints[rightPointID_].y() << endl;

        return true;
    }

    return false;
}


bool Foam::sloshingHistory::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const vectorField& fsPoints =
        mesh.boundaryMesh()[freeSurfacePatchID_].localPoints();

    if (Pstream::master() && time_.outputTime())
    {
        historyFilePtr_->precision(12);

        (*historyFilePtr_) << time_.value() << tab
            << fsPoints[leftPointID_].y() << tab
            << fsPoints[rightPointID_].y() << endl;

        return true;
    }

    return false;
}


bool Foam::sloshingHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
