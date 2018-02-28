/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "movingOversetRegion.H"
#include "transformField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::movingOversetRegion::calcMotionMask() const
{
    if (motionMaskPtr_)
    {
        FatalErrorIn("void movingOversetRegion::calcMotionMask() const")
            << "Motion mask for movingOversetRegion " << name()
            << " already calculated"
            << abort(FatalError);
    }

    motionMaskPtr_ = new scalarField(mesh().allPoints().size(), 0);
    scalarField& mask = *motionMaskPtr_;

    const labelListList& cp = mesh().cellPoints();

    forAll (movingZoneNames_, mzI)
    {
        // Find moving cell zone
        const label zoneID =
            mesh().cellZones().findZoneID(movingZoneNames_[mzI]);

        if (zoneID < 0)
        {
            // Zone not found.  Consider if this is valid in parallel?
            FatalErrorIn("void movingOversetRegion::calcMotionMask() const")
                << "Cannot find moving zone " << movingZoneNames_[mzI]
                << " for movingOversetRegion " << name()
                << abort(FatalError);
        }

        // Get zone cells and mark vertices
        const cellZone& mz = mesh().cellZones()[zoneID];

        forAll (mz, cellI)
        {
            const labelList& curCp = cp[mz[cellI]];

            forAll (curCp, pointI)
            {
                mask[curCp[pointI]] = 1;
            }
        }
    }

    // if (debug)
    {
        // Count moving points
        label nMovingPoints = 0;

        forAll (mask, i)
        {
            if (mask[i] > 0)
            {
                nMovingPoints++;
            }
        }

        Info<< "movingOversetRegion " << name() << ": "
            << returnReduce(nMovingPoints, sumOp<scalar>()) << " moving points"
            << endl;
    }
}


void Foam::movingOversetRegion::clearOut()
{
    deleteDemandDrivenData(motionMaskPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingOversetRegion::movingOversetRegion
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    sbmfPtr_(solidBodyMotionFunction::New(dict, mesh_.time())),
    movingZoneNames_(dict.lookup("movingZones")),
    motionMaskPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingOversetRegion::~movingOversetRegion()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::scalarField& Foam::movingOversetRegion::motionMask() const
{
    if (!motionMaskPtr_)
    {
        calcMotionMask();
    }

    return *motionMaskPtr_;
}


const Foam::tmp<Foam::pointField>
Foam::movingOversetRegion::motionIncrement
(
    const pointField& undisplacedPoints
) const
{
    return motionMask()*
        (
            transform(sbmfPtr_->transformation(), undisplacedPoints)
          - undisplacedPoints
        );
}


// ************************************************************************* //
