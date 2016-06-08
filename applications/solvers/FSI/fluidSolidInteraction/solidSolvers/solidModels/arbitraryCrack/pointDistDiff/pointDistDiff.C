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

Description
    Linear distance-based motion diffusion.

\*---------------------------------------------------------------------------*/

#include "pointDistDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "elementFields.H"
#include "patchWave.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointDistDiff, 0);
    addToRunTimeSelectionTable(motionDiff, pointDistDiff, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pointDistDiff::pointDistDiff(const tetMotionSolver& mSolver)
:
    motionDiff(mSolver),
    patchNames_(mSolver.lookup("distancePatches")),
    motionGamma_
    (
        IOobject
        (
            "motionGamma",
            tetMesh().time().timeName(),
            tetMesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tetMesh(),
        dimensionedScalar("1.0", dimless, 1.0)
    )
{
    motionGamma_.internalField() = pow(1.0/L(), 3.0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointDistDiff::~pointDistDiff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::pointDistDiff::L() const
{
    const polyMesh& m = tetMesh()();
    const vectorField& p = tetMesh().points();
    const vectorField& C = m.cellCentres();


    const DynamicList<label>& fp = mSolver().fixedPoints();

    tmp<scalarField> tDist(new scalarField(C.size(), GREAT));
    scalarField& dist = tDist();

    if (fp.size())
    {
        forAll(fp, pointI)
        {
            vector origin = p[fp[pointI]];
            if (m.nGeometricD() != 3)
            {
                const Vector<label>& directions = m.geometricD();
                for (direction dir = 0; dir < directions.nComponents; dir++)
                {
                    if (directions[dir] == -1)
                    {
                        origin.component(dir) = 0;
//                         break;
                    }
                }
            }

            forAll(C, cellI)
            {
                scalar curDist = mag(C[cellI] - origin);

                if (curDist < dist[cellI])
                {
                    dist[cellI] = curDist;
                }
            }
        }

//         return tDist;

//         vector origin = vector::zero;
//         forAll(fp, pointI)
//         {
//             if (fp[pointI] < m.points().size())
//             {
//                 origin = m.points()[fp[pointI]];
//                 break;
//             }
//         }
//         if (m.nGeometricD() == 2)
//         {
//             const Vector<label>& directions = m.geometricD();
//             for (direction dir = 0; dir < directions.nComponents; dir++)
//             {
//                 if (directions[dir] == -1)
//                 {
//                     origin.component(dir) = 0;
//                     break;
//                 }
//             }
//         }
//         Info << "origin: " << origin << endl;

//         return mag(m.cellCentres()-origin);
    }
    else
    {
        dist = 1.0;
//         return tmp<scalarField>(new scalarField(m.nCells(), 1.0));
    }

    return tDist;
}


void Foam::pointDistDiff::correct()
{
    motionGamma().internalField() = pow(1.0/L(), 3.0);
}


// ************************************************************************* //
