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

#include "maxEdgeLengthDelta.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(maxEdgeLengthDelta, 0);
addToRunTimeSelectionTable(LESdelta, maxEdgeLengthDelta, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void maxEdgeLengthDelta::calcDelta()
{
    label nD = mesh().nGeometricD();

    const cellList& cells = mesh().cells();
    const faceList& faces = mesh().faces();
    const pointField& points = mesh().points();

    scalarField& d = delta_.internalField();

    // Initialise d to a large number
    d = 0;


    if (nD == 3)
    {
        forAll (cells, cellI)
        {
            edgeList ce = cells[cellI].edges(faces);

            forAll (ce, ceI)
            {
                d[cellI] = Foam::max(d[cellI], ce[ceI].mag(points));
            }
        }
    }
    else if (nD == 2)
    {
        WarningIn("maxEdgeLengthDelta::calcDelta()")
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;

        // Find the dead direction
        const Vector<label>& directions = mesh().geometricD();

        vector deadDir = vector::zero;

        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                deadDir[dir] = 1;
                break;
            }
        }

        forAll (cells, cellI)
        {
            edgeList ce = cells[cellI].edges(faces);

            forAll (ce, ceI)
            {
                // Calculate edge magnitude and unit vector
                scalar magE = ce[ceI].mag(points);
                vector eN = ce[ceI].vec(points)/magE;

                // Filter the edges in the dead direction: the dot-product
                // of eN and deadDir should be 1, but in 2-D checking
                // can be relaxed.  HJ, 8/Dec/2010
                if (mag(eN & deadDir) < 0.9)
                {
                    d[cellI] = Foam::max(d[cellI], magE);
                }
            }
        }

    }
    else
    {
        FatalErrorIn("maxEdgeLengthDelta::calcDelta()")
            << "Case is not 3D or 2D, LES is not applicable"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

maxEdgeLengthDelta::maxEdgeLengthDelta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dd
)
:
    LESdelta(name, mesh),
    deltaCoeff_(readScalar(dd.subDict(type() + "Coeffs").lookup("deltaCoeff")))
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void maxEdgeLengthDelta::read(const dictionary& dd)
{
    dd.subDict(type() + "Coeffs").lookup("deltaCoeff") >> deltaCoeff_;
    calcDelta();
}


void maxEdgeLengthDelta::correct()
{
    if (mesh_.changing())
    {
        calcDelta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
