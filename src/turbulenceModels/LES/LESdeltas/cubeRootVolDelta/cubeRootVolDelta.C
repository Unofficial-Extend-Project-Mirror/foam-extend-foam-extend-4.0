/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "cubeRootVolDelta.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cubeRootVolDelta, 0);
addToRunTimeSelectionTable(LESdelta, cubeRootVolDelta, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void cubeRootVolDelta::calcDelta()
{
    const Vector<label>& directions = mesh().directions();
    label nD = (directions.nComponents + cmptSum(directions))/2;

    if (nD == 3)
    {
        delta_.internalField() = deltaCoeff_*pow(mesh().V(), 1.0/3.0);
    }
    else if (nD == 2)
    {
        WarningIn("cubeRootVolDelta::calcDelta()")
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;

        scalar thickness = 0.0;
        for (direction dir=0; dir<directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                boundBox bb(mesh().points(), false);

                thickness = (bb.max() - bb.min())[dir];
            }
        }

        delta_.internalField() = deltaCoeff_*sqrt(mesh().V()/thickness);
    }
    else
    {
        delta_.internalField() = deltaCoeff_*pow(mesh().V(), 1.0/3.0);

        WarningIn("cubeRootVolDelta::calcDelta()")
            << "Case is not 3D or 2D, LES is not applicable"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cubeRootVolDelta::cubeRootVolDelta
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

void cubeRootVolDelta::read(const dictionary& dd)
{
    dd.subDict(type() + "Coeffs").lookup("deltaCoeff") >> deltaCoeff_;
    calcDelta();
}


void cubeRootVolDelta::correct()
{
    if (mesh_.changing())
    {
        calcDelta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
