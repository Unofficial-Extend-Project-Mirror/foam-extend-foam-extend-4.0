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

\*----------------------------------------------------------------------------*/

#include "moleculeCloud.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineParticleTypeNameAndDebug(molecule, 0);
    defineTemplateTypeNameAndDebug(Cloud<molecule>, 0);
};


template<>
const char* Foam::NamedEnum<Foam::moleculeCloud::
integrationMethods, 2>::names[] =
{
    "verletLeapfrog",
    "predictorCorrector"
};

const Foam::NamedEnum<Foam::moleculeCloud::integrationMethods, 2>
        Foam::moleculeCloud::integrationMethodNames_;


Foam::scalar Foam::moleculeCloud::transTol = 1e-12;

Foam::scalar Foam::moleculeCloud::kb = 1.380650277e-23;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moleculeCloud::moleculeCloud(const polyMesh& mesh)
:
    Cloud<molecule>(mesh, "moleculeCloud", false),
    mesh_(mesh),
    referredInteractionList_(*this)
{
    molecule::readFields(*this);

#   include "moleculeCloudReadMDParameters.H"

    buildCellInteractionLists();

    buildCellReferralLists();

    buildCellOccupancy();
}


Foam::moleculeCloud::moleculeCloud
(
    const polyMesh& mesh,
    label nMol,
    const labelField& id,
    const scalarField& mass,
    const vectorField& positions,
    const labelField& cells,
    const vectorField& U,
    const vectorField& A,
    const labelField& tethered,
    const vectorField& tetherPositions
)
:
    Cloud<molecule>(mesh, "moleculeCloud", false),
    mesh_(mesh),
    referredInteractionList_(*this)
{
    molecule::readFields(*this);

    clear();

    // This clear ()is here for the moment to stop existing files
    // being appended to, this would be better accomplished by getting
    // mesh.removeFiles(mesh.instance()); (or equivalent) to work.

    int i;

    const Cloud<molecule>& cloud = *this;

    for (i=0; i<nMol; i++)
    {
        addParticle
        (
            new molecule
            (
                cloud,
                positions[i],
                cells[i],
                mass[i],
                U[i],
                A[i],
                tetherPositions[i],
                tethered[i],
                id[i]
            )
        );
    }
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::moleculeCloud::writeFields() const
{
    molecule::writeFields(*this);
}


// ************************************************************************* //
