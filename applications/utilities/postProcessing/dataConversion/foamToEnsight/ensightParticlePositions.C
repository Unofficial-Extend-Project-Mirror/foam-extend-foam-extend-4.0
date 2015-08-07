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

#include "ensightParticlePositions.H"
#include "fvMesh.H"
#include "passiveParticle.H"
#include "CloudTemplate.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "itoa.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void ensightParticlePositions
(
    const Foam::fvMesh& mesh,
    const Foam::fileName& postProcPath,
    const Foam::word& timeFile,
    const Foam::word& cloudName,
    const bool dataExists
)
{
    if (dataExists)
    {
        Info<< "Converting cloud " << cloudName << " positions" <<  endl;
    }
    else
    {
        Info<< "Creating empty cloud " << cloudName << " positions" << endl;
    }

    const Time& runTime = mesh.time();

    fileName ensightFileName(timeFile + "." + cloudName);
    OFstream ensightFile
    (
        postProcPath/ensightFileName,
        ios_base::out|ios_base::trunc,
        runTime.writeFormat(),
        runTime.writeVersion(),
        runTime.writeCompression()
    );

    // Output header
    ensightFile
        << cloudName.c_str() << nl
        << "particle coordinates" << nl;

    if (dataExists)
    {
        Cloud<passiveParticle> parcels(mesh, cloudName, false);

        // Set Format
        ensightFile.setf(ios_base::scientific, ios_base::floatfield);
        ensightFile.precision(5);

        ensightFile<< setw(8) << parcels.size() << nl;

        label nParcels = 0;

        // Output positions
        forAllConstIter(Cloud<passiveParticle>, parcels, elmnt)
        {
            const vector& p = elmnt().position();

            ensightFile
                << setw(8) << ++nParcels
                << setw(12) << p.x() << setw(12) << p.y() << setw(12) << p.z()
                << nl;
        }
    }
    else
    {
        label nParcels = 0;
        ensightFile<< setw(8) << nParcels << nl;
    }
}


// ************************************************************************* //
