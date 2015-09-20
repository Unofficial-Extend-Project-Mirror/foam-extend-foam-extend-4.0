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

#include "lagrangianWriter.H"
#include "writeFuns.H"
#include "CloudTemplate.H"
#include "passiveParticle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::lagrangianWriter::lagrangianWriter
(
    const vtkMesh& vMesh,
    const bool binary,
    const fileName& fName,
    const word& cloudName,
    const bool dummyCloud
)
:
    vMesh_(vMesh),
    binary_(binary),
    fName_(fName),
    cloudName_(cloudName),
    os_(fName.c_str())
{
    const fvMesh& mesh = vMesh_.mesh();

    // Write header
    writeFuns::writeHeader(os_, binary_, mesh.time().caseName());
    os_ << "DATASET POLYDATA" << std::endl;

    if (dummyCloud)
    {
        nParcels_ = 0;

        os_ << "POINTS " << nParcels_ << " float" << std::endl;
    }
    else
    {
        Cloud<passiveParticle> parcels(mesh, cloudName_, false);

        nParcels_ = parcels.size();

        os_ << "POINTS " << nParcels_ << " float" << std::endl;

        DynamicList<floatScalar> partField(3*parcels.size());

        forAllConstIter(Cloud<passiveParticle>, parcels, elmnt)
        {
            writeFuns::insert(elmnt().position(), partField);
        }
        writeFuns::write(os_, binary_, partField);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lagrangianWriter::writeParcelHeader(const label nFields)
{
    os_ << "POINT_DATA " << nParcels_ << std::endl
        << "FIELD attributes " << nFields
        << std::endl;
}


// ************************************************************************* //
