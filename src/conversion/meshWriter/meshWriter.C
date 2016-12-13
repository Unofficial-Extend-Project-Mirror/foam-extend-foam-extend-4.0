/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

\*----------------------------------------------------------------------------*/

#include "meshWriter.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::cellModel* Foam::meshWriter::unknownModel = Foam::cellModeller::
lookup
(
    "unknown"
);


const Foam::cellModel* Foam::meshWriter::tetModel = Foam::cellModeller::
lookup
(
    "tet"
);


const Foam::cellModel* Foam::meshWriter::pyrModel = Foam::cellModeller::
lookup
(
    "pyr"
);


const Foam::cellModel* Foam::meshWriter::prismModel = Foam::cellModeller::
lookup
(
    "prism"
);


const Foam::cellModel* Foam::meshWriter::hexModel = Foam::cellModeller::
lookup
(
    "hex"
);


Foam::string Foam::meshWriter::defaultMeshName = "meshExport";
Foam::string Foam::meshWriter::defaultSurfaceName = "surfExport";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshWriter::meshWriter(const polyMesh& mesh, const scalar scaleFactor)
:
    mesh_(mesh),
    scaleFactor_(scaleFactor),
    writeBoundary_(true),
    boundaryRegion_(),
    cellTable_(),
    cellTableId_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshWriter::~meshWriter()
{}


// ************************************************************************* //
