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

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "X3DsurfaceFormatCore.H"
#include "clock.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::X3DsurfaceFormatCore::writeHeader
(
    Ostream& os
)
{
    os  <<
        "<?xml version='1.0' encoding='UTF-8'?>\n"
        "<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.0//EN\" \"http://www.web3d.org/specifications/x3d-3.0.dtd\">\n"
        "<X3D\n"
        "  version='3.0'\n"
        "  profile='Immersive'\n"
        "  xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance'\n"
        "  xsd:noNamespaceSchemaLocation='http://www.web3d.org/specifications/x3d-3.0.xsd'\n"
        "  >\n";
}


void Foam::fileFormats::X3DsurfaceFormatCore::writeAppearance
(
    Ostream& os
)
{
    os  <<
        "  <Appearance>\n"
        "   <Material"
        " diffuseColor='0.8 0.8 0.8'"
        " specularColor='1.0 1.0 1.0'"
        " shininess='0.5'"
        " transparency='0.0'"
        " />\n"           // end material
        "  </Appearance>\n";
}


// ************************************************************************* //
