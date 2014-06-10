/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description
    Declaration of IOLongList ClassNames for IOLists that do not have .C files.

\*---------------------------------------------------------------------------*/

#include "IOLongListInstances.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(IOLongList<label>, labelIOListPMG);
    defineCompoundTypeName(IOLongList<point>, pointIOFieldPMG);
    //defineCompoundTypeName(IOLongList<face>, faceIOListPMG);
    //defineCompoundTypeName(IOLongList<cell>, cellIOListPMG);
    //addCompoundToRunTimeSelectionTable(IOLongList<label>, labelIOLongList);

    defineTemplateTypeNameAndDebugWithName(labelIOListPMG, "labelList", 0);
    defineTemplateTypeNameAndDebugWithName(pointIOFieldPMG, "vectorField", 0);
    //defineTemplateTypeNameAndDebugWithName(faceIOListPMG, "faceList", 0);
    //defineTemplateTypeNameAndDebugWithName(cellIOListPMG, "cellList", 0);
}

// ************************************************************************* //
