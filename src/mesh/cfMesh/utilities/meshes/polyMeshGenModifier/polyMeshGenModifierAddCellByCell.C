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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenModifierAddCellByCell.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenModifierAddCellByCell::polyMeshGenModifierAddCellByCell
(
    polyMeshGen& mesh
)
:
    polyMeshGenModifier(mesh),
    nFaces_(mesh.faces().size()),
    newFaces_(nFaces_),
    nCells_(mesh.cells().size()),
    newCells_(nCells_)
{
    this->pointFaces();
    faceListPMG& faces = this->facesAccess();
    forAll(faces, faceI)
        newFaces_[faceI].transfer(faces[faceI]);
    
    cellListPMG& cells = this->cellsAccess();
    forAll(cells, cellI)
        newCells_[cellI].transfer(cells[cellI]);
};
            
// Destructor
polyMeshGenModifierAddCellByCell::~polyMeshGenModifierAddCellByCell()
{
    faceListPMG& faces = this->facesAccess();
    faces.setSize(nFaces_);
    forAll(faces, faceI)
        faces[faceI].transfer(newFaces_[faceI]);
    
    cellListPMG& cells = this->cellsAccess();
    cells.setSize(newCells_.size());
    forAll(cells, cellI)
        cells[cellI].transfer(newCells_[cellI]);
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifierAddCellByCell::addCell(const faceList& cellFaces)
{
    cell c(cellFaces.size());
    
    VRWGraph& pointFaces = this->pointFaces();
    
    forAll(cellFaces, faceI)
    {
        const face& f = cellFaces[faceI];
        
        const label pointI = f[0];
        
        label fLabel(-1);
        forAllRow(pointFaces, pointI, pfI)
        {
            const label faceI = pointFaces(pointI, pfI);
            
            if( newFaces_[faceI] == f )
            {
                fLabel = faceI;
                break;
            }
        }
            
        if( fLabel == -1 )
        {
            newFaces_.append(f);
            c[faceI] = nFaces_;
            forAll(f, pI)
                pointFaces.append(f[pI], nFaces_);
            
            ++nFaces_;
        }
        else
        {
            c[faceI] = fLabel;
        }
    }
    
    newCells_.append(c);
    ++nCells_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
