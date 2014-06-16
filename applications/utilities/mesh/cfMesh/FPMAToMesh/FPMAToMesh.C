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
    Writes the mesh in fpma format readable by AVL's CfdWM

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMeshGenModifier.H"
#include "IFstream.H"

#include "Map.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{ 

#   include "setRootCase.H"
#   include "createTime.H"

    fileName inFileName;
    
    Info << "Reading mesh from file " << endl;
    cin >> inFileName;
    
    IFstream file(inFileName);

    polyMeshGen pmg(runTime);
    polyMeshGenModifier meshModifier(pmg);
    
    label counter;
    
    //- read the number of vertices
    pointFieldPMG& points = meshModifier.pointsAccess();
    file >> counter;
    
    //- read points from file
    points.setSize(counter);
    forAll(points, pointI)
    {
    point p;
    file >> p.x();
    file >> p.y();
    file >> p.z();
    
    points[pointI] = p;
    }
    
    //- read the number of faces
    file >> counter;
    
    faceListPMG& faces = meshModifier.facesAccess();
    
    //- read faces from file
    faces.setSize(counter);
    forAll(faces, faceI)
    {
        file >> counter;
    
    face f;
    f.setSize(counter);
    
    forAll(f, pI)
      file >> f[pI];
    
    faces[faceI] = f.reverseFace();
    }
    
    //- read the number of cells
    file >> counter;
    
    //- read cells from file
    cellListPMG& cells = meshModifier.cellsAccess();
    cells.setSize(counter);
    
    forAll(cells, cellI)
    {
    file >> counter;
    
    cell& c = cells[cellI];
    
    c.setSize(counter);
    
    forAll(c, fI)
        file >> c[fI];
    }
    
    //- read selections
    file >> counter;
    
    wordList patchNames;
    Map<label> subsetToPatch;
    
    for(label setI=0;setI<counter;++setI)
    {
    word sName;
    file >> sName;
    
    label type;
    file >> type;
    
    label nEntries;
    file >> nEntries;
    
    switch( type )
    {
        case 3:
        {
        //- face selection
        const label id = pmg.addFaceSubset(sName);
        
        patchNames.setSize(patchNames.size()+1);
        patchNames[patchNames.size()-1] = sName;
        subsetToPatch.insert(id, patchNames.size()-1);
        
        Info << "Reading face selection " << sName << endl;
        
        for(label i=0;i<nEntries;++i)
        {
            label entryI;
            file >> entryI;
            pmg.addFaceToSubset(id, entryI);
        }
        } break;
        case 2:
        {
        //- cell selection
        const label id = pmg.addCellSubset(sName);
        
        for(label i=0;i<nEntries;++i)
        {
            label entryI;
            file >> entryI;
            pmg.addCellToSubset(id, entryI);
        }
        } break;
        case 1:
        {
        //- node selection
        const label id = pmg.addPointSubset(sName);
        
        for(label i=0;i<nEntries;++i)
        {
            label entryI;
            file >> entryI;
            pmg.addPointToSubset(id, entryI);
        }
        } break;
    };
    }
    
    //- create patches from face selections
    VRWGraph boundaryFaces;
    labelLongList boundaryOwner;
    labelLongList boundaryPatches;
    
    const labelList& owner = pmg.owner();
    DynList<label> faceSubsets;
    pmg.faceSubsetIndices(faceSubsets);
    
    forAll(faceSubsets, setI)
    {
    labelLongList setFaces;
    pmg.facesInSubset(faceSubsets[setI], setFaces);
    
    forAll(setFaces, i)
    {
        boundaryFaces.appendList(faces[setFaces[i]]);
        boundaryOwner.append(owner[setFaces[i]]);
        boundaryPatches.append(subsetToPatch[faceSubsets[setI]]);
    }
    
    pmg.removeFaceSubset(faceSubsets[setI]);
    }
    
    meshModifier.reorderBoundaryFaces();
    meshModifier.replaceBoundary
    (
    patchNames,
    boundaryFaces,
    boundaryOwner,
    boundaryPatches
    );
    
    pmg.write();
    
    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
