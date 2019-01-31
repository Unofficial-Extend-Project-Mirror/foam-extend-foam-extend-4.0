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

\*---------------------------------------------------------------------------*/

#include "fpmaMesh.H"
#include "IOmanip.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyMeshGen
Foam::fpmaMesh::fpmaMesh(const polyMeshGen& mesh)
:
    mesh_(mesh)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fpmaMesh::~fpmaMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fpmaMesh::writePoints(Foam::OFstream& fpmaGeometryFile) const
{
    fpmaGeometryFile << mesh_.points().size() << nl;
    const pointFieldPMG& points = mesh_.points();
    forAll(points, pointI)
    {
        const point& p = points[pointI];
        fpmaGeometryFile << p.x() << ' ' << p.y() << ' ' << p.z() << ' ';
    }
    
    fpmaGeometryFile << nl;
}

void fpmaMesh::writeCells(OFstream& fpmaGeometryFile) const
{
    const cellListPMG& cells = mesh_.cells();
    
    fpmaGeometryFile << cells.size() << nl;
    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];
        
        fpmaGeometryFile << c.size();
        forAll(c, fI)
            fpmaGeometryFile << ' ' << c[fI];
        fpmaGeometryFile << nl;
    }
}

void Foam::fpmaMesh::writeFaces(OFstream& fpmaGeometryFile) const
{
    const faceListPMG& faces = mesh_.faces();
    fpmaGeometryFile << faces.size() << nl;
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];
        
        fpmaGeometryFile << f.size();
        forAllReverse(f, pI)
            fpmaGeometryFile << ' ' << f[pI];
        fpmaGeometryFile << nl;
    }
}

void Foam::fpmaMesh::writeSubsets(Foam::OFstream& fpmaGeometryFile) const
{
    //- write patches as face selections
    const PtrList<boundaryPatch>& patches = mesh_.boundaries();
    
    label nSubsets(0);
    
    nSubsets += patches.size();
    DynList<label> indices;
    mesh_.pointSubsetIndices(indices);
    nSubsets += indices.size();
    Info << "Mesh has " << indices.size() << " point subsets" << endl;
    mesh_.faceSubsetIndices(indices);
    nSubsets += indices.size();
    Info << "Mesh has " << indices.size() << " face subsets" << endl;
    mesh_.cellSubsetIndices(indices);
    nSubsets += indices.size();
    Info << "Mesh has " << indices.size() << " cell subsets" << endl;
    
    fpmaGeometryFile << nSubsets << nl;
    
    //- write patches as face selections
    forAll(patches, patchI)
    {
        label start = patches[patchI].patchStart();
        const label size = patches[patchI].patchSize();
        
        fpmaGeometryFile << patches[patchI].patchName() << nl;
        fpmaGeometryFile << 3 << nl;
        fpmaGeometryFile << size << nl;
        for(label i=0;i<size;++i)
            fpmaGeometryFile << start++ << ' ';
        fpmaGeometryFile << nl;
    }
    
    //- write node selections
    mesh_.pointSubsetIndices(indices);
    forAll(indices, indexI)
    {
        labelLongList nodesInSubset;
        mesh_.pointsInSubset(indices[indexI], nodesInSubset);
        
        fpmaGeometryFile << mesh_.pointSubsetName(indices[indexI]) << nl;
        fpmaGeometryFile << 1 << nl;
        fpmaGeometryFile << nodesInSubset.size() << nl;
        forAll(nodesInSubset, i)
            fpmaGeometryFile << nodesInSubset[i] << ' ';
        fpmaGeometryFile << nl;
    }
    
    //- write face selections
    mesh_.faceSubsetIndices(indices);
    forAll(indices, indexI)
    {
        labelLongList facesInSubset;
        mesh_.facesInSubset(indices[indexI], facesInSubset);
        
        fpmaGeometryFile << mesh_.faceSubsetName(indices[indexI]) << nl;
        fpmaGeometryFile << 3 << nl;
        fpmaGeometryFile << facesInSubset.size() << nl;
        forAll(facesInSubset, i)
            fpmaGeometryFile << facesInSubset[i] << ' ';
        fpmaGeometryFile << nl;
    }
    
    //- write cell selections
    mesh_.cellSubsetIndices(indices);
    forAll(indices, indexI)
    {
        labelLongList cellsInSubset;
        mesh_.cellsInSubset(indices[indexI], cellsInSubset);
        
        fpmaGeometryFile << mesh_.cellSubsetName(indices[indexI]) << nl;
        fpmaGeometryFile << 2 << nl;
        fpmaGeometryFile << cellsInSubset.size() << nl;
        forAll(cellsInSubset, i)
            fpmaGeometryFile << cellsInSubset[i] << ' ';
        fpmaGeometryFile << nl;
    }
}


void fpmaMesh::write(OFstream& fpmaGeometryFile) const
{
    writePoints(fpmaGeometryFile);
    
    writeFaces(fpmaGeometryFile);

    writeCells(fpmaGeometryFile);

    writeSubsets(fpmaGeometryFile);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
