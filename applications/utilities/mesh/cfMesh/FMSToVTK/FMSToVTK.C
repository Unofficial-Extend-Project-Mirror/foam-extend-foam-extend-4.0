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
    cfMesh utility to convert a surface file to VTK multiblock dataset
    format, including the patches, feature edges and surface features.

Author
    Ivor Clifford <ivor.clifford@psi.ch>

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "autoPtr.H"
#include "OFstream.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "xmlTag.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write the supplied pointList in serial vtkPolyData format
void writePointsToVTK
(
    const fileName& fn,
    const string& title,
    const UList<point>& points
)
{
    xmlTag xmlRoot("VTKFile");
    xmlRoot.addAttribute("type", "PolyData");
    
    xmlTag& xmlPolyData = xmlRoot.addChild("PolyData");
    
    xmlTag& xmlPiece = xmlPolyData.addChild("Piece");
    xmlPiece.addAttribute("NumberOfPoints", points.size());
    
    xmlTag& xmlPoints = xmlPiece.addChild("Points");
    
    xmlTag& xmlPointData = xmlPoints.addChild("DataArray");
    xmlPointData.addAttribute("type", "Float32");
    xmlPointData.addAttribute("NumberOfComponents", 3);
    xmlPointData.addAttribute("format", "ascii");
    xmlPointData << points;
    
    OFstream os(fn);
    os << xmlRoot << endl;
    
    Info << "Created " << fn << endl;
}


//- Write the supplied addressed pointList in serial vtkPolyData format
void writePointsToVTK
(
    const fileName& fn,
    const string& title,
    const UList<point>& points,
    unallocLabelList& addr
)
{
    // Create subaddressed points
    pointField newPoints(addr.size());
    
    forAll(addr, i)
    {
        newPoints[i] = points[addr[i]];
    }
    
    writePointsToVTK
    (
        fn,
        title,
        newPoints
    );
}


//- Write the supplied edgeList in serial vtkPolyData format
void writeEdgesToVTK
(
    const fileName& fn,
    const string& title,
    const UList<point>& points,
    const LongList<edge>& edges
)
{
    labelList connectivity(edges.size());
    
    forAll(edges, edgeI)
    {
        connectivity[edgeI] = 2*(edgeI+1);
    }
    
    xmlTag xmlRoot("VTKFile");
    xmlRoot.addAttribute("type", "PolyData");
    
    xmlTag& xmlPolyData = xmlRoot.addChild("PolyData");
    
    xmlTag& xmlPiece = xmlPolyData.addChild("Piece");
    xmlPiece.addAttribute("NumberOfPoints", points.size());
    xmlPiece.addAttribute("NumberOfLines", edges.size());
    
    xmlTag& xmlPoints = xmlPiece.addChild("Points");
    
    xmlTag& xmlPointData = xmlPoints.addChild("DataArray");
    xmlPointData.addAttribute("type", "Float32");
    xmlPointData.addAttribute("NumberOfComponents", 3);
    xmlPointData.addAttribute("format", "ascii");
    xmlPointData << points;

    xmlTag& xmlEdges = xmlPiece.addChild("Lines");

    xmlTag& xmlEdgeData = xmlEdges.addChild("DataArray");
    xmlEdgeData.addAttribute("type", "Int32");
    xmlEdgeData.addAttribute("Name", "connectivity");
    xmlEdgeData.addAttribute("format", "ascii");
    xmlEdgeData << edges;

    xmlTag& xmlConnectData = xmlEdges.addChild("DataArray");
    xmlConnectData.addAttribute("type", "Int32");
    xmlConnectData.addAttribute("Name", "offsets");
    xmlConnectData.addAttribute("format", "ascii");
    xmlConnectData << connectivity;
    
    OFstream os(fn);
    os << xmlRoot << endl;

    Info << "Created " << fn << endl;
}

//- Write the supplied edgeList subset in serial vtkPolyData format
void writeEdgesToVTK
(
    const fileName& fn,
    const string& title,
    const UList<point>& points,
    const LongList<edge>& edges,
    const unallocLabelList& addr
)
{
    // Remove unused points and create subaddressed edges
    DynamicList<point> newPoints;
    labelList newPointAddr(points.size(), -1);
    LongList<edge> newEdges(addr.size());
    
    forAll(addr, addrI)
    {
        label edgeI = addr[addrI];
        
        const edge& curEdge = edges[edgeI];
        edge& newEdge = newEdges[addrI];
        
        forAll(curEdge, i)
        {
            label pointId = curEdge[i];
            
            if (newPointAddr[pointId] == -1)
            {
                newPoints.append(points[pointId]);
                newPointAddr[pointId] = newPoints.size()-1;
            }
            
            newEdge[i] = newPointAddr[pointId];
        }
    }
    
    writeEdgesToVTK
    (
        fn,
        title,
        newPoints,
        newEdges
    );
}


//- Write the supplied facet list in serial vtkPolyData format
void writeFacetsToVTK
(
    const fileName& fn,
    const string& title,
    const UList<point>& points,
    const LongList<labelledTri>& facets
)
{
    labelList connectivity(facets.size());
    
    forAll(facets, faceI)
    {
        connectivity[faceI] = 3*(faceI+1);
    }
    
    labelList regionData(facets.size());
    
    forAll(facets, faceI)
    {
        regionData[faceI] =  facets[faceI].region();
    }
    
    xmlTag xmlRoot("VTKFile");
    xmlRoot.addAttribute("type", "PolyData");
    
    xmlTag& xmlPolyData = xmlRoot.addChild("PolyData");
    
    xmlTag& xmlPiece = xmlPolyData.addChild("Piece");
    xmlPiece.addAttribute("NumberOfPoints", points.size());
    xmlPiece.addAttribute("NumberOfPolys", facets.size());
    
    xmlTag& xmlPoints = xmlPiece.addChild("Points");
    
    xmlTag& xmlPointData = xmlPoints.addChild("DataArray");
    xmlPointData.addAttribute("type", "Float32");
    xmlPointData.addAttribute("NumberOfComponents", 3);
    xmlPointData.addAttribute("format", "ascii");
    xmlPointData << points;

    xmlTag& xmlPolys = xmlPiece.addChild("Polys");

    xmlTag& xmlPolyDataArray = xmlPolys.addChild("DataArray");
    xmlPolyDataArray.addAttribute("type", "Int32");
    xmlPolyDataArray.addAttribute("Name", "connectivity");
    xmlPolyDataArray.addAttribute("format", "ascii");
    xmlPolyDataArray << facets;

    xmlTag& xmlConnectData = xmlPolys.addChild("DataArray");
    xmlConnectData.addAttribute("type", "Int32");
    xmlConnectData.addAttribute("Name", "offsets");
    xmlConnectData.addAttribute("format", "ascii");
    xmlConnectData << connectivity;
    
    xmlTag& xmlCellData = xmlPiece.addChild("CellData");
    
    xmlTag& xmlCellDataArray = xmlCellData.addChild("DataArray");
    xmlCellDataArray.addAttribute("type", "Int32");
    xmlCellDataArray.addAttribute("Name", "region");
    xmlCellDataArray.addAttribute("format", "ascii");
    xmlCellDataArray << regionData;
    
    OFstream os(fn);
    os << xmlRoot << endl;
    
    Info << "Created " << fn << endl;
}


//- Write an addressed subset of the supplied facet list
//-  in serial vtkPolyData format
void writeFacetsToVTK
(
    const fileName& fn,
    const string& title,
    const pointField& points,
    const LongList<labelledTri>& facets,
    const unallocLabelList& addr
)
{        
    // Remove unused points and create subaddressed facets
    DynamicList<point> newPoints;
    labelList newPointAddr(points.size(), -1);
    LongList<labelledTri> newFacets(addr.size());
    
    forAll(addr, addrI)
    {
        label faceI = addr[addrI];
        
        const labelledTri& facet = facets[faceI];
        const FixedList<label, 3>& pointIds = facet;
        FixedList<label, 3> newPointIds;
        
        forAll(pointIds, i)
        {
            label pointId = pointIds[i];
            
            if (newPointAddr[pointId] == -1)
            {
                newPoints.append(points[pointId]);
                newPointAddr[pointId] = newPoints.size()-1;
            }
            
            newPointIds[i] = newPointAddr[pointId];
        }
        
        newFacets[addrI] = labelledTri
        (
            newPointIds[0],
            newPointIds[1],
            newPointIds[2],
            facet.region()
        );
    }
    
    writeFacetsToVTK
    (
        fn,
        title,
        newPoints,
        newFacets
    );
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("output prefix");
    argList args(argc, argv);

    // Process commandline arguments
    fileName inFileName(args.args()[1]);
    fileName outPrefix(args.args()[2]);
    
    // Read original surface
    triSurf origSurf(inFileName);
    
    const pointField& points = origSurf.points();
    const LongList<labelledTri>& facets = origSurf.facets();
    const LongList<edge>& edges = origSurf.featureEdges();
    const geometricSurfacePatchList& patches = origSurf.patches();
    
    label index = 0;
    
    // Create file structure for multiblock dataset
    mkDir(outPrefix);
    mkDir(outPrefix + "/patches");
    mkDir(outPrefix + "/pointSubsets");
    mkDir(outPrefix + "/edgeSubsets");
    mkDir(outPrefix + "/faceSubsets");
    
    // Create VTK multiblock dataset file
    xmlTag xmlRoot("VTKFile");
    xmlRoot.addAttribute("type", "vtkMultiBlockDataSet");
    xmlRoot.addAttribute("version", "1.0");
    xmlRoot.addAttribute("byte_order", "LittleEndian");
    
    xmlTag& xmlDataSet = xmlRoot.addChild("vtkMultiBlockDataSet");
        
    // Write faces and feature edges
    {
        fileName fn = outPrefix / "facets.vtp";
        
        writeFacetsToVTK
        (
            outPrefix / "facets.vtp",
            outPrefix,
            points,
            facets
        );
        
        xmlTag& tag = xmlDataSet.addChild("DataSet");
        tag.addAttribute("index", Foam::name(index++));
        tag.addAttribute("name", "facets");
        tag.addAttribute("file", fn);
    }
    
    {
        fileName fn = outPrefix / "featureEdges.vtp";
        
        writeEdgesToVTK
        (
            outPrefix / "featureEdges.vtp",
            "featureEdges",
            points,
            edges
        );
    
        xmlTag& tag = xmlDataSet.addChild("DataSet");
        tag.addAttribute("index", Foam::name(index++));
        tag.addAttribute("name", "featureEdges");
        tag.addAttribute("file", fn);
    }
    
    // Write patches
    // Create patch addressing
    List<DynamicList<label> > patchAddr(patches.size());
    
    forAll(facets, faceI)
    {
        patchAddr[facets[faceI].region()].append(faceI);
    }
    
    {
        xmlTag& xmlBlock = xmlDataSet.addChild("Block");
        xmlBlock.addAttribute("index", Foam::name(index++));
        xmlBlock.addAttribute("name", "patches");
        
        forAll(patches, patchI)
        {            
            word patchName = patches[patchI].name();
            
            fileName fn = outPrefix / "patches" / patchName + ".vtp";
            
            writeFacetsToVTK
            (
                fn,
                patchName,
                points,
                facets,
                patchAddr[patchI]
            );
            
            xmlTag& tag = xmlBlock.addChild("DataSet");
            tag.addAttribute("index", Foam::name(patchI));
            tag.addAttribute("name", patchName);
            tag.addAttribute("file", fn);
        }
    }
    
    // Write point subsets
    {
        xmlTag& xmlBlock = xmlDataSet.addChild("Block");
        xmlBlock.addAttribute("index", Foam::name(index++));
        xmlBlock.addAttribute("name", "pointSubsets");
        
        DynList<label> subsetIndices;
        labelList subsetAddr;
        
        origSurf.pointSubsetIndices(subsetIndices);
        
        forAll(subsetIndices, id)
        {
            word subsetName = origSurf.pointSubsetName(id);
            origSurf.pointsInSubset(id, subsetAddr);
            
            fileName fn = outPrefix / "pointSubsets" / subsetName + ".vtp";
            
            writePointsToVTK
            (
                fn,
                subsetName,
                points,
                subsetAddr
            );
            
            xmlTag& tag = xmlBlock.addChild("DataSet");
            tag.addAttribute("index", Foam::name(id));
            tag.addAttribute("name", subsetName);
            tag.addAttribute("file", fn);
        }
    }
    
    // Write edge subsets
    {
        xmlTag& xmlBlock = xmlDataSet.addChild("Block");
        xmlBlock.addAttribute("index", Foam::name(index++));
        xmlBlock.addAttribute("name", "edgeSubsets");
        
        DynList<label> subsetIndices;
        labelList subsetAddr;
        
        origSurf.edgeSubsetIndices(subsetIndices);
        
        forAll(subsetIndices, id)
        {
            word subsetName = origSurf.edgeSubsetName(id);
            origSurf.edgesInSubset(id, subsetAddr);
            
            fileName fn = outPrefix / "edgeSubsets" / subsetName + ".vtp";
            
            writeEdgesToVTK
            (
                fn,
                subsetName,
                points,
                edges,
                subsetAddr
            );
            
            xmlTag& tag = xmlBlock.addChild("DataSet");
            tag.addAttribute("index", Foam::name(id));
            tag.addAttribute("name", subsetName);
            tag.addAttribute("file", fn);
        }
    }
    
    // Write facet subsets
    {
        xmlTag& xmlBlock = xmlDataSet.addChild("Block");
        xmlBlock.addAttribute("index", Foam::name(index++));
        xmlBlock.addAttribute("name", "faceSubsets");
        
        DynList<label> subsetIndices;
        labelList subsetAddr;
                
        origSurf.facetSubsetIndices(subsetIndices);
        
        forAll(subsetIndices, id)
        {            
            word subsetName = origSurf.facetSubsetName(id);
            origSurf.facetsInSubset(id, subsetAddr);
            
            fileName fn = outPrefix / "faceSubsets"  / subsetName + ".vtp";
            
            writeFacetsToVTK
            (
                fn,
                subsetName,
                points,
                facets,
                subsetAddr
            );
            
            xmlTag& tag = xmlBlock.addChild("DataSet");
            tag.addAttribute("index", Foam::name(id));
            tag.addAttribute("name", subsetName);
            tag.addAttribute("file", fn);
        }
    }    
        
    OFstream os(outPrefix + ".vtm");
    os << xmlRoot << endl;
    
    Info << "Created " << outPrefix + ".vtm" << endl;
    
    Info << "End\n" << endl;
    
    return 0;
}

// ************************************************************************* //
