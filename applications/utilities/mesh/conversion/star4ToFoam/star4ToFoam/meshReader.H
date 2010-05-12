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

Class
    meshReader

Description
    This class supports creating polyMeshes with baffles.
    The boundary definitions are given as cell/face.
    The derived classes are responsible for providing the protected data.
    This implementation is somewhat messy, but could/should be restructured
    to provide a more generalized reader (at the moment it has been written
    for converting pro-STAR data).

    The meshReader supports celltable information (see new user's guide entry).

SourceFiles
    calcPointCells.C
    createPolyBoundary.C
    createPolyCells.C
    meshReader.C
    meshReaderAux.C

\*---------------------------------------------------------------------------*/

#ifndef meshReader_H
#define meshReader_H

#include "fvMesh.H"
// #include "polyPatchList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
			Class meshReader Declaration
\*---------------------------------------------------------------------------*/

class meshReader
{
    // Private data

    //- Database
    const Time& runTime_;

    //- Number of internal faces for polyMesh
    label nInternalFaces_;

    //- Polyhedral mesh boundary patch start indices and dimensions
    labelList patchStarts_;
    labelList patchSizes_;

    //- association between faces
    List<labelPair> interfaces_;

    //- cells and their faces associated with each baffle
    labelListList baffleCellIds_;
    labelListList baffleFaceIds_;

    //- Global face list for polyMesh
    faceList meshFaces_;

    //- Cells as polyhedra for polyMesh
    cellList cellPolys_;

protected:
    // Protected data

    // Pointers to cell shape models
    static const cellModel* unknownModel;
    static const cellModel* tetModel;
    static const cellModel* pyrModel;
    static const cellModel* prismModel;
    static const cellModel* hexModel;

    //- lookup original Cell number for a given cell
    labelList origCellId_;

    //- raw boundary faces identified by cell and face number
    //  for each patch, parallel lists of cells/faces
    labelListList boundaryCells_;
    labelListList boundaryFaces_;

    //- Boundary patch types
    wordList patchTypes_;

    //- Boundary patch names
    wordList patchNames_;

    //- Boundary patch physical types
    wordList patchPhysicalTypes_;

    //- List of faces for every cell
    faceListList cellFaces_;

    //- List of each baffle face
    faceList baffleFaces_;

    // cell table id for each cell
    labelList cellTableId_;

    // cell table persistent data saved as a dictionary
    Map<dictionary> cellTable_;

    // subclasses are required to supply this information
    virtual void readMesh(const fileName &, const scalar scaleFactor) = 0;

private:
    // Private data

    //- Point-cell addressing. Used for topological analysis
    // Warning. This point cell addressing list potentially contains
    // duplicate cell entries. Use additional checking
    mutable labelListList* pointCellsPtr_;

    // Private static data members

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshReader(const meshReader&);

        //- Disallow default bitwise assignment
        void operator=(const meshReader&);

        void writeCellTable(const polyMesh&);
        void writeInterfaces(const polyMesh&);

        void writeAux(const polyMesh& mesh)
        {
	    writeCellTable(mesh);
	    writeInterfaces(mesh);
	}

        //- Calculate pointCells
        void calcPointCells() const;

        const labelListList& pointCells() const;

        //- Make polyhedral cells and global faces if the mesh is polyhedral
        void createPolyCells();

        //- add in boundary face
        void addPolyBoundaryFace
        (
	    const label cellId,
	    const label cellFaceId,
	    const label nCreatedFaces
	);

        //- Make polyhedral boundary from shape boundary
        // (adds more faces to the face list)
        void createPolyBoundary();

        //- Add polyhedral boundary
        List<polyPatch*> polyBoundaryPatches(const polyMesh&);

        //- Clear extra storage before creation of the mesh to remove
        //  a memory peak
        void clearExtraStorage();

        void volFieldCellTableId(fvMesh &);

public:

    // Constructors
    meshReader
    (
	const fileName& file_or_prefix,
	const Time& runtime,
	const scalar scaleFactor = 1.0
    );

    // Destructor
    virtual ~meshReader();

    // Member Functions

    //- Return mesh points
    virtual const pointField& points() const = 0;

    //- Return list of faces for every cell
    const faceListList & cellFaces() const 
    {
	return cellFaces_;
    }

    //- Write mesh
    void writeMesh();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
