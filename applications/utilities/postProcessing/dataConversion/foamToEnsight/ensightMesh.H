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

Class
    Foam::ensightMesh

Description

SourceFiles
    ensightMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ensightMesh_H
#define ensightMesh_H

#include "cellSets.H"
#include "faceSets.H"
#include "HashTable.H"
#include "HashSet.H"
#include "fvMesh.H"
#include "OFstream.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class argList;

/*---------------------------------------------------------------------------*\
                           Class ensightMesh Declaration
\*---------------------------------------------------------------------------*/

class ensightMesh
{
public:
        class nFacePrimitives
        {
        public:

            label nPoints;
            label nTris;
            label nQuads;
            label nPolys;

            nFacePrimitives()
            :
                nPoints(0),
                nTris(0),
                nQuads(0),
                nPolys(0)
            {}
        };

private:

    // Private data

        //- Reference to the FOAM mesh
        const fvMesh& mesh_;

        //- Set binary file output
        bool binary_;

        //- The ensight part id for the first patch
        label patchPartOffset_;

        cellSets meshCellSets_;

        List<faceSets> boundaryFaceSets_;

        wordList allPatchNames_;

        List<labelList> allPatchProcs_;

        wordHashSet patchNames_;

        HashTable<nFacePrimitives> nPatchPrims_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ensightMesh(const ensightMesh&);

        //- Disallow default bitwise assignment
        void operator=(const ensightMesh&);

        void writePoints
        (
            const scalarField& pointsComponent,
            OFstream& ensightGeometryFile
        ) const;

        cellShapeList ensMap
        (
            const cellShapeList& cellShapes,
            const labelList& prims
        ) const;

        cellShapeList ensMap
        (
            const cellShapeList& cellShapes,
            const labelList& hexes,
            const labelList& wedges
        ) const;

        void writePrims
        (
            const cellShapeList& cellShapes,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        void writePolysNFaces
        (
            const labelList& polys,
            const cellList& cellFaces,
            OFstream& ensightGeometryFile
        ) const;

        void writePolysNPointsPerFace
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            OFstream& ensightGeometryFile
        ) const;

        void writePolysPoints
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        void writeAllPolys
        (
            const labelList& pointOffsets,
            OFstream& ensightGeometryFile
        ) const;

        void writeAllPrims
        (
            const char* key,
            const label nPrims,
            const cellShapeList& cellShapes,
            const labelList& pointOffsets,
            OFstream& ensightGeometryFile
        ) const;

        void writeFacePrims
        (
            const faceList& patchFaces,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        faceList ensMap
        (
            const faceList& patchFaces,
            const labelList& prims
        ) const;

        void writeAllFacePrims
        (
            const char* key,
            const labelList& prims,
            const label nPrims,
            const faceList& patchFaces,
            const labelList& pointOffsets,
            const labelList& patchProcessors,
            OFstream& ensightGeometryFile
        ) const;

        void writeNSidedNPointsPerFace
        (
            const faceList& patchFaces,
            OFstream& ensightGeometryFile
        ) const;

        void writeNSidedPoints
        (
            const faceList& patchFaces,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        void writeAllNSided
        (
            const labelList& prims,
            const label nPrims,
            const faceList& patchFaces,
            const labelList& pointOffsets,
            const labelList& patchProcessors,
            OFstream& ensightGeometryFile
        ) const;

        void writeAscii
        (
            const fileName& postProcPath,
            const word& prepend,
            const label timeIndex,
            Ostream& ensightCaseFile
        ) const;

        void writeBinary
        (
            const fileName& postProcPath,
            const word& prepend,
            const label timeIndex,
            Ostream& ensightCaseFile
        ) const;

        void writePrimsBinary
        (
            const cellShapeList& cellShapes,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeAllPrimsBinary
        (
            const char* key,
            const label nPrims,
            const cellShapeList& cellShapes,
            const labelList& pointOffsets,
            std::ofstream& ensightGeometryFile
        ) const;

        void writePolysNFacesBinary
        (
            const labelList& polys,
            const cellList& cellFaces,
            std::ofstream& ensightGeometryFile
        ) const;

        void writePolysNPointsPerFaceBinary
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            std::ofstream& ensightGeometryFile
        ) const;

        void writePolysPointsBinary
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeAllPolysBinary
        (
            const labelList& pointOffsets,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeAllFacePrimsBinary
        (
            const char* key,
            const labelList& prims,
            const label nPrims,
            const faceList& patchFaces,
            const labelList& pointOffsets,
            const labelList& patchProcessors,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeFacePrimsBinary
        (
            const faceList& patchFaces,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeNSidedPointsBinary
        (
            const faceList& patchFaces,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeNSidedNPointsPerFaceBinary
        (
            const faceList& patchFaces,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeAllNSidedBinary
        (
            const labelList& prims,
            const label nPrims,
            const faceList& patchFaces,
            const labelList& pointOffsets,
            const labelList& patchProcessors,
            std::ofstream& ensightGeometryFile
        ) const;

public:

    // Constructors

        //- Construct from fvMesh
        ensightMesh(const fvMesh&, const argList& args, const bool binary);


    //- Destructor
    ~ensightMesh();


    // Member Functions

        // Access

            const fvMesh& mesh() const
            {
                return mesh_;
            }

            const cellSets& meshCellSets() const
            {
                return meshCellSets_;
            }

            const List<faceSets>& boundaryFaceSets() const
            {
                return boundaryFaceSets_;
            }

            const wordList& allPatchNames() const
            {
                return allPatchNames_;
            }

            const List<labelList>& allPatchProcs() const
            {
                return allPatchProcs_;
            }

            const wordHashSet& patchNames() const
            {
                return patchNames_;
            }

            const HashTable<nFacePrimitives>& nPatchPrims() const
            {
                return nPatchPrims_;
            }

            //- The ensight part id for the first patch
            label patchPartOffset() const
            {
                return patchPartOffset_;
            }


    // I-O

        void write
        (
            const fileName& postProcPath,
            const word& prepend,
            const label timeIndex,
            Ostream& ensightCaseFile
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
