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

Class
    Foam::meshWriters::Elmer

Description
    Writes polyMesh in CSC/Elmer format: generates the files

       elmerMesh.header
       elmerMesh.nodes
       elmerMesh.elements
       elmerMesh.boundary
       elmerMesh.names

    The cellTableId and cellTable information are used (if available).
    Otherwise the cellZones are used (if available).

SourceFiles
    ElmerMeshWriter.C

\*---------------------------------------------------------------------------*/

#ifndef ElmerMeshWriter_H
#define ElmerMeshWriter_H

#include "meshWriter.H"
#include "wordReList.H"

namespace Foam
{

namespace meshWriters
{

class Elmer
:
    public meshWriter
{

private:

    // Elmer element type labels
    #define ELMER_ETYPE_BAD    0   // Not a valid element type
    #define ELMER_ETYPE_TRIA   303
    #define ELMER_ETYPE_QUAD   404
    #define ELMER_ETYPE_TET    504
    #define ELMER_ETYPE_PRISM  706
    #define ELMER_ETYPE_PYRAM  605
    #define ELMER_ETYPE_HEX    808

    const wordRe& faceZoneExcludePattern;

    //- Disallow default bitwise copy construct
    Elmer(const Elmer&);

    //- Disallow default bitwise assignment
    void operator=(const Elmer&);

    //- Get the Elmer face element type ID based on the number of vertices of a face
    label getFaceType(const label nvert, const word &zname) const;

    //- Save the mesh.header file containing the global counters
    bool writeHeader() const;

    //- Save the mesh.names file containing the names of all zones & patches
    void writeNames() const;

    //- Save the mesh.nodes file
    void writeNodes() const;

    //- Save the mesh.elements file. This file contains only the volume
    //- elements, not the faces on the boundaries
    void writeElements () const;

    //- Save the mesh.boundary file. This file contains only the faces elements
    //- on boundaries.
    void writeBoundary () const;

    //- Build the cell table information
    void getCellTable ();


public:

    // Constructors

        //- Construct from polyMesh and option arguments.
        Elmer
        (
            const polyMesh&,
            const wordRe& excludePattern,
            const scalar scaleFactor = 1.0
        );

    // Destructor

        virtual ~Elmer();

    // Member Functions

        //- Write mesh in elmer format
        virtual bool write(const fileName& dirName = fileName::null) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // End namespace meshWriters


#endif

// ************************************************************************* //
