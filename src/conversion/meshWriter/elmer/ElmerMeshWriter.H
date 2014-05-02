/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
class Elmer: public meshWriter
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

        const wordRe* faceZoneExcludePattern;

        //- Disallow default bitwise copy construct
        Elmer(const Elmer&);

        //- Disallow default bitwise assignment
        void operator=(const Elmer&);

        label getFaceType(const label nvert, const word &zname) const;
        bool  writeHeader   () const;
        void  writeNames    () const;
        void  writeNodes    () const;
        void  writeElements () const;
        void  writeBoundary () const;
        void  getCellTable  ();

public:
        Elmer(const polyMesh&, const scalar scaleFactor = 1.0, const wordRe* excludePattern = 0);
        virtual ~Elmer();
        virtual bool write(const fileName& dirName = fileName::null) const;
};
}
}

#endif

// ************************************************************************* //
