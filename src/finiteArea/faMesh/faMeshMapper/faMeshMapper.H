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
    Foam::faMeshMapper

Description
    Class holds all the necessary information for mapping fields associated
    with faMesh

Note
    In order to capture all necessary mesh sizes and mapping data, mapper
    is created with the OLD mesh, and provides new mesh data.
    In the process, field mapping information is assembled from the old faMesh
    and the mapping data

SourceFiles
    faMeshMapper.C

\*---------------------------------------------------------------------------*/

#ifndef faMeshMapper_H
#define faMeshMapper_H

#include "faceMapper.H"
#include "faAreaMapper.H"
#include "faEdgeMapper.H"
#include "faBoundaryMeshMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class faMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                        Class faMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class faMeshMapper
{
    // Private data

        //- Reference to mesh
        const faMesh& mesh_;


        // Old mesh data

            //- Number of old points
            label nOldPoints_;

            //- Number of old edges
            label nOldEdges_;

            //- Number of old internal edges
            label nOldInternalEdges_;

            //- Number of old faces
            label nOldFaces_;

            //- Old patch sizes
            labelList oldPatchSizes_;

            //- Old patch starts
            labelList oldPatchStarts_;

            //- Old patch edgeFaces
            labelListList oldPatchEdgeFaces_;


        // Mappers

            //- Area mapper
            faAreaMapper areaMap_;

            //- Edge mapper
            faEdgeMapper edgeMap_;

            //- Boundary mapper
            faBoundaryMeshMapper boundaryMap_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faMeshMapper(const faMeshMapper&);

        //- Disallow default bitwise assignment
        void operator=(const faMeshMapper&);


public:

    // Constructors

        //- Construct from components
        faMeshMapper(const faMesh& mesh, const mapPolyMesh& mpm);


    // Member Functions

        //- Return reference to mesh
        const faMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to objectRegistry storing fields. Can be
        //  removed once fields stored on pointMesh.
        const objectRegistry& thisDb() const
        {
            return mesh_.thisDb();
        }


        // Basic sizing information

            //- Return number of old points
            label nOldPoints() const
            {
                return nOldPoints_;
            }

            //- Return number of old edges
            label nOldEdges() const
            {
                return nOldEdges_;
            };

            //- Return number of old internal edges
            label nOldInternalEdges() const
            {
                return nOldInternalEdges_;
            };

            //- Return number of old faces
            label nOldFaces() const
            {
                return nOldFaces_;
            };

            //- Return old patch sizes
            const labelList& oldPatchSizes() const
            {
                return oldPatchSizes_;
            };

            //- Return old patch starts
            const labelList& oldPatchStarts() const
            {
                return oldPatchStarts_;
            };

            //- Return old patch edgeFaces
            const labelListList& oldPatchEdgeFaces() const
            {
                return oldPatchEdgeFaces_;
            };


        // Mappers

            //- Return surface mapper
            const faAreaMapper& areaMap() const
            {
                return areaMap_;
            }

            //- Return edge mapper
            const faEdgeMapper& edgeMap() const
            {
                return edgeMap_;
            }

            //- Return boundary mapper
            const faBoundaryMeshMapper& boundaryMap() const
            {
                return boundaryMap_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
