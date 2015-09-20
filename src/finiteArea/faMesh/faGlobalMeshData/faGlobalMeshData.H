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
    faGlobalMeshData

Description
    Various mesh related information for a parallel run.  Incomplete.
    HJ, 3/Jul/2007

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef faGlobalMeshData_H
#define faGlobalMeshData_H

#include "faProcessorTopology.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class faMesh;

/*---------------------------------------------------------------------------*\
                        Class faGlobalMeshData Declaration
\*---------------------------------------------------------------------------*/

class faGlobalMeshData
:
    public faProcessorTopology
{
    // Private data

        //- Reference to mesh
        const faMesh& mesh_;

        // Globally shared point addressing

            //- Total number of global points
            label nGlobalPoints_;

            //- Indices of local points that are globally shared
            labelList sharedPointLabels_;

            //- Indices of globally shared points in the master list
            //  This list contains all the shared points in the mesh
            labelList sharedPointAddr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        faGlobalMeshData(const faGlobalMeshData&);

        //- Disallow default bitwise assignment
        void operator=(const faGlobalMeshData&);


public:

    //- Runtime type information
    ClassName("faGlobalMeshData");


    // Constructors

        //- Construct from mesh
        faGlobalMeshData(const faMesh& mesh);


    // Destructor

        ~faGlobalMeshData();


    // Member Functions

        // Access

            //- Return mesh reference
            const faMesh& mesh() const;

        // Globally shared point addressing

            //- Return number of globally shared points
            label nGlobalPoints() const
            {
                return nGlobalPoints_;
            }

            //- Return indices of local points that are globally shared
            const labelList& sharedPointLabels() const
            {
                return sharedPointLabels_;
            }

            //- Return addressing into the complete globally shared points
            //  list
            const labelList& sharedPointAddr() const
            {
                return sharedPointAddr_;
            }

            //- Change global mesh data given a topological change.
            void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
