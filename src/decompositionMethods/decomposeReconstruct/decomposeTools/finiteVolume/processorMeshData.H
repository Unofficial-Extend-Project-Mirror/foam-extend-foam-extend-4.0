/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::processorMeshData

Description
    The class holds mesh data required to assemble a single processor mesh
    from the decomposition.  The data is raw but can be assembled into
    a complete mesh

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    processorMeshData.C

\*---------------------------------------------------------------------------*/

#ifndef processorMeshData_H
#define processorMeshData_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class processorMeshData Declaration
\*---------------------------------------------------------------------------*/

class processorMeshData
{
    // Private data

        // Internal mesh data

            //- Points
            pointField procPoints_;

            //- Faces
            faceList procFaces_;

            //- Face owner
            labelList procOwner_;

            //- Face neighbour
            labelList procNeighbour_;


        // Boundary data, separated into existing and new patches

            //- Boundary addressing
            labelList boundaryAddressing_;

            //- Regular patch sizes
            labelList patchSizes_;

            //- Regular patch starts
            labelList patchStarts_;

            //- Processor patch sizes
            labelList procPatchSizes_;

            //- Processor patch starts
            labelList procPatchSizes_;


        // Zones

            //- Point zones
            labelListList pointZones_;

            //- Face zones
            labelListList faceZones_;

            //- Face zone flips
            labelListList faceZoneFlips_;

            //- Cell zones
            labelListList cellZones_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        processorMeshData(const processorMeshData&);

        //- Disallow default bitwise assignment
        void operator=(const processorMeshData&);


public:

    // Static data members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct from components
        processorMeshData(const dataType& data);

        //- Construct from Istream
        processorMeshData(Istream&);

        //- Construct as copy
        processorMeshData(const processorMeshData&);


    // Selectors

        //- Select null constructed
        static autoPtr<processorMeshData> New();


    //- Destructor
    ~processorMeshData();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const processorMeshData&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, processorMeshData&);
        friend Ostream& operator<<(Ostream&, const processorMeshData&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "processorMeshDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
