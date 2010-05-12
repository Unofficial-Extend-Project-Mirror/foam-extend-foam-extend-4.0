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
    Foam::cellAddressing

Description
    Additional addressing for a cellModel. (cellModel themselves only
    contain cell-points). Used in cellFeatures class.

SourceFiles
    cellAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef cellAddressing_H
#define cellAddressing_H

#include "cellModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class cellModel;

/*---------------------------------------------------------------------------*\
                           Class cellAddressing Declaration
\*---------------------------------------------------------------------------*/

class cellAddressing
{
    // Private data

        faceList modelFaces_;

        edgeList edges_;

        labelListList edgeFaces_;

        labelListList faceEdges_;

        labelListList pointEdges_;

    // Private static functions        

        static labelList makeIdentity(const label size);

    // Private static functions        

        label findEdge(const edge& e); 

public:

    // Constructors

        //- Construct from components
        cellAddressing(const cellModel& model);

    // Destructor

    // Member Functions

        // Access

            inline const faceList& modelFaces() const
            {
                return modelFaces_;
            }

            inline const edgeList& edges() const
            {
                return edges_;
            }

            inline const labelListList& edgeFaces() const
            {
                return edgeFaces_;
            }

            inline const labelListList& faceEdges() const
            {
                return faceEdges_;
            }

            inline const labelListList& pointEdges() const
            {
                return pointEdges_;
            }


        // Check

        // Edit

        // Write


    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "cellAddressingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
