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
    Foam::centreEdge

Description
    Implicit description of faceCentre to cellCentre edge
    (from tet decomposition).

    Stores
      - face label
      - whether in owner or neighbour of face

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centreEdge_H
#define centreEdge_H

#include "label.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class centreEdge Declaration
\*---------------------------------------------------------------------------*/

class centreEdge
{
    // Private data

        //- face label
        label faceLabel_;

        //- is in owner or neighbour of face
        bool inOwner_;

public:

    // Public classes

        //- Hash function
        class centreEdgeHash
        {

        public:

            centreEdgeHash()
            {}

            //- simple hashing function of labels
            label operator()(const centreEdge& ce, const label tableSize) const
            {
                return (ce.faceLabel() + ce.inOwner()) % tableSize;
            }
        };


    // Constructors

        //- Construct null
        inline centreEdge()
        :
            faceLabel_(-1),
            inOwner_(-1)
        {}

        //- Construct from components
        inline centreEdge(const label faceLabel, const bool inOwner)
        :
            faceLabel_(faceLabel),
            inOwner_(inOwner)
        {}

    // Member Functions

        label faceLabel() const
        {   
            return faceLabel_;
        }

        bool inOwner() const
        {   
            return inOwner_;
        }

        template <class T>
        T interpolate
        (
            const primitiveMesh& mesh,
            const Field<T>& cellField,
            const Field<T>& faceField,
            const scalar weight
        ) const
        {
            label cellI;
            if (inOwner_)
            {
                cellI = mesh.faceOwner()[faceLabel_];
            }
            else
            {
                cellI = mesh.faceNeighbour()[faceLabel_];
            }

            return
                (1-weight)*faceField[faceLabel_]
              + weight*cellField[cellI];
        }

        point coord(const primitiveMesh& mesh, const scalar weight) const
        {
            return interpolate
            (
                mesh,
                mesh.cellCentres(),
                mesh.faceCentres(),
                weight
            );
        }


    // Member Operators

        bool operator==(const centreEdge& ce) const
        {
            return
                (faceLabel() == ce.faceLabel())
             && (inOwner() == ce.inOwner());
        }


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream& os, const centreEdge& ce)
        {
            os  << token::BEGIN_LIST
                << ce.faceLabel_ << token::SPACE
                << ce.inOwner_
                << token::END_LIST;

            // Check state of Ostream
            os.check("Ostream& operator<<(Ostream&, const centreEdge&)");

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
