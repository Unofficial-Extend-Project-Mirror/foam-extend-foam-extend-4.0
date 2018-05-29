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
    oversetFringe

Description
    Virtual base class for overset fringe selection.

    oversetFringe
    - is owned by a particular oversetRegion
    - defines the hole cells for its region
    - defines acceptor cells for its region
    - controls the iterative process for finding suitable overlap (as needed by
      overlapFringe algorithm)
    - defines donorAcceptorList used by oversetRegion to create combined data

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Contributor
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

SourceFiles
    oversetFringe.C
    newOversetFringe.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFringe_H
#define oversetFringe_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "donorAcceptorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class oversetRegion;

/*---------------------------------------------------------------------------*\
                         Class oversetFringe Declaration
\*---------------------------------------------------------------------------*/

class oversetFringe
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to oversetRegion
        const oversetRegion& region_;

        //- Flag indicating whether a suitable overlap has been found
        mutable bool foundSuitableOverlap_;


    // Private Member Functions

        //- Disallow copy construct
        oversetFringe(const oversetFringe&);

        //- Disallow default bitwise assignment
        void operator=(const oversetFringe&);


protected:

    // Protected Member Functions

        //- Set whether a suitable overlap has been found (called by derived
        //  classes during calls to updateIteration())
        void updateSuitableOverlapFlag(const bool isFound) const
        {
            foundSuitableOverlap_ = isFound;
        }


public:

    //- Runtime type information
    TypeName("oversetFringe");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            oversetFringe,
            dictionary,
            (
                const fvMesh& mesh,
                const oversetRegion& region,
                const dictionary& dict
            ),
            (mesh, region, dict)
        );


    // Selectors

        //- Return a reference to the selected thermal model
        static autoPtr<oversetFringe> New
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        oversetFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    // Destructor

        virtual ~oversetFringe()
        {}


    // Member Functions

        // Access

            //- Return reference to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return reference to region
            const oversetRegion& region() const
            {
                return region_;
            }


        // Iteration control

            //- Has a suitable overlap been found for this fringe (region)?
            bool foundSuitableOverlap() const
            {
                return foundSuitableOverlap_;
            }

            //- Update iteration, returns whether a suitable fringe has been
            //  found
            virtual bool updateIteration
            (
                donorAcceptorList& donorAcceptorRegionData
            ) const = 0;


        // Access to fringe data

            //- Return list of deactivated (hole) cells
            //  Fringe hole cells are collected in addition to geometric hole
            //  cells, which fall outside of all donor regions
            virtual const labelList& fringeHoles() const = 0;

            //- Return list of candidate acceptor cells
            virtual const labelList& candidateAcceptors() const = 0;

            //- Return reference to final donor acceptor pairs
            virtual donorAcceptorList& finalDonorAcceptors() const = 0;


        //- Update the fringe
        virtual void update() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
