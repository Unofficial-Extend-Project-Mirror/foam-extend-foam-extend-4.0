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
    donorBasedLayeredOverlapFringe

Description
    Fringe algorithm based on donors from different regions. The algorithm waits
    until the final donor/acceptor assembly has been performed for all regions
    whose donor region is this region. Then, all donors are collected and the
    acceptors are cells neighbouring the donors nLayers towards the interior.
    Interior is defined either by user-specified points for each region or as a
    centre of volume of the donor cells in particular region.

    This fringe algorithm is intended to be used along with the faceCells fringe
    on the other side, where the cell sizes are not significantly different from
    each other (e.g. 10:1 cell ratio, where the finer cells are found on the
    background mesh will probably be problematic to correctly set-up).

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

SourceFiles
    donorBasedLayeredOverlapFringe.C

\*---------------------------------------------------------------------------*/

#ifndef donorBasedLayeredOverlapFringe_H
#define donorBasedLayeredOverlapFringe_H

#include "oversetFringe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class donorBasedLayeredOverlapFringe Declaration
\*---------------------------------------------------------------------------*/

class donorBasedLayeredOverlapFringe
:
    public oversetFringe
{
    // Private data

        //- Regions IDs from which the donors will be collected as a starting
        //  point
        labelList connectedRegionIDs_;

        //- Optional list of points representing a rough estimate of the centre
        //  for each underlying connected region. If these are not provided, the
        //  centre is calculated as the centre of all donors for a given
        //  connected region
        List<point> regionCentrePoints_;

        //- How many layers to move away from connected region donors to define
        //  acceptor (and holes)
        label nLayers_;

        //- Fringe hole cells
        mutable labelList* fringeHolesPtr_;

        //- Acceptor cells
        mutable labelList* acceptorsPtr_;

        //- Final donor/acceptor pairs for this region (fringe)
        mutable donorAcceptorList* finalDonorAcceptorsPtr_;


    // Private static data

        //- Distance tolerance to determine propagation direction. Note:
        //  absolute value, default = 0. Might be useful is some strange cases
        static const debug::tolerancesSwitch distTol_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        donorBasedLayeredOverlapFringe(const donorBasedLayeredOverlapFringe&);

        //- Disallow default bitwise assignment
        void operator=(const donorBasedLayeredOverlapFringe&);


        // Calculate hole and acceptor addressing
        void calcAddressing() const;

        // Clear addressing
        void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("donorBasedLayeredOverlap");


    // Constructors

        //- Construct from dictionary
        donorBasedLayeredOverlapFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    // Destructor

        virtual ~donorBasedLayeredOverlapFringe();


    // Member Functions

        //- Update iteration. Note: invalidates parameter
        virtual bool updateIteration
        (
            donorAcceptorList& donorAcceptorRegionData
        ) const;

        //- Return list of deactivated (hole) cells
        //  Fringe hole cells are collected in addition to geometric hole
        //  cells, which fall outside of all donor regions
        virtual const labelList& fringeHoles() const;

        //- Return list of acceptor cells
        virtual const labelList& candidateAcceptors() const;

        //- Return list of final donor acceptor pairs. Note: caller may
        //  invalidate finalDonorAcceptorsPtr_ for optimisation purposes
        virtual donorAcceptorList& finalDonorAcceptors() const;

        //- Update the fringe
        virtual void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
