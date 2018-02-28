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
    Foam::oversetRegion

Description
    Overset mesh region

    Overset region:
    - collects its cells by cell zone
    - contains list of donor regions
    - contains list of patch names used to define a hole (walls)
    - contains the fringe selection algorithm

    Overset region provides data on
    - cells within region (local processor)
    - local acceptor, donor and hole cells
    - eligible local donors: cells that are neither acceptors nor holes

    Overset region controls the iterative process of fringe assembly with the
    help of fringe selection algorithm. The flowchart of the algorithm is:
    1. Someone (usually oversetMesh) asks region for donors and acceptors, which
       triggers the lazy evaluation in calcDonorAcceptorCells(). During the
       donor/acceptor assembly, all regions are interdependent because eligible
       donors and holes could change during the iterative process. Hence, asking
       a certain region for its donors/acceptors actually triggers the fringe
       assembly in all regions,
    2. The iterative process starts by going through each region, updating
       current donor/acceptors (see updateDonorAcceptors) which does the
       following:
       - Asks fringe handler a set of "fresh" acceptors,
       - Finds donors for these acceptors,
       - Ships of donor/acceptor pairs to fringe handler, which decided whether
         a suitable overlap has been found
    3. If the suitable overlap has been found for ALL regions, then we finalise
       the donor/acceptor pairs by calling finaliseDonorAcceptors(), which
       actually allocates requested fields (acceptorCells and donorCells)

Ideas to further improve efficiency
    1. During the possibly iterative process of finding suitable overlap with
       the help of fringe selection algorithm (see calcDonorAcceptorCells()), we
       are sending/receiving full donorAcceptor objects. A single donor/acceptor
       object holds: donor and acceptor cell numbers, processor IDs and cell
       centres. On top of that, it holds additional "extended donors" for a
       given acceptor, which is a list containing eligible neighbouring donors
       for a given "master" donor (see donorAcceptor). Something might be gained
       if we do the following:
           - Create a bareDonorAcceptor class which will hold only the data
             required by the iterative process (e.g. acceptor cell centre, donor
             cell number, processor ID and cell centre). Then when the iterative
             process is finished, simply finalise the data by extending bare
             data and communicate it once again.
           - As a first step, it would be sufficient just to remove the search
             for extended donors during the iterative process.
    2. I've written the iterative process in a way that it works both for serial
       and parallel runs (enabled by mapDistribute class design), which means
       that we might have significant overhead when doing serial runs due to
       unnecessary allocations, indirections and access. I've written it this
       way since I currently do not care about loss of efficiency for serial
       runs. VV, 2/Feb/2017.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Contributor
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

SourceFiles
    oversetRegion.C

\*---------------------------------------------------------------------------*/

#ifndef oversetRegion_H
#define oversetRegion_H

#include "fvMesh.H"
#include "donorAcceptorList.H"
#include "boundBox.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "optimisationSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class oversetMesh;
class oversetFringe;

/*---------------------------------------------------------------------------*\
                        Class oversetRegion Declaration
\*---------------------------------------------------------------------------*/

class oversetRegion
{
    // Private data

        //- Overset region name and cell zone name
        const word name_;

        //- Overset region index
        const label index_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Overset mesh reference
        const oversetMesh& oversetMesh_;

        //- Cell zone index
        const label zoneIndex_;

        //- List of region names that provide donors
        //  Note: ordering is important as it defines precedence
        wordList donorRegionNames_;

        //- Fringe selection algorithm
        autoPtr<oversetFringe> fringePtr_;


    // Demand-driven data

        //- List of region indices that provide donors
        //  Note: ordering is important as it defines precedence
        mutable labelList* donorRegionsPtr_;

        //- Acceptor region labels
        mutable labelList* acceptorRegionsPtr_;

        //- Acceptor cell data
        mutable donorAcceptorList* acceptorCellsPtr_;

        //- Donor cell data.  Donor cell is defined as the nearest cell in
        //  donor regions to the acceptor cell
        mutable donorAcceptorList* donorCellsPtr_;

        //- Cut hole cell labels
        mutable labelList* cutHoleCellsPtr_;

        //- Complete hole cell labels (cut + fringe holes)
        mutable labelList* holeCellsPtr_;

        //- Eligible donor cells from current region
        //  Eligible donors are region cells that are not holes or acceptors
        mutable labelList* eligibleDonorCellsPtr_;


        // Hole cutting tools

            //- Hole patches as a triangular surface
            mutable triSurface* holeTriMeshPtr_;

            //- Triangular hole surface search algorithm.  Contains complete
            //  cutting surface
            //  Used in hole cutting
            mutable triSurfaceSearch* holeSearchPtr_;


        // Donor search tools

            //- Overset region local (processor) bounding box
            mutable boundBox* localBoundsPtr_;

            //- Overset region global (region) bounding box
            mutable boundBox* globalBoundsPtr_;

            //- Cell octree for search in local region
            //  Used in donor cell search
            mutable indexedOctree<treeDataCell>* cellSearchPtr_;


            // Data for optimizing parallel comms

                //- List of bounding boxes for each region on each processor,
                //  e.g. procBoundBoxes[3][2] contains a bounding box of region
                //  2 on processor 3
                mutable List<List<boundBox> >* procBoundBoxesPtr_;

                //- Whether to use processor bounding box to determine where to
                //  send candidate acceptors. By default, this is switched on
                //  for parallel efficiency but it might fail for overset
                //  regions where one of the region bounding boxes is not fully
                //  within the other bounding box. Controls the behaviour of
                //  calcProcBoundBoxes()
                Switch useLocalBoundBoxes_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        oversetRegion(const oversetRegion&);

        //- Disallow default bitwise assignment
        void operator=(const oversetRegion&);


        //- Prepare hole triangular mesh
        void calcHoleTriMesh() const;

        //- Calculate donor regions
        void calcDonorRegions() const;

        //- Calculate acceptor regions
        void calcAcceptorRegions() const;

        //- Calculate donor and acceptor cells
        void calcDonorAcceptorCells() const;

        //- Calculate cut hole cells
        void calcCutHoleCells() const;

        //- Calculate complete hole cells (cut + fringe)
        void calcHoleCells() const;

        //- Calculate eligible donor cells
        void calcEligibleDonorCells() const;

        //- Calculate bounds
        void calcBounds() const;

        //- Calculate cell octree
        void calcCellSearch() const;

        //- Calculate bounding boxes on all processors
        void calcProcBoundBoxes() const;

        //- Clear all data
        void clearOut() const;

        //- Update donor/acceptors for global overset fringe algorithm control
        //  (called within calcDonorAcceptorCells(), even from other regions).
        //  Returns true if this region found suitable set of donors
        bool updateDonorAcceptors() const;

        //- Finalise donor/acceptor, called within calcDonorAcceptorCells()
        //  after a suitable overlap has been found for all regions
        void finaliseDonorAcceptors() const;


public:

    // Constructors

        //- Construct from dictionary
        oversetRegion
        (
            const word& name,
            const label index,
            const fvMesh& mesh,
            const oversetMesh& oversetMesh,
            const dictionary& dict
        );


    //- Destructor
    ~oversetRegion();


    // Member Functions

        // Access

            //- Return region name and cell zone name
            const word& name() const
            {
                return name_;
            }

            //- Return region index
            label index() const
            {
                return index_;
            }

            //- Return mesh reference
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return overset mesh reference
            const oversetMesh& overset() const
            {
                return oversetMesh_;
            }

            //- Return region cell zone index
            label zoneIndex() const
            {
                return zoneIndex_;
            }

            //- Return region cell zone
            const cellZone& zone() const
            {
                return mesh_.cellZones()[zoneIndex()];
            }

            //- Return region cell labels (implemented in terms of zone())
            const labelList& regionCells() const
            {
                return zone();
            }

            //- Return list of donor region indices
            const labelList& donorRegions() const;

            //- Return list of acceptor region indices
            const labelList& acceptorRegions() const;


        // Donor-acceptor addressing

            //- Return acceptor cell list
            const donorAcceptorList& acceptors() const;

            //- Return donor cell list
            const donorAcceptorList& donors() const;

            //- Return cut holes only
            const labelList& cutHoles() const;

            //- Return complete hole cell list (cut + fringe holes)
            const labelList& holes() const;

            //- Return eligible donor cell list
            const labelList& eligibleDonors() const;


        // Hole cutting search functionality

            //- Are hole patches present?
            bool holePatchesPresent() const;

            //- Hole triangular surface
            const triSurface& holeTriMesh() const;

            //- Return triangular surface search object for hole patches
            const triSurfaceSearch& holeSearch() const;


        // Donor cell search functionality

            //- Return local (processor) bounding box for region
            const boundBox& localBounds() const;

            //- Return global (region) bounding box for region
            const boundBox& globalBounds() const;

            //- Return cell octree for search in local region
            //  Only eligible donors are added to the list
            const indexedOctree<treeDataCell>& cellSearch() const;


            // Function for avoiding unnecessary parallel comms overheads

                //- Return list of region bounding boxes for each processor.
                //  Used to optimize out parallel comms since there is no need
                //  to send the acceptor data to a certain processor if the
                //  bounding box of the donor region on a given processor does
                //  not contain acceptor point
                const List<List<boundBox> >& procBoundBoxes() const;


        //- Update for topo change or mesh motion
        void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
