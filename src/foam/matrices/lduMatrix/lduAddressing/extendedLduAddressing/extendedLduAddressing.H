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
    Foam::extendedLduAddressing

Description
    The class contains data and calculates demand driven extended ldu
    addressing of arbitrary level. Level 0 is not allowed as it is the same
    as ordinary lduAddressing. Level 1 considers neighbours of neighbours
    for a given cell, etc.

    The addressing is completely defined with extended variants of:
    owner/neighbour, losort, ownerStart and losortStart: all being
    demand-driven data. In order to define extended owner/neighbour lists,
    additional "virtual" faces are assumed between cells not sharing a face.

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved

SourceFiles
    extendedLduAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef extendedLduAddressing_H
#define extendedLduAddressing_H

#include "lduAddressing.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class extendedLduAddressing Declaration
\*---------------------------------------------------------------------------*/

class extendedLduAddressing
{
    // Private data

        //- Reference to basic lduAddressing
        const lduAddressing& lduAddr_;

        //- Extension level
        const label p_;


    // Demand-driven data

        //- Extended lower (owner) addressing
        mutable labelList* extendedLowerPtr_;

        //- Extended upper (neighbour) addressing
        mutable labelList* extendedUpperPtr_;

        //- Mapping between ldu and extendedLdu faces for lower/upper
        mutable labelList* faceMapPtr_;

        //- Extended losort addressing
        mutable labelList* extendedLosortPtr_;

        //- Extended owner start addressing
        mutable labelList* extendedOwnerStartPtr_;

        //- Extended losort start addressing
        mutable labelList* extendedLosortStartPtr_;


    // Private Member Functions

        // Copy control

            //- Disallow default bitwise copy construct
            extendedLduAddressing(const extendedLduAddressing&);

            //- Disallow default bitwise assignment
            void operator=(const extendedLduAddressing&);


        // Helper functions

            //- Clear all demand driven data - helper function
            void clearOut() const;

            //- Recursive function that inserts cell neighbours in a hash set
            //  Parameters:
            //  masterCellI - searching for extended neighbours of this cell
            //  nbrCells - neighbouring cells of a cell (based on level)
            //  cellCells - passing as a parameter for minor optimisation
            //  extendedNbrs - HashSet containing labels of extended neighbours
            //  curLevel - current level needed for recursion control
            void markNeighbours
            (
                const label& masterCellI,
                const labelList& nbrCells,
                const labelListList& cellCells,
                HashSet<label>& extendedNbrs,
                label curLevel
            ) const;


        // Calculation of demand driven data

            //- Calculate cells cells using lower/upper addressing
            //  into provided storage
            void calcCellCells(labelListList& cc) const;

            //- Calculate extended lower/upper
            void calcExtendedLowerUpper() const;

            //- Calculate face map
            void calcFaceMap() const;

            //- Calculate extended losort
            void calcExtendedLosort() const;

            //- Calculate extended owner start
            void calcExtendedOwnerStart() const;

            //- Calculate extended losort start
            void calcExtendedLosortStart() const;


public:

    // Declare name of the class and its debug switch
    TypeName("extendedLduAddressing");


    // Constructors

        //- Construct given lduAddressing
        extendedLduAddressing
        (
            const lduAddressing& lduAddr,
            const label extensionLevel
        );


    // Destructor

        virtual ~extendedLduAddressing();


    // Member Functions

        // Access

            //- Return basic lduAddressing
            const lduAddressing& lduAddr() const
            {
                return lduAddr_;
            }

            //- Return lower addressing (i.e. lower label = upper triangle)
            const unallocLabelList& lowerAddr() const
            {
                return lduAddr_.lowerAddr();
            }

            //- Return upper addressing (i.e. upper label)
            const unallocLabelList& upperAddr() const
            {
                return lduAddr_.upperAddr();
            }

            //- Return extended lower addressing
            const unallocLabelList& extendedLowerAddr() const;

            //- Return extended upper addressing
            const unallocLabelList& extendedUpperAddr() const;

            //- Return basic to extended ldu addressing mapping list
            const unallocLabelList& faceMap() const;

            //- Return extended losort addressing
            const unallocLabelList& extendedLosortAddr() const;

            //- Return extended owner start addressing
            const unallocLabelList& extendedOwnerStartAddr() const;

            //- Return extended losort start addressing
            const unallocLabelList& extendedLosortStartAddr() const;

            //- Return off-diagonal index given extended owner and neighbour
            //  label
            label extendedTriIndex(const label a, const label b) const;


        // Mesh object base class functions

            //- Update after mesh motion:
            //  Needed for mesh object - does nothing
            virtual bool movePoints() const
            {
                return true;
            }

            //- Update after topology change:
            //  Deletes all demand-driven data when mesh changes
            virtual bool updateMesh(const mapPolyMesh& mpm) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
