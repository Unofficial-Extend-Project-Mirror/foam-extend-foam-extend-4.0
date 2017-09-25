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
    Foam::polyRefinementHistory

Description
    All refinement history. Used in unrefinement.

    - visibleCells: valid for the current mesh and contains per cell -1
      (cell unrefined) or an index into splitCells_.
    - splitCells: for every split contains the parent (also index into
      splitCells) and optionally a subsplit as n indices into splitCells, where
      n is the number of points of a cell. Note that the numbers in splitCells
      are not cell labels, they are purely indices into splitCells.

    E.g. 2 cells, cell 1 gets refined so end up with 9 cells:
    @verbatim
        // splitCells
        9
        (
        -1 (1 2 3 4 5 6 7 8)
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        )

        // visibleCells
        9(-1 1 2 3 4 5 6 7 8)
    @endverbatim


    So cell0 (visibleCells=-1) is unrefined.
    Cells 1-8 have all valid splitCells entries which are:
      - parent:0
      - subsplits:0()

    The parent 0 refers back to the splitcell entries.


SourceFiles
    polyRefinementHistory.C

\*---------------------------------------------------------------------------*/

#ifndef polyRefinementHistory_H
#define polyRefinementHistory_H

#include "DynamicList.H"
#include "dynamicLabelList.H"
#include "labelList.H"
#include "FixedList.H"
#include "SLList.H"
#include "autoPtr.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapPolyMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                     Class polyRefinementHistory Declaration
\*---------------------------------------------------------------------------*/

class polyRefinementHistory
:
    public regIOobject
{
public:

    class splitPolyCell
    {
    public:

        // Index to original splitCell this cell was refined off from
        // -1: top level cell
        // -2: free splitCell (so should also be in freeSplitCells_)
        label parent_;

        //- cells this cell was refined into
        autoPtr<labelList> addedCellsPtr_;

        //- Construct null (parent = -1)
        splitPolyCell();

        //- Construct from parent
        splitPolyCell(const label parent);

        //- Construct from Istream
        splitPolyCell(Istream& is);

        //- Construct as deep copy
        splitPolyCell(const splitPolyCell&);

        //- Copy operator since autoPtr otherwise 'steals' storage.
        void operator=(const splitPolyCell& s)
        {
            // Check for assignment to self
            if (this == &s)
            {
                FatalErrorIn
                (
                    "splitPolyCell::operator=(const Foam::splitPolyCell&)"
                )   << "Attempted assignment to self"
                    << abort(FatalError);
            }

            parent_ = s.parent_;

            addedCellsPtr_.reset
            (
                s.addedCellsPtr_.valid()
              ? new labelList(s.addedCellsPtr_())
              : NULL
            );
        }

        bool operator==(const splitPolyCell& s) const
        {
            if (addedCellsPtr_.valid() != s.addedCellsPtr_.valid())
            {
                return false;
            }
            else if (parent_ != s.parent_)
            {
                return false;
            }
            else if (addedCellsPtr_.valid())
            {
                return addedCellsPtr_() == s.addedCellsPtr_();
            }
            else
            {
                return true;
            }
        }

        bool operator!=(const splitPolyCell& s) const
        {
            return !operator==(s);
        }

        friend Istream& operator>>(Istream&, splitPolyCell&);
        friend Ostream& operator<<(Ostream&, const splitPolyCell&);
    };


private:

    TypeName("polyRefinementHistory");

    // Private data

        //- Storage for splitCells
        DynamicList<splitPolyCell> splitCells_;

        //- Unused indices in splitCells
        dynamicLabelList freeSplitCells_;

        //- Currently visible cells. Indices into splitCells.
        labelList visibleCells_;


    // Private member functions

        //- Debug write
        static void writeEntry
        (
            const List<splitPolyCell>&,
            const splitPolyCell&
        );

        //- Debug write
        static void writeDebug
        (
            const labelList&,
            const List<splitPolyCell>&
        );

        //- Check consistency of structure, i.e. indices into splitCells_.
        void checkIndices() const;

        //- Allocate a splitCell. Return index in splitCells_.
        label allocateSplitCell
        (
            const label parent,
            const label i,
            const label nCells
        );

        //- Free a splitCell.
        void freeSplitCell(const label index);

        //- Mark entry in splitCells. Recursively mark its parent and subs.
        void markSplit
        (
            const label,
            labelList& oldToNew,
            DynamicList<splitPolyCell>&
        ) const;

        void countProc
        (
            const label index,
            const label newProcNo,
            labelList& splitCellProc,
            labelList& splitCellNum
        ) const;

public:

    // Constructors

        //- Construct (read) given an IOobject
        polyRefinementHistory(const IOobject&);

        //- Construct (read) or construct null
        polyRefinementHistory
        (
            const IOobject&,
            const List<splitPolyCell>& splitCells,
            const labelList& visibleCells
        );

        //- Construct (read) or construct from initial number of cells
        //  (all visible)
        polyRefinementHistory(const IOobject&, const label nCells);

        //- Construct as copy
        polyRefinementHistory(const IOobject&, const polyRefinementHistory&);

        //- Construct from Istream
        polyRefinementHistory(const IOobject&, Istream&);


    // Member Functions

        //- Per cell in the current mesh (i.e. visible) either -1 (unrefined)
        //  or an index into splitCells.
        const labelList& visibleCells() const
        {
            return visibleCells_;
        }

        //- Storage for splitPolyCells.
        const DynamicList<splitPolyCell>& splitCells() const
        {
            return splitCells_;
        }

        //- Cache of unused indices in splitCells
        const dynamicLabelList& freeSplitCells() const
        {
            return freeSplitCells_;
        }

        //- Is there unrefinement history. Note that this will fall over if
        //  there are 0 cells in the mesh. But this gives problems with
        //  lots of other programs anyway.
        bool active() const
        {
            return visibleCells_.size() > 0;
        }

        //- Get parent of cell
        label parentIndex(const label cellI) const
        {
            label index = visibleCells_[cellI];

            if (index < 0)
            {
                FatalErrorIn("polyRefinementHistory::parentIndex(const label)")
                    << "Cell " << cellI << " is not visible"
                    << abort(FatalError);
            }
            return splitCells_[index].parent_;
        }

        //- Store splitting of cell into n (number of points for a cell)
        void storeSplit
        (
            const label cellI,
            const labelList& addedCells
        );

        //- Store combining n cells into master
        void combineCells
        (
            const label masterCellI,
            const labelList& combinedCells
        );

        //- Update numbering for mesh changes
        void updateMesh(const mapPolyMesh&);

        //- Update numbering for subsetting
        void subset
        (
            const labelList& pointMap,
            const labelList& faceMap,
            const labelList& cellMap
        );

        //- Update local numbering for mesh redistribution.
        //  Can only distribute clusters sent across in one go; cannot
        //  handle parts recombined in multiple passes.
        void distribute(const mapDistributePolyMesh&);

        //- Compact splitCells_. Removes all freeSplitCells_ elements.
        void compact();

        //- Extend/shrink storage. additional visibleCells_ elements get
        //  set to -1.
        void resize(const label nCells);

        //- Debug write
        void writeDebug() const;


        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, polyRefinementHistory&);
        friend Ostream& operator<<(Ostream&, const polyRefinementHistory&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
