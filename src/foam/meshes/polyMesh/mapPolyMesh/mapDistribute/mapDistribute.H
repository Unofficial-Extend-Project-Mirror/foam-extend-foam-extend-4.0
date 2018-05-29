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
    Foam::mapDistribute

Description
    Class containing processor-to-processor mapping information.

    We store mapping from the bits-to-send to the complete starting list
    (subXXXMap) and from the received bits to their location in the new
    list (constructXXXMap).

Note:
    Schedule is a list of processor pairs (one send, one receive. One of
    them will be myself) which forms a scheduled (i.e. non-buffered) exchange.
    See distribute on how to use it.
    Note2: number of items send on one processor have to equal the number
    of items received on the other processor.

SourceFiles
    mapDistribute.C

\*---------------------------------------------------------------------------*/

#ifndef mapDistribute_H
#define mapDistribute_H

#include "labelList.H"
#include "labelPair.H"
#include "Pstream.H"
#include "boolList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;
class globalIndex;

/*---------------------------------------------------------------------------*\
                        Class mapDistribute Declaration
\*---------------------------------------------------------------------------*/

class mapDistribute
{
    // Private data

        //- Size of reconstructed data
        label constructSize_;

        //- Maps from subsetted data back to original data
        labelListList subMap_;

        //- Maps from subsetted data to new reconstructed data
        labelListList constructMap_;

        //- Schedule
        mutable autoPtr<List<labelPair> > schedulePtr_;


   // Private Member Functions

        static void checkReceivedSize
        (
            const label procI,
            const label expectedSize,
            const label receivedSize
        );

        void calcCompactAddressing
        (
            const globalIndex& globalNumbering,
            const labelList& elements,
            List<Map<label> >& compactMap
        ) const;

        void calcCompactAddressing
        (
            const globalIndex& globalNumbering,
            const labelListList& elements,
            List<Map<label> >& compactMap
        ) const;

        void exchangeAddressing
        (
            const int tag,
            const globalIndex& globalNumbering,
            labelList& elements,
            List<Map<label> >& compactMap,
            labelList& compactStart
        );

        void exchangeAddressing
        (
            const int tag,
            const globalIndex& globalNumbering,
            labelListList& elements,
            List<Map<label> >& compactMap,
            labelList& compactStart
        );


public:

    // Declare name of the class and its debug switch
    ClassName("mapDistribute");


    // Constructors

        //- Construct from components
        mapDistribute
        (
            const label constructSize,
            const labelListList& subMap,
            const labelListList& constructMap
        );

        //- (optionally destructively) construct from components
        mapDistribute
        (
            const label constructSize,
            labelListList& subMap,
            labelListList& constructMap,
            const bool reUse                // clone or reuse
        );

        //- Construct from reverse addressing: per data item the send
        //  processor and the receive processor. All processors get same data.
        mapDistribute
        (
            const labelList& sendProcs,
            const labelList& recvProcs
        );


        //- Construct from list of (possibly) remote elements in globalIndex
        //  numbering (or -1). Determines compact numbering (see above) and
        //  distribute map to get data into this ordering and renumbers the
        //  elements to be in compact numbering.
        mapDistribute
        (
            const globalIndex& globalNumbering,
            labelList& elements,
            List<Map<label> >& compactMap,
            const int tag = Pstream::msgType()
        );

        //- Special variant that works with the info sorted into bins
        //  according to local indices. E.g. think cellCells where
        //  cellCells[localCellI] is a list of global cells
        mapDistribute
        (
            const globalIndex&,
            labelListList& cellCells,
            List<Map<label> >& compactMap,
            const int tag = Pstream::msgType()
        );

        //- Construct copy
        mapDistribute(const mapDistribute&);


    // Member Functions

        // Access

            //- Constructed data size
            label constructSize() const
            {
                return constructSize_;
            }

            //- Constructed data size
            label& constructSize()
            {
                return constructSize_;
            }

            //- From subsetted data back to original data
            const labelListList& subMap() const
            {
                return subMap_;
            }

            //- From subsetted data back to original data
            labelListList& subMap()
            {
                return subMap_;
            }

            //- From subsetted data to new reconstructed data
            const labelListList& constructMap() const
            {
                return constructMap_;
            }

            //- From subsetted data to new reconstructed data
            labelListList& constructMap()
            {
                return constructMap_;
            }

            //- Calculate a schedule. See above.
            static List<labelPair> schedule
            (
                const labelListList& subMap,
                const labelListList& constructMap,
                const int tag
            );

            //- Return a schedule. Demand driven. See above.
            const List<labelPair>& schedule() const;


        // Other

            //- Helper for construct from globalIndex. Renumbers element
            //  (in globalIndex numbering) into compact indices.
            static label renumber
            (
                const globalIndex&,
                const List<Map<label> >& compactMap,
                const label globalElement
            );

            //- Compact maps. Gets per field a bool whether it is used locally
            //  and works out itself what this side and sender side can remove
            //  from maps.
            void compact
            (
                const boolList& elemIsUsed,
                const int tag = Pstream::msgType()
            );

            //- Distribute data. Note:schedule only used for Pstream::scheduled
            //  for now, all others just use send-to-all, receive-from-all.
            template<class T>
            static void distribute
            (
                const Pstream::commsTypes commsType,
                const List<labelPair>& schedule,
                const label constructSize,
                const labelListList& subMap,
                const labelListList& constructMap,
                List<T>&,
                const int tag = Pstream::msgType()
            );

            //- Distribute data. If multiple processors writing to same
            //  position adds contributions using cop.
            template<class T, class CombineOp>
            static void distribute
            (
                const Pstream::commsTypes commsType,
                const List<labelPair>& schedule,
                const label constructSize,
                const labelListList& subMap,
                const labelListList& constructMap,
                List<T>&,
                const CombineOp& cop,
                const T& nullValue,
                const int tag = Pstream::msgType()
            );

            //- Distribute data using default commsType.
            template<class T>
            void distribute
            (
                 List<T>& fld,
                const int tag = Pstream::msgType()
            ) const
            {
                if
                (
                    Pstream::defaultCommsType() == Pstream::nonBlocking
                 && contiguous<T>()
                )
                {
                    distribute
                    (
                        Pstream::nonBlocking,
                        List<labelPair>(),
                        constructSize_,
                        subMap_,
                        constructMap_,
                        fld
                    );
                }
                else if (Pstream::defaultCommsType() == Pstream::scheduled)
                {
                    distribute
                    (
                        Pstream::scheduled,
                        schedule(),
                        constructSize_,
                        subMap_,
                        constructMap_,
                        fld
                    );
                }
                else
                {
                    distribute
                    (
                        Pstream::blocking,
                        List<labelPair>(),
                        constructSize_,
                        subMap_,
                        constructMap_,
                        fld
                    );
                }
            }

            //- Debug: print layout
            void printLayout(Ostream& os) const;

            //- Correct for topo change
            void updateMesh(const mapPolyMesh&)
            {
                notImplemented
                (
                    "mapDistribute::updateMesh(const mapPolyMesh&)"
                );
            }


    // Member Operators

        void operator=(const mapDistribute&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mapDistributeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
