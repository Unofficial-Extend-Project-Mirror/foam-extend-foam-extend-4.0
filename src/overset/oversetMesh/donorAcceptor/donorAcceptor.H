/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::donorAcceptor

Description
    Class holds donor and acceptor data.  Used for searching and communications

SourceFiles

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Contributor
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef donorAcceptor_H
#define donorAcceptor_H

#include "label.H"
#include "point.H"
#include "DynamicList.H"
#include "boolList.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class donorAcceptor;

inline bool operator==(const donorAcceptor&, const donorAcceptor&);
inline bool operator!=(const donorAcceptor&, const donorAcceptor&);

Istream& operator>>(Istream&, donorAcceptor&);
Ostream& operator<<(Ostream&, const donorAcceptor&);


/*---------------------------------------------------------------------------*\
                        Class donorAcceptor Declaration
\*---------------------------------------------------------------------------*/

class donorAcceptor
{
public:

    // Public data types
    typedef DynamicList<label, 10> DynamicLabelList;
    typedef DynamicList<point, 10> DynamicPointList;


private:

    // Private data

        // Acceptor side

            //- Acceptor cell number
            label acceptorCell_;

            //- Acceptor processor number
            label acceptorProcNo_;

            //- Acceptor location
            point acceptorPoint_;


        // Donor side

            //- Acceptor cell number
            label donorCell_;

            //- Donor processor number
            label donorProcNo_;

            //- Donor location
            point donorPoint_;

            //- Extended donor cell numbers
            DynamicLabelList extendedDonorCells_;

            //- Extended donor points, used to construct weights in different
            //  overset interpolation schemes
            DynamicPointList extendedDonorPoints_;

            // Note that the processor number is the same for all extended
            // donors


public:

    // Constructors

        //- Construct null
        inline donorAcceptor()
        {}

        //- Construct from acceptor data
        inline donorAcceptor
        (
            const label acceptorCell,
            const label acceptorProcNo,
            const point& acceptorPoint
        )
        :
            acceptorCell_(acceptorCell),
            acceptorProcNo_(acceptorProcNo),
            acceptorPoint_(acceptorPoint),
            donorCell_(-1),
            donorProcNo_(-1),
            donorPoint_(vector::zero),
            extendedDonorCells_(),
            extendedDonorPoints_()
        {}

        //- Construct from Istream
        inline donorAcceptor(Istream& is)
        :
            acceptorCell_(readLabel(is)),
            acceptorProcNo_(readLabel(is)),
            acceptorPoint_(is),
            donorCell_(readLabel(is)),
            donorProcNo_(readLabel(is)),
            donorPoint_(is),
            extendedDonorCells_(is),
            extendedDonorPoints_(is)
        {}

        //- Copy constructor - default


    //- Destructor - default


    // Member Functions

        // Acceptor side

            //- Return acceptor cell number
            label acceptorCell() const
            {
                return acceptorCell_;
            }

            //- Return access to acceptor cell number
            label& acceptorCell()
            {
                return acceptorCell_;
            }

            //- Return acceptor processor number
            label acceptorProcNo() const
            {
                return acceptorProcNo_;
            }

            //- Return access to acceptor processor number
            label& acceptorProcNo()
            {
                return acceptorProcNo_;
            }

            //- Return acceptor location
            const point& acceptorPoint() const
            {
                return acceptorPoint_;
            }


        // Donor side

            //- Has a donor been found?
            bool donorFound() const
            {
                return donorCell_ > -1;
            }

            //- Return donor cell number
            label donorCell() const
            {
                return donorCell_;
            }

            //- Return access to donor cell number
            label& donorCell()
            {
                return donorCell_;
            }

            //- Return donor processor number
            label donorProcNo() const
            {
                return donorProcNo_;
            }

            //- Return access to donor processor number
            label& donorProcNo()
            {
                return donorProcNo_;
            }

            //- Return donor point
            const point& donorPoint() const
            {
                if (donorFound())
                {
                    return donorPoint_;
                }
                else
                {
                    return vector::max;
                }
            }

            //- Return extended donor cell numbers
            const DynamicLabelList& extendedDonorCells() const
            {
                return extendedDonorCells_;
            }

            //- Return extended donor cell centres
            const DynamicPointList& extendedDonorPoints() const
            {
                return extendedDonorPoints_;
            }

            //- Return donor to acceptor distance
            scalar distance() const
            {
                if (donorFound())
                {
                    return mag(acceptorPoint_ - donorPoint_);
                }
                else
                {
                    return GREAT;
                }
            }

        // Edit

            //- Set hit: donor found
            void setDonor
            (
                const label& donorCell,
                const label& donorProcNo,
                const point& donorPoint
            )
            {
                donorCell_ = donorCell;
                donorProcNo_ = donorProcNo;
                donorPoint_ = donorPoint;
            }

            //- Set extended donors by going through neighbours of currently set
            //  "best" donor, taking into account whether the cell is eligible
            //  donor (not acceptor or hole).
            //  Note 1: we do not check whether "best" donor is actually set for
            //          better performance.
            //  Note 2: setting extended donors is not necessary for injection
            //          interpolation, consider optimising for this case.
            void setExtendedDonors
            (
                const boolList& eligibleDonors,
                const labelListList& cCells,
                const vectorField& cellCentres
            )
            {
                // Get neighbours of this donor
                const labelList& nbrs = cCells[donorCell_];

                // Loop through neighbours and mark eligible ones
                forAll (nbrs, nbrI)
                {
                    const label& nbrCellI = nbrs[nbrI];

                    if (eligibleDonors[nbrCellI])
                    {
                        extendedDonorCells_.append(nbrCellI);
                        extendedDonorPoints_.append(cellCentres[nbrCellI]);
                    }
                }

                // Note: no need to shrink DynamicLists for parallel efficiency
                // as IOstream operators first cast into List, thus sending only
                // necessary data (not the capacity of the DynamicList).
            }

            //- Copy extended donors from another candidate donor
            void setExtendedDonors(const donorAcceptor& rd)
            {
                extendedDonorCells_ = rd.extendedDonorCells_;
                extendedDonorPoints_ = rd.extendedDonorPoints_;
            }

            //- Reject donor: used in overlap minimisation algorithms via Donor
            //  Suitability Function
            void rejectDonor()
            {
                donorCell_ = -1;
                donorProcNo_ = -1;
                donorPoint_ = vector::zero;
            }


    // Member Operators

        void operator=(const donorAcceptor& rd)
        {
            acceptorCell_ = rd.acceptorCell_;
            acceptorProcNo_ = rd.acceptorProcNo_;
            acceptorPoint_ = rd.acceptorPoint_;

            donorCell_ = rd.donorCell_;
            donorProcNo_ = rd.donorProcNo_;
            donorPoint_ = rd.donorPoint_;

            extendedDonorCells_ = rd.extendedDonorCells_;
            extendedDonorPoints_ = rd.extendedDonorPoints_;
        }


    // Friend Operators

        friend bool operator==
        (
            const donorAcceptor& a,
            const donorAcceptor& b
        )
        {
            return
                a.acceptorCell_ == b.acceptorCell_
             && a.acceptorProcNo_ == b.acceptorProcNo_
             && a.acceptorPoint_ == b.acceptorPoint_
             && a.donorCell_ == b.donorCell_
             && a.donorProcNo_ == b.donorProcNo_;

             // Note: do not check whether extended neighbours are the same, we
             // assume they will be if donor data is the same
        }

        friend bool operator!=
        (
            const donorAcceptor& a,
            const donorAcceptor& b
        )
        {
            return !(a == b);
        }


    // IOstream Operators

        //- Istream operator
        friend Istream& operator>>(Istream& is, donorAcceptor& rd)
        {
            is  >> rd.acceptorCell_
                >> rd.acceptorProcNo_
                >> rd.acceptorPoint_
                >> rd.donorCell_
                >> rd.donorProcNo_
                >> rd.donorPoint_
                >> rd.extendedDonorCells_
                >> rd.extendedDonorPoints_;

            return is;
        }

        //- Ostream operator
        friend Ostream& operator<<(Ostream& os, const donorAcceptor& rd)
        {
            os  << rd.acceptorCell_ << token::SPACE
                << rd.acceptorProcNo_ << token::SPACE
                << rd.acceptorPoint_ << token::SPACE
                << rd.donorCell_ << token::SPACE
                << rd.donorProcNo_ << token::SPACE
                << rd.donorPoint_ << token::SPACE
                << rd.extendedDonorCells_ << token::SPACE
                << rd.extendedDonorPoints_;

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
