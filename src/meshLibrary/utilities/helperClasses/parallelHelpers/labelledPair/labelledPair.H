/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    labelledPair

Description
    A class containing point label and meshOctreeCubeCoordinates.
    It is used for exchanging data over processors

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef labelledPair_H
#define labelledPair_H

#include "labelPair.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class labelledPair Declaration
\*---------------------------------------------------------------------------*/

class labelledPair
{
    // Private data
        //- label
        label pLabel_;

        //- pair data
        labelPair pair_;

    public:

    // Constructors
        //- Null construct
        labelledPair()
        :
            pLabel_(-1),
            pair_()
        {}

        //- Construct from label and pair
        labelledPair
        (
            const label pl,
            const labelPair& lp
        )
        :
            pLabel_(pl),
            pair_(lp)
        {}

    // Destructor
        ~labelledPair()
        {}

    // Member functions
        //- return pair label
        inline label pairLabel() const
        {
            return pLabel_;
        }

        //- return the pair
        inline const labelPair& pair() const
        {
            return pair_;
        }

    // Member operators

        inline void operator=(const labelledPair& lpp)
        {
            pLabel_ = lpp.pLabel_;
            pair_ = lpp.pair_;
        }

        inline bool operator==
        (
            const labelledPair& lpp
        ) const
        {
            if( pLabel_ != lpp.pLabel_ )
                return false;

            const labelPair& op = lpp.pair_;

            if( (pair_.first() == op.first()) && (pair_.second() == op.second()) )
                return true;
            if( (pair_.first() == op.second()) && (pair_.second() == op.first()) )
                return true;

            return false;
        }

        inline bool operator!=
        (
            const labelledPair& lcc
        ) const
        {
            return !this->operator==(lcc);
        }

        inline bool operator<(const labelledPair& lpp) const
        {
            if( pLabel_ < lpp.pLabel_ )
            {
                return true;
            }
            else if( pLabel_ > lpp.pLabel_ )
            {
                return false;
            }

            if(
                (pair_.first() + pair_.second()) <
                (lpp.pair().first() + lpp.pair().second())
            )
                return true;

            if(
                Foam::min(pair_.first(), pair_.second()) <
                Foam::min(lpp.pair().first(), lpp.pair().second())
            )
                return true;

            return false;
        }

    // Friend operators
        friend Ostream& operator<<
        (
            Ostream& os,
            const labelledPair& lpp
        )
        {
            os << token::BEGIN_LIST;
            os << lpp.pLabel_ << token::SPACE;
            os << lpp.pair_ << token::END_LIST;

            // Check state of Ostream
            os.check
            (
                "operator<<(Ostream&, const labelledPair&"
            );

            return os;
        }

        friend Istream& operator>>
        (
            Istream& is,
            labelledPair& lpp
        )
        {
            // Read beginning of labelledPair
            is.readBegin("labelledPair");

            is >> lpp.pLabel_;
            is >> lpp.pair_;

            // Read end of labelledPair
            is.readEnd("labelledPair");

            // Check state of Istream
            is.check("operator>>(Istream&, labelledPair");

            return is;
        }
};

//- Specify data associated with labelledPair
//- type is contiguous
template<>
inline bool contiguous<labelledPair>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
