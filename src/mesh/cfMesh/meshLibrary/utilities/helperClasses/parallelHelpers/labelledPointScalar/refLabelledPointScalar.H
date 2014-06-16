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
    refLabelledPointScalar

Description
    A class containing a label and labelledPointScalar. It is used for
    exchanging data over processors

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef refLabelledPointScalar_H
#define refLabelledPointScalar_H

#include "labelledPointScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class refLabelledPointScalar Declaration
\*---------------------------------------------------------------------------*/

class refLabelledPointScalar
{
    // Private data
        //- point label
        label pLabel_;

        //- labelledPointScalar
        labelledPointScalar lps_;

    public:

    // Constructors
        //- Null construct
        refLabelledPointScalar()
        :
            pLabel_(-1),
            lps_()
        {}

        //- Construct from label and labelledPointScalar
        refLabelledPointScalar(const label pl, const labelledPointScalar& lps)
        :
            pLabel_(pl),
            lps_(lps)
        {}

    // Destructor
        ~refLabelledPointScalar()
        {}

    // Member functions
        //- return object label
        inline label objectLabel() const
        {
            return pLabel_;
        }

        //- return labelledPointScalar
        inline const labelledPointScalar& lps() const
        {
            return lps_;
        }

    // Member operators

        inline void operator=(const refLabelledPointScalar& lps)
        {
            pLabel_ = lps.pLabel_;
            lps_ = lps.lps_;
        }

        inline bool operator==(const refLabelledPointScalar& lps) const
        {
            if( pLabel_ == lps.pLabel_ )
                return true;

            return false;
        }

        inline bool operator!=(const refLabelledPointScalar& lps) const
        {
            return !this->operator==(lps);
        }

    // Friend operators
        friend Ostream& operator<<(Ostream& os, const refLabelledPointScalar& lps)
        {
            os << token::BEGIN_LIST;
            os << lps.pLabel_ << token::SPACE;
            os << lps.lps_ << token::END_LIST;

            // Check state of Ostream
            os.check("operator<<(Ostream&, const labelledPointScalarS&");

            return os;
        }

        friend Istream& operator>>(Istream& is, refLabelledPointScalar& lps)
        {
            // Read beginning of refLabelledPointScalar
            is.readBegin("refLabelledPointScalar");

            is >> lps.pLabel_;
            is >> lps.lps_;

            // Read end of refLabelledPointScalar
            is.readEnd("refLabelledPointScalar");

            // Check state of Istream
            is.check("operator>>(Istream&, refLabelledPointScalar");

            return is;
        }
};

//- Specify data associated with refLabelledPointScalar type is contiguous
template<>
inline bool contiguous<refLabelledPointScalar>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
