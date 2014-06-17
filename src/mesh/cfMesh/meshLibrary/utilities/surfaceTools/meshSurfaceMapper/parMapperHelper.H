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
    parMapperHelper

Description
    A class containing point coordinates, moving distance,
    global label of the mapped point and projected patch. It is used for
    exchanging data over processors

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef parMapperHelper_H
#define parMapperHelper_H

#include "label.H"
#include "point.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class parMapperHelper Declaration
\*---------------------------------------------------------------------------*/

class parMapperHelper
{
    // Private data
        point p_;
        scalar movingDist_;
        label globalLabel_;
        label patch_;
    
    public:
        
    // Constructors
        explicit inline parMapperHelper()
        :
            p_(),
            movingDist_(),
            globalLabel_(-1),
            patch_(-1)
        {};
        
        explicit inline parMapperHelper
        (
            const point& p,
            const scalar& movingDst,
            const label globalLabel,
            const label patchLabel
        )
        :
            p_(p),
            movingDist_(movingDst),
            globalLabel_(globalLabel),
            patch_(patchLabel)
        {}
    
    // Destructor
        ~parMapperHelper()
        {};

    // Member functions
        //- return point coordinates
        inline const point& coordinates() const
        {
            return p_;
        }
        
        //- return moving distance
        inline const scalar& movingDistance() const
        {
            return movingDist_;
        }
        
        //- return global label
        inline label globalLabel() const
        {
            return globalLabel_;
        }
        
        //- return patch the point is mapped to
        inline label pointPatch() const
        {
            return patch_;
        }
        
    // Operators
        
        inline void operator=(const parMapperHelper& h)
        {
            p_ = h.p_;
            movingDist_ = h.movingDist_;
            globalLabel_ = h.globalLabel_;
            patch_ = h.patch_;
        }
        
        inline bool operator!=(const parMapperHelper& h) const
        {
            if( globalLabel_ != h.globalLabel_ )
                return true;
            
            return false;
        }
        
    // Friend operators
        friend Ostream& operator<<(Ostream& os, const parMapperHelper& h)
        {
            os << token::BEGIN_LIST;
            os << h.p_ << token::SPACE;
            os << h.movingDist_ << token::SPACE;
            os << h.globalLabel_ << token::SPACE;
            os << h.patch_ << token::END_LIST;

            // Check state of Ostream
            os.check("operator<<(Ostream&, const parMapperHelper&");
            return os;
        }
        
        friend Istream& operator>>(Istream& is, parMapperHelper& h)
        {
            // Read beginning of parMapperHelper
            is.readBegin("parMapperHelper");
        
            is >> h.p_;
            is >> h.movingDist_;
            is >> h.globalLabel_;
            is >> h.patch_;
            
            // Read end of parMapperHelper
            is.readEnd("parMapperHelper");
        
            // Check state of Istream
            is.check("operator>>(Istream&, parMapperHelper");
            
            return is;
        }
};

//- Specify data associated with parMapperHelper type is contiguous
template<>
inline bool contiguous<parMapperHelper>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
