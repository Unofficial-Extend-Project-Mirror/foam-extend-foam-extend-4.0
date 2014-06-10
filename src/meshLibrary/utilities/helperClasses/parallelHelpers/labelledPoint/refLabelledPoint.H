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
    refLabelledPoint

Description
    A class containing the label of the object it is associated to and
    a labelledPoint

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef refLabelledPoint_H
#define refLabelledPoint_H

#include "labelledPoint.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class refLabelledPoint Declaration
\*---------------------------------------------------------------------------*/
    
class refLabelledPoint
{
    // Private data
        //- label of the object it is associated to
        label objectLabel_;
    
        //- point to be transferred
        labelledPoint p_;
    
    public:
        
    // Constructors
        //- Null construct
        refLabelledPoint()
        :
            objectLabel_(-1),
            p_()
        {}
        
        //- Construct from label and labelledPoint
        refLabelledPoint(const label pl, const labelledPoint& p)
        :
            objectLabel_(pl),
            p_(p)
        {}
            
    // Destructor
        ~refLabelledPoint()
        {}
            
    // Member functions
        //- return label of the object it is associated to
        inline label objectLabel() const
        {
            return objectLabel_;
        }
        
        //- return labelledPoint
        inline const labelledPoint& lPoint() const
        {
            return p_;
        }
        
    // Member operators
        
        inline void operator=(const refLabelledPoint& lp)
        {
            objectLabel_ = lp.objectLabel_;
            p_ = lp.p_;
        }
        
        inline bool operator==(const refLabelledPoint& lp) const
        {
            if( objectLabel_ == lp.objectLabel_ )
                return true;
            
            return false;
        }
        
        inline bool operator!=(const refLabelledPoint& lp) const
        {
            return !this->operator==(lp);
        }
        
    // Friend operators
        friend Ostream& operator<<(Ostream& os, const refLabelledPoint& lp)
        {
            os << token::BEGIN_LIST;
            os << lp.objectLabel_ << token::SPACE;
            os << lp.p_ << token::END_LIST;

            // Check state of Ostream
            os.check("operator<<(Ostream&, const refLabelledPoint&");

            return os;
        }
        
        friend Istream& operator>>(Istream& is, refLabelledPoint& lp)
        {
            // Read beginning of refLabelledPoint
            is.readBegin("refLabelledPoint");
        
            is >> lp.objectLabel_;
            is >> lp.p_;
            
            // Read end of refLabelledPoint
            is.readEnd("refLabelledPoint");
        
            // Check state of Istream
            is.check("operator>>(Istream&, refLabelledPoint");
            
            return is;
        }
};

//- Specify data associated with refLabelledPoint type is contiguous
template<>
inline bool contiguous<refLabelledPoint>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
