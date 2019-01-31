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
    labelledPointScalar

Description
    A class containing a label, coordinates and scalar. It is used for
    exchanging data over processors

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef labelledPointScalar_H
#define labelledPointScalar_H

#include "label.H"
#include "point.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class labelledPointScalar Declaration
\*---------------------------------------------------------------------------*/
    
class labelledPointScalar
{
    // Private data
        //- point label
        label pLabel_;
    
        //- point coordinates
        point coords_;
    
        //- scalar data
        scalar weight_;
    
    public:
        
    // Constructors
        //- Null construct
        labelledPointScalar()
        :
            pLabel_(-1),
            coords_(vector::zero),
            weight_(0.0)
        {}
        
        //- Construct from point and label
        labelledPointScalar(const label pl, const point& p, const scalar s)
        :
            pLabel_(pl),
            coords_(p),
            weight_(s)
        {}
            
    // Destructor
        ~labelledPointScalar()
        {}
            
    // Member functions
        //- return point label
        inline label pointLabel() const
        {
            return pLabel_;
        }
        
        inline label& pointLabel()
        {
            return pLabel_;
        }
        
        //- return point coordinates
        inline const point& coordinates() const
        {
            return coords_;
        }
        
        inline point& coordinates()
        {
            return coords_;
        }
        
        //- return scalar value
        inline const scalar& scalarValue() const
        {
            return weight_;
        }
        
        inline scalar& scalarValue()
        {
            return weight_;
        }
        
    // Member operators
        
        inline void operator=(const labelledPointScalar& lps)
        {
            pLabel_ = lps.pLabel_;
            coords_ = lps.coords_;
            weight_ = lps.weight_;
        }
        
        inline bool operator==(const labelledPointScalar& lps) const
        {
            if( pLabel_ == lps.pLabel_ )
                return true;
            
            return false;
        }
        
        inline bool operator!=(const labelledPointScalar& lps) const
        {
            return !this->operator==(lps);
        }
        
    // Friend operators
        friend Ostream& operator<<(Ostream& os, const labelledPointScalar& lps)
        {
            os << token::BEGIN_LIST;
            os << lps.pLabel_ << token::SPACE;
            os << lps.coords_ << token::SPACE;
            os << lps.weight_ << token::END_LIST;

            // Check state of Ostream
            os.check("operator<<(Ostream&, const labelledPointScalarS&");

            return os;
        }
        
        friend Istream& operator>>(Istream& is, labelledPointScalar& lps)
        {
            // Read beginning of labelledPointScalar
            is.readBegin("labelledPointScalar");
        
            is >> lps.pLabel_;
            is >> lps.coords_;
            is >> lps.weight_;
            
            // Read end of labelledPointScalar
            is.readEnd("labelledPointScalar");
        
            // Check state of Istream
            is.check("operator>>(Istream&, labelledPointScalar");
            
            return is;
        }
};

//- Specify data associated with labelledPointScalar type is contiguous
template<>
inline bool contiguous<labelledPointScalar>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
