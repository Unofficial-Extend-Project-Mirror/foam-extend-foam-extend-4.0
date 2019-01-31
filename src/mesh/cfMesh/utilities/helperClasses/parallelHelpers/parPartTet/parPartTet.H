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
    parPartTet

Description
    Holds labels and coordinates of points making a tet.
    It is used for exchanging data over processors

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef parPartTet_H
#define parPartTet_H

#include "labelledPoint.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class parPartTet Declaration
\*---------------------------------------------------------------------------*/
        
class parPartTet
{
    // Private data
        labelledPoint pts_[4];
    
public:
        
    // Constructors
    
        inline parPartTet()
        {}
        
        explicit inline parPartTet
        (
            const labelledPoint& p0,
            const labelledPoint& p1,
            const labelledPoint& p2,
            const labelledPoint& p3
        )
        {
            pts_[0] = p0;
            pts_[1] = p1;
            pts_[2] = p2;
            pts_[3] = p3;
        }
        
    // Destructor
        
        ~parPartTet()
        {}
            
    // Member functions
        
    // Member operators
            
        inline const labelledPoint& operator[](const label i) const
        {
            return pts_[i];
        }
        
        inline bool operator !=(const parPartTet& /*ptf*/) const
        {
            Serr << "Not implemented" << endl;
            ::exit(1);
            
            return true;
        }
        
    // Friend operators
        
        inline friend Ostream& operator<<(Ostream& os, const parPartTet& ppt)
        {
            os << token::BEGIN_LIST;
            os << ppt.pts_[0] << token::SPACE;
            os << ppt.pts_[1] << token::SPACE;
            os << ppt.pts_[2] << token::SPACE;
            os << ppt.pts_[3];
            os << token::END_LIST;

            // Check state of Ostream
            os.check("operator<<(Ostream&, const parPartTet&");
            return os;
        }
        
        inline friend Istream& operator>>(Istream& is, parPartTet& ppt)
        {
            // Read beginning of parPartTet
            is.readBegin("parPartTet");
        
            for(label i=0;i<4;++i)
                is >> ppt.pts_[i];
            
            // Read end of parHelper
            is.readEnd("parPartTet");
        
            // Check state of Istream
            is.check("operator>>(Istream&, parPartTet");
            
            return is;
        }
};

//- Specify data associated with parPartTet type is contiguous
template<>
inline bool contiguous<parPartTet>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
