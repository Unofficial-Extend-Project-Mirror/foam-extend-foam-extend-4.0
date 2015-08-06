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

Description


\*---------------------------------------------------------------------------*/

#include "meshOctreeCubeCoordinates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline meshOctreeCubeCoordinates::meshOctreeCubeCoordinates()
:
    posX_(-1),
    posY_(-1),
    posZ_(-1),
    level_(0)
{}

inline meshOctreeCubeCoordinates::meshOctreeCubeCoordinates
(
    const label posX,
    const label posY,
    const label posZ,
    const direction level
)
:
    posX_(posX),
    posY_(posY),
    posZ_(posZ),
    level_(level)
{}

inline meshOctreeCubeCoordinates::meshOctreeCubeCoordinates
(
    const meshOctreeCubeCoordinates& cc
)
:
    posX_(cc.posX_),
    posY_(cc.posY_),
    posZ_(cc.posZ_),
    level_(cc.level_)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member functions

inline direction meshOctreeCubeCoordinates::level() const
{
    return level_;
}

inline label meshOctreeCubeCoordinates::posX() const
{
    return posX_;
}

inline label meshOctreeCubeCoordinates::posY() const
{
    return posY_;
}

inline label meshOctreeCubeCoordinates::posZ() const
{
    return posZ_;
}

inline meshOctreeCubeCoordinates meshOctreeCubeCoordinates::refineForPosition
(
    const label i
) const
{
    //- create new boxes in z-order fashion
    FixedList<label, 3> addPos(0);
    if( i & 1 )
        addPos[0] = 1;
    if( i & 2 )
        addPos[1] = 1;
    if( i & 4 )
        addPos[2] = 1;

    return meshOctreeCubeCoordinates
    (
        2*posX_+addPos[0],
        2*posY_+addPos[1],
        posZ_<0?posZ_:2*posZ_+addPos[2],
        level_+1
    );
}

inline meshOctreeCubeCoordinates meshOctreeCubeCoordinates::reduceLevelBy
(
    const direction diff
) const
{
    const label divider = (1 << diff);

    return meshOctreeCubeCoordinates
    (
        posX_ / divider,
        posY_ / divider,
        posZ_<0?posZ_:posZ_ / divider,
        level_ - diff
    );
}

inline meshOctreeCubeCoordinates meshOctreeCubeCoordinates::reduceToLevel
(
    const direction l
) const
{
    const direction diff = level_ - l;

    return reduceLevelBy(diff);
}

inline meshOctreeCubeCoordinates meshOctreeCubeCoordinates::increaseToLevelMin
(
    const direction l
) const
{
    const direction diff = l - level_;
    return meshOctreeCubeCoordinates
    (
        posX_ << diff,
        posY_ << diff,
        posZ_<0?posZ_:posZ_ << diff,
        l
    );
}

inline meshOctreeCubeCoordinates meshOctreeCubeCoordinates::increaseToLevelMax
(
    const direction l
) const
{
    const direction diff = l - level_;
    return meshOctreeCubeCoordinates
    (
        ((posX_ + 1) << diff) - 1,
        ((posY_ + 1) << diff) - 1,
        posZ_<0?posZ_:((posZ_ + 1) << diff) - 1,
        l
    );
}

inline void meshOctreeCubeCoordinates::cubeBox
(
    const boundBox& rootBox,
    point& min,
    point& max
) const
{
    vector dc = (rootBox.max() - rootBox.min());
    const label shift = 1 << level_;

    dc.x() /= shift;
    dc.y() /= shift;
    if( posZ_>=0 )
        dc.z() /= shift;

    min.x() = rootBox.min().x() + dc.x() * this->posX();
    min.y() = rootBox.min().y() + dc.y() * this->posY();
    if( posZ_ >= 0 )
    {
        min.z() = rootBox.min().z() + dc.z() * this->posZ();
    }
    else
    {
        min.z() = rootBox.min().z();
    }

    max = min + dc;
}

inline point meshOctreeCubeCoordinates::centre
(
    const boundBox& rootBox
) const
{
    point min, max;
    cubeBox(rootBox, min, max);
    return (max + min) / 2.0;
}

inline scalar meshOctreeCubeCoordinates::size
(
    const boundBox& rootBox
) const
{
    point min, max;
    cubeBox(rootBox, min, max);
    return (max.x() - min.x());
}

inline void meshOctreeCubeCoordinates::neighbourRange
(
    meshOctreeCubeCoordinates& minCoord,
    meshOctreeCubeCoordinates& maxCoord
) const
{
    minCoord = this->refineForPosition(0);
    if( minCoord.posX_ > 0 )
        minCoord.posX_ -= 1;
    if( minCoord.posY_ > 0 )
        minCoord.posY_ -= 1;
    if( minCoord.posZ_ > 0 )
        minCoord.posZ_ -= 1;

    maxCoord = this->refineForPosition(7);
    const label maxc = (1 << maxCoord.level_) - 1;
    if( maxCoord.posX_ < maxc )
        maxCoord.posX_ += 1;
    if( maxCoord.posY_ < maxc )
        maxCoord.posY_ += 1;
    if( maxCoord.posZ_ < maxc )
        maxCoord.posZ_ += 1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Operators

inline meshOctreeCubeCoordinates meshOctreeCubeCoordinates::operator+
(
    const meshOctreeCubeCoordinates& cc
) const
{
    meshOctreeCubeCoordinates ccc
    (
        posX_ + cc.posX_,
        posY_ + cc.posY_,
        posZ_ + cc.posZ_,
        level_ + cc.level_
    );

    return ccc;
}

inline void meshOctreeCubeCoordinates::operator=
(
    const meshOctreeCubeCoordinates& cc
)
{
    posX_ = cc.posX_;
    posY_ = cc.posY_;
    posZ_ = cc.posZ_;
    level_ = cc.level_;
}

inline bool meshOctreeCubeCoordinates::operator==
(
    const meshOctreeCubeCoordinates& cc
) const
{
    if(
        (level_ != cc.level_) ||
        (posX_ != cc.posX_) ||
        (posY_ != cc.posY_) ||
        (posZ_ != cc.posZ_)
    )
        return false;

    return true;
}

inline bool meshOctreeCubeCoordinates::operator!=
(
    const meshOctreeCubeCoordinates& cc
) const
{
    if(
        (level_ != cc.level_) ||
        (posX_ != cc.posX_) ||
        (posY_ != cc.posY_) ||
        (posZ_ != cc.posZ_)
    )
        return true;

    return false;
}

inline bool meshOctreeCubeCoordinates::operator<=
(
    const meshOctreeCubeCoordinates& cc
) const
{
    const direction minL = Foam::min(level_, cc.level_);
    const meshOctreeCubeCoordinates curr = this->reduceToLevel(minL);
    const meshOctreeCubeCoordinates cctr = cc.reduceToLevel(minL);

    const label resx = curr.posX_ ^ cctr.posX_;
    const label resy = curr.posY_ ^ cctr.posY_;
    const label resz = curr.posZ_ ^ cctr.posZ_;

    label max(resx), dir(0);
    if( (max <= resy) || ((max ^ resy) < resy) )
    {
        max = resy;
        dir = 1;
    }
    if( (max <= resz) || ((max ^ resz) < resz) )
    {
        max = resz;
        dir = 2;
    }

    switch( dir )
    {
        case 0:
        {
            if( curr.posX_ <= cctr.posX_ )
                return true;

            return false;
        } break;
        case 1:
        {
            if( curr.posY_ <= cctr.posY_ )
                return true;

            return false;
        } break;
        case 2:
        {
            if( curr.posZ_ <= cctr.posZ_ )
                return true;

            return false;
        } break;
    };

    return false;
}

inline bool meshOctreeCubeCoordinates::operator>=
(
    const meshOctreeCubeCoordinates& cc
) const
{
    const direction minL = Foam::min(level_, cc.level_);
    const meshOctreeCubeCoordinates curr = this->reduceToLevel(minL);
    const meshOctreeCubeCoordinates cctr = cc.reduceToLevel(minL);

    const label resx = curr.posX_ ^ cctr.posX_;
    const label resy = curr.posY_ ^ cctr.posY_;
    const label resz = curr.posZ_ ^ cctr.posZ_;

    label max(resx), dir(0);
    if( (max <= resy) || ((max ^ resy) < resy) )
    {
        max = resy;
        dir = 1;
    }
    if( (max <= resz) || ((max ^ resz) < resz) )
    {
        max = resz;
        dir = 2;
    }

    switch( dir )
    {
        case 0:
        {
            if( curr.posX_ >= cctr.posX_ )
                return true;

            return false;
        } break;
        case 1:
        {
            if( curr.posY_ >= cctr.posY_ )
                return true;

            return false;
        } break;
        case 2:
        {
            if( curr.posZ_ >= cctr.posZ_ )
                return true;

            return false;
        } break;
    };

    return false;
}

inline bool meshOctreeCubeCoordinates::operator<
(
    const meshOctreeCubeCoordinates& cc
) const
{
    const direction minL = Foam::min(level_, cc.level_);
    const meshOctreeCubeCoordinates curr = this->reduceToLevel(minL);
    const meshOctreeCubeCoordinates cctr = cc.reduceToLevel(minL);

    const label resx = curr.posX_ ^ cctr.posX_;
    const label resy = curr.posY_ ^ cctr.posY_;
    const label resz = curr.posZ_ ^ cctr.posZ_;

    label max(resx), dir(0);
    if( (max <= resy) || ((max ^ resy) < resy) )
    {
        max = resy;
        dir = 1;
    }
    if( (max <= resz) || ((max ^ resz) < resz) )
    {
        max = resz;
        dir = 2;
    }

    switch( dir )
    {
        case 0:
        {
            if( curr.posX_ < cctr.posX_ )
                return true;

            return false;
        } break;
        case 1:
        {
            if( curr.posY_ < cctr.posY_ )
                return true;

            return false;
        } break;
        case 2:
        {
            if( curr.posZ_ < cctr.posZ_ )
                return true;

            return false;
        } break;
    };

    return false;
}

inline bool meshOctreeCubeCoordinates::operator>
(
    const meshOctreeCubeCoordinates& cc
) const
{
    const direction minL = Foam::min(level_, cc.level_);
    const meshOctreeCubeCoordinates curr = this->reduceToLevel(minL);
    const meshOctreeCubeCoordinates cctr = cc.reduceToLevel(minL);

    const label resx = curr.posX_ ^ cctr.posX_;
    const label resy = curr.posY_ ^ cctr.posY_;
    const label resz = curr.posZ_ ^ cctr.posZ_;

    label max(resx), dir(0);
    if( (max <= resy) || ((max ^ resy) < resy) )
    {
        max = resy;
        dir = 1;
    }
    if( (max <= resz) || ((max ^ resz) < resz) )
    {
        max = resz;
        dir = 2;
    }

    switch( dir )
    {
        case 0:
        {
            if( curr.posX_ > cctr.posX_ )
                return true;

            return false;
        } break;
        case 1:
        {
            if( curr.posY_ > cctr.posY_ )
                return true;

            return false;
        } break;
        case 2:
        {
            if( curr.posZ_ > cctr.posZ_ )
                return true;

            return false;
        } break;
    };

    return false;
}

inline Istream& operator>>
(
    Istream& is,
    meshOctreeCubeCoordinates& cc
)
{
    // Read beginning of meshOctreeCubeCoordinates
    is.readBegin("meshOctreeCubeCoordinates");

    label l;
    is >> l;
    cc.level_ = l;
    is >> cc.posX_;
    is >> cc.posY_;
    is >> cc.posZ_;

    // Read end of meshOctreeCubeCoordinates
    is.readEnd("meshOctreeCubeCoordinates");

    // Check state of Istream
    is.check("operator>>(Istream&, meshOctreeCubeCoordinates");

    return is;
}

inline Ostream& operator<<
(
    Ostream& os,
    const meshOctreeCubeCoordinates& cc
)
{
    os << token::BEGIN_LIST;

    os << label(cc.level_) << token::SPACE;
    os << cc.posX_ << token::SPACE;
    os << cc.posY_ << token::SPACE;
    os << cc.posZ_ << token::END_LIST;

    // Check state of Ostream
    os.check("operator<<(Ostream&, const meshOctreeCubeCoordinates");

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
