/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    polyMeshGen2DEngine

Description
    A simple engine which provides topological information of a 2D mesh
    and allows for maintaining consistency

SourceFiles
    polyMeshGen2DEngine.C

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGen2DEngine.H"
#include "boolList.H"
#include "labelList.H"
#include "boundBox.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const boolList& polyMeshGen2DEngine::activeFace() const
{
    if( !activeFacePtr_ )
        findActiveFaces();

    return *activeFacePtr_;
}

inline const labelList& polyMeshGen2DEngine::activeFaceLabels() const
{
    if( !activeFaceLabelsPtr_ )
        findActiveFaceLabels();

    return *activeFaceLabelsPtr_;
}

inline const boolList& polyMeshGen2DEngine::zMinPoints() const
{
    if( !zMinPointPtr_ )
        findZMinPoints();

    return *zMinPointPtr_;
}

inline const labelList& polyMeshGen2DEngine::zMinPointLabels() const
{
    if( !zMinPointLabelsPtr_ )
        findZMinPointLabels();

    return *zMinPointLabelsPtr_;
}

inline const labelList& polyMeshGen2DEngine::zMinToZMax() const
{
    if( !zMinToZMaxPtr_ )
        findZMinOffsetPoints();

    return *zMinToZMaxPtr_;
}

inline const boolList& polyMeshGen2DEngine::zMaxPoints() const
{
    if( !zMaxPointPtr_ )
        findZMaxPoints();

    return *zMaxPointPtr_;
}

inline const labelList& polyMeshGen2DEngine::zMaxPointLabels() const
{
    if( !zMaxPointLabelsPtr_ )
        findZMaxPointLabels();

    return *zMaxPointLabelsPtr_;
}

inline const labelList& polyMeshGen2DEngine::zMaxToZMin() const
{
    if( !zMaxToZMinPtr_ )
        findZMaxOffsetPoints();

    return *zMaxToZMinPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
