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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline label extrudeLayer::addressingCalculator::origFaceLabel
(
    const label extrudedI
) const
{
    return extrudedFaces_[extrudedI].second();
}

inline label extrudeLayer::addressingCalculator::positionInFace
(
    const label extrudedI,
    const label pointI
) const
{
    const face& f = faces_[extrudedFaces_[extrudedI].first()];
    
    return f.which(pointI);
}
                
inline label extrudeLayer::addressingCalculator::origPointLabel
(
    const label extrudedI,
    const label pos
) const
{
    const face& of = faces_[extrudedFaces_[extrudedI].second()];
    
    if( pairOrientation_[extrudedI] )
    {
        return of[pos];
    }
    else
    {
        return of[(of.size()-pos)%of.size()];
    }
    
    FatalErrorIn
    (
        "label extrudeLayer::addressingCalculator::origPointLabel"
        "(const label, const label) const"
    ) << "Cannot find point for the given position" << abort(FatalError);
    
    return -1;
}

inline label extrudeLayer::addressingCalculator::origPoint
(
    const label extrudedI,
    const label pointI
) const
{
    const face& f = faces_[extrudedFaces_[extrudedI].first()];
    const face& of = faces_[extrudedFaces_[extrudedI].second()];
    const label pos = f.which(pointI);
    
    if( pairOrientation_[extrudedI] )
    {
        return of[pos];
    }
    else
    {
        return of[(of.size()-pos)%of.size()];
    }
    
    FatalErrorIn
    (
        "label extrudeLayer::addressingCalculator::origPoint"
        "(const label, const label) const"
    ) << "Cannot find point for the given position" << abort(FatalError);
    
    return -1;
}

inline label extrudeLayer::addressingCalculator::faceSharingEdge
(
    const label extrudedI,
    const label eI
) const
{
    const face& f = faces_[extrudedFaces_[extrudedI].first()];
    
    const label pointI = f[eI];
    const label nextI = f.nextLabel(eI);
    
    label otherFace(-1);
    forAllRow(pointExtruded_, pointI, pfI)
    {
        const label currFaceI = pointExtruded_(pointI, pfI);
        
        if( currFaceI == extrudedI )
            continue;
        
        if( pointExtruded_.contains(nextI, currFaceI) )
        {
            if( otherFace != -1 )
                FatalErrorIn
                (
                    "label extrudeLayer::addressingCalculator::faceSharingEdge"
                    "(const label, const label) const"
                ) << "Expected only one such face"
                    << abort(FatalError);
            
            otherFace = currFaceI;
        }
    }
    
    return otherFace;
}

inline void extrudeLayer::addressingCalculator::facesSharingEdge
(
    const label start,
    const label end,
    DynList<label>& edgeFaces
) const
{
    edgeFaces.clear();
    
    forAllRow(pointExtruded_, start, pfI)
    {
        const label currFaceI = pointExtruded_(start, pfI);
        
        if( pointExtruded_.contains(end, currFaceI) )
            edgeFaces.append(currFaceI);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
