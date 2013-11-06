/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Interpolate point field
template<class Type>
tmp<Field<Type> >
coupleMap::pointInterpolate
(
    const Map<label>& mPointMap,
    const Map<label>& sPointMap,
    const Field<Type>& pf,
    bool reverse
) const
{
    const Map<label>& pointMap = entityMap(coupleMap::POINT);

    if (pf.size() != pointMap.size())
    {
        FatalErrorIn
        (
            "coupleMap::pointInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << pointMap.size() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            pointMap.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    if (reverse)
    {
        forAllConstIter(Map<label>, pointMap, pIter)
        {
            result[mPointMap[pIter.key()]] = pf[sPointMap[pIter()]];
        }
    }
    else
    {
        forAllConstIter(Map<label>, pointMap, pIter)
        {
            result[sPointMap[pIter()]] = pf[mPointMap[pIter.key()]];
        }
    }

    return tresult;
}


template<class Type>
tmp<Field<Type> >
coupleMap::pointInterpolate
(
    const Map<label>& mPointMap,
    const Map<label>& sPointMap,
    const tmp<Field<Type> >& tpf,
    bool reverse
) const
{
    tmp<Field<Type> > tint =
    (
        pointInterpolate<Type>(mPointMap, sPointMap, tpf(), reverse)
    );
    tpf.clear();
    return tint;
}


//- Interpolate point field
template<class Type>
tmp<Field<Type> >
coupleMap::faceInterpolate
(
    const label mStart,
    const label sStart,
    const Field<Type>& pf,
    bool reverse
) const
{
    const Map<label>& faceMap = entityMap(coupleMap::FACE);

    if (pf.size() != faceMap.size())
    {
        FatalErrorIn
        (
            "coupleMap::faceInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << faceMap.size() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            faceMap.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    if (reverse)
    {
        forAllConstIter(Map<label>, faceMap, fIter)
        {
            result[fIter.key() - mStart] = pf[fIter() - sStart];
        }
    }
    else
    {
        forAllConstIter(Map<label>, faceMap, fIter)
        {
            result[fIter() - sStart] = pf[fIter.key() - mStart];
        }
    }

    return tresult;
}


template<class Type>
tmp<Field<Type> >
coupleMap::faceInterpolate
(
    const label mStart,
    const label sStart,
    const tmp<Field<Type> >& tpf,
    bool reverse
) const
{
    tmp<Field<Type> > tint =
    (
        faceInterpolate<Type>(mStart, sStart, tpf(), reverse)
    );
    tpf.clear();
    return tint;
}


inline label& coupleMap::patchIndex() const
{
    return patchIndex_;
}


inline label coupleMap::masterIndex() const
{
    return masterIndex_;
}


inline label coupleMap::slaveIndex() const
{
    return slaveIndex_;
}


inline bool coupleMap::isTwoDMesh() const
{
    return twoDMesh_;
}


inline bool coupleMap::isLocal() const
{
    return isLocal_;
}


inline bool coupleMap::isProcessor() const
{
    return !isLocal_;
}


inline bool coupleMap::isSend() const
{
    return isSend_;
}


inline bool coupleMap::isRecv() const
{
    return !isSend_;
}


inline FixedList<label,10>& coupleMap::nEntities() const
{
    return nEntities_;
}


inline label& coupleMap::nEntities(const label eType) const
{
    return nEntities_[eType];
}


inline Map<label>& coupleMap::entityMap(const label eType) const
{
    return entityMap_[eType];
}


inline Map<label>& coupleMap::reverseEntityMap(const label eType) const
{
    return reverseEntityMap_[eType];
}


inline FixedList<labelList,12>& coupleMap::entityBuffer() const
{
    return entityBuffer_;
}


inline labelList& coupleMap::entityBuffer(const label eType) const
{
    return entityBuffer_[eType];
}


inline labelList& coupleMap::entityIndices() const
{
    return entityIndices_;
}


inline List<coupleMap::opType>& coupleMap::entityOperations() const
{
    return entityOperations_;
}


inline pointField& coupleMap::moveNewPoints() const
{
    return moveNewPoints_;
}


inline pointField& coupleMap::moveOldPoints() const
{
    return moveOldPoints_;
}


} // End namespace Foam

// ************************************************************************* //
