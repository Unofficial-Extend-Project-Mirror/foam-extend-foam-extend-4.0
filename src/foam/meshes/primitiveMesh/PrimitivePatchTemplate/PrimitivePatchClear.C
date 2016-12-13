/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "PrimitivePatchTemplate.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
clearGeom()
{
    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            << "clearGeom() : clearing geometric data"
            << endl;
    }

    deleteDemandDrivenData(localPointsPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(faceNormalsPtr_);
    deleteDemandDrivenData(pointNormalsPtr_);
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
clearTopology()
{
    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            << "clearTopology() : clearing patch addressing"
            << endl;
    }

    // group created and destroyed together
    if (edgesPtr_ && faceFacesPtr_ && edgeFacesPtr_ && faceEdgesPtr_)
    {
        delete edgesPtr_;
        edgesPtr_ = NULL;

        delete faceFacesPtr_;
        faceFacesPtr_ = NULL;

        delete edgeFacesPtr_;
        edgeFacesPtr_ = NULL;

        delete faceEdgesPtr_;
        faceEdgesPtr_ = NULL;
    }

    deleteDemandDrivenData(boundaryPointsPtr_);
    deleteDemandDrivenData(pointEdgesPtr_);
    deleteDemandDrivenData(pointFacesPtr_);
    deleteDemandDrivenData(edgeLoopsPtr_);
    deleteDemandDrivenData(localPointOrderPtr_);
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
clearPatchMeshAddr()
{
    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            << "clearPatchMeshAddr() : "
            << "clearing patch-mesh addressing"
            << endl;
    }

    deleteDemandDrivenData(meshPointsPtr_);
    deleteDemandDrivenData(meshPointMapPtr_);
    deleteDemandDrivenData(localFacesPtr_);
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
clearOut()
{
    clearGeom();
    clearTopology();
    clearPatchMeshAddr();
}


// ************************************************************************* //
