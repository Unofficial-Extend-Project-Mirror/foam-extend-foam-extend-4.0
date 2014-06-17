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

#include "triSurf.H"
#include "demandDrivenData.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"

#include "gzstream.h"

#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurf::readFromFTR(const fileName& fName)
{
    IFstream fStream(fName);

    fStream >> triSurfFacets::patches_;

    fStream >> triSurfPoints::points_;

    fStream >> triSurfFacets::triangles_;
}

void triSurf::writeToFTR(const fileName& fName) const
{
    OFstream fStream(fName);

    fStream << triSurfFacets::patches_;

    fStream << nl;

    fStream << triSurfPoints::points_;

    fStream << nl;

    fStream << triSurfFacets::triangles_;
}

void triSurf::readFromFMS(const fileName& fName)
{
    IFstream fStream(fName);

    //- read the list of patches defined on the surface mesh
    fStream >> triSurfFacets::patches_;

    //- read points
    fStream >> triSurfPoints::points_;

    //- read surface triangles
    fStream >> triSurfFacets::triangles_;

    //- read feature edges
    fStream >> triSurfFeatureEdges::featureEdges_;

    List<meshSubset> subsets;

    //- read point subsets
    fStream >> subsets;
    forAll(subsets, subsetI)
        triSurfPoints::pointSubsets_.insert(subsetI, subsets[subsetI]);

    subsets.clear();

    //- read facet subsets
    fStream >> subsets;
    forAll(subsets, subsetI)
        triSurfFacets::facetSubsets_.insert(subsetI, subsets[subsetI]);

    subsets.clear();

    //- read subsets on feature edges
    fStream >> subsets;
    forAll(subsets, subsetI)
        triSurfFeatureEdges::featureEdgeSubsets_.insert
        (
            subsetI,
            subsets[subsetI]
        );
}

void triSurf::writeToFMS(const fileName& fName) const
{
    OFstream fStream(fName);

    //- write patches
    fStream << triSurfFacets::patches_;

    fStream << nl;

    //- write points
    fStream << triSurfPoints::points_;

    fStream << nl;

    //- write triangles
    fStream << triSurfFacets::triangles_;

    fStream << nl;

    //- write feature edges
    fStream << triSurfFeatureEdges::featureEdges_;

    fStream << nl;

    //- write point subsets
    List<meshSubset> subsets;
    label i(0);
    subsets.setSize(pointSubsets_.size());
    forAllConstIter(Map<meshSubset>, pointSubsets_, it)
        subsets[i++] = it();
    fStream << subsets;

    fStream << nl;

    //- write subsets of facets
    subsets.setSize(triSurfFacets::facetSubsets_.size());
    i = 0;
    forAllConstIter(Map<meshSubset>, triSurfFacets::facetSubsets_, it)
        subsets[i++] = it();
    fStream << subsets;

    fStream << nl;

    //- write subets of feature edges
    subsets.setSize(triSurfFeatureEdges::featureEdgeSubsets_.size());
    i = 0;
    forAllConstIter
    (
        Map<meshSubset>,
        triSurfFeatureEdges::featureEdgeSubsets_,
        it
    )
        subsets[i++] = it();
    fStream << subsets;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurf::triSurf()
:
    triSurfPoints(),
    triSurfFacets(),
    triSurfFeatureEdges(),
    triSurfAddressing(triSurfPoints::points_, triSurfFacets::triangles_)
{}

//- Construct from parts
triSurf::triSurf
(
    const LongList<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    const edgeLongList& featureEdges,
    const pointField& points
)
:
    triSurfPoints(points),
    triSurfFacets(triangles, patches),
    triSurfFeatureEdges(featureEdges),
    triSurfAddressing(triSurfPoints::points_, triSurfFacets::triangles_)
{}

//- Read from file
triSurf::triSurf(const fileName& fName)
:
    triSurfPoints(),
    triSurfFacets(),
    triSurfFeatureEdges(),
    triSurfAddressing(triSurfPoints::points_, triSurfFacets::triangles_)
{
    readSurface(fName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triSurf::~triSurf()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurf::readSurface(const fileName& fName)
{
    if( fName.ext() == "fms" || fName.ext() == "FMS" )
    {
        readFromFMS(fName);
    }
    else if( fName.ext() == "ftr" || fName.ext() == "FTR" )
    {
        readFromFTR(fName);
    }
    else
    {
        triSurface copySurface(fName);

        //- copy the points
        triSurfPoints::points_.setSize(copySurface.points().size());
        forAll(copySurface.points(), pI)
            triSurfPoints::points_[pI] = copySurface.points()[pI];

        //- copy the triangles
        triSurfFacets::triangles_.setSize(copySurface.size());
        forAll(copySurface, tI)
            triSurfFacets::triangles_[tI] = copySurface[tI];

        //- copy patches
        triSurfFacets::patches_ = copySurface.patches();
    }
}

void triSurf::writeSurface(const fileName& fName) const
{
    if( fName.ext() == "fms" || fName.ext() == "FMS" )
    {
        writeToFMS(fName);
    }
    else if( fName.ext() == "ftr" || fName.ext() == "FTR" )
    {
        writeToFTR(fName);
    }
    else
    {
        const pointField& pts = this->points();
        const LongList<labelledTri>& facets = this->facets();
        const geometricSurfacePatchList& patches = this->patches();

        List<labelledTri> newTrias(facets.size());
        forAll(facets, tI)
            newTrias[tI] = facets[tI];

        triSurface newSurf(newTrias, patches, pts);
        newSurf.write(fName);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
