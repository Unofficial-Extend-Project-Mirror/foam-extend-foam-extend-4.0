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
#include "gzstream.h"
#include "triSurface.H"

#include "helperFunctions.H"

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

void triSurf::topologyCheck()
{
    const pointField& pts = this->points();
    const LongList<labelledTri>& trias = this->facets();

    //- check for inf and nan points
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(pts, pointI)
    {
        const point& p = pts[pointI];

        if( help::isnan(p) || help::isinf(p) )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            {
                FatalErrorIn
                (
                    "void triSurf::topologyCheck()"
                ) << "Point " << pointI << " has invalid coordinates "
                  << p << exit(FatalError);
            }
        }
    }

    //- check whether the nodes are within the scope
    //- report duplicate nodes in the same triangle
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(trias, triI)
    {
        const labelledTri& ltri = trias[triI];

        forAll(ltri, pI)
        {
            if( ltri[pI] < 0 || (ltri[pI] >= pts.size()) )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                FatalErrorIn
                (
                    "void triSurf::topologyCheck()"
                ) << "Point " << ltri[pI] << " in triangle " << ltri
                  << " is out of scope 0 " << pts.size() << exit(FatalError);
            }

            if( ltri[pI] == ltri[(pI+1)%3] || ltri[pI] == ltri[(pI+2)%3] )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                WarningIn
                (
                    "void triSurf::topologyCheck()"
                ) << "Triangle " << ltri << " has duplicated points. "
                  << "This may cause problems in the meshing process!" << endl;
            }
        }
    }

    //- check feature edges
    const edgeLongList& featureEdges = this->featureEdges();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(featureEdges, eI)
    {
        const edge& fe = featureEdges[eI];

        forAll(fe, pI)
        {
            if( fe[pI] < 0 || (fe[pI] >= pts.size()) )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                FatalErrorIn
                (
                    "void triSurf::topologyCheck()"
                ) << "Feature edge " << fe << " point " << fe[pI]
                  << " is out of scope 0 " << pts.size() << exit(FatalError);
            }
        }

        if( fe.start() == fe.end() )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            WarningIn
            (
                "void triSurf::topologyCheck()"
            ) << "Feature edge " << fe << " has duplicated points. "
              << "This may cause problems in the meshing process!" << endl;
        }
    }

    //- check point subsets
    DynList<label> subsetIds;
    this->pointSubsetIndices(subsetIds);
    forAll(subsetIds, i)
    {
        labelLongList elmts;
        this->pointsInSubset(subsetIds[i], elmts);

        forAll(elmts, elmtI)
        {
            const label elI = elmts[elmtI];

            if( elI < 0 || elI >= pts.size() )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                FatalErrorIn
                (
                    "void triSurf::topologyCheck()"
                ) << "Point " << elI << " in point subset "
                  << this->pointSubsetName(subsetIds[i])
                  << " is out of scope 0 " << pts.size() << exit(FatalError);
            }
        }
    }

    //- check face subsets
    subsetIds.clear();
    this->facetSubsetIndices(subsetIds);
    forAll(subsetIds, i)
    {
        labelLongList elmts;
        this->facetsInSubset(subsetIds[i], elmts);

        forAll(elmts, elmtI)
        {
            const label elI = elmts[elmtI];

            if( elI < 0 || elI >= trias.size() )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                FatalErrorIn
                (
                    "void triSurf::topologyCheck()"
                ) << "Triangle " << elI << " in facet subset "
                  << this->facetSubsetName(subsetIds[i])
                  << " is out of scope 0 " << trias.size() << exit(FatalError);
            }
        }
    }

    //- check feature edge subsets
    subsetIds.clear();
    this->edgeSubsetIndices(subsetIds);
    forAll(subsetIds, i)
    {
        labelLongList elmts;
        this->edgesInSubset(subsetIds[i], elmts);

        forAll(elmts, elmtI)
        {
            const label elI = elmts[elmtI];

            if( elI < 0 || elI >= featureEdges.size() )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                FatalErrorIn
                (
                    "void triSurf::topologyCheck()"
                ) << "Feature edge " << elI << " in edge subset "
                  << this->edgeSubsetName(subsetIds[i])
                  << " is out of scope 0 " << featureEdges.size()
                  << exit(FatalError);
            }
        }
    }
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
{
    topologyCheck();
}

//- Read from file
triSurf::triSurf(const fileName& fName)
:
    triSurfPoints(),
    triSurfFacets(),
    triSurfFeatureEdges(),
    triSurfAddressing(triSurfPoints::points_, triSurfFacets::triangles_)
{
    readSurface(fName);

    topologyCheck();
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
