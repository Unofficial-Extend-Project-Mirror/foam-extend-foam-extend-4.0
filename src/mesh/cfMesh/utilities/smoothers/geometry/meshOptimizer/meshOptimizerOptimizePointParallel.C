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

#include "demandDrivenData.H"
#include "meshOptimizer.H"
#include "polyMeshGenAddressing.H"

#include "refLabelledPoint.H"
#include "labelledPointScalar.H"
#include "helperFunctionsPar.H"

#include <map>

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOptimizer::laplaceSmoother::laplacianParallel
(
    const labelLongList& procPoints,
    const bool smoothOnlySurfaceNodes
)
{
    if( !Pstream::parRun() )
        return;

    if( returnReduce(procPoints.size(), sumOp<label>()) == 0 )
        return;

    const polyMeshGenAddressing& addressing = mesh_.addressingData();
    const VRWGraph& pPoints = mesh_.addressingData().pointPoints();
    pointFieldPMG& points = mesh_.points();

    //- exchange data between processors
    const labelLongList& globalPointLabel = addressing.globalPointLabel();
    const VRWGraph& pointAtProcs = addressing.pointAtProcs();
    const Map<label>& globalToLocal = addressing.globalToLocalPointAddressing();

    //- create storage for the data
    std::map<label, labelledPoint> localData;

    //- create data which shall be exchanged with other processors
    std::map<label, LongList<refLabelledPoint> > exchangeData;
    forAll(procPoints, pI)
    {
        const label pointI = procPoints[pI];

        localData.insert(std::make_pair(pointI, labelledPoint(0,vector::zero)));
        labelledPoint& lpd = localData[pointI];

        forAllRow(pPoints, pointI, ppI)
        {
            const label nei = pPoints(pointI, ppI);

            if( smoothOnlySurfaceNodes && (vertexLocation_[nei] & INSIDE) )
                continue;

            if( pointAtProcs.sizeOfRow(nei) != 0 )
            {
                label pMin(Pstream::nProcs());
                forAllRow(pointAtProcs, nei, procI)
                {
                    const label procJ = pointAtProcs(nei, procI);
                    if( (procJ < pMin) && pointAtProcs.contains(pointI, procJ) )
                        pMin = procJ;
                }

                if( pMin != Pstream::myProcNo() )
                    continue;
            }

            ++lpd.pointLabel();
            lpd.coordinates() += points[nei];
        }

        forAllRow(pointAtProcs, pointI, procI)
        {
            const label neiProc = pointAtProcs(pointI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            if( exchangeData.find(neiProc) == exchangeData.end() )
            {
                exchangeData.insert
                (
                    std::make_pair(neiProc, LongList<refLabelledPoint>())
                );
            }

            //- add data to the list which will be sent to other processor
            LongList<refLabelledPoint>& dts = exchangeData[neiProc];
            dts.append(refLabelledPoint(globalPointLabel[pointI], lpd));
        }
    }

    //- exchange data with other processors
    LongList<refLabelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const refLabelledPoint& lp = receivedData[i];
        const label pointI = globalToLocal[lp.objectLabel()];

        labelledPoint& lpd = localData[pointI];

        lpd.pointLabel() += lp.lPoint().pointLabel();
        lpd.coordinates() += lp.lPoint().coordinates();
    }

    //- create new positions of nodes
    forAll(procPoints, pI)
    {
        const label pointI = procPoints[pI];

        const labelledPoint& lp = localData[pointI];

        if( lp.pointLabel() != 0 )
        {
            const point newP = lp.coordinates() / lp.pointLabel();

            points[pointI] = newP;
        }
    }

    # ifdef DEBUGSmooth
    //- check
    std::map<label, LongList<labelledPoint> > check;
    forAll(procPoints, pI)
    {
        const label pointI = procPoints[pI];

        forAllRow(pointAtProcs, pointI, i)
        {
            const label procI = pointAtProcs(pointI, i);
            if( procI == Pstream::myProcNo() )
                continue;
            if( check.find(procI) == check.end() )
            {
                check.insert(std::make_pair(procI, LongList<labelledPoint>()));
            }

            LongList<labelledPoint>& data = check[procI];
            data.append(labelledPoint(globalPointLabel[pointI],points[pointI]));
        }
    }

    LongList<labelledPoint> data;
    help::exchangeMap(check, data);

    forAll(data, i)
    {
        const label pointI = globalToLocal[data[i].pointLabel()];

        if( mag(points[pointI] - data[i].coordinates()) > SMALL )
            Pout << "Crap " << globalPointLabel[pointI] << " coordinates "
                << points[pointI] << " point there " << data[i] << endl;
    }
    # endif
}

void meshOptimizer::laplaceSmoother::laplacianPCParallel
(
    const labelLongList& procPoints
)
{
    if( !Pstream::parRun() )
        return;

    if( returnReduce(procPoints.size(), sumOp<label>()) == 0 )
        return;

    const polyMeshGenAddressing& addressing = mesh_.addressingData();
    const VRWGraph& pCells = mesh_.addressingData().pointCells();
    const vectorField& centres = mesh_.addressingData().cellCentres();
    pointFieldPMG& points = mesh_.points();

    //- exchange data between processors
    const labelLongList& globalPointLabel = addressing.globalPointLabel();
    const VRWGraph& pointAtProcs = addressing.pointAtProcs();
    const Map<label>& globalToLocal = addressing.globalToLocalPointAddressing();

    //- create storage for the data
    std::map<label, labelledPoint> localData;

    //- create data which shall be exchanged with other processors
    std::map<label, LongList<refLabelledPoint> > exchangeData;
    forAll(procPoints, pI)
    {
        const label pointI = procPoints[pI];

        if( vertexLocation_[pointI] & LOCKED )
            continue;

        localData.insert(std::make_pair(pointI, labelledPoint(0,vector::zero)));
        labelledPoint& lpd = localData[pointI];

        forAllRow(pCells, pointI, pcI)
        {
            const label cellI = pCells(pointI, pcI);

            ++lpd.pointLabel();
            lpd.coordinates() += centres[cellI];
        }

        forAllRow(pointAtProcs, pointI, procI)
        {
            const label neiProc = pointAtProcs(pointI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            if( exchangeData.find(neiProc) == exchangeData.end() )
            {
                exchangeData.insert
                (
                    std::make_pair(neiProc, LongList<refLabelledPoint>())
                );
            }

            //- add data to the list which will be sent to other processor
            LongList<refLabelledPoint>& dts = exchangeData[neiProc];
            dts.append(refLabelledPoint(globalPointLabel[pointI], lpd));
        }
    }

    //- exchange data with other processors
    LongList<refLabelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const refLabelledPoint& lp = receivedData[i];
        const label pointI = globalToLocal[lp.objectLabel()];

        labelledPoint& lpd = localData[pointI];

        lpd.pointLabel() += lp.lPoint().pointLabel();
        lpd.coordinates() += lp.lPoint().coordinates();
    }

    //- create new positions of nodes
    forAll(procPoints, pI)
    {
        const label pointI = procPoints[pI];

        const labelledPoint& lp = localData[pointI];

        if( lp.pointLabel() != 0 )
        {
            const point newP = lp.coordinates() / lp.pointLabel();

            points[pointI] = newP;
        }
    }

    # ifdef DEBUGSmooth
    //- check
    std::map<label, LongList<labelledPoint> > check;
    forAll(procPoints, pI)
    {
        const label pointI = procPoints[pI];

        forAllRow(pointAtProcs, pointI, i)
        {
            const label procI = pointAtProcs(pointI, i);
            if( procI == Pstream::myProcNo() )
                continue;
            if( check.find(procI) == check.end() )
            {
                check.insert(std::make_pair(procI, LongList<labelledPoint>()));
            }

            LongList<labelledPoint>& data = check[procI];
            data.append(labelledPoint(globalPointLabel[pointI],points[pointI]));
        }
    }

    LongList<labelledPoint> data;
    help::exchangeMap(check, data);

    forAll(data, i)
    {
        const label pointI = globalToLocal[data[i].pointLabel()];

        if( mag(points[pointI] - data[i].coordinates()) > SMALL )
            Pout << "Crap " << globalPointLabel[pointI] << " coordinates "
                << points[pointI] << " point there " << data[i] << endl;
    }
    # endif
}

void meshOptimizer::laplaceSmoother::laplacianWPCParallel
(
    const labelLongList& procPoints
)
{
    if( !Pstream::parRun() )
        return;

    if( returnReduce(procPoints.size(), sumOp<label>()) == 0 )
        return;

    const polyMeshGenAddressing& addressing = mesh_.addressingData();
    const VRWGraph& pCells = mesh_.addressingData().pointCells();
    const vectorField& centres = mesh_.addressingData().cellCentres();
    const scalarField& volumes = mesh_.addressingData().cellVolumes();
    pointFieldPMG& points = mesh_.points();

    //- exchange data between processors
    const labelLongList& globalPointLabel = addressing.globalPointLabel();
    const VRWGraph& pointAtProcs = addressing.pointAtProcs();
    const Map<label>& globalToLocal = addressing.globalToLocalPointAddressing();

    //- create storage for the data
    std::map<label, labelledPointScalar> localData;

    //- create data which shall be exchanged with other processors
    std::map<label, LongList<labelledPointScalar> > exchangeData;
    forAll(procPoints, pI)
    {
        const label pointI = procPoints[pI];

        if( vertexLocation_[pointI] & LOCKED )
            continue;

        localData.insert
        (
            std::make_pair
            (
                pointI,
                labelledPointScalar(globalPointLabel[pointI], vector::zero, 0.)
            )
        );
        labelledPointScalar& lps = localData[pointI];

        forAllRow(pCells, pointI, pcI)
        {
            const label cellI = pCells(pointI, pcI);

            const scalar w = Foam::max(volumes[cellI], VSMALL);
            lps.coordinates() += w * centres[cellI];
            lps.scalarValue() += w;
        }

        forAllRow(pointAtProcs, pointI, procI)
        {
            const label neiProc = pointAtProcs(pointI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            if( exchangeData.find(neiProc) == exchangeData.end() )
            {
                exchangeData.insert
                (
                    std::make_pair(neiProc, LongList<labelledPointScalar>())
                );
            }

            //- add data to the list which will be sent to other processor
            LongList<labelledPointScalar>& dts = exchangeData[neiProc];
            dts.append(lps);
        }
    }

    //- exchange data with other processors
    LongList<labelledPointScalar> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const labelledPointScalar& lps = receivedData[i];
        const label pointI = globalToLocal[lps.pointLabel()];

        labelledPointScalar& lp = localData[pointI];

        lp.scalarValue() += lps.scalarValue();
        lp.coordinates() += lps.coordinates();
    }

    //- create new positions of nodes
    forAll(procPoints, pI)
    {
        const label pointI = procPoints[pI];

        if( vertexLocation_[pointI] & LOCKED )
            continue;

        const labelledPointScalar& lp = localData[pointI];

        if( lp.pointLabel() != 0 )
        {
            const point newP = lp.coordinates() / lp.scalarValue();

            points[pointI] = newP;
        }
    }

    # ifdef DEBUGSmooth
    //- check
    std::map<label, LongList<labelledPoint> > check;
    forAll(procPoints, pI)
    {
        const label pointI = procPoints[pI];

        forAllRow(pointAtProcs, pointI, i)
        {
            const label procI = pointAtProcs(pointI, i);
            if( procI == Pstream::myProcNo() )
                continue;
            if( check.find(procI) == check.end() )
            {
                check.insert(std::make_pair(procI, LongList<labelledPoint>()));
            }

            LongList<labelledPoint>& data = check[procI];
            data.append(labelledPoint(globalPointLabel[pointI],points[pointI]));
        }
    }

    LongList<labelledPoint> data;
    help::exchangeMap(check, data);

    forAll(data, i)
    {
        const label pointI = globalToLocal[data[i].pointLabel()];

        if( mag(points[pointI] - data[i].coordinates()) > SMALL )
            Pout << "Crap " << globalPointLabel[pointI] << " coordinates "
                << points[pointI] << " point there " << data[i] << endl;
    }
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
