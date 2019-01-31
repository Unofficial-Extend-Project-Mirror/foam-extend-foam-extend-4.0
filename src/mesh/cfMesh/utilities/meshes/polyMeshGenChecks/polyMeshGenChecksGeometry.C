/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the origina author
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenChecks.H"
#include "polyMeshGenAddressing.H"
#include "pyramidPointFaceRef.H"
#include "helperFunctions.H"
#include "tetrahedron.H"
#include "syncTools.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace polyMeshGenChecks
{

bool checkClosedBoundary(const polyMeshGen& mesh, const bool report)
{
    // Loop through all boundary faces and sum up the face area vectors.
    // For a closed boundary, this should be zero in all vector components
    vector sumClosed(vector::zero);
    scalar sumMagClosedBoundary = 0;

    const vectorField& areas = mesh.addressingData().faceAreas();

    for(label faceI = mesh.nInternalFaces();faceI<areas.size();++faceI)
    {
        sumClosed += areas[faceI];
        sumMagClosedBoundary += mag(areas[faceI]);
    }

    // Check the openness in the maximum direction (this is APPROXIMATE!)
    scalar maxOpen = max
    (
        sumClosed.component(vector::X),
        max
        (
            sumClosed.component(vector::Y),
            sumClosed.component(vector::Z)
        )
    );

    reduce(sumClosed, sumOp<vector>());
    reduce(maxOpen, maxOp<scalar>());

    if( maxOpen > SMALL*max(1.0, sumMagClosedBoundary) )
    {
        SeriousErrorIn
        (
            "bool checkClosedBoundary(const polyMeshGen&, const bool report)"
        )   << "Possible hole in boundary description" << endl;

        Info<< "Boundary openness in x-direction = "
            << sumClosed.component(vector::X) << endl;

        Info<< "Boundary openness in y-direction = "
            << sumClosed.component(vector::Y) << endl;

        Info<< "Boundary openness in z-direction = "
            << sumClosed.component(vector::Z) << endl;

        return true;
    }
    else
    {
        if( report )
        {
            Info<< "Boundary openness in x-direction = "
                << sumClosed.component(vector::X) << endl;

            Info<< "Boundary openness in y-direction = "
                << sumClosed.component(vector::Y) << endl;

            Info<< "Boundary openness in z-direction = "
                << sumClosed.component(vector::Z) << endl;

            Info<< "Boundary closed (OK)." << endl;
        }

        return false;
    }
}

bool checkClosedCells
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar aspectWarn,
    labelHashSet* setPtr
)
{
    // Check that all cells labels are valid
    const cellListPMG& cells = mesh.cells();
    const label nFaces = mesh.faces().size();

    label nErrorClosed = 0;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided) reduction(+ : nErrorClosed)
    # endif
    forAll(cells, cI)
    {
        const cell& curCell = cells[cI];

        if( min(curCell) < 0 || max(curCell) > nFaces )
        {
            WarningIn
            (
                "bool checkClosedCells("
                "const polyMeshGen&, const bool, const scalar, labelHashSet*)"
            )   << "Cell " << cI << " contains face labels out of range: "
                << curCell << " Max face index = " << nFaces << endl;

            if( setPtr )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                setPtr->insert(cI);
            }

            ++nErrorClosed;
        }
    }

    if( nErrorClosed > 0 )
    {
        SeriousErrorIn
        (
            "bool checkClosedCells("
            "const polyMeshGen&, const bool, const scalar, labelHashSet*)"
        )  << nErrorClosed << " cells with invalid face labels found"
            << endl;

        return true;
    }

    // Loop through cell faces and sum up the face area vectors for each cell.
    // This should be zero in all vector components
    vectorField sumClosed(cells.size(), vector::zero);

    scalarField sumMagClosed(cells.size(), 0.0);

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const vectorField& areas = mesh.addressingData().faceAreas();

    forAll(own, faceI)
    {
        // Add to owner
        sumClosed[own[faceI]] += areas[faceI];
        sumMagClosed[own[faceI]] += mag(areas[faceI]);

        if( nei[faceI] == -1 )
            continue;

        // Subtract from neighbour
        sumClosed[nei[faceI]] -= areas[faceI];
        sumMagClosed[nei[faceI]] += mag(areas[faceI]);
    }

    label nOpen = 0;
    scalar maxOpenCell = 0;

    label nAspect = 0;
    scalar maxAspectRatio = 0;

    const scalarField& vols = mesh.addressingData().cellVolumes();

    // Check the sums
    forAll(sumClosed, cellI)
    {
        scalar maxOpen = max
        (
            sumClosed[cellI].component(vector::X),
            max
            (
                sumClosed[cellI].component(vector::Y),
                sumClosed[cellI].component(vector::Z)
            )
        );

        maxOpenCell = max(maxOpenCell, maxOpen);

        if( maxOpen > SMALL*max(1.0, sumMagClosed[cellI]) )
        {
            if( report )
            {
                Pout<< "Cell " << cellI << " is not closed. "
                    << "Face area vectors sum up to " << mag(sumClosed[cellI])
                    << " directionwise " << sumClosed[cellI] << " or "
                    << mag(sumClosed[cellI])
                       /(mag(sumMagClosed[cellI]) + VSMALL)
                    << " of the area of all faces of the cell. " << endl
                    << "    Area magnitudes sum up to "
                    << sumMagClosed[cellI] << endl;
            }

            if( setPtr )
            {
                setPtr->insert(cellI);
            }

            ++nOpen;
        }

        scalar aspectRatio =
            1.0/6.0*sumMagClosed[cellI]/pow(vols[cellI], 2.0/3.0);

        maxAspectRatio = max(maxAspectRatio, aspectRatio);

        if( aspectRatio > aspectWarn )
        {
            if( report )
            {
                Pout<< "High aspect ratio for cell " << cellI
                    << ": " << aspectRatio << endl;
            }

            ++nAspect;
        }
    }

    reduce(nOpen, sumOp<label>());
    reduce(maxOpenCell, maxOp<scalar>());

    reduce(nAspect, sumOp<label>());
    reduce(maxAspectRatio, maxOp<scalar>());

    if( nOpen > 0 )
    {
        SeriousErrorIn
        (
            "bool checkClosedCells("
            "const polyMeshGen&, const bool, const scalar, labelHashSet*)"
        )   << nOpen << " open cells found. Max cell openness: "
            << maxOpenCell << endl;

        return true;
    }

    if( nAspect > 0 )
    {
        SeriousErrorIn
        (
            "bool checkClosedCells("
            "const polyMeshGen&, const bool, const scalar, labelHashSet*)"
        )   << nAspect << " high aspect ratio cells found.  "
            << "Max aspect ratio: " << maxAspectRatio
            << endl;

        return true;
    }

    if( report )
    {
        Info<< "Max cell openness = " << maxOpenCell
            << "  Max aspect ratio = " << maxAspectRatio
            << ".  All cells OK.\n" << endl;
    }

    return false;
}

bool checkCellVolumes
(
    const polyMeshGen& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    const scalarField& vols = mesh.addressingData().cellVolumes();

    scalar minVolume = GREAT;
    scalar maxVolume = -GREAT;

    label nNegVolCells = 0;

    forAll(vols, cellI)
    {
        if( vols[cellI] < VSMALL )
        {
            if( report )
                SeriousErrorIn
                (
                    "bool checkCellVolumes("
                    "const polyMeshGen&, const bool, labelHashSet*)"
                )   << "Zero or negative cell volume detected for cell "
                    << cellI << ".  Volume = " << vols[cellI] << endl;

            if( setPtr )
                setPtr->insert(cellI);

            ++nNegVolCells;
        }

        minVolume = min(minVolume, vols[cellI]);
        maxVolume = max(maxVolume, vols[cellI]);
    }

    reduce(minVolume, minOp<scalar>());
    reduce(maxVolume, maxOp<scalar>());
    reduce(nNegVolCells, sumOp<label>());

    if( minVolume < VSMALL )
    {
        SeriousErrorIn
        (
            "bool checkCellVolumes("
            "const polyMeshGen&, const bool, labelHashSet*)"
        )   << "Zero or negative cell volume detected.  "
            << "Minimum negative volume: "
            << minVolume << ".\nNumber of negative volume cells: "
            << nNegVolCells << ".  This mesh is invalid"
            << endl;

        return true;
    }
    else
    {
        if( report )
        {
            Info<< "Min volume = " << minVolume
                << ". Max volume = " << maxVolume
                << ".  Total volume = " << sum(vols)
                << ".  Cell volumes OK.\n" << endl;
        }

        return false;
    }
}

bool checkFaceAreas
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar minFaceArea,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    const vectorField& faceAreas = mesh.addressingData().faceAreas();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    scalar minArea = VGREAT;
    scalar maxArea = -VGREAT;

    # ifdef USE_OMP
    # pragma omp parallel if( own.size() > 100 )
    # endif
    {
        scalar localMaxArea(-VGREAT), localMinArea(VGREAT);

        # ifdef USE_OMP
        # pragma omp for schedule(guided)
        # endif
        forAll(faceAreas, faceI)
        {
            if( changedFacePtr && !changedFacePtr->operator[](faceI) )
                continue;

            const scalar magFaceArea = mag(faceAreas[faceI]);

            if( magFaceArea < minFaceArea )
            {
                if( report )
                {
                    if( nei[faceI] != -1 )
                    {
                        Pout<< "Zero or negative face area detected for "
                            << "internal face " << faceI << " between cells "
                            << own[faceI] << " and " << nei[faceI]
                            << ".  Face area magnitude = "
                            << magFaceArea << endl;
                    }
                    else
                    {
                        Pout<< "Zero or negative face area detected for "
                            << "boundary face " << faceI << " next to cell "
                            << own[faceI] << ".  Face area magnitude = "
                            << magFaceArea << endl;
                    }
                }

                if( setPtr )
                {
                    # ifdef USE_OMP
                    # pragma omp critical
                    # endif
                    setPtr->insert(faceI);
                }
            }

            localMinArea = Foam::min(localMinArea, magFaceArea);
            localMaxArea = Foam::max(localMaxArea, magFaceArea);
        }

        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            minArea = Foam::min(minArea, localMinArea);
            maxArea = Foam::max(maxArea, localMaxArea);
        }
    }

    reduce(minArea, minOp<scalar>());
    reduce(maxArea, maxOp<scalar>());

    if( minArea < VSMALL )
    {
        SeriousErrorIn
        (
            "bool checkFaceAreas("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*)"
        )   << "Zero or negative face area detected.  Minimum negative area: "
            << minArea << ". This mesh is invalid"
            << endl;

        return true;
    }
    else
    {
        if( report )
        {
            Info<< "Minumum face area = " << minArea
                << ". Maximum face area = " << maxArea
                << ".  Face area magnitudes OK.\n" << endl;
        }

        return false;
    }
}

bool checkTetQuality
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar minTetQuality,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    const pointFieldPMG& points = mesh.points();
    const faceListPMG& faces = mesh.faces();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const vectorField& cCentres = mesh.addressingData().cellCentres();

    label nBadFaces(0);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100) reduction(+:nBadFaces)
    # endif
    forAll(owner, faceI)
    {
        if( changedFacePtr && !(*changedFacePtr)[faceI] )
            continue;

        const face& f = faces[faceI];

        const label nDecomposed = f.size() - 1;

        forAll(f, pI)
        {
            bool badQualityFace(false);

            for(label j=1;j<nDecomposed;++j)
            {
                const label fpJ = f[(pI+j) % f.size()];
                const label nfpJ = f[(pI+j+1) % f.size()];

                //- check the tet on the neighbour side
                tetrahedron<point, point> tetOwn
                (
                    cCentres[owner[faceI]],
                    points[f[pI]],
                    points[fpJ],
                    points[nfpJ]
                );

                const scalar tetQualityOwn = help::tetQuality(tetOwn);

                if( tetQualityOwn < minTetQuality )
                {
                    ++nBadFaces;

                    if( report )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(output)
                        # endif
                        Pout<< "Face " << faceI
                            << " has a triangle that points the wrong way."
                            << endl
                            << "Tet quality: " << tetQualityOwn
                            << " Face " << faceI
                            << endl;
                    }

                    if( setPtr )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(insertingBadFaces)
                        # endif
                        setPtr->insert(faceI);
                    }

                    //- found a problematic face. Do not search further
                    badQualityFace = true;
                    continue;
                }

                if( neighbour[faceI] < 0 )
                    continue;

                //- check the tet on the neighbour side
                tetrahedron<point, point> tetNei
                (
                    cCentres[neighbour[faceI]],
                    points[f[pI]],
                    points[nfpJ],
                    points[fpJ]
                );

                const scalar tetQualityNei = help::tetQuality(tetNei);

                if( tetQualityNei < minTetQuality )
                {
                    ++nBadFaces;

                    if( report )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(output)
                        # endif
                        Pout<< "Face " << faceI
                            << " has a triangle that points the wrong way."
                            << endl
                            << "Tet quality: " << tetQualityNei
                            << " Face " << faceI
                            << endl;
                    }

                    if( setPtr )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(insertingBadFaces)
                        # endif
                        setPtr->insert(faceI);
                    }

                    //- found a problematic face. Do not search further
                    badQualityFace = true;
                    continue;
                }
            }

            if( badQualityFace )
                continue;
        }
    }

    reduce(nBadFaces, sumOp<label>());

    if( Pstream::parRun() && nBadFaces && setPtr )
    {
        //- make sure that processor faces are marked on both sides
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        //- send and receive data where needed
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            labelLongList markedFaces;
            for(label faceI=0;faceI<size;++faceI)
            {
                if( setPtr->found(start+faceI) )
                    markedFaces.append(faceI);
            }

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                markedFaces.byteSize()
            );

            toOtherProc << markedFaces;
        }

        forAll(procBoundaries, patchI)
        {
            labelList receivedData;
            IPstream fromOtheProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );
            fromOtheProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();
            forAll(receivedData, i)
                setPtr->insert(start+receivedData[i]);
        }
    }

    if( nBadFaces != 0 )
    {
        WarningIn
        (
            "bool checkTetQuality("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*)"
        )   << "Found " << nBadFaces
            << " faces with negative tet decomposition (minTetQuality < "
            << minTetQuality << ")." << endl;

        return true;
    }

    return false;
}


//- check twist of faces
bool checkMinTwist
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar minTwist,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)

{
    if( (minTwist < (-1.-SMALL)) || (minTwist > (1.+SMALL)) )
    {
        FatalErrorIn
        (
            "bool checkMinTwist("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*"
        )   << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }

    const faceListPMG& faces = mesh.faces();

    label nWarped = 0;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const vectorField& fCentres = mesh.addressingData().faceCentres();
    const vectorField& centres = mesh.addressingData().cellCentres();
    const pointFieldPMG& points = mesh.points();

    const label nInternalFaces = mesh.nInternalFaces();

    # ifdef USE_OMP
    # pragma omp parallel if ( nInternalFaces > 1000 ) reduction(+:nWarped)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(guided, 100)
        # endif
        for(label faceI=0; faceI<nInternalFaces;++faceI)
        {
            if( changedFacePtr && !(*changedFacePtr)[faceI] )
                continue;

            const face& f = faces[faceI];

            if( f.size() > 3 )
            {
                vector nf = centres[nei[faceI]] - centres[own[faceI]];
                const scalar magNf = mag(nf) + VSMALL;
                nf /= magNf;

                forAll(f, fpI)
                {
                    const triangle<point, point> triangle
                    (
                        points[f[fpI]],
                        points[f.nextLabel(fpI)],
                        fCentres[faceI]
                    );

                    const vector triArea = triangle.normal();

                    const scalar magTri = (mag(triArea) + VSMALL);

                    if( magTri > VSMALL &&
                            ((nf & (triArea/magTri)) < minTwist) )
                    {
                        ++nWarped;

                        if( setPtr )
                        {
                            # ifdef USE_OMP
                            # pragma omp critical(badFace)
                            # endif
                            setPtr->insert(faceI);
                        }

                        break;
                    }
                }
            }
        }

        //- boundary faces
        const label start = nInternalFaces;
        label end = faces.size();
        if( Pstream::parRun() )
        {
            const PtrList<processorBoundaryPatch>& procBoundaries =
                    mesh.procBoundaries();
            end = procBoundaries[0].patchStart();
        }

        # ifdef USE_OMP
        # pragma omp parallel for schedule(static, 1) reduction(+:nWarped)
        # endif
        for(label faceI = start; faceI < end;++faceI)
        {
            const face& f = faces[faceI];

            vector nf(vector::zero);

            if( f.size() > 3 )
            {
                nf = fCentres[faceI] - centres[own[faceI]];
                nf /= mag(nf) + VSMALL;
            }

            forAll(f, fpI)
            {
                const triangle<point, point> triangle
                (
                    points[f[fpI]],
                    points[f.nextLabel(fpI)],
                    fCentres[faceI]
                );

                const vector triArea = triangle.normal();

                const scalar magTri = mag(triArea);

                if( magTri > VSMALL &&
                        ((nf & triArea/magTri) < minTwist) )
                {
                    nWarped++;

                    if( setPtr )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(badFace)
                        # endif
                        setPtr->insert(faceI);
                    }
                    break;
                }
            }
        }
    }

    //- if running in parallel
    if( Pstream::parRun() )
    {
        //- check parallel boundaries
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            vectorField cCentres(procBoundaries[patchI].patchSize());

            forAll(cCentres, faceI)
                cCentres[faceI] = centres[own[start+faceI]];

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                cCentres.byteSize()
            );

            toOtherProc << cCentres;
        }

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            vectorField otherCentres;

            IPstream fromOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> otherCentres;

            forAll(otherCentres, fI)
            {
                const label faceI = start + fI;

                if( changedFacePtr && !(*changedFacePtr)[faceI] )
                    continue;

                const point& cOwn = centres[own[faceI]];
                const point& cNei = otherCentres[fI];

                const face& f = faces[faceI];

                if( f.size() > 3 )
                {
                    vector nf = cNei - cOwn;
                    nf /= (mag(nf) + VSMALL);

                    forAll(f, fpI)
                    {
                        const triangle<point, point> triangle
                        (
                            points[f[fpI]],
                            points[f.nextLabel(fpI)],
                            fCentres[faceI]
                        );

                        const vector triArea = triangle.normal();

                        const scalar magTri = (mag(triArea) + VSMALL);

                        if( magTri > VSMALL &&
                                ((nf & triArea/magTri) < minTwist) )
                        {
                            ++nWarped;

                            if( setPtr )
                            {
                                setPtr->insert(faceI);
                            }

                            break;
                        }
                    }
                }
            }
        }
    }

    reduce(nWarped,sumOp<label>());
    if( report )
    {
        if( nWarped > 0 )
        {
            Info << "There are " << nWarped
                 << " faces with cosine of the angle"
                 << " between triangle normal and face normal less than "
                 << minTwist << nl << endl;
        }
        else
        {
            Info << "All faces are flat in that the cosine of the angle"
                 << " between triangle normal and face normal less than "
                 << minTwist << nl << endl;
        }
    }

    if( nWarped > 0 )
    {
        if( report )
        {
            WarningIn
            (
                "bool checkMinTwist("
                "const polyMeshGen&, const bool, const scalar,"
                " labelHashSet*, const boolList*"
            )   << nWarped << " faces with severe warpage "
                << "(cosine of the angle between triangle normal and "
                << "face normal < " << minTwist << ") found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }

}


bool checkCellDeterminant
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar warnDet,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    const cellListPMG& cells = mesh.cells();
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const label nInternalFaces = mesh.nInternalFaces();

    scalar minDet = VGREAT;
    scalar sumDet = 0.0;
    label nSumDet = 0;
    label nWarnDet = 0;
    boolList affectedCells(cells.size(), false);
    if( changedFacePtr )
    {
        const boolList& changedFaces = *changedFacePtr;
        forAll(changedFaces, fI)
        {
            if( changedFaces[fI] )
            {
                affectedCells[own[fI]] = true;

                if( fI >= 0 && fI < nInternalFaces )
                {
                    affectedCells[nei[fI]] = true;
                }
            }
        }
    }

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        scalar localMinDet(VGREAT);

        # ifdef USE_OMP
        # pragma omp for schedule(guided, 50) \
        reduction(+: nSumDet, sumDet, nWarnDet)
        # endif
        forAll(cells, cI)
        {
            if( affectedCells[cI] )
            {
                const cell& c = cells[cI];
                const vectorField& areas = mesh.addressingData().faceAreas();

                tensor areaSum(tensor::zero);
                scalar magAreaSum = 0.0;

                forAll(c, fI)
                {
                    const label faceI = c[fI];

                    const scalar magArea = mag(areas[faceI]) + VSMALL;

                    magAreaSum += magArea;
                    areaSum += areas[faceI]*(areas[faceI]/magArea);
                }

                const scalar scaledDet =
                        det(areaSum/magAreaSum)/0.037037037037037;

                localMinDet = min(localMinDet, scaledDet);
                sumDet += scaledDet;
                ++nSumDet;

                if( scaledDet < warnDet )
                {
                    if( setPtr )
                    {
                        //Insert all faces of the cell
                        forAll(c, fI)
                        {
                            const label faceI = c[fI];

                            # ifdef USE_OMP
                            # pragma omp critical(badFace)
                            # endif
                            setPtr->insert(faceI);
                        }
                    }
                    ++nWarnDet;
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp critical(minDet)
        # endif
        minDet = min(minDet, localMinDet);
    }

    reduce(nSumDet, sumOp<label>());
    reduce(nWarnDet, sumOp<label>());
    reduce(minDet, minOp<scalar>());
    reduce(sumDet, sumOp<scalar>());

    if( report )
    {
        if( nSumDet > 0 )
        {
            Info << "Cell determinant (1 = uniform cube) : average = "
                 << sumDet / nSumDet << " min = " << minDet << endl;
        }

        if( nWarnDet > 0 )
        {
            Info << "There are " << nWarnDet
                 << " cells with determinant < " << warnDet << '.' << nl
                 << endl;
        }

        else
        {
            Info << "All faces have determinant > " << warnDet << '.' << nl
                 << endl;
        }
    }

    if( nWarnDet > 0 )
    {
        if( report )
        {
            WarningIn
            (
                "bool checkCellDeterminant("
                " const polyMeshGen&, const bool, const scalar,"
                " labelHashSet*, const boolList*"
            )   << nWarnDet << " cells with determinant < " << warnDet
                << " found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}

void checkMinVolRatio
(
    const polyMeshGen& mesh,
    scalarField& volRatio,
    const boolList* changedFacePtr
)
{
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const label nInternalFaces = mesh.nInternalFaces();
    const scalarField& vols = mesh.addressingData().cellVolumes();

    volRatio.setSize(own.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule (dynamic, 100)
    # endif
    for(label faceI=0; faceI<nInternalFaces;++faceI)
    {
        volRatio[faceI] = 1.0;

        if( changedFacePtr && !changedFacePtr->operator[](faceI) )
            continue;

        const scalar ownVol
        (
            Foam::mag(vols[own[faceI]])
        );

        const scalar neiVol
        (
            Foam::mag(vols[nei[faceI]])
        );

        volRatio[faceI] = Foam::min(neiVol, ownVol)/max(neiVol, ownVol);
    }

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
                mesh.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            scalarField vol(procBoundaries[patchI].patchSize());
            const label start = procBoundaries[patchI].patchStart();

            forAll(vol, fI)
            {
                const label faceI = start + fI;
                vol[fI] = Foam::mag(vols[own[faceI]]);
            }

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                vol.byteSize()
            );

            toOtherProc << vol;

        }

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            scalarField otherVol;

            IPstream fromOtheProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtheProc >> otherVol;

            forAll(otherVol, fI)
            {
                const label faceI = start + fI;

                if( changedFacePtr && !changedFacePtr->operator[](faceI) )
                    continue;

                const scalar ownVol =
                        Foam::mag(vols[own[faceI]]);
                const scalar neiVol = otherVol[fI];
                volRatio[faceI] =
                        Foam::min(neiVol, ownVol)/max(neiVol, ownVol);
            }
        }
    }
}

bool checkMinVolRatio
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar warnVolRatio,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    scalarField volRatio;
    checkMinVolRatio(mesh, volRatio, changedFacePtr);

    //- for all internal faces check the uniformity of the mesh at faces
    const label nInternalFaces = mesh.nInternalFaces();

    scalar maxVolRatio = 0.0;
    scalar minVolRatio = VGREAT;

    scalar sumVolRatio = 0.0;

    label nWarnRatio = 0;

    # ifdef USE_OMP
    # pragma omp parallel \
    reduction(+ : sumVolRatio, nWarnRatio)
    # endif
    {
        scalar localMinVolRatio(VGREAT);
        scalar localMaxVolRatio(0.0);

        # ifdef USE_OMP
        # pragma omp for schedule(guided, 50)
        # endif
        for(label faceI=0;faceI<nInternalFaces;++faceI)
        {
            if( changedFacePtr && !changedFacePtr->operator[](faceI) )
                continue;

            const scalar ratio = volRatio[faceI];

            if( ratio < warnVolRatio )
            {
                if( setPtr )
                {
                    # ifdef USE_OMP
                    # pragma omp critical(badFace)
                    # endif
                    setPtr->insert(faceI);
                }

                ++nWarnRatio;
            }

            localMaxVolRatio = Foam::max(localMaxVolRatio, ratio);
            localMinVolRatio = Foam::min(localMinVolRatio, ratio);
            sumVolRatio += ratio;
        }

        # ifdef USE_OMP
        # pragma omp critical(minVolRatio)
        # endif
        {
            maxVolRatio = Foam::max(maxVolRatio , localMaxVolRatio);
            minVolRatio = Foam::min(minVolRatio , localMinVolRatio);
        }
    }

    label counter = nInternalFaces;

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            for(label fI=0;fI<size;++fI)
            {
                const label faceI = start + fI;

                const scalar ratio = volRatio[faceI];

                if(ratio < warnVolRatio)
                {
                    if( setPtr )
                        setPtr->insert(faceI);

                    ++nWarnRatio;
                }

                maxVolRatio = Foam::max(maxVolRatio , ratio);
                minVolRatio = Foam::min(minVolRatio , ratio);
                sumVolRatio += ratio;
            }

            if( procBoundaries[patchI].owner() )
                counter += size;
        }
    }

    reduce(maxVolRatio, maxOp<scalar>());
    reduce(minVolRatio, minOp<scalar>());
    reduce(sumVolRatio, sumOp<scalar>());
    reduce(nWarnRatio, sumOp<label>());
    reduce(counter, sumOp<label>());

    if( counter > 0 )
    {
        if( report && minVolRatio < warnVolRatio )
            Info << "Small ratio on " << nWarnRatio << "faces."
                 << endl;
    }

    if( report )
    {
        if( counter > 0 )
            Info << "Mesh volume ratio Max: "
                 << maxVolRatio
                 << " Min: " << minVolRatio
                 << " average: " << sumVolRatio/counter
                 << endl;
    }

    return false;
}


bool checkTriangleTwist
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar minTwist,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    if( (minTwist < (-1.-SMALL)) || (minTwist > (1.+SMALL)) )
    {
        FatalErrorIn
        (
            "bool checkTriangleTwist("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*"
        )   << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }

    const faceListPMG& faces = mesh.faces();

    label nWarped = 0;

    const vectorField& fCentres = mesh.addressingData().faceCentres();
    const pointFieldPMG& points = mesh.points();

    # ifdef USE_OMP
    # pragma omp parallel if ( faces.size() > 1000 ) reduction(+:nWarped)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(faces, faceI)
        {
            if( changedFacePtr && !changedFacePtr->operator[](faceI) )
                continue;

            const face& f = faces[faceI];

            if( f.size() > 3 )
            {
                //- find starting triangle with non-zero area
                bool isTwisted(false);
                forAll(f, pI)
                {
                    if( isTwisted )
                        break;

                    const triangle<point, point> tri
                    (
                        points[f[pI]],
                        points[f.nextLabel(pI)],
                        fCentres[faceI]
                    );

                    vector prevN = tri.normal();

                    const scalar magTri = mag(prevN) + VSMALL;

                    prevN /= magTri;

                    //- compare with other possible triangles
                    for(label pJ=pI+1;pJ<f.size();++pJ)
                    {
                        const triangle<point, point> tri
                        (
                            points[f[pJ]],
                            points[f.nextLabel(pJ)],
                            fCentres[faceI]
                        );

                        vector triN = tri.normal();

                        scalar magTriN = mag(triN) + VSMALL;

                        triN /= magTriN;

                        if( (prevN & triN) < minTwist )
                        {
                            isTwisted = true;
                            break;
                        }
                    }
                }

                if( isTwisted )
                {
                    ++nWarped;

                    if( setPtr )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(badFace)
                        # endif
                        setPtr->insert(faceI);
                    }
                }
            }
        }
    }

    reduce(nWarped, sumOp<label>());

    if( report )
    {
        if( nWarped > 0 )
        {
            Info << "There are " << nWarped
                 << " faces with cosine of the angle "
                 << " between consecutive triangle normals less than "
                 << minTwist << nl << endl;
        }
        else
        {
            Info << "All faces are flat in that the cosine of the angle"
                 << " between consecutive triangle normals is less than "
                 << minTwist << nl << endl;
        }
    }

    if( nWarped > 0 )
    {
        if( report )
        {
            WarningIn
            (
                "bool checkTriangleTwist("
                "const polyMeshGen&, const bool, const scalar,"
                " labelHashSet*, const boolList*"
            )   << nWarped << " faces with severe warpage "
                << "(cosine of the angle between consecutive triangle normals"
                << " < " << minTwist << ")found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool checkCellPartTetrahedra
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar minPartTet,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    const pointFieldPMG& points = mesh.points();
    const faceListPMG& faces = mesh.faces();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const vectorField& fCentres = mesh.addressingData().faceCentres();
    const vectorField& cCentres = mesh.addressingData().cellCentres();

    label nNegVolCells = 0;

    # ifdef USE_OMP
    # pragma omp parallel for if( owner.size() > 100 ) \
    schedule(guided) reduction(+ : nNegVolCells)
    # endif
    forAll(owner, faceI)
    {
        if( changedFacePtr && !(*changedFacePtr)[faceI] )
            continue;

        const face& f = faces[faceI];

        bool badFace(false);

        forAll(f, eI)
        {
            const tetrahedron<point, point> tetOwn
            (
                fCentres[faceI],
                points[f.nextLabel(eI)],
                points[f[eI]],
                cCentres[owner[faceI]]
            );

            if( tetOwn.mag() < minPartTet )
            {
                if( report )
                {
                    # ifdef USE_OMP
                    # pragma omp critical(report)
                    # endif
                    Pout<< "Zero or negative cell volume detected for cell "
                        << owner[faceI] << "." << endl;
                }

                badFace = true;
            }

            if( neighbour[faceI] < 0 )
                continue;

            const tetrahedron<point, point> tetNei
            (
                fCentres[faceI],
                points[f[eI]],
                points[f.nextLabel(eI)],
                cCentres[neighbour[faceI]]
            );

            if( tetNei.mag() < minPartTet )
            {
                if( report )
                {
                    # ifdef USE_OMP
                    # pragma omp critical(report)
                    # endif
                    Pout<< "Zero or negative cell volume detected for cell "
                        << neighbour[faceI] << "." << endl;
                }

                badFace = true;
            }
        }

        if( badFace )
        {
            if( setPtr )
            {
                # ifdef USE_OMP
                # pragma omp critical(badFace)
                # endif
                setPtr->insert(faceI);
            }

            ++nNegVolCells;
        }
    }

    if( setPtr )
    {
        //- ensure that faces are selected at both sides
        const PtrList<processorBoundaryPatch>& procBnd = mesh.procBoundaries();
        forAll(procBnd, patchI)
        {
            const label start = procBnd[patchI].patchStart();
            const label size = procBnd[patchI].patchSize();

            labelLongList sendData;
            for(label faceI=0;faceI<size;++faceI)
            {
                if( setPtr->found(faceI+start) )
                    sendData.append(faceI);
            }

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBnd[patchI].neiProcNo(),
                sendData.byteSize()
            );

            toOtherProc << sendData;
        }

        forAll(procBnd, patchI)
        {
            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::commsTypes::blocking,
                procBnd[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            const label start = procBnd[patchI].patchStart();
            forAll(receivedData, i)
                setPtr->insert(start+receivedData[i]);
        }
    }

    reduce(nNegVolCells, sumOp<label>());

    if( nNegVolCells != 0 )
    {
        WarningIn
        (
            "bool checkCellPartTetrahedra("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*)"
        )   << nNegVolCells << " zero or negative part tetrahedra detected."
            << endl;

        return true;
    }
    else
    {
        if( report )
            Info<< "Part cell tetrahedra OK.\n" << endl;

        return false;
    }
}

void checkFaceDotProduct
(
    const polyMeshGen& mesh,
    scalarField& faceDotProduct,
    const boolList* changedFacePtr
)
{
    //- for all internal faces check theat the d dot S product is positive
    const vectorField& centres = mesh.addressingData().cellCentres();
    const vectorField& areas = mesh.addressingData().faceAreas();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const label nInternalFaces = mesh.nInternalFaces();

    faceDotProduct.setSize(own.size());
    faceDotProduct = 1.0;

    # ifdef USE_OMP
    # pragma omp parallel if( nInternalFaces > 1000 )
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        for(label faceI=0;faceI<nInternalFaces;++faceI)
        {
            if( changedFacePtr && !(*changedFacePtr)[faceI] )
                continue;

            const vector d = centres[nei[faceI]] - centres[own[faceI]];
            const vector& s = areas[faceI];

            faceDotProduct[faceI] = (d & s)/(mag(d)*mag(s) + VSMALL);
        }
    }

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            vectorField cCentres(procBoundaries[patchI].patchSize());
            forAll(cCentres, faceI)
                cCentres[faceI] = centres[own[start+faceI]];

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                cCentres.byteSize()
            );

            toOtherProc << cCentres;
        }

        forAll(procBoundaries, patchI)
        {
            vectorField otherCentres;
            IPstream fromOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> otherCentres;

            //- calculate skewness at processor faces
            const label start = procBoundaries[patchI].patchStart();
            # ifdef USE_OMP
            # pragma omp parallel
            # endif
            {
                # ifdef USE_OMP
                # pragma omp for schedule(dynamic, 100)
                # endif
                forAll(otherCentres, fI)
                {
                    const label faceI = start + fI;

                    if( changedFacePtr && !(*changedFacePtr)[faceI] )
                        continue;

                    const point& cOwn = centres[own[faceI]];
                    const point& cNei = otherCentres[fI];
                    const vector d = cNei - cOwn;
                    const vector& s = areas[faceI];

                    faceDotProduct[faceI] = (d & s)/(mag(d)*mag(s) + VSMALL);
                }
            }
        }
    }
}

bool checkFaceDotProduct
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar nonOrthWarn,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    //- for all internal faces check theat the d dot S product is positive
    scalarField faceDotProduct;
    checkFaceDotProduct(mesh, faceDotProduct, changedFacePtr);

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const label nInternalFaces = mesh.nInternalFaces();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(nonOrthWarn/180.0*M_PI);

    scalar minDDotS = VGREAT;
    scalar sumDDotS = 0.0;

    label severeNonOrth = 0;
    label errorNonOrth = 0;
    label counter = 0;

    # ifdef USE_OMP
    # pragma omp parallel if( nInternalFaces > 1000 ) \
    reduction(+ : severeNonOrth, errorNonOrth, sumDDotS)
    # endif
    {
        scalar localMinDDotS(VGREAT);
        # ifdef USE_OMP
        # pragma omp for schedule(guided)
        # endif
        for(label faceI=0;faceI<nInternalFaces;++faceI)
        {
            if( changedFacePtr && !(*changedFacePtr)[faceI] )
                continue;

            const scalar dDotS = faceDotProduct[faceI];

            if( dDotS < severeNonorthogonalityThreshold )
            {
                if( dDotS > SMALL )
                {
                    if( report )
                    {
                        // Severe non-orthogonality but mesh still OK
                        # ifdef USE_OMP
                        # pragma omp critical(report)
                        # endif
                        Pout<< "Severe non-orthogonality for face " << faceI
                            << " between cells " << own[faceI]
                            << " and " << nei[faceI]
                            << ": Angle = "
                            << ::acos(dDotS)/M_PI*180.0
                            << " deg." << endl;
                    }

                    if( setPtr )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(badFace)
                        # endif
                        setPtr->insert(faceI);
                    }

                    ++severeNonOrth;
                }
                else
                {
                    ++errorNonOrth;

                    if( setPtr )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(badFace)
                        # endif
                        setPtr->insert(faceI);
                    }
                }
            }

            localMinDDotS = Foam::min(dDotS, localMinDDotS);
            sumDDotS += dDotS;
        }

        # ifdef USE_OMP
        # pragma omp critical(minDDotS)
        # endif
        minDDotS = Foam::min(minDDotS, localMinDDotS);
    }

    counter += nInternalFaces;

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        const label start = procBoundaries[0].patchStart();

        # ifdef USE_OMP
        # pragma omp parallel \
        reduction(+ : severeNonOrth, errorNonOrth, sumDDotS, counter)
        # endif
        {
            scalar localMinDDotS(VGREAT);

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            for(label faceI=start;faceI<own.size();++faceI)
            {
                if( changedFacePtr && !(*changedFacePtr)[start+faceI] )
                            continue;

                const scalar dDotS = faceDotProduct[faceI];

                if( dDotS < severeNonorthogonalityThreshold )
                {
                    if( dDotS > SMALL )
                    {
                        if( report )
                        {
                            // Severe non-orthogonality but mesh still OK
                            # ifdef USE_OMP
                            # pragma omp critical(report)
                            # endif
                            {
                                const scalar angle
                                (
                                    Foam::acos(dDotS) /
                                    M_PI * 180.0
                                );
                                Pout<< "Severe non-orthogonality for face "
                                    << start+faceI
                                    << ": Angle = "
                                    << angle
                                    << " deg." << endl;
                            }
                        }

                        if( setPtr )
                        {
                            # ifdef USE_OMP
                            # pragma omp critical(badFace)
                            # endif
                            setPtr->insert(start+faceI);
                        }

                        ++severeNonOrth;
                    }
                    else
                    {
                        ++errorNonOrth;

                        if( setPtr )
                        {
                            # ifdef USE_OMP
                            # pragma omp critical(badFace)
                            # endif
                            setPtr->insert(start+faceI);
                        }
                    }
                }

                localMinDDotS = Foam::min(dDotS, localMinDDotS);
                sumDDotS += 0.5 * dDotS;

                ++counter;
            }

            # ifdef USE_OMP
            # pragma omp critical(minDDotS)
            # endif
            minDDotS = Foam::min(minDDotS, localMinDDotS);
        }
    }

    reduce(minDDotS, minOp<scalar>());
    reduce(sumDDotS, sumOp<scalar>());
    reduce(severeNonOrth, sumOp<label>());
    reduce(errorNonOrth, sumOp<label>());

    reduce(counter, sumOp<label>());

    // Only report if there are some internal faces
    if( counter > 0 )
    {
        if( report && (minDDotS < severeNonorthogonalityThreshold) )
        {
            Info<< "Number of non-orthogonality errors: " << errorNonOrth
                << ". Number of severely non-orthogonal faces: "
                << severeNonOrth  << "." << endl;
        }
    }

    if( report )
    {
        if( counter > 0 )
        {
            Info<< "Mesh non-orthogonality Max: "
                << ::acos(minDDotS)/M_PI*180.0
                << " average: " <<
                   ::acos(sumDDotS/counter)/M_PI*180.0
                << endl;
        }
    }

    if( errorNonOrth > 0 )
    {
        WarningIn
        (
            "checkFaceDotProduct("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*)"
        )   << "Error in non-orthogonality detected" << endl;

        return true;
    }
    else
    {
        if( report )
            Info<< "Non-orthogonality check OK.\n" << endl;

        return false;
    }
}

bool checkFacePyramids
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar minPyrVol,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    // check whether face area vector points to the cell with higher label
    const vectorField& ctrs = mesh.addressingData().cellCentres();

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const faceListPMG& faces = mesh.faces();
    const pointFieldPMG& points = mesh.points();

    label nErrorPyrs = 0;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided) reduction(+ : nErrorPyrs)
    # endif
    forAll(faces, faceI)
    {
        if( changedFacePtr && !(*changedFacePtr)[faceI] )
            continue;

        // Create the owner pyramid - it will have negative volume
        const scalar pyrVol = pyramidPointFaceRef
        (
            faces[faceI],
            ctrs[owner[faceI]]
        ).mag(points);

        bool badFace(false);

        if( pyrVol > -minPyrVol )
        {
            if( report )
            {
                # ifdef USE_OMP
                # pragma omp critical(report)
                # endif
                Pout<< "bool checkFacePyramids("
                    << "const bool, const scalar, labelHashSet*) : "
                    << "face " << faceI << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVol
                    << " Face " << faces[faceI] << " area: "
                    << faces[faceI].mag(points)
                    << " Owner cell: " << owner[faceI] << endl
                    << "Owner cell vertex labels: "
                    << mesh.cells()[owner[faceI]].labels(faces)
                    << endl;
            }

            badFace = true;
        }

        if( neighbour[faceI] != -1 )
        {
            // Create the neighbour pyramid - it will have positive volume
            const scalar pyrVol =
                pyramidPointFaceRef
                (
                    faces[faceI],
                    ctrs[neighbour[faceI]]
                ).mag(points);

            if( pyrVol < minPyrVol )
            {
                if( report )
                {
                    # ifdef USE_OMP
                    # pragma omp critical(report)
                    # endif
                    Pout<< "bool checkFacePyramids("
                        << "const bool, const scalar, labelHashSet*) : "
                        << "face " << faceI << " points the wrong way. " << endl
                        << "Pyramid volume: " << -pyrVol
                        << " Face " << faces[faceI] << " area: "
                        << faces[faceI].mag(points)
                        << " Neighbour cell: " << neighbour[faceI] << endl
                        << "Neighbour cell vertex labels: "
                        << mesh.cells()[neighbour[faceI]].labels(faces)
                        << endl;
                }

                badFace = true;
            }
        }

        if( badFace )
        {
            if( setPtr )
            {
                # ifdef USE_OMP
                # pragma omp critical(badFace)
                # endif
                setPtr->insert(faceI);
            }

            ++nErrorPyrs;
        }
    }

    reduce(nErrorPyrs, sumOp<label>());

    if( setPtr )
    {
        //- make sure that processor faces are marked on both sides
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        //- send and receive data where needed
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            labelLongList markedFaces;
            for(label faceI=0;faceI<size;++faceI)
            {
                if( setPtr->found(start+faceI) )
                    markedFaces.append(faceI);
            }

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                markedFaces.byteSize()
            );

            toOtherProc << markedFaces;
        }

        forAll(procBoundaries, patchI)
        {
            labelList receivedData;
            IPstream fromOtheProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );
            fromOtheProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();
            forAll(receivedData, i)
                setPtr->insert(start+receivedData[i]);
        }
    }

    if( nErrorPyrs > 0 )
    {
        if( Pstream::master() )
            WarningIn
            (
                "bool checkFacePyramids("
                "const polyMeshGen&, const bool, const scalar,"
                " labelHashSet*, const boolList*)"
            )   << "Error in face pyramids: " << nErrorPyrs
                << " faces pointing the wrong way!"
                << endl;

        return true;
    }
    else
    {
        if( report )
            Info<< "Face pyramids OK.\n" << endl;

        return false;
    }
}

void checkFaceSkewness
(
    const polyMeshGen& mesh,
    scalarField& faceSkewness,
    const boolList* changedFacePtr
)
{
    //- Warn if the skew correction vector is more than skewWarning times
    //- larger than the face area vector
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const label nInternalFaces = mesh.nInternalFaces();

    const vectorField& centres = mesh.addressingData().cellCentres();
    const vectorField& fCentres = mesh.addressingData().faceCentres();

    faceSkewness.setSize(own.size());
    faceSkewness = 0.0;

    //- check internal faces
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    for(label faceI=0;faceI<nInternalFaces;++faceI)
    {
        if( changedFacePtr && !changedFacePtr->operator[](faceI) )
            continue;

        const scalar dOwn = mag(fCentres[faceI] - centres[own[faceI]]);
        const scalar dNei = mag(fCentres[faceI] - centres[nei[faceI]]);

        const point faceIntersection =
            centres[own[faceI]]*dNei/(dOwn+dNei)
          + centres[nei[faceI]]*dOwn/(dOwn+dNei);

        faceSkewness[faceI] =
            mag(fCentres[faceI] - faceIntersection)
            /(mag(centres[nei[faceI]] - centres[own[faceI]]) + VSMALL);
    }

    if( Pstream::parRun() )
    {
        //- check parallel boundaries
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            vectorField cCentres(procBoundaries[patchI].patchSize());
            forAll(cCentres, faceI)
                cCentres[faceI] = centres[own[start+faceI]];

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                cCentres.byteSize()
            );

            toOtherProc << cCentres;
        }

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            vectorField otherCentres;
            IPstream fromOtheProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtheProc >> otherCentres;

            //- calculate skewness at processor faces
            # ifdef USE_OMP
            # pragma omp parallel
            # endif
            {
                # ifdef USE_OMP
                # pragma omp for schedule(dynamic, 100)
                # endif
                forAll(otherCentres, fI)
                {
                    const label faceI = start + fI;

                    if(
                        changedFacePtr &&
                        !changedFacePtr->operator[](faceI)
                    )
                        continue;

                    const point& cOwn = centres[own[faceI]];
                    const point& cNei = otherCentres[fI];
                    const scalar dOwn = mag(fCentres[faceI] - cOwn);
                    const scalar dNei = mag(fCentres[faceI] - cNei);

                    const point faceIntersection =
                        cOwn*dNei/(dOwn+dNei)
                      + cNei*dOwn/(dOwn+dNei);

                    faceSkewness[faceI] =
                        mag(fCentres[faceI] - faceIntersection)
                        /(mag(cOwn - cNei) + VSMALL);
                }
            }
        }
    }

    //- boundary faces
    const faceListPMG& faces = mesh.faces();
    const pointFieldPMG& points = mesh.points();
    const PtrList<boundaryPatch>& boundaries = mesh.boundaries();
    forAll(boundaries, patchI)
    {
        label faceI = boundaries[patchI].patchStart();
        const label end = faceI + boundaries[patchI].patchSize();

        for(;faceI<end;++faceI)
        {
            const vector d = fCentres[faceI] - centres[own[faceI]];

            vector n = faces[faceI].normal(points);
            const scalar magn = mag(n);
            if( magn > VSMALL )
            {
                n /= magn;
            }
            else
            {
                continue;
            }

            const vector dn = (n & d) * n;

            faceSkewness[faceI] = mag(d - dn) / (mag(d) + VSMALL);
        }
    }
}

bool checkFaceSkewness
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar warnSkew,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    scalarField faceSkewness;
    checkFaceSkewness(mesh, faceSkewness, changedFacePtr);

    //- Warn if the skew correction vector is more than skewWarning times
    //- larger than the face area vector
    scalar maxSkew = 0.0;
    scalar sumSkew = 0.0;

    label nWarnSkew = 0;

    //- check faces
    # ifdef USE_OMP
    # pragma omp parallel reduction(+ : sumSkew, nWarnSkew)
    # endif
    {
        scalar localMaxSkew(0.0);

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(faceSkewness, faceI)
        {
            if( changedFacePtr && !changedFacePtr->operator[](faceI) )
                continue;

            const scalar skewness = faceSkewness[faceI];

            // Check if the skewness vector is greater than the PN vector.
            if( skewness > warnSkew )
            {
                if( report )
                {
                    # ifdef USE_OMP
                    # pragma omp critical(report)
                    # endif
                    Pout<< " Severe skewness for face " << faceI
                        << " skewness = " << skewness << endl;
                }

                if( setPtr )
                {
                    # ifdef USE_OMP
                    # pragma omp critical(badFace)
                    # endif
                    setPtr->insert(faceI);
                }

                ++nWarnSkew;
            }

            localMaxSkew = Foam::max(localMaxSkew, skewness);
            sumSkew += skewness;
        }

        # ifdef USE_OMP
        # pragma omp critical(maxSkew)
        # endif
        maxSkew = Foam::max(maxSkew, localMaxSkew);
    }

    reduce(maxSkew, maxOp<scalar>());
    reduce(sumSkew, sumOp<scalar>());
    reduce(nWarnSkew, sumOp<label>());

    if( nWarnSkew > 0 )
    {
        WarningIn
        (
            "checkFaceSkewness("
            "const polyMeshGen&, const bool, const scalar,"
            "labelHashSet*, const boolList*)"
        )   << "Large face skewness detected.  Max skewness = " << maxSkew
            << " Average skewness = " << sumSkew/faceSkewness.size()
            << ".\nThis may impair the quality of the result." << nl
            << nWarnSkew << " highly skew faces detected."
            << endl;

        return true;
    }
    else
    {
        if( report )
            Info<< "Max skewness = " << maxSkew
                << " Average skewness = " << sumSkew/faceSkewness.size()
                << ".  Face skewness OK.\n" << endl;

        return false;
    }
}

void checkFaceUniformity
(
    const polyMeshGen& mesh,
    scalarField& faceUniformity,
    const boolList* changedFacePtr
)
{
    //- for all internal faces check the uniformity of the mesh at faces
    const vectorField& centres = mesh.addressingData().cellCentres();
    const vectorField& fCentres = mesh.addressingData().faceCentres();

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const label nInternalFaces = mesh.nInternalFaces();

    faceUniformity.setSize(owner.size());
    faceUniformity = 0.5;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    for(label faceI=0;faceI<nInternalFaces;++faceI)
    {
        if( changedFacePtr && !changedFacePtr->operator[](faceI) )
            continue;

        const scalar dOwn
        (
            Foam::mag(centres[owner[faceI]] - fCentres[faceI])
        );
        const scalar dNei
        (
            Foam::mag(centres[neighbour[faceI]] - fCentres[faceI])
        );

        faceUniformity[faceI] = Foam::min(dOwn, dNei) / (dOwn + dNei);
    }

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            scalarField dst(procBoundaries[patchI].patchSize());
            const label start = procBoundaries[patchI].patchStart();

            forAll(dst, faceI)
            {
                const label fI = start + faceI;
                dst[faceI] = Foam::mag(centres[owner[fI]] - fCentres[fI]);
            }

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                dst.byteSize()
            );

            toOtherProc << dst;
        }

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            scalarField otherDst;
            IPstream fromOtheProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtheProc >> otherDst;

            forAll(otherDst, faceI)
            {
                const label fI = start + faceI;
                const scalar dOwn =
                    Foam::mag(centres[owner[fI]] - fCentres[fI]);
                const scalar dNei = otherDst[faceI];
                faceUniformity[fI] = Foam::min(dOwn, dNei) / (dOwn + dNei);
            }
        }
    }
}

bool checkFaceUniformity
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar warnUniform,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    scalarField faceUniformity;
    checkFaceUniformity(mesh, faceUniformity, changedFacePtr);

    //- for all internal faces check the uniformity of the mesh at faces
    const label nInternalFaces = mesh.nInternalFaces();

    scalar maxUniformity = 0.0;
    scalar minUniformity = VGREAT;

    scalar sumUniformity = 0.0;

    label severeNonUniform = 0;

    # ifdef USE_OMP
    # pragma omp parallel \
    reduction(+ : sumUniformity, severeNonUniform)
    # endif
    {
        scalar localMinUniformity(VGREAT);
        scalar localMaxUniformity(0.0);

        # ifdef USE_OMP
        # pragma omp for schedule(guided)
        # endif
        for(label faceI=0;faceI<nInternalFaces;++faceI)
        {
            if( changedFacePtr && !changedFacePtr->operator[](faceI) )
                continue;

            const scalar uniformity = faceUniformity[faceI];

            if( uniformity < warnUniform )
            {
                if( setPtr )
                {
                    # ifdef USE_OMP
                    # pragma omp critical(badFace)
                    # endif
                    setPtr->insert(faceI);
                }

                ++severeNonUniform;
            }

            localMaxUniformity = Foam::max(localMaxUniformity, uniformity);
            localMinUniformity = Foam::min(localMinUniformity, uniformity);
            sumUniformity += uniformity;
        }

        # ifdef USE_OMP
        # pragma omp critical(maxUniformity)
        # endif
        {
            maxUniformity = Foam::max(maxUniformity, localMaxUniformity);
            minUniformity = Foam::min(minUniformity, localMinUniformity);
        }
    }

    label counter = nInternalFaces;

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            for(label fI=0;fI<size;++fI)
            {
                const label faceI = start + fI;

                const scalar uniformity = faceUniformity[faceI];

                if( uniformity < warnUniform )
                {
                    if( setPtr)
                        setPtr->insert(faceI);

                    ++severeNonUniform;
                }

                maxUniformity = Foam::max(maxUniformity, uniformity);
                minUniformity = Foam::min(minUniformity, uniformity);
                sumUniformity += 0.5 * uniformity;
            }

            if( procBoundaries[patchI].owner() )
                counter += size;
        }
    }

    reduce(maxUniformity, maxOp<scalar>());
    reduce(minUniformity, minOp<scalar>());
    reduce(sumUniformity, sumOp<scalar>());
    reduce(severeNonUniform, sumOp<label>());
    reduce(counter, sumOp<label>());

    // Only report if there are some internal faces
    if( counter > 0 )
    {
        if( report && (minUniformity < warnUniform) )
            Info<< "Number of severely non-uniform faces: "
                << severeNonUniform << "." << endl;
    }

    if( report )
    {
        if( counter > 0 )
            Info<< "Mesh non-uniformity Max: "
                << maxUniformity
                << " Min: " << minUniformity
                << " average: " << sumUniformity/counter
                << endl;
    }

    return false;
}

void checkVolumeUniformity
(
    const polyMeshGen& mesh,
    scalarField& faceUniformity,
    const boolList* changedFacePtr
);

bool checkVolumeUniformity
(
    const polyMeshGen&,
    const bool /*report*/,
    const scalar /*warnUniform*/,
    labelHashSet* /*setPtr*/,
    const boolList* /*changedFacePtr*/
)
{

    return false;
}

bool checkFaceAngles
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar maxDeg,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    if( maxDeg < -SMALL || maxDeg > 180+SMALL )
    {
        FatalErrorIn
        (
            "bool checkFaceAngles("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*)"
        )   << "maxDeg should be [0..180] but is now " << maxDeg
            << abort(FatalError);
    }

    const scalar maxSin = Foam::sin(maxDeg/180.0*M_PI);

    const pointFieldPMG& points = mesh.points();
    const faceListPMG& faces = mesh.faces();
    vectorField faceNormals(mesh.addressingData().faceAreas());
    faceNormals /= mag(faceNormals) + VSMALL;

    scalar maxEdgeSin = 0.0;

    label nConcave = 0;

    # ifdef USE_OMP
    # pragma omp parallel reduction(+ : nConcave)
    # endif
    {
        scalar localMaxEdgeSin(0.0);
        label errorFaceI(-1);

        # ifdef USE_OMP
        # pragma omp for schedule(guided)
        # endif
        forAll(faces, faceI)
        {
            if( changedFacePtr && !changedFacePtr->operator[](faceI) )
                continue;

            const face& f = faces[faceI];

            // Get edge from f[0] to f[size-1];
            vector ePrev(points[f[0]] - points[f[f.size()-1]]);
            scalar magEPrev = mag(ePrev);
            ePrev /= magEPrev + VSMALL;

            forAll(f, fp0)
            {
                // Get vertex after fp
                label fp1 = f.fcIndex(fp0);

                // Normalized vector between two consecutive points
                vector e10(points[f[fp1]] - points[f[fp0]]);
                scalar magE10 = mag(e10);
                e10 /= magE10 + VSMALL;

                if( magEPrev > SMALL && magE10 > SMALL )
                {
                    vector edgeNormal = ePrev ^ e10;
                    scalar magEdgeNormal = mag(edgeNormal);

                    if( magEdgeNormal < maxSin )
                    {
                        // Edges (almost) aligned -> face is ok.
                    }
                    else
                    {
                        // Check normal
                        edgeNormal /= magEdgeNormal;

                        if( (edgeNormal & faceNormals[faceI]) < SMALL )
                        {
                            if( faceI != errorFaceI )
                            {
                                // Count only one error per face.
                                errorFaceI = faceI;
                                ++nConcave;
                            }

                            if( setPtr )
                            {
                                # ifdef USE_OMP
                                # pragma omp critical
                                # endif
                                {
                                    setPtr->insert(faceI);
                                }
                            }

                            localMaxEdgeSin =
                                Foam::max(localMaxEdgeSin, magEdgeNormal);
                        }
                    }
                }

                ePrev = e10;
                magEPrev = magE10;
            }
        }

        # ifdef USE_OMP
        # pragma omp critical
        # endif
        maxEdgeSin = Foam::max(maxEdgeSin, localMaxEdgeSin);
    }

    reduce(nConcave, sumOp<label>());
    reduce(maxEdgeSin, maxOp<scalar>());

    if( report )
    {
        if( maxEdgeSin > SMALL )
        {
            scalar maxConcaveDegr =
                Foam::asin(Foam::min(1.0, maxEdgeSin))
             * 180.0/M_PI;

            Warning<< "There are " << nConcave
                << " faces with concave angles between consecutive"
                << " edges. Max concave angle = "
                << maxConcaveDegr
                << " degrees.\n" << endl;
        }
        else
        {
            Info<< "All angles in faces are convex or less than "  << maxDeg
                << " degrees concave.\n" << endl;
        }
    }

    if( nConcave > 0 )
    {
        WarningIn
        (
            "bool checkFaceAngles("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*)"
        )   << nConcave  << " face points with severe concave angle (> "
            << maxDeg << " deg) found.\n"
            << endl;

        return true;
    }
    else
    {
        return false;
    }
}

// Check warpage of faces. Is calculated as the difference between areas of
// individual triangles and the overall area of the face (which ifself is
// is the average of the areas of the individual triangles).
bool checkFaceFlatness
(
    const polyMeshGen& mesh,
    const bool report,
    const scalar warnFlatness,
    labelHashSet* setPtr,
    const boolList* changedFacePtr
)
{
    if( warnFlatness < 0 || warnFlatness > 1 )
    {
        FatalErrorIn
        (
            "bool checkFaceFlatness("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*)"
        )   << "warnFlatness should be [0..1] but is now " << warnFlatness
            << abort(FatalError);
    }

    const pointFieldPMG& points = mesh.points();
    const faceListPMG& faces = mesh.faces();
    const vectorField& fctrs = mesh.addressingData().faceCentres();

    // Areas are calculated as the sum of areas. (see
    // polyMeshGenAddressingFaceCentresAndAreas.C)
    scalarField magAreas(mag(mesh.addressingData().faceAreas()));

    label nWarped = 0;

    scalar minFlatness = GREAT;
    scalar sumFlatness = 0;
    label nSummed = 0;

    # ifdef USE_OMP
    # pragma omp parallel if( faces.size() > 1000 ) \
    reduction(+ : nSummed, nWarped) reduction(+ : sumFlatness)
    # endif
    {
        scalar minFlatnessProc = VGREAT;

        # ifdef USE_OMP
        # pragma omp for schedule(guided)
        # endif
        forAll(faces, faceI)
        {
            if( changedFacePtr && !(*changedFacePtr)[faceI] )
                continue;

            const face& f = faces[faceI];

            if( f.size() > 3 && magAreas[faceI] > VSMALL )
            {
                const point& fc = fctrs[faceI];

                //- Calculate the sum of magnitude of areas and compare
                //- the magnitude of sum of areas.

                scalar sumA = 0.0;

                forAll(f, fp)
                {
                    const point& thisPoint = points[f[fp]];
                    const point& nextPoint = points[f.nextLabel(fp)];

                    // Triangle around fc.
                    vector n = 0.5*((nextPoint - thisPoint)^(fc - thisPoint));
                    sumA += mag(n);
                }

                scalar flatness = magAreas[faceI] / (sumA+VSMALL);

                sumFlatness += flatness;
                ++nSummed;

                minFlatnessProc = Foam::min(minFlatnessProc, flatness);

                if( flatness < warnFlatness )
                {
                    ++nWarped;

                    if( setPtr )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        setPtr->insert(faceI);
                    }
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            minFlatness = Foam::min(minFlatness, minFlatnessProc);
        }
    }

    if( Pstream::parRun() && setPtr )
    {
        //- make sure that processor faces are marked on both sides
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        List<DynList<label> > markedFaces(procBoundaries.size());
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            for(label i=0;i<size;++i)
                if( setPtr->found(start+i) )
                    markedFaces[patchI].append(i);
        }

        //- exchange list sizes
        forAll(procBoundaries, patchI)
        {
            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                sizeof(label)
            );

            toOtherProc << markedFaces[patchI].size();
        }

        labelList nMarkedOnOtherProcs(procBoundaries.size());
        forAll(procBoundaries, patchI)
        {
            IPstream fromOtheProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                sizeof(label)
            );

            fromOtheProc >> nMarkedOnOtherProcs[patchI];
        }

        //- exchange data
        forAll(procBoundaries, patchI)
        {
            if( markedFaces[patchI].size() == 0 )
                continue;

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                markedFaces[patchI].byteSize()
            );

            toOtherProc << markedFaces[patchI];
        }

        forAll(procBoundaries, patchI)
        {
            if( nMarkedOnOtherProcs[patchI] == 0 )
                continue;

            labelList receivedData;
            IPstream fromOtheProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                nMarkedOnOtherProcs[patchI]*sizeof(label)
            );
            fromOtheProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();
            forAll(receivedData, i)
                if( !setPtr->found(start+receivedData[i]) )
                    setPtr->insert(start+receivedData[i]);
        }
    }

    reduce(nWarped, sumOp<label>());
    reduce(minFlatness, minOp<scalar>());

    reduce(nSummed, sumOp<label>());
    reduce(sumFlatness, sumOp<scalar>());

    if( report )
    {
        if( nSummed > 0 )
        {
            Info<< "Face flatness (1 = flat, 0 = butterfly) : average = "
                << sumFlatness / nSummed << "  min = " << minFlatness << endl;
        }

        if( nWarped> 0 )
        {
            Info<< "There are " << nWarped
                << " faces with ratio between projected and actual area < "
                << warnFlatness
                << ".\nMinimum ratio (minimum flatness, maximum warpage) = "
                << minFlatness << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the ratio between projected"
                << " and actual area is > " << warnFlatness << nl << endl;
        }
    }

    if( nWarped > 0 )
    {
        WarningIn
        (
            "bool checkFaceFlatness("
            "const polyMeshGen&, const bool, const scalar,"
            " labelHashSet*, const boolList*)"
        )   << nWarped  << " faces with severe warpage (flatness < "
            << warnFlatness << ") found.\n"
            << endl;

        return true;
    }
    else
    {
        return false;
    }
}

label findBadFacesAdditionalChecks
(
    const polyMeshGen& mesh,
    const bool report,
    labelHashSet& badFaces,
    const boolList* activeFacePtr

)
{
    //Default maximum non-orthogonality
    scalar maxNonOrtho = 65.0;

    //Default maximum skewness
    scalar maxSkewness = 2.0;

    //Default face uniformity
    scalar fcUniform = 0.1;

    //Default volume uniformity
    scalar volUniform = 0.1;

    //Default maximum face angle
    scalar maxAngle = 10;

    //Default minimum quality of tetrahedra
    scalar minTetQuality = VSMALL;

    //Default minimum face twist
    scalar minFaceTwist = VSMALL;

    //Default minimum cell determinant
    scalar minCellDeterminant = 1e-3;

    //Default minimum volume ratio
    scalar minVolRatio = 0.01;

    //Default minimum face triangle twist
    scalar minFaceTriangleTwist = VSMALL;

    //Default minimum area of a face
    scalar minimumFaceArea = VSMALL;

    //Default minimum volume of the face pyramid
    scalar minPyramidVolume = VSMALL;

    //Default face flatness
    scalar faceFlatness = 0.8;

    //Default minimum tetrahedral part of a cell
    scalar minCellPartTetrahedra = VSMALL;

    //Check whether quality criteria is specified by user
    //if( mesh_.returnTime().foundObject<IOdictionary>("meshDict") )
    if( mesh.returnTime().foundObject<IOdictionary>("meshDict") )
    {
        const dictionary& meshDict =
            mesh.returnTime().lookupObject<IOdictionary>("meshDict");

        if( meshDict.found("meshQualitySettings") )
        {
            const dictionary& qualityDict =
                meshDict.subDict("meshQualitySettings");

            //Reading maximum non-orthogonality defined by the user
            if( qualityDict.found("maxNonOrthogonality") )
            {
                maxNonOrtho =
                    readScalar
                    (
                        qualityDict.lookup("maxNonOrthogonality")
                    );

                polyMeshGenChecks::checkFaceDotProduct
                (
                    mesh,
                    report,
                    maxNonOrtho,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading maximum skewness defined by the user
            if( qualityDict.found("maxSkewness") )
            {
                maxSkewness =
                    readScalar
                    (
                        qualityDict.lookup("maxSkewness")
                    );

                polyMeshGenChecks::checkFaceSkewness
                (
                    mesh,
                    report,
                    maxSkewness,
                    &badFaces,
                    activeFacePtr
                );

            }

            //Reading face uniformity specified by the user
            if( qualityDict.found("fcUniform") )
            {
                fcUniform =
                    readScalar
                    (
                        qualityDict.lookup("fcUniform")
                    );

                polyMeshGenChecks::checkFaceUniformity
                (
                    mesh,
                    report,
                    fcUniform,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading volume uniformity specified by the user
            if( qualityDict.found("volUniform") )
            {
                volUniform = readScalar(qualityDict.lookup("volUniform"));

                polyMeshGenChecks::checkVolumeUniformity
                (
                    mesh,
                    report,
                    volUniform,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading maximum face angle specified by the user
            if( qualityDict.found("maxAngle") )
            {
                maxAngle = readScalar(qualityDict.lookup("maxAngle"));

                polyMeshGenChecks::checkFaceAngles
                (
                    mesh,
                    report,
                    maxAngle,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading minimum quality of tetrahedra defined by the user
            if( qualityDict.found("minTetQuality") )
            {
                minTetQuality = readScalar(qualityDict.lookup("minTetQuality"));

                polyMeshGenChecks::checkTetQuality
                (
                    mesh,
                    report,
                    minTetQuality,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading minimum face twist defined by the user
            if( qualityDict.found("minFaceTwist") )
            {
                minFaceTwist = readScalar(qualityDict.lookup("minFaceTwist"));

                polyMeshGenChecks::checkMinTwist
                (
                    mesh,
                    report,
                    minFaceTwist,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading minimum cell determinant defined by the user
            if( qualityDict.found("minCellDeterminant") )
            {
                minCellDeterminant =
                    readScalar(qualityDict.lookup("minCellDeterminant"));

                polyMeshGenChecks::checkCellDeterminant
                (
                    mesh,
                    report,
                    minCellDeterminant,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading minimum cell volume ratio defined by the user
            if( qualityDict.found("minVolRatio") )
            {
                minVolRatio = readScalar(qualityDict.lookup("minVolRatio"));

                polyMeshGenChecks::checkMinVolRatio
                (
                    mesh,
                    report,
                    minVolRatio,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading minimum face triangle twist defined by the user
            if( qualityDict.found("minFaceTriangleTwist") )
            {
                minFaceTriangleTwist =
                    readScalar(qualityDict.lookup("minFaceTriangleTwist"));

                polyMeshGenChecks::checkTriangleTwist
                (
                    mesh,
                    report,
                    minFaceTriangleTwist,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading minimum volume of the face pyramid defined by the user
            if( qualityDict.found("minPyramidVolume") )
            {
                minPyramidVolume =
                    readScalar(qualityDict.lookup("minPyramidVolume"));

                polyMeshGenChecks::checkFacePyramids
                (
                    mesh,
                    report,
                    minPyramidVolume,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading minimum area of a face defined by the user
            if( qualityDict.found("minimumFaceArea") )
            {
                minimumFaceArea =
                    readScalar(qualityDict.lookup("minimumFaceArea"));

                polyMeshGenChecks::checkFaceAreas
                (
                    mesh,
                    report,
                    minimumFaceArea,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading face flatness defined by the user
            if( qualityDict.found("faceFlatness") )
            {
                faceFlatness = readScalar(qualityDict.lookup("faceFlatness"));

                polyMeshGenChecks::checkFaceFlatness
                (
                    mesh,
                    report,
                    faceFlatness,
                    &badFaces,
                    activeFacePtr
                );
            }

            //Reading minimum tetrahedral part of a cell defined by the user
            if( qualityDict.found("minCellPartTetrahedra") )
            {
                minCellPartTetrahedra =
                    readScalar(qualityDict.lookup("minCellPartTetrahedra"));

                polyMeshGenChecks::checkCellPartTetrahedra
                (
                    mesh,
                    report,
                    minCellPartTetrahedra,
                    &badFaces,
                    activeFacePtr
                );

            }
        }
    }

    label nBadFaces = badFaces.size();
    reduce(nBadFaces, sumOp<label>());
    return nBadFaces;
}

label findBadFacesRelaxed
(
    const polyMeshGen& mesh,
    labelHashSet& badFaces,
    const bool report,
    const boolList* activeFacePtr
)
{
    badFaces.clear();

    polyMeshGenChecks::findBadFacesAdditionalChecks
    (
        mesh,
        report,
        badFaces,
        activeFacePtr
    );

    polyMeshGenChecks::checkFacePyramids
    (
        mesh,
        report,
        VSMALL,
        &badFaces,
        activeFacePtr
    );

    polyMeshGenChecks::checkFaceAreas
    (
        mesh,
        report,
        VSMALL,
        &badFaces,
        activeFacePtr
    );

    const label nBadFaces = returnReduce(badFaces.size(), sumOp<label>());

    return nBadFaces;
}

label findBadFaces
(
    const polyMeshGen& mesh,
    labelHashSet& badFaces,
    const bool report,
    const boolList* activeFacePtr
)
{
    badFaces.clear();

    polyMeshGenChecks::findBadFacesAdditionalChecks
    (
        mesh,
        report,
        badFaces,
        activeFacePtr
    );

    polyMeshGenChecks::checkFacePyramids
    (
        mesh,
        report,
        VSMALL,
        &badFaces,
        activeFacePtr
    );

    polyMeshGenChecks::checkFaceFlatness
    (
        mesh,
        report,
        0.8,
        &badFaces,
        activeFacePtr
    );

    polyMeshGenChecks::checkCellPartTetrahedra
    (
        mesh,
        report,
        VSMALL,
        &badFaces,
        activeFacePtr
    );

    polyMeshGenChecks::checkFaceAreas
    (
        mesh,
        report,
        VSMALL,
        &badFaces,
        activeFacePtr
    );

    const label nBadFaces = returnReduce(badFaces.size(), sumOp<label>());

    return nBadFaces;
}

label findLowQualityFaces
(
    const polyMeshGen& mesh,
    labelHashSet& badFaces,
    const bool report,
    const boolList* activeFacePtr
)
{
    badFaces.clear();

    polyMeshGenChecks::findBadFacesAdditionalChecks
    (
        mesh,
        report,
        badFaces,
        activeFacePtr
    );

    polyMeshGenChecks::checkFaceDotProduct
    (
        mesh,
        report,
        65.0,
        &badFaces,
        activeFacePtr
    );

    polyMeshGenChecks::checkFaceSkewness
    (
        mesh,
        report,
        2.0,
        &badFaces,
        activeFacePtr
    );

    const label nBadFaces = returnReduce(badFaces.size(), sumOp<label>());

    return nBadFaces;
}

label findWorstQualityFaces
(
    const polyMeshGen& mesh,
    labelHashSet& badFaces,
    const bool /*report*/,
    const boolList* activeFacePtr,
    const scalar relativeThreshold
)
{
    badFaces.clear();

    scalarField checkValues;
    polyMeshGenChecks::checkFaceDotProduct
    (
        mesh,
        checkValues,
        activeFacePtr
    );

    scalar minNonOrtho = returnReduce(min(checkValues), minOp<scalar>());
    const scalar warnNonOrtho =
        minNonOrtho + relativeThreshold * (1.0 - minNonOrtho);

    Info << "Worst non-orthogonality " << Foam::acos(minNonOrtho) * 180.0 / M_PI
         << " selecting faces with non-orthogonality greater than "
         << (Foam::acos(warnNonOrtho) * 180.0 / M_PI) << endl;

    forAll(checkValues, faceI)
    {
        if
        (
            activeFacePtr && activeFacePtr->operator[](faceI) &&
            checkValues[faceI] < warnNonOrtho
        )
            badFaces.insert(faceI);
    }

    polyMeshGenChecks::checkFaceSkewness
    (
        mesh,
        checkValues,
        activeFacePtr
    );

    const scalar maxSkew = returnReduce(max(checkValues), maxOp<scalar>());
    const scalar warnSkew = (1.0 - relativeThreshold) * maxSkew;
    forAll(checkValues, faceI)
        if
        (
            activeFacePtr && activeFacePtr->operator[](faceI) &&
            checkValues[faceI] > warnSkew
        )
            badFaces.insert(faceI);

    Info << "Maximum skewness in the mesh is " << maxSkew
         << " selecting faces with skewness greater than " << warnSkew << endl;

    const label nBadFaces = returnReduce(badFaces.size(), sumOp<label>());

    Info << "Selected " << nBadFaces
         << " out of " << returnReduce(checkValues.size(), sumOp<label>())
         << " faces" << endl;

    return nBadFaces;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace polyMeshGenChecks

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


