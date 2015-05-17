/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "immersedBoundaryFvPatch.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::FieldField<Foam::Field, Type> >
Foam::immersedBoundaryFvPatch::sendAndReceive
(
    const Field<Type>& psi
) const
{
    tmp<FieldField<Field, Type> > tprocPsi
    (
        new FieldField<Field, Type>(Pstream::nProcs())
    );
    FieldField<Field, Type>& procPsi = tprocPsi();

    forAll (procPsi, procI)
    {
        procPsi.set
        (
            procI,
            new Field<Type>
            (
                ibProcCentres()[procI].size(),
                pTraits<Type>::zero
            )
        );
    }

    if (Pstream::parRun())
    {
        // Send
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Do not send empty lists
                if (!ibProcCells()[procI].empty())
                {
                    Field<Type> curPsi(psi, ibProcCells()[procI]);

                    // Parallel data exchange
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        curPsi.size()*sizeof(Type)
                    );

                    toProc << curPsi;
                }
            }
        }

        // Receive
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Do not receive empty lists
                if (!procPsi[procI].empty())
                {
                    // Parallel data exchange
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        procPsi[procI].size()*sizeof(Type)
                    );

                    fromProc >> procPsi[procI];
                }
            }
        }
    }

    return tprocPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toIbPoints
(
    const Field<Type>& triValues
) const
{
    if (triValues.size() != ibMesh().size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatch::toIbPoints\n"
            "(\n"
            "    const Field<Type>& triValues\n"
            ") const"
        )   << "Field size does not correspond to size of immersed boundary "
            << "triangulated surface for patch " << name() << nl
            << "Field size = " << triValues.size()
            << " surface size = " << ibMesh().size()
            << abort(FatalError);
    }

    const labelList& ibc = ibCells();

    tmp<Field<Type> > tIbPsi
    (
        new Field<Type>(ibc.size(), pTraits<Type>::zero)
    );
    Field<Type>& ibPsi = tIbPsi();

    const labelList& hf = hitFaces();

    // Assuming triSurface data is on triangles
    forAll (ibPsi, cellI)
    {
        ibPsi[cellI] = triValues[hf[cellI]];
    }

//     const vectorField& p = ibPoints();
//     const List<labelledTri>& faces = ibMesh();
//     const vectorField& triPoints = ibMesh().points();

//     // Assuming triSurface data is on vertices
//     forAll (ibPsi, cellI)
//     {
//         const labelledTri& tri = faces[hf[cellI]];
//         triPointRef triPt = faces[hf[cellI]].tri(triPoints);

//         ibPsi[cellI] =
//             triValues[tri[0]]*triPt.Ni(0, p[cellI])
//           + triValues[tri[1]]*triPt.Ni(1, p[cellI])
//           + triValues[tri[2]]*triPt.Ni(2, p[cellI]);
//     }

    return tIbPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toIbPoints
(
    const tmp<Field<Type> >& ttriValues
) const
{
    tmp<Field<Type> > tint = toIbPoints(ttriValues());
    ttriValues.clear();
    return tint;

}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toTriFaces
(
    const Field<Type>& ibValues
) const
{
    if (ibValues.size() != ibCells().size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatch::toTriFaces\n"
            "(\n"
            "    const Field<Type>& ibValues\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "triangulated surface for patch " << name() << nl
            << "Field size = " << ibValues.size()
            << " IB points size = " << ibCells().size()
            << abort(FatalError);
    }

    const labelListList& ctfAddr = cellsToTriAddr();
    const scalarListList& ctfWeights = cellsToTriWeights();

    tmp<Field<Type> > tIbPsi
    (
        new Field<Type>(ctfAddr.size(), pTraits<Type>::zero)
    );
    Field<Type>& ibPsi = tIbPsi();

    // Do interpolation
    forAll (ctfAddr, triI)
    {
        const labelList& curAddr = ctfAddr[triI];
        const scalarList& curWeights = ctfWeights[triI];

        forAll (curAddr, cellI)
        {
            ibPsi[triI] += curWeights[cellI]*ibValues[curAddr[cellI]];
        }
    }

    return tIbPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toTriFaces
(
    const tmp<Field<Type> >& tibValues
) const
{
    tmp<Field<Type> > tint = toTriFaces(tibValues());
    tibValues.clear();
    return tint;

}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toSamplingPoints
(
    const Field<Type>& cellValues
) const
{
    if (cellValues.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatch::toSamplingPoints\n"
            "(\n"
            "    const Field<Type>& cellValues\n"
            ") const"
        )   << "Field size does not correspond to cell centres "
            << "for patch " << name() << nl
            << "Field size = " << cellValues.size()
            << " nCells = " << mesh_.nCells()
            << abort(FatalError);
    }

    // Get addressing
    const labelList& ibc = ibCells();
    const labelListList& ibcc = ibCellCells();
    const List<List<labelPair> >& ibcProcC = ibCellProcCells();

    // Get weights
    const scalarListList& cellWeights = ibSamplingWeights();
    const scalarListList& cellProcWeights = ibSamplingProcWeights();

    tmp<Field<Type> > tIbPsi
    (
        new Field<Type>(ibc.size(), pTraits<Type>::zero)
    );
    Field<Type>& ibPsi = tIbPsi();

    // Do interpolation, local cell data
    forAll (ibc, cellI)
    {
        const labelList& curAddr = ibcc[cellI];
        const scalarList& curWeights = cellWeights[cellI];

        forAll (curAddr, ccI)
        {
            ibPsi[cellI] += curWeights[ccI]*cellValues[curAddr[ccI]];
        }
    }

    // Parallel communication for psi
    FieldField<Field, Type> procCellValues = sendAndReceive(cellValues);

    // Do interpolation, cell data from other processors
    forAll (ibc, cellI)
    {
        const List<labelPair>& curProcCells = ibcProcC[cellI];
        const scalarList& curProcWeights = cellProcWeights[cellI];

        forAll (curProcCells, cpcI)
        {
            ibPsi[cellI] +=
                curProcWeights[cpcI]*
                procCellValues
                [
                    curProcCells[cpcI].first()
                ]
                [
                    curProcCells[cpcI].second()
                ];
        }
    }

    return tIbPsi;
}


// ************************************************************************ //
