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

Class
    meshOps

Description
    Various utility functions that perform mesh-related operations.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "Pair.H"
#include "meshOps.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace meshOps
{

// Utility method to find the common edge between two faces.
inline bool findCommonEdge
(
    const label first,
    const label second,
    const UList<labelList>& faceEdges,
    label& common
)
{
    bool found = false;

    const labelList& fEi = faceEdges[first];
    const labelList& fEj = faceEdges[second];

    forAll(fEi, edgeI)
    {
        forAll(fEj, edgeJ)
        {
            if (fEi[edgeI] == fEj[edgeJ])
            {
                common = fEi[edgeI];

                found = true;
                break;
            }
        }

        if (found)
        {
            break;
        }
    }

    return found;
}


// Utility method to find the interior (quad) / boundary (tri) faces
// for an input quad-face and adjacent triangle-prism cell.
inline void findPrismFaces
(
    const label fIndex,
    const label cIndex,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& neighbour,
    FixedList<face,2>& bdyf,
    FixedList<label,2>& bidx,
    FixedList<face,2>& intf,
    FixedList<label,2>& iidx
)
{
    label indexO = 0, indexI = 0;

    const cell& c = cells[cIndex];

    forAll(c, i)
    {
        label faceIndex = c[i];

        // Don't count the face under consideration
        if (faceIndex != fIndex)
        {
            const face& fi = faces[faceIndex];

            if (neighbour[faceIndex] == -1)
            {
                if (fi.size() == 3)
                {
                    // Triangular face on the boundary
                    bidx[indexO] = faceIndex;
                    bdyf[indexO++] = fi;
                }
                else
                {
                    // This seems to be a non-triangular face on the boundary
                    // Consider this as "interior" and move on
                    iidx[indexI] = faceIndex;
                    intf[indexI++] = fi;
                }
            }
            else
            {
                // Face on the interior
                iidx[indexI] = faceIndex;
                intf[indexI++] = fi;
            }
        }
    }
}


// Utility method to find the isolated point given two triangular faces.
//  - Returns the point on checkFace that does not belong to baseFace.
inline label findIsolatedPoint
(
    const face& baseFace,
    const face& checkFace
)
{
    // Get the fourth point
    forAll(checkFace, pointI)
    {
        if
        (
            checkFace[pointI] != baseFace[0] &&
            checkFace[pointI] != baseFace[1] &&
            checkFace[pointI] != baseFace[2]
        )
        {
            return checkFace[pointI];
        }
    }

    return -1;
}


// Utility method to find the isolated point on a triangular face
// that doesn't lie on the specified edge. Also returns the point next to it.
inline void findIsolatedPoint
(
    const face& f,
    const edge& e,
    label& ptIndex,
    label& nextPtIndex
)
{
    // Check the first point
    if ( f[0] != e.start() && f[0] != e.end() )
    {
        ptIndex = f[0];
        nextPtIndex = f[1];
        return;
    }

    // Check the second point
    if ( f[1] != e.start() && f[1] != e.end() )
    {
        ptIndex = f[1];
        nextPtIndex = f[2];
        return;
    }

    // Check the third point
    if ( f[2] != e.start() && f[2] != e.end() )
    {
        ptIndex = f[2];
        nextPtIndex = f[0];
        return;
    }

    // This bit should never happen.
    FatalErrorIn
    (
        "inline void meshOps::findIsolatedPoint"
        "(const face&, const edge&, label&, label&)"
    )
        << "Cannot find isolated point in face " << f << endl
        << " Using edge: " << e
        << abort(FatalError);
}


// Given a face and cell index, find the apex point for a tet cell.
inline label tetApexPoint
(
    const label cIndex,
    const label fIndex,
    const UList<face>& faces,
    const UList<cell>& cells
)
{
    label apexPoint = -1;
    bool foundApex = false;

    const cell& cellToCheck = cells[cIndex];
    const face& baseFace = faces[fIndex];

    forAll(cellToCheck, faceI)
    {
        const face& faceToCheck = faces[cellToCheck[faceI]];

        apexPoint = findIsolatedPoint(baseFace, faceToCheck);

        if (apexPoint > -1)
        {
            foundApex = true;
            break;
        }
    }

    if (!foundApex)
    {
        Pout<< "fIndex: " << fIndex << ":: " << faces[fIndex] << endl;
        Pout<< "cIndex: " << cIndex << ":: " << cellToCheck << endl;

        forAll(cellToCheck, faceI)
        {
            Pout<< '\t' << cellToCheck[faceI] << ":: "
                << faces[cellToCheck[faceI]] << endl;
        }

        FatalErrorIn
        (
            "\n\n"
            "inline label meshOps::tetApexPoint\n"
            "(\n"
            "    const label cIndex,\n"
            "    const label fIndex,\n"
            "    const UList<face>& faces,\n"
            "    const UList<cell>& cells\n"
            ")\n"
        )
            << "Could not find an apex point in cell: " << cIndex
            << " given face: " << fIndex
            << abort(FatalError);
    }

    return apexPoint;
}


// Given a cell index, find the centroid / volume of a cell.
//  - If checking is enabled, return whether the cell is closed
inline void cellCentreAndVolume
(
    const label cIndex,
    const vectorField& points,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& owner,
    vector& centre,
    scalar& volume
)
{
    // Reset inputs
    volume = 0.0;
    centre = vector::zero;

    const cell& cellToCheck = cells[cIndex];

    // Average face-centres to get an estimate centroid
    vector cEst(vector::zero), fCentre(vector::zero);

    forAll(cellToCheck, faceI)
    {
        const face& checkFace = faces[cellToCheck[faceI]];

        cEst += checkFace.centre(points);
    }

    cEst /= cellToCheck.size();

    forAll(cellToCheck, faceI)
    {
        const face& checkFace = faces[cellToCheck[faceI]];

        if (checkFace.size() == 3)
        {
            tetPointRef tpr
            (
                points[checkFace[2]],
                points[checkFace[1]],
                points[checkFace[0]],
                cEst
            );

            scalar tetVol = tpr.mag();

            // Flip if necessary
            if (owner[cellToCheck[faceI]] != cIndex)
            {
                tetVol *= -1.0;
            }

            // Accumulate volume-weighted tet centre
            centre += tetVol*tpr.centre();

            // Accumulate tet volume
            volume += tetVol;
        }
        else
        {
            forAll(checkFace, pI)
            {
                tetPointRef tpr
                (
                    points[checkFace[pI]],
                    points[checkFace.prevLabel(pI)],
                    checkFace.centre(points),
                    cEst
                );

                scalar tetVol = tpr.mag();

                // Flip if necessary
                if (owner[cellToCheck[faceI]] != cIndex)
                {
                    tetVol *= -1.0;
                }

                // Accumulate volume-weighted tet centre
                centre += tetVol*tpr.centre();

                // Accumulate tet volume
                volume += tetVol;
            }
        }
    }

    centre /= volume + VSMALL;
}


// Determine whether a segment intersects a triangular face
inline bool segmentTriFaceIntersection
(
    const triPointRef& faceToCheck,
    const linePointRef& edgeToCheck,
    vector& intPoint
)
{
    // Fetch references
    const vector& p1 = edgeToCheck.start();
    const vector& p2 = edgeToCheck.end();

    // Define face-normal
    vector n = faceToCheck.normal();
    n /= mag(n) + VSMALL;

    // Compute uValue
    scalar numerator = n & (faceToCheck.a() - p1);
    scalar denominator = n & (p2 - p1);

    // Check if the edge is parallel to the face
    if (mag(denominator) < VSMALL)
    {
        return false;
    }

    scalar u = (numerator / denominator);

    // Check for intersection along line.
    if ((u >= 0.0) && (u <= 1.0))
    {
        // Compute point of intersection
        intPoint = p1 + u*(p2 - p1);

        // Also make sure that intPoint lies within the face
        if (pointInTriFace(faceToCheck, intPoint, false))
        {
            return true;
        }
    }

    // Failed to fall within edge-bounds, or within face
    return false;
}


// Determine whether the particular point lies
// inside the given triangular face
inline bool pointInTriFace
(
    const triPointRef& faceToCheck,
    const vector& cP,
    bool testCoplanarity
)
{
    // Copy inputs
    const vector& a = faceToCheck.a();
    const vector& b = faceToCheck.b();
    const vector& c = faceToCheck.c();

    // Compute the normal
    vector nf = faceToCheck.normal();

    if ( ((0.5 * ((b - a)^(cP - a))) & nf) < 0.0)
    {
        return false;
    }

    if ( ((0.5 * ((c - b)^(cP - b))) & nf) < 0.0)
    {
        return false;
    }

    if ( ((0.5 * ((a - c)^(cP - c))) & nf) < 0.0)
    {
        return false;
    }

    if (testCoplanarity)
    {
        vector ftp(a - cP);

        // Normalize vectors
        nf /= mag(nf) + VSMALL;
        ftp /= mag(ftp) + VSMALL;

        if (mag(ftp & nf) > VSMALL)
        {
            return false;
        }
    }

    // Passed test with all edges
    return true;
}


// Parallel blocking send
inline void pWrite
(
    const label toID,
    const label& data,
    bool blocking
)
{
    OPstream::write
    (
        (blocking ? Pstream::blocking : Pstream::nonBlocking),
        toID,
        reinterpret_cast<const char*>(&data),
        sizeof(label)
    );
}


// Parallel blocking receive
inline void pRead
(
    const label fromID,
    label& data,
    bool blocking
)
{
    IPstream::read
    (
        (blocking ? Pstream::blocking : Pstream::nonBlocking),
        fromID,
        reinterpret_cast<char*>(&data),
        sizeof(label)
    );
}


// Parallel non-blocking send for fixed lists
template <class Type, unsigned Size>
inline void pWrite
(
    const label toID,
    const FixedList<Type, Size>& data
)
{
    OPstream::write
    (
        Pstream::blocking,
        toID,
        reinterpret_cast<const char*>(&data[0]),
        Size*sizeof(Type)
    );
}


// Parallel non-blocking receive for fixed lists
template <class Type, unsigned Size>
inline void pRead
(
    const label fromID,
    FixedList<Type, Size>& data
)
{
    IPstream::read
    (
        Pstream::blocking,
        fromID,
        reinterpret_cast<char*>(data.begin()),
        Size*sizeof(Type)
    );
}


// Parallel non-blocking send for lists
template <class Type>
inline void pWrite
(
    const label toID,
    const UList<Type>& data
)
{
    OPstream::write
    (
        Pstream::nonBlocking,
        toID,
        reinterpret_cast<const char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Parallel non-blocking receive for lists
template <class Type>
inline void pRead
(
    const label fromID,
    UList<Type>& data
)
{
    IPstream::read
    (
        Pstream::nonBlocking,
        fromID,
        reinterpret_cast<char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Wait for buffer transfer completion.
inline void waitForBuffers()
{
    if (Pstream::parRun())
    {
        OPstream::waitRequests();
        IPstream::waitRequests();
    }
}


// Method to insert a label between two labels in a list
// Assumes that all labels are unique.
inline void insertLabel
(
    const label newLabel,
    const label labelA,
    const label labelB,
    labelList& list
)
{
    // Create a new list
    bool found = false;
    label origSize = list.size();
    labelList newList(origSize + 1);

    label index = 0, nextI = -1;

    // Start a linear search
    forAll(list, itemI)
    {
        newList[index++] = list[itemI];

        nextI = list.fcIndex(itemI);

        if
        (
            (
                (list[itemI] == labelA && list[nextI] == labelB) ||
                (list[itemI] == labelB && list[nextI] == labelA)
            ) &&
           !found
        )
        {
            found = true;
            newList[index++] = newLabel;
        }
    }

    if (!found)
    {
        FatalErrorIn
        (
            "inline void meshOps::insertLabel"
            "(const label, const label, const label, labelList&)"
        )   << nl << "Cannot insert " << newLabel
            << " in list: " << list << nl
            << " Labels: "
            << labelA << " and " << labelB
            << " were not found in sequence."
            << abort(FatalError);
    }

    // Transfer the list
    list.transfer(newList);
}


// Utility method to replace a label in a given list
inline void replaceLabel
(
     const label original,
     const label replacement,
     labelList& list
)
{
    label index = -1;

    if ((index = findIndex(list, original)) > -1)
    {
        list[index] = replacement;
    }
    else
    {
        FatalErrorIn
        (
            "inline void label meshOps::replaceLabel"
            "(const label, const label, labelList&)"
        )   << nl << "Cannot find " << original
            << " in list: " << list << nl
            << " Label: " << replacement
            << " was not used in replacement."
            << abort(FatalError);
    }
}


// Utility method to size-up the list to include an item
template <class Type>
inline void sizeUpList
(
    const Type& item,
    List<Type>& list
)
{
    list.setSize(list.size() + 1, item);
}


// Utility method to size-down the list to remove an item
template <class Type>
inline void sizeDownList
(
    const Type& item,
    List<Type>& list
)
{
    label index = -1;

    if ((index = findIndex(list, item)) > -1)
    {
        meshOps::removeIndex(index, list);
    }
    else
    {
        FatalErrorIn
        (
            "inline void meshOps::sizeDownList"
            "(const Type& item, List<Type>& list)"
        )
            << nl << "Item: " << item
            << " was not found in list. " << nl
            << " List: " << nl << list
            << abort(FatalError);
    }
}


// Remove an item at a particular index in the list
template <class Type>
inline void removeIndex
(
    const label index,
    List<Type>& list
)
{
    // Create a new list
    List<Type> newList(list.size() - 1);

    // Copy individual items
    label n = 0;

    forAll(list, itemI)
    {
        if (itemI == index)
        {
            continue;
        }

        newList[n++] = list[itemI];
    }

    // Overwrite
    list.transfer(newList);
}


} // End namespace meshOps


} // End namespace Foam

// ************************************************************************* //
