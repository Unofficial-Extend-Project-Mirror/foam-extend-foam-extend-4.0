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

#include "indexedOctree.H"
#include "linePointRef.H"
#include "triSurface.H"
#include "meshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class Type>
Foam::scalar Foam::indexedOctree<Type>::perturbTol_ = 10*SMALL;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Does bb intersect a sphere around sample? Or is any corner point of bb
// closer than nearestDistSqr to sample.
template <class Type>
bool indexedOctree<Type>::intersects
(
    const point& p0,
    const point& p1,
    const scalar nearestDistSqr,
    const point& sample
)
{
    // Find out where sample is in relation to bb.
    // Find nearest point on bb.
    scalar distSqr = 0;

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        scalar d0 = p0[dir] - sample[dir];
        scalar d1 = p1[dir] - sample[dir];

        if ((d0 > 0) != (d1 > 0))
        {
            // sample inside both extrema. This component does not add any
            // distance.
        }
        else if (mag(d0) < mag(d1))
        {
            distSqr += d0*d0;
        }
        else
        {
            distSqr += d1*d1;
        }

        if (distSqr > nearestDistSqr)
        {
            return false;
        }
    }

    return true;
}


// Does bb intersect a sphere around sample? Or is any corner point of bb
// closer than nearestDistSqr to sample.
template <class Type>
bool indexedOctree<Type>::intersects
(
    const treeBoundBox& parentBb,
    const direction octant,
    const scalar nearestDistSqr,
    const point& sample
)
{
    //- Speeded up version of
    //     treeBoundBox subBb(parentBb.subBbox(mid, octant))
    //     intersects
    //     (
    //          subBb.min(),
    //          subBb.max(),
    //          nearestDistSqr,
    //          sample
    //     )

    const point& min = parentBb.min();
    const point& max = parentBb.max();

    point other;

    if (octant & treeBoundBox::RIGHTHALF)
    {
        other.x() = max.x();
    }
    else
    {
        other.x() = min.x();
    }

    if (octant & treeBoundBox::TOPHALF)
    {
        other.y() = max.y();
    }
    else
    {
        other.y() = min.y();
    }

    if (octant & treeBoundBox::FRONTHALF)
    {
        other.z() = max.z();
    }
    else
    {
        other.z() = min.z();
    }

    const point mid(0.5*(min+max));

    return intersects(mid, other, nearestDistSqr, sample);
}


//
// Construction helper routines
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

// Split list of indices into 8 bins
template <class Type>
void indexedOctree<Type>::divide
(
    const labelList& indices,
    const treeBoundBox& bb,
    labelListList& result
) const
{
    List<DynamicList<label> > subIndices(8);
    for (direction octant = 0; octant < subIndices.size(); octant++)
    {
        subIndices[octant].setSize(indices.size()/8);
    }

    // Precalculate bounding boxes.
    FixedList<treeBoundBox, 8> subBbs;
    for (direction octant = 0; octant < subBbs.size(); octant++)
    {
        subBbs[octant] = bb.subBbox(octant);
    }

    forAll(indices, i)
    {
        label shapeI = indices[i];

        for (direction octant = 0; octant < 8; octant++)
        {
            if (shapes_.overlaps(shapeI, subBbs[octant]))
            {
                subIndices[octant].append(shapeI);
            }
        }
    }

    result.setSize(8);
    for (direction octant = 0; octant < subIndices.size(); octant++)
    {
        subIndices[octant].shrink();
        result[octant].transfer(subIndices[octant]);
        subIndices[octant].clear();
    }
}


// Subdivide the (content) node.
template <class Type>
typename indexedOctree<Type>::node indexedOctree<Type>::divide
(
    const treeBoundBox& bb,
    DynamicList<labelList>& contents,
    const label contentI
) const
{
    const labelList& indices = contents[contentI];

    node nod;

    if
    (
        bb.min()[0] >= bb.max()[0]
     || bb.min()[1] >= bb.max()[1]
     || bb.min()[2] >= bb.max()[2]
    )
    {
        FatalErrorIn("indexedOctree<Type>::divide")
            << "Badly formed bounding box:" << bb
            << abort(FatalError);
    }

    nod.bb_ = bb;
    nod.parent_ = -1;

    labelListList dividedIndices(8);
    divide(indices, bb, dividedIndices);

    // Have now divided the indices into 8 (possibly empty) subsets.
    // Replace current contentI with the first (non-empty) subset.
    // Append the rest.
    bool replaced = false;

    for (direction octant = 0; octant < dividedIndices.size(); octant++)
    {
        labelList& subIndices = dividedIndices[octant];

        if (subIndices.size() > 0)
        {
            if (!replaced)
            {
                contents[contentI].transfer(subIndices);
                nod.subNodes_[octant] = contentPlusOctant(contentI, octant);
                replaced = true;
            }
            else
            {
                // Store at end of contents.
                // note dummy append + transfer trick
                label sz = contents.size();
                contents.append(labelList(0));
                contents[sz].transfer(subIndices);
                nod.subNodes_[octant] = contentPlusOctant(sz, octant);
            }
        }
        else
        {
            // Mark node as empty
            nod.subNodes_[octant] = emptyPlusOctant(octant);
        }
    }

    return nod;
}


// Split any contents node with more than minSize elements.
template <class Type>
void indexedOctree<Type>::splitNodes
(
    const label minSize,
    DynamicList<indexedOctree<Type>::node>& nodes,
    DynamicList<labelList>& contents
) const
{
    label currentSize = nodes.size();

    // Take care to loop only over old nodes.
    // Also we loop over the same DynamicList which gets modified and
    // moved so make sure not to keep any references!
    for (label nodeI = 0; nodeI < currentSize; nodeI++)
    {
        for
        (
            direction octant = 0;
            octant < nodes[nodeI].subNodes_.size();
            octant++
        )
        {
            labelBits index = nodes[nodeI].subNodes_[octant];

            if (isNode(index))
            {
                // tree node. Leave intact.
            }
            else if (isContent(index))
            {
                label contentI = getContent(index);

                if (contents[contentI].size() > minSize)
                {
                    // Create node for content.

                    // Find the bounding box for the subnode
                    const node& nod = nodes[nodeI];
                    const treeBoundBox bb(nod.bb_.subBbox(octant));

                    node subNode(divide(bb, contents, contentI));
                    subNode.parent_ = nodeI;
                    label sz = nodes.size();
                    nodes.append(subNode);
                    nodes[nodeI].subNodes_[octant] = nodePlusOctant(sz, octant);
                }
            }
        }
    }
}


// Reorder contents to be in same order as nodes. Returns number of nodes on
// the compactLevel.
template <class Type>
label indexedOctree<Type>::compactContents
(
    DynamicList<node>& nodes,
    DynamicList<labelList>& contents,
    const label compactLevel,
    const label nodeI,
    const label level,

    List<labelList>& compactedContents,
    label& compactI
)
{
    const node& nod = nodes[nodeI];

    label nNodes = 0;

    if (level < compactLevel)
    {
        for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
        {
            labelBits index = nod.subNodes_[octant];

            if (isNode(index))
            {
                nNodes += compactContents
                (
                    nodes,
                    contents,
                    compactLevel,
                    getNode(index),
                    level+1,
                    compactedContents,
                    compactI
                );
            }
        }
    }
    else if (level == compactLevel)
    {
        // Compact all content on this level
        for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
        {
            labelBits index = nod.subNodes_[octant];

            if (isContent(index))
            {
                label contentI = getContent(index);

                compactedContents[compactI].transfer(contents[contentI]);

                // Subnode is at compactI. Adapt nodeI to point to it
                nodes[nodeI].subNodes_[octant] =
                    contentPlusOctant(compactI, octant);

                compactI++;
            }
            else if (isNode(index))
            {
                nNodes++;
            }
        }
    }
    return nNodes;
}


// Pre-calculates wherever possible the volume status per node/subnode.
// Recurses to determine status of lowest level boxes. Level above is
// combination of octants below.
template <class Type>
typename indexedOctree<Type>::volumeType indexedOctree<Type>::calcVolumeType
(
    const label nodeI
) const
{
    const node& nod = nodes_[nodeI];

    volumeType myType = UNKNOWN;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        volumeType subType;

        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            // tree node. Recurse.
            subType = calcVolumeType(getNode(index));
        }
        else if (isContent(index))
        {
            // Contents. Depending on position in box might be on either
            // side.
            subType = MIXED;
        }
        else
        {
            // No data in this octant. Set type for octant acc. to the mid
            // of its bounding box.
            const treeBoundBox subBb = nod.bb_.subBbox(octant);

            subType = volumeType(shapes_.getVolumeType(*this, subBb.mid()));
        }

        // Store octant type
        nodeTypes_.set((nodeI<<3)+octant, subType);

        // Combine sub node types into type for treeNode. Result is 'mixed' if
        // types differ among subnodes.
        if (myType == UNKNOWN)
        {
            myType = subType;
        }
        else if (subType != myType)
        {
            myType = MIXED;
        }
    }
    return myType;
}


template <class Type>
typename indexedOctree<Type>::volumeType indexedOctree<Type>::getVolumeType
(
    const label nodeI,
    const point& sample
) const
{
    const node& nod = nodes_[nodeI];

    direction octant = nod.bb_.subOctant(sample);

    volumeType octantType = volumeType(nodeTypes_.get((nodeI<<3)+octant));

    if (octantType == INSIDE)
    {
        return octantType;
    }
    else if (octantType == OUTSIDE)
    {
        return octantType;
    }
    else if (octantType == UNKNOWN)
    {
        // Can happen for e.g. non-manifold surfaces.
        return octantType;
    }
    else if (octantType == MIXED)
    {
        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            // Recurse
            volumeType subType = getVolumeType(getNode(index), sample);

            return subType;
        }
        else if (isContent(index))
        {
            // Content. Defer to shapes.
            return volumeType(shapes_.getVolumeType(*this, sample));
        }
        else
        {
            // Empty node. Cannot have 'mixed' as its type since not divided
            // up and has no items inside it.
            FatalErrorIn
            (
                "indexedOctree<Type>::getVolumeType"
                "(const label, const point&)"
            )   << "Sample:" << sample << " node:" << nodeI
                << " with bb:" << nodes_[nodeI].bb_ << nl
                << "Empty subnode has invalid volume type MIXED."
                << abort(FatalError);

            return UNKNOWN;
        }
    }
    else
    {
        FatalErrorIn
        (
            "indexedOctree<Type>::getVolumeType"
            "(const label, const point&)"
        )   << "Sample:" << sample << " at node:" << nodeI
            << " octant:" << octant
            << " with bb:" << nod.bb_.subBbox(octant) << nl
            << "Node has invalid volume type " << octantType
            << abort(FatalError);

        return UNKNOWN;
    }
}


template <class Type>
typename indexedOctree<Type>::volumeType indexedOctree<Type>::getSide
(
    const vector& outsideNormal,
    const vector& vec
)
{
    if ((outsideNormal&vec) >= 0)
    {
        return OUTSIDE;
    }
    else
    {
        return INSIDE;
    }
}


//
// Query routines
// ~~~~~~~~~~~~~~
//

// Find nearest point starting from nodeI
template <class Type>
void indexedOctree<Type>::findNearest
(
    const label nodeI,
    const point& sample,

    scalar& nearestDistSqr,
    label& nearestShapeI,
    point& nearestPoint
) const
{
    const node& nod = nodes_[nodeI];

    // Determine order to walk through octants
    FixedList<direction, 8> octantOrder;
    nod.bb_.searchOrder(sample, octantOrder);

    // Go into all suboctants (one containing sample first) and update nearest.
    for (direction i = 0; i < 8; i++)
    {
        direction octant = octantOrder[i];

        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            label subNodeI = getNode(index);

            const treeBoundBox& subBb = nodes_[subNodeI].bb_;

            if (intersects(subBb.min(), subBb.max(), nearestDistSqr, sample))
            {
                findNearest
                (
                    subNodeI,
                    sample,

                    nearestDistSqr,
                    nearestShapeI,
                    nearestPoint
                );
            }
        }
        else if (isContent(index))
        {
            if
            (
                intersects
                (
                    nod.bb_,
                    octant,
                    nearestDistSqr,
                    sample
                )
            )
            {
                shapes_.findNearest
                (
                    contents_[getContent(index)],
                    sample,

                    nearestDistSqr,
                    nearestShapeI,
                    nearestPoint
                );
            }
        }
    }
}


// Find nearest point to line.
template <class Type>
void indexedOctree<Type>::findNearest
(
    const label nodeI,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& nearestShapeI,
    point& linePoint,
    point& nearestPoint
) const
{
    const node& nod = nodes_[nodeI];
    const treeBoundBox& nodeBb = nod.bb_;

    // Determine order to walk through octants
    FixedList<direction, 8> octantOrder;
    nod.bb_.searchOrder(ln.centre(), octantOrder);

    // Go into all suboctants (one containing sample first) and update nearest.
    for (direction i = 0; i < 8; i++)
    {
        direction octant = octantOrder[i];

        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            const treeBoundBox& subBb = nodes_[getNode(index)].bb_;

            if (subBb.intersects(tightest))
            {
                findNearest
                (
                    getNode(index),
                    ln,

                    tightest,
                    nearestShapeI,
                    linePoint,
                    nearestPoint
                );
            }
        }
        else if (isContent(index))
        {
            const treeBoundBox subBb(nodeBb.subBbox(octant));

            if (subBb.intersects(tightest))
            {
                shapes_.findNearest
                (
                    contents_[getContent(index)],
                    ln,

                    tightest,
                    nearestShapeI,
                    linePoint,
                    nearestPoint
                );
            }
        }
    }
}


// Walk tree to neighbouring node. Gets current position as
// node and octant in this node and walks in the direction given by
// the faceID (one of treeBoundBox::LEFTBIT, RIGHTBIT etc.)
// Returns false if edge of tree hit.
template <class Type>
bool indexedOctree<Type>::walkToNeighbour
(
    const point& facePoint,
    const direction faceID,         // direction to walk in
    label& nodeI,
    direction& octant
) const
{
    // Find out how to test for octant. Say if we want to go left we need
    // to traverse up the tree until we hit a node where our octant is
    // on the right.

    direction octantMask = 0;
    direction valueMask = 0;

    if ((faceID & treeBoundBox::LEFTBIT) != 0)
    {
        // We want to go left so check if in right octant.
        octantMask |= treeBoundBox::RIGHTHALF;
        valueMask |= treeBoundBox::RIGHTHALF;
    }
    else if ((faceID & treeBoundBox::RIGHTBIT) != 0)
    {
        octantMask |= treeBoundBox::RIGHTHALF;  // valueMask already 0
    }

    if ((faceID & treeBoundBox::BOTTOMBIT) != 0)
    {
        octantMask |= treeBoundBox::TOPHALF;
        valueMask |= treeBoundBox::TOPHALF;
    }
    else if ((faceID & treeBoundBox::TOPBIT) != 0)
    {
        octantMask |= treeBoundBox::TOPHALF;
    }

    if ((faceID & treeBoundBox::BACKBIT) != 0)
    {
        octantMask |= treeBoundBox::FRONTHALF;
        valueMask |= treeBoundBox::FRONTHALF;
    }
    else if ((faceID & treeBoundBox::FRONTBIT) != 0)
    {
        octantMask |= treeBoundBox::FRONTHALF;
    }

    // Go up until we have chance to cross to the wanted direction
    while (valueMask != (octant & octantMask))
    {
        label parentNodeI = nodes_[nodeI].parent_;

        if (parentNodeI == -1)
        {
            // Reached edge of tree
            return false;
        }

        const node& parentNode = nodes_[parentNodeI];

        // Find octant nodeI is in.
        direction parentOctant = 255;

        for (direction i = 0; i < parentNode.subNodes_.size(); i++)
        {
            labelBits index = parentNode.subNodes_[i];

            if (isNode(index) && getNode(index) == nodeI)
            {
                parentOctant = i;
                break;
            }
        }

        if (parentOctant == 255)
        {
            FatalErrorIn("walkToNeighbour")
                << abort(FatalError);
        }
        nodeI = parentNodeI;
        octant = parentOctant;
    }

    // So now we hit the correct parent node. Switch to the correct octant
    octant ^= octantMask;

    // See if we need to travel down. Note that we already go into the
    // the first level ourselves (instead of having findNode decide) since
    // the facePoint is now exactly on the mid of the node so there could
    // be truncation problems.
    labelBits index = nodes_[nodeI].subNodes_[octant];

    if (isNode(index))
    {
        labelBits node = findNode(getNode(index), facePoint);

        nodeI = getNode(node);
        octant = getOctant(node);
    }

    return true;
}


// Return unique number for the face of a bounding box a point is on.
// (number is single bit but not really nessecary)
// Return 0 if point not on any face of bb.
template <class Type>
direction indexedOctree<Type>::getFace(const treeBoundBox& bb, const point& pt)
{
    direction faceID = 0;

    if (pt.x() <= bb.min().x())
    {
        faceID |= treeBoundBox::LEFTBIT;
    }
    if (pt.x() >= bb.max().x())
    {
        faceID |= treeBoundBox::RIGHTBIT;
    }

    if (pt.y() <= bb.min().y())
    {
        faceID |= treeBoundBox::BOTTOMBIT;
    }
    if (pt.y() >= bb.max().y())
    {
        faceID |= treeBoundBox::TOPBIT;
    }

    if (pt.z() <= bb.min().z())
    {
        faceID |= treeBoundBox::BACKBIT;
    }
    if (pt.z() >= bb.max().z())
    {
        faceID |= treeBoundBox::FRONTBIT;
    }
    return faceID;
}


// Traverse a node. If intersects a triangle return first intersection point.
// Else return the bounxing box face hit:
//  hitInfo.point = coordinate of intersection of ray with bounding box
//  faceID  = index of bounding box face
template <class Type>
void indexedOctree<Type>::traverseNode
(
    const bool findAny,
    const point& start,
    const point& end,
    const vector& dir,
    const label nodeI,
    const direction octant,

    pointIndexHit& hitInfo,
    direction& faceID
) const
{
    const node& nod = nodes_[nodeI];

    labelBits index = nod.subNodes_[octant];

    if (isContent(index))
    {
        const labelList& indices = contents_[getContent(index)];

        if (findAny)
        {
            // Find any intersection

            forAll(indices, elemI)
            {
                label shapeI = indices[elemI];

                point pt;
                bool hit = shapes_.intersects(shapeI, start, end, pt);

                if (hit)
                {
                    // Hit so pt is nearer than nearestPoint.
                    // Update hit info
                    hitInfo.setHit();
                    hitInfo.setIndex(shapeI);
                    hitInfo.setPoint(pt);
                    return;
                }
            }
        }
        else
        {
            // Find nearest intersection.

            point nearestPoint(end);

            forAll(indices, elemI)
            {
                label shapeI = indices[elemI];

                point pt;
                bool hit = shapes_.intersects(shapeI, start, nearestPoint, pt);

                if (hit)
                {
                    // Hit so pt is nearer than nearestPoint.
                    nearestPoint = pt;
                    // Update hit info
                    hitInfo.setHit();
                    hitInfo.setIndex(shapeI);
                    hitInfo.setPoint(pt);
                }
            }

            if (hitInfo.hit())
            {
                // Found intersected shape.
                return;
            }
        }
    }

    // Nothing intersected in this node
    // Traverse to end of node. Do by ray tracing back from end and finding
    // intersection with bounding box of node.

    point pt;

    if (isNode(index))
    {
        const treeBoundBox& subBb = nodes_[getNode(index)].bb_;

        if (!subBb.intersects(end, start, pt))
        {
            faceID = 0;

            WarningIn("indexedOctree<Type>::traverseNode")
                << "Did not hit side of box " << subBb
                << " with ray from " << start << " to " << end
                << endl;
        }
        else
        {
            faceID = getFace(subBb, pt);
        }
    }
    else
    {
        const treeBoundBox subBb(nod.bb_.subBbox(octant));

        if (!subBb.intersects(end, start, pt))
        {
            faceID = 0;

            WarningIn("indexedOctree<Type>::traverseNode")
                << "Did not hit side of box " << subBb
                << " with ray from " << start << " to " << end
                << endl;
        }
        else
        {
            faceID = getFace(subBb, pt);
        }
    }


    // Return miss. Set misspoint to face.
    hitInfo.setPoint(pt);
}


// Find first intersection
template <class Type>
pointIndexHit indexedOctree<Type>::findLine
(
    const bool findAny,
    const point& treeStart,
    const point& treeEnd,
    const label startNodeI,
    const direction startOctant
) const
{
    const vector dir(treeEnd - treeStart);

    // Current node as parent+octant
    label nodeI = startNodeI;
    direction octant = startOctant;

    // Current position. Initialize to miss
    pointIndexHit hitInfo(false, treeStart, -1);

    // Side of current bounding box current point is on
    direction faceID = 0;

    while (true)
    {
        // Ray-trace to end of current node. Updates point (either on triangle
        // in case of hit or on node bounding box in case of miss)

        point startPoint(hitInfo.rawPoint());

        traverseNode
        (
            findAny,
            startPoint,     // Note: pass in copy since hitInfo
                            // also used in return value.
            treeEnd,
            dir,
            nodeI,
            octant,

            hitInfo,
            faceID
        );

        // Did we hit a triangle?
        if (hitInfo.hit())
        {
            break;
        }

        if (faceID == 0)
        {
            // Was never inside the tree. Return miss.
            break;
        }

        //Pout<< "findLine : treeStart:" << treeStart
        //    << " treeEnd:" << treeEnd
        //    << " tracked from " << startPoint << " to edge of nodeBb:"
        //    << hitInfo.missPoint()
        //    << " node:" << nodeI << " octant:" << octant
        //    << " subBb:"
        //    << nodes_[nodeI].bb_.subBbox(octant)
        //    << endl;


        // Nothing hit so we are on face of bounding box (given as node+octant+
        // faceID). Traverse to neighbouring node.

        bool ok = walkToNeighbour
        (
            hitInfo.rawPoint(), // point on face
            faceID,             // direction to walk in
            nodeI,
            octant
        );

        if (!ok)
        {
            // Hit the edge of the tree. Return miss.
            break;
        }
    }
    return hitInfo;
}


// Find first intersection
template <class Type>
pointIndexHit indexedOctree<Type>::findLine
(
    const bool findAny,
    const point& start,
    const point& end
) const
{
    pointIndexHit hitInfo;

    if (nodes_.size() > 0)
    {
        const treeBoundBox& treeBb = nodes_[0].bb_;

        direction startBit = treeBb.posBits(start);
        direction endBit = treeBb.posBits(end);

        if ((startBit & endBit) != 0)
        {
            // Both start and end outside domain and in same block.
            return pointIndexHit(false, vector::zero, -1);
        }

        point trackStart(start);
        point trackEnd(end);

        if (startBit != 0)
        {
            // Track start to inside domain.
            if (!treeBb.intersects(start, end, trackStart))
            {
                return pointIndexHit(false, vector::zero, -1);
            }
        }

        if (endBit != 0)
        {
            // Track end to inside domain.
            if (!treeBb.intersects(end, trackStart, trackEnd))
            {
                return pointIndexHit(false, vector::zero, -1);
            }
        }

        // Find lowest level tree node that start is in.
        labelBits index = findNode(0, trackStart);

        label parentNodeI = getNode(index);
        direction octant = getOctant(index);

        hitInfo = findLine
        (
            findAny,
            trackStart,
            trackEnd,
            parentNodeI,
            octant
        );
    }

    return hitInfo;
}


template <class Type>
void indexedOctree<Type>::findBox
(
    const label nodeI,
    const treeBoundBox& searchBox,
    labelHashSet& elements
) const
{
    const node& nod = nodes_[nodeI];
    const treeBoundBox& nodeBb = nod.bb_;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            const treeBoundBox& subBb = nodes_[getNode(index)].bb_;

            if (subBb.intersects(searchBox))
            {
                findBox(getNode(index), searchBox, elements);
            }
        }
        else if (isContent(index))
        {
            const treeBoundBox subBb(nodeBb.subBbox(octant));

            if (subBb.intersects(searchBox))
            {
                const labelList& indices = contents_[getContent(index)];

                forAll(indices, i)
                {
                    label shapeI = indices[i];

                    if (shapes_.overlaps(shapeI, searchBox))
                    {
                        elements.insert(shapeI);
                    }
                }
            }
        }
    }
}


// Number of elements in node.
template <class Type>
label indexedOctree<Type>::countElements(const labelBits index) const
{
    label nElems = 0;

    if (isNode(index))
    {
        // tree node.
        label nodeI = getNode(index);

        const node& nod = nodes_[nodeI];

        for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
        {
            nElems += countElements(nod.subNodes_[octant]);
        }
    }
    else if (isContent(index))
    {
        nElems += contents_[getContent(index)].size();
    }
    else
    {
        // empty node
    }

    return nElems;
}


template <class Type>
void indexedOctree<Type>::writeOBJ
(
    const label nodeI,
    const direction octant
) const
{
    OFstream str
    (
        "node" + Foam::name(nodeI) + "_octant" + Foam::name(octant) + ".obj"
    );

    labelBits index = nodes_[nodeI].subNodes_[octant];

    treeBoundBox subBb;

    if (isNode(index))
    {
        subBb = nodes_[getNode(index)].bb_;
    }
    else if (isContent(index))
    {
        subBb = nodes_[nodeI].bb_.subBbox(octant);
    }

    Pout<< "dumpContentNode : writing node:" << nodeI << " octant:" << octant
        << " to " << str.name() << endl;

    label vertI = 0;

    // Dump bounding box
    pointField bbPoints(subBb.points());

    label pointVertI = vertI;
    forAll(bbPoints, i)
    {
        meshTools::writeOBJ(str, bbPoints[i]);
        vertI++;
    }

    forAll(treeBoundBox::edges, i)
    {
        const edge& e = treeBoundBox::edges[i];

        str<< "l " << e[0]+pointVertI+1 << ' ' << e[1]+pointVertI+1 << nl;
    }


    //// Dump triangles
    //if (isContent(index))
    //{
    //    const labelList& indices = contents_[getContent(index)];
    //    const triSurface& surf = shapes_.surface();
    //    const pointField& points = surf.points();
    //
    //    forAll(indices, i)
    //    {
    //        label shapeI = indices[i];
    //
    //        const labelledTri& f = surf[shapeI];
    //
    //        meshTools::writeOBJ(str, points[f[0]]);
    //        vertI++;
    //        meshTools::writeOBJ(str, points[f[1]]);
    //        vertI++;
    //        meshTools::writeOBJ(str, points[f[2]]);
    //        vertI++;
    //
    //        str<< "l " << vertI-2 << ' ' << vertI-1 << ' ' << vertI << ' '
    //            << vertI-2 << nl;
    //    }
    //}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Type>
indexedOctree<Type>::indexedOctree(const Type& shapes)
:
    shapes_(shapes),
    nodes_(0),
    contents_(0),
    nodeTypes_(0)
{}


template <class Type>
indexedOctree<Type>::indexedOctree
(
    const Type& shapes,
    const List<node>& nodes,
    const labelListList& contents
)
:
    shapes_(shapes),
    nodes_(nodes),
    contents_(contents),
    nodeTypes_(0)
{}


template <class Type>
indexedOctree<Type>::indexedOctree
(
    const Type& shapes,
    const treeBoundBox& bb,
    const label maxLevels,          // maximum number of levels
    const scalar maxLeafRatio,
    const scalar maxDuplicity
)
:
    shapes_(shapes),
    nodes_(0),
    contents_(0),
    nodeTypes_(0)
{
    if (shapes.size() == 0)
    {
        return;
    }

    // Start off with one node with all shapes in it.
    DynamicList<node> nodes(label(shapes.size() / maxLeafRatio));
    DynamicList<labelList> contents(label(shapes.size() / maxLeafRatio));
    contents.append(identity(shapes.size()));

    // Create topnode.
    node topNode(divide(bb, contents, 0));
    nodes.append(topNode);


    // Now have all contents at level 1. Create levels by splitting levels
    // above.

    label nLevels = 1;

    for (; nLevels < maxLevels; nLevels++)
    {
        // Count number of references into shapes (i.e. contents)
        label nEntries = 0;
        forAll(contents, i)
        {
            nEntries += contents[i].size();
        }

        if (debug)
        {
            Pout<< "indexedOctree<Type>::indexedOctree:" << nl
                << "    nLevels:" << nLevels << nl
                << "    nEntries per treeLeaf:" << nEntries/contents.size()
                << nl
                << "    nEntries per shape (duplicity):"
                << nEntries/shapes.size()
                << nl
                << endl;
        }

        if
        (
            //nEntries < maxLeafRatio*contents.size()
         // ||
            nEntries > maxDuplicity*shapes.size()
        )
        {
            break;
        }


        // Split nodes with more than maxLeafRatio elements
        label nOldNodes = nodes.size();
        splitNodes
        (
            label(maxLeafRatio),
            nodes,
            contents
        );

        if (nOldNodes == nodes.size())
        {
            break;
        }
    }

    // Shrink
    nodes.shrink();
    contents.shrink();


    // Compact such that deeper level contents are always after the
    // ones for a shallower level. This way we can slice a coarser level
    // off the tree.
    contents_.setSize(contents.size());
    label compactI = 0;

    label level = 0;

    while (true)
    {
        compactContents
        (
            nodes,
            contents,
            level,
            0,
            0,
            contents_,
            compactI
        );

        if (compactI == contents_.size())
        {
            // Transferred all contents to contents_ (in order breadth first)
            break;
        }

        level++;
    }
    nodes_.transfer(nodes);
    nodes.clear();


    if (debug)
    {
        label nEntries = 0;
        forAll(contents_, i)
        {
            nEntries += contents_[i].size();
        }

        Pout<< "indexedOctree<Type>::indexedOctree : finished construction:"
            << nl
            << "    shapes:" << shapes.size() << nl
            << "    nLevels:" << nLevels << nl
            << "    treeNodes:" << nodes_.size() << nl
            << "    nEntries:" << nEntries << nl
            << "        per treeLeaf:" << nEntries/contents.size() << nl
            << "        per shape (duplicity):" << nEntries/shapes.size() << nl
            << endl;
    }
}


template <class Type>
indexedOctree<Type>::indexedOctree
(
    const Type& shapes,
    Istream& is
)
:
    shapes_(shapes),
    nodes_(is),
    contents_(is),
    nodeTypes_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Type>
scalar& Foam::indexedOctree<Type>::perturbTol()
{
    return perturbTol_;
}

template <class Type>
pointIndexHit indexedOctree<Type>::findNearest
(
    const point& sample,
    const scalar startDistSqr
) const
{
    scalar nearestDistSqr = startDistSqr;
    label nearestShapeI = -1;
    point nearestPoint;

    if (nodes_.size() == 0)
    {
        nearestPoint = vector::zero;
    }
    else
    {
        findNearest
        (
            0,
            sample,

            nearestDistSqr,
            nearestShapeI,
            nearestPoint
        );
    }

    return pointIndexHit(nearestShapeI != -1, nearestPoint, nearestShapeI);
}


template <class Type>
pointIndexHit indexedOctree<Type>::findNearest
(
    const linePointRef& ln,
    treeBoundBox& tightest,
    point& linePoint
) const
{
    label nearestShapeI = -1;
    point nearestPoint;

    if (nodes_.size() == 0)
    {
        nearestPoint = vector::zero;
    }
    else
    {
        findNearest
        (
            0,
            ln,

            tightest,
            nearestShapeI,
            linePoint,
            nearestPoint
        );
    }

    return pointIndexHit(nearestShapeI != -1, nearestPoint, nearestShapeI);
}


// Find nearest intersection
template <class Type>
pointIndexHit indexedOctree<Type>::findLine
(
    const point& start,
    const point& end
) const
{
    return findLine(false, start, end);
}


// Find nearest intersection
template <class Type>
pointIndexHit indexedOctree<Type>::findLineAny
(
    const point& start,
    const point& end
) const
{
    return findLine(true, start, end);
}


template <class Type>
labelList indexedOctree<Type>::findBox(const boundBox& searchBox) const
{
    // Storage for labels of shapes inside bb. Size estimate.
    labelHashSet elements(shapes_.size() / 100);

    if (nodes_.size() > 0)
    {
        findBox(0, searchBox, elements);
    }

    return elements.toc();
}


// Find node (as parent+octant) containing point
template <class Type>
labelBits indexedOctree<Type>::findNode
(
    const label nodeI,
    const point& sample
) const
{
    if (nodes_.size() == 0)
    {
        // Empty tree. Return what?
        return nodePlusOctant(nodeI, 0);
    }

    const node& nod = nodes_[nodeI];

    direction octant = nod.bb_.subOctant(sample);

    labelBits index = nod.subNodes_[octant];

    if (isNode(index))
    {
        // Recurse
        return findNode(getNode(index), sample);
    }
    else if (isContent(index))
    {
        // Content. Return treenode+octant
        return nodePlusOctant(nodeI, octant);
    }
    else
    {
        // Empty. Return treenode+octant
        return nodePlusOctant(nodeI, octant);
    }
}


// Determine type (inside/outside/mixed) per node.
template <class Type>
typename indexedOctree<Type>::volumeType indexedOctree<Type>::getVolumeType
(
    const point& sample
) const
{
    if (nodes_.size() == 0)
    {
        return UNKNOWN;
    }

    if (nodeTypes_.size() != 8*nodes_.size())
    {
        // Calculate type for every octant of node.

        nodeTypes_.setSize(8*nodes_.size());
        nodeTypes_ = UNKNOWN;

        calcVolumeType(0);

        if (debug)
        {
            label nUNKNOWN = 0;
            label nMIXED = 0;
            label nINSIDE = 0;
            label nOUTSIDE = 0;

            forAll(nodeTypes_, i)
            {
                volumeType type = volumeType(nodeTypes_.get(i));

                if (type == UNKNOWN)
                {
                    nUNKNOWN++;
                }
                else if (type == MIXED)
                {
                    nMIXED++;
                }
                else if (type == INSIDE)
                {
                    nINSIDE++;
                }
                else if (type == OUTSIDE)
                {
                    nOUTSIDE++;
                }
                else
                {
                    FatalErrorIn("getVolumeType") << abort(FatalError);
                }
            }

            Pout<< "indexedOctree<Type>::getVolumeType : "
                << " bb:" << bb()
                << " nodes_:" << nodes_.size()
                << " nodeTypes_:" << nodeTypes_.size()
                << " nUNKNOWN:" << nUNKNOWN
                << " nMIXED:" << nMIXED
                << " nINSIDE:" << nINSIDE
                << " nOUTSIDE:" << nOUTSIDE
                << endl;
        }
    }

    return getVolumeType(0, sample);
}


// Print contents of nodeI
template <class Type>
void indexedOctree<Type>::print
(
    prefixOSstream& os,
    const bool printContents,
    const label nodeI
) const
{
    const node& nod = nodes_[nodeI];
    const treeBoundBox& bb = nod.bb_;

    os  << "nodeI:" << nodeI << " bb:" << bb << nl
        << "parent:" << nod.parent_ << nl
        << "n:" << countElements(nodePlusOctant(nodeI, 0)) << nl;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        const treeBoundBox subBb(bb.subBbox(octant));

        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            // tree node.
            label subNodeI = getNode(index);

            os  << "octant:" << octant
                << " node: n:" << countElements(index)
                << " bb:" << subBb << endl;

            string oldPrefix = os.prefix();
            os.prefix() = "  " + oldPrefix;

            print(os, printContents, subNodeI);

            os.prefix() = oldPrefix;
        }
        else if (isContent(index))
        {
            const labelList& indices = contents_[getContent(index)];

            os  << "octant:" << octant
                << " content: n:" << indices.size()
                << " bb:" << subBb;

            if (printContents)
            {
                os << " contents:";
                forAll(indices, j)
                {
                    os  << ' ' << indices[j];
                }
            }
            os  << endl;
        }
        else
        {
            os  << "octant:" << octant << " empty:" << subBb << endl;
        }
    }
}


// Print contents of nodeI
template <class Type>
bool indexedOctree<Type>::write(Ostream& os) const
{
    os << *this;

    return os.good();
}


template <class Type>
Ostream& operator<<(Ostream& os, const indexedOctree<Type>& t)
{
    return
        os  << t.bb() << token::SPACE << t.nodes()
            << token::SPACE << t.contents();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
