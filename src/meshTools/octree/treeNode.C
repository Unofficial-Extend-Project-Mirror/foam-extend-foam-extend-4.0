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

Description

\*---------------------------------------------------------------------------*/

#include "treeNode.H"
#include "octree.H"
#include "treeLeaf.H"
#include "treeBoundBox.H"
#include "long.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class Type>
const label treeNode<Type>::leafOffset = 100;

template <class Type>
const labelList treeNode<Type>::dummy(1);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
void treeNode<Type>::setAsNode(const label octant)
{
    subNodeTypes_ |= (0x1 << octant);
}


template <class Type>
void treeNode<Type>::setAsLeaf(const label octant)
{
    subNodeTypes_ &= ~(0x1 << octant);
}


// Set pointer to sub node
template <class Type>
void treeNode<Type>::setNodePtr(const label octant, treeElem<Type>* treeNodePtr)
{
    setAsNode(octant);
    subNodes_[octant] = treeNodePtr;
}


// Set pointer to sub leaf
template <class Type>
void treeNode<Type>::setLeafPtr(const label octant, treeElem<Type>* treeLeafPtr)
{
    setAsLeaf(octant);
    subNodes_[octant] = treeLeafPtr;
}


template <class Type>
void treeNode<Type>::setVolType(const label octant, const label type)
{
    if ((type < 0) || (type > 3))
    {
        FatalErrorIn("treeNode<Type>::setVolType(const label, const label)")
            << "Type " << type << " not within range 0..3" << endl;
    }

    // Clear out two bits at position 2*octant
    volType_ &= ~(0x3 << 2*octant);

    // Add the two bits of type
    volType_ |= (type << 2*octant);
}


template <class Type>
void treeNode<Type>::space(Ostream& os, const label n)
{
    for(label i=0; i<n; i++)
    {
        os<< ' ';
    }
}


// look in single octant starting from <start>
template <class Type>
const treeLeaf<Type>* treeNode<Type>::findLeafLineOctant
(
    const int level,
    const Type& shapes,
    const label octant,
    const vector& direction,
    point& start,
    const point& end
) const
{
    static const char* functionName =
        "treeNode<Type>::findLeafLineOctant"
        "(const int, const Type&, const label, const vector&,"
        " point&, const point&)";

    if (debug & 2)
    {
        space(Pout, 2*level);
        Pout<< "findLeafLineOctant : bb:" << this->bb()
            << "  start:" << start
            << "  end:" << end
            << "  mid:" << mid()
            << " Searching octant:" << octant
            << endl;
    }

    if (subNodes()[octant])
    {
        if (isNode(octant))
        {
            // Node: recurse into subnodes
            const treeNode<Type>* subNodePtr = getNodePtr(octant);

            if (subNodePtr->bb().contains(direction, start))
            {
                // Search on lower level
                const treeLeaf<Type>* subLeafPtr =
                    subNodePtr->findLeafLine
                    (
                        level + 1,
                        shapes,
                        start,
                        end
                    );

                if (debug & 2)
                {
                    space(Pout, 2*level);
                    Pout<< "findLeafLineOctant : bb:" << this->bb()
                        << " returning from sub treeNode"
                        << " with start:" << start << "  subLeaf:"
                        << long(subLeafPtr) << endl;
                }

                return subLeafPtr;
            }
            else
            {
                FatalErrorIn(functionName)
                    << "Sub node " << subNodePtr->bb()
                    << " at octant " << octant
                    << " does not contain start " << start
                    << abort(FatalError);
            }
        }
        else
        {
            // Leaf
            const treeLeaf<Type>* subLeafPtr = getLeafPtr(octant);

            if (subLeafPtr->bb().contains(direction, start))
            {
                // Step to end of subleaf bb
                point tmp;
                if 
                (
                   !subLeafPtr->bb().intersects
                    (
                        end,
                        start,
                        tmp
                    )
                )
                {
                    FatalErrorIn(functionName)
                        << "Sub leaf contains start " << start
                        << " but line does not intersect its bb "
                        << subLeafPtr->bb()
                        << abort(FatalError);
                }
                start = tmp;

                if (debug & 2)
                {
                    space(Pout, 2*level);
                    Pout<< "findLeafLineOctant : returning from intersecting"
                        << " treeLeaf " << subLeafPtr->bb()
                        << " with start:" << start << "  subLeaf:"
                        << long(subLeafPtr) << endl;
                }

                return subLeafPtr;
            }
            else
            {
                FatalErrorIn(functionName)
                    << "Sub leaf " << subLeafPtr->bb()
                    << " at octant " << octant
                    << " does not contain start " << start
                    << abort(FatalError);
            }
        }
    }
    else
    {
        // Empty subNode. Transfer across.
        const treeBoundBox emptyBb = this->bb().subBbox(mid(), octant);

        if (emptyBb.contains(direction, start))
        {
            if (debug & 2)
            {
                space(Pout, 2*level);
                Pout<< "findLeafLineOctant : Empty node. Octant:" << octant
                    << "  start:" << start
                    << "  bb:" << this->bb()
                    << "  emptyBb:" << emptyBb << endl;
            }

            // Update start by clipping to emptyBb
            point tmp;
            if 
            (
               !emptyBb.intersects
                (
                    end,
                    start,
                    tmp
                )
            )
            {
                FatalErrorIn(functionName)
                    << "Empty node contains start " << start
                    << " but line does not intersect its (calculated)"
                    << " bb " << emptyBb
                    << endl << "This might be due to truncation error"
                    << abort(FatalError);
            }
            start = tmp;

            if (debug & 2)
            {
                space(Pout, 2*level);
                Pout<< "findLeafLineOctant : returning from intersecting with"
                    << " empty " << emptyBb
                    << " with start:" << start << "  subLeaf:" << 0 << endl;
            }

            return NULL;
        }
        else
        {
            FatalErrorIn(functionName)
                << "Empty node " << emptyBb
                << " at octant " << octant
                << " does not contain start " << start
                << abort(FatalError);
        }
    }

    FatalErrorIn(functionName)
        << "Octant " << octant << " of cube " << this->bb()
        << " does not contain start " << start
        << abort(FatalError);

    return NULL;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template <class Type>
treeNode<Type>::treeNode(const treeBoundBox& bb)
:
    treeElem<Type>(bb),
    treeNodeName(),
    mid_(bb.mid()),
    subNodeTypes_(0),
    volType_(0)
{
    for(label octant=0; octant<8; octant++)
    {
        subNodes_[octant] = NULL;
        setVolType(octant, octree<Type>::UNKNOWN);
    }
}


// Construct from Istream
template <class Type>
treeNode<Type>::treeNode
(
    Istream& is
)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class Type>
treeNode<Type>::~treeNode()
{
    for(int octant=0; octant<8; octant++)
    {
        if (subNodes()[octant])
        {
            if (isNode(octant))
            {
                delete getNodePtr(octant);
            }
            else
            {
                delete getLeafPtr(octant);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Distributes cells to subLeaves
template <class Type>
void treeNode<Type>::distribute
(
    const label level,
    octree<Type>& top,
    const Type& shapes,
    const labelList& indices
)
{
    if (debug & 1)
    {
        space(Pout, level);
        Pout<< "treeNode::distributing " << indices.size() << endl;
    }

    // Create subLeaves if nessecary
    for(label octant=0; octant<8; octant++)
    {
        if (subNodes()[octant])
        {
            printNode(Pout, level);
            FatalErrorIn
            (
                "treeNode<Type>::distribute(const label, octree<Type>&, "
                "const Type&, const labelList&)"
            )   << "subNode already available at octant:" << octant
                << abort(FatalError);                
        }
        else
        {
            treeLeaf<Type>* subLeafPtr =
                new treeLeaf<Type>
                (
                    this->bb().subBbox(mid(), octant),
                    indices.size()
                );

            top.setLeaves(top.nLeaves() + 1);
            setLeafPtr(octant, subLeafPtr);
        }
    }


    // add cells to correct sub leaf
    forAll(indices, i)
    {
        const label shapei = indices[i];

        for(label octant=0; octant<8; octant++)
        {
            treeLeaf<Type>* leafPtr = getLeafPtr(octant);

            if (shapes.overlaps(shapei, leafPtr->bb()))
            {
                if (debug == 1)
                {
                    space(Pout, level);
                    Pout<< "inserting " << shapei;
                    shapes.write(Pout, shapei);
                    Pout<< " into " << leafPtr->bb() << endl;
                }
                leafPtr->insert(shapei);
                top.setEntries(top.nEntries() + 1);
            }
        }
    }

    // Trim size of subLeaves
    for(label octant=0; octant<8; octant++)
    {
        treeLeaf<Type>* subLeafPtr = getLeafPtr(octant);

        if (subLeafPtr->size() == 0)
        {
            // Contains no data. Delete.
            setLeafPtr(octant, NULL);
            delete subLeafPtr;
            top.setLeaves(top.nLeaves() - 1);
        }
        else
        {
            // Trim to actual size.
            subLeafPtr->trim();
        }
    }

    if (debug & 1)
    {
        space(Pout, level);
        Pout<< "end of treeNode::distribute" << endl;
    }
}


// Descends to refineLevel and checks the subLeaves for redistribution
template <class Type>
void treeNode<Type>::redistribute
(
    const label level,
    octree<Type>& top,
    const Type& shapes,
    const label refineLevel
)
{
    if (debug & 1)
    {
        space(Pout, level);
        Pout<< "treeNode::redistribute with level:" << level
            << "  refineLevel:" << refineLevel << endl;
    }

    // Descend to correct level
    if (level < refineLevel)
    {
        for(label octant=0; octant<8; octant++)
        {
            if (subNodes()[octant])
            {
                if (isNode(octant))
                {
                    getNodePtr(octant)->redistribute
                    (
                        level+1,
                        top,
                        shapes,
                        refineLevel
                    );
                }
            }
        }
    }
    else
    {
        // Reached correct (should also be deepest) level of treeNode
        if (debug & 1)
        {
            space(Pout, level);
            Pout<< "treeNode::redistribute : now at correct level" << endl;
        }

        // handle redistribution of sub leaves
        for(label octant=0; octant<8; octant++)
        {
            if (subNodes()[octant])
            {
                if (isNode(octant))
                {
                    FatalErrorIn
                    (
                        "treeNode<Type>::redistribute(const int, octree& top,"
                        "const int, const treeBoundBox&)"
                    )   << "found treeNode instead of treeLeaf" << endl
                        << abort(FatalError);
                }
                else
                {
                    treeLeaf<Type>* leafPtr = getLeafPtr(octant);

                    treeLeaf<Type>* newSubPtr = leafPtr->redistribute
                    (
                        level,
                        top,
                        shapes
                    );

                    if (newSubPtr && (newSubPtr != leafPtr))
                    {
                        // redistribute has created nodePtr
                        // so delete no longer used subPtr and update info
                        if (debug & 1)
                        {
                            Pout<< "deleting "
                                << top.nEntries() - leafPtr->size()
                                << " entries" << endl;
                        }
                        top.setEntries(top.nEntries() - leafPtr->size());

                        delete leafPtr;

                        top.setLeaves(top.nLeaves() - 1);

                        setNodePtr(octant, newSubPtr);
                    }
                }
            }
        }
        if (debug & 1)
        {
            space(Pout, level);
            Pout<< "end of treeNode::redistribute for correct level" << endl;
        }
    }

    if (debug & 1)
    {
        space(Pout, level);
        Pout<< "return from treeNode::redistribute with bb:" << this->bb()
            << endl;
    }
}


// Set type of node.
template <class Type>
label treeNode<Type>::setSubNodeType
(
    const label level,
    octree<Type>& top,
    const Type& shapes
)
{
    if (debug & 4)
    {
        space(Pout, level);
        Pout<< "treeNode::setSubNodeType with level:" << level
            << "   bb:" << this->bb() << endl;
    }

    label myType = -1;

    for(label octant=0; octant<8; octant++)
    {
        label subType = -1;

        if (subNodes()[octant])
        {
            if (isNode(octant))
            {
                subType = getNodePtr(octant)->setSubNodeType
                (
                    level+1,
                    top,
                    shapes
                );
            }
            else
            {
                subType = getLeafPtr(octant)->setSubNodeType
                (
                    level+1,
                    top,
                    shapes
                );
            }
        }
        else
        {
            // No data in this one. Set type for octant acc. to its bounding
            // box.
            const treeBoundBox subBb = this->bb().subBbox(mid(), octant);

            subType = shapes.getSampleType(top, subBb.mid());
        }

        if (debug & 4)
        {
            space(Pout, level);
            Pout<< "treeNode::setSubNodeType : setting octant with bb:"
                << this->bb().subBbox(mid(), octant)
                << "  to type:" << octree<Type>::volType(subType) << endl;
        }
        setVolType(octant, subType);

        // Combine sub node types into type for treeNode. Result is 'mixed' if
        // types differ among subnodes.
        if (myType == -1)
        {
            myType = subType;
        }
        else if (subType != myType)
        {
            myType = octree<Type>::MIXED;
        }
    }

    if (debug & 4)
    {
        space(Pout, level);
        Pout<< "return from treeNode::setSubNodeType with type:"
            << octree<Type>::volType(myType)
            << "  bb:" << this->bb() << endl;
    }

    return myType;
}


// Get type of node.
template <class Type>
label treeNode<Type>::getSampleType
(
    const label level,
    const octree<Type>& top,
    const Type& shapes,
    const point& sample
) const
{
    if (debug & 4)
    {
        space(Pout, level);
        Pout<< "treeNode::getSampleType with level:" << level
            << " bb:" << this->bb() << "  sample:" << sample << endl;
    }

    // Determine octant of bb. If on edge just use whichever octant.
    bool onEdge = false;

    label octant = this->bb().subOctant(mid(), sample, onEdge);

    label type = getVolType(octant);

    if (type == octree<Type>::MIXED)
    {
        // At this level multiple sub types. Recurse to resolve.
        if (subNodes()[octant])
        {
            if (isNode(octant))
            {
                // Node: recurse into subnodes
                type = getNodePtr(octant)->getSampleType
                (
                    level + 1,
                    top,
                    shapes,
                    sample
                );
            }
            else
            {
                // Leaf
                type = getLeafPtr(octant)->getSampleType
                (
                    level + 1,
                    top,
                    shapes,
                    sample
                );
            }
        }
        else
        {
            // Problem: empty subnode should have a type
            FatalErrorIn
            (
                "treeNode<Type>::getSampleType"
                "(const label, octree<Type>&, const Type&, const point&)"
            )   << "Empty node bb:" << this->bb().subBbox(mid(), octant)
                << " has non-mixed type:"
                << octree<Type>::volType(type)
                << abort(FatalError);
        }
    }

    if (type == octree<Type>::MIXED)
    {
        FatalErrorIn
        (
            "treeNode<Type>::getSampleType"
            "(const label, octree<Type>&, const Type&, const point&)"
        )   << "Type is MIXED when searching for " << sample
            << " at level " << this->bb() << endl
            << "This probably is because the octree has not been constructed"
            << " with search facility." << exit(FatalError);
    }

    if (debug & 4)
    {
        space(Pout, level);
        Pout<< "return from treeNode::getSampleType with type:"
            << octree<Type>::volType(type)
            << "  bb:" << this->bb()
            << "  sample:" << sample << endl;
    }
    return type;
}


template <class Type>
label treeNode<Type>::find
(
    const Type& shapes,
    const point& sample
) const
{
    // Find octant of sample. Don't care if on edge (since any item on edge
    // will have been inserted in both subcubes)
    bool onEdge = false;

    label octant = this->bb().subOctant(mid(), sample, onEdge);

    if (subNodes()[octant])
    {
        if (isNode(octant))
        {
            // Node: recurse into subnodes
            return getNodePtr(octant)->find
            (
                shapes,
                sample
            );
        }
        else
        {
            // Leaf: let leaf::find handle this
            return getLeafPtr(octant)->find(shapes, sample);
        }
    }
    return -1;
}


template <class Type>
bool treeNode<Type>::findTightest
(
    const Type& shapes,
    const point& sample,
    treeBoundBox& tightest
) const
{
    bool changed = false;

    // Get best guess for starting octant
    bool onEdge = false;

    label sampleOctant = this->bb().subOctant(mid(), sample, onEdge);

    // Go into all suboctants (one containing sample first) and update tightest.
    // Order of visiting is if e.g. sampleOctant = 5:
    //  5 1 2 3 4 0 6 7
    for(label octanti=0; octanti<8; octanti++)
    {
        label octant;
        if (octanti == 0)
        {
            // Use sampleOctant first
            octant = sampleOctant;
        }
        else if (octanti == sampleOctant)
        {
            octant = 0;
        }
        else
        {
            octant = octanti;
        }

        if (subNodes()[octant])
        {
            if (isNode(octant))
            {
                // Node: recurse into subnodes
                const treeNode<Type>* subNodePtr = getNodePtr(octant);

                if (subNodePtr->bb().intersects(tightest))
                {
                    // there might be a better fit inside this subNode
                    changed |=
                        subNodePtr->findTightest
                        (
                            shapes,
                            sample,
                            tightest
                        );
                }
            }
            else
            {
                // Leaf: let leaf::find handle this
                const treeLeaf<Type>* subLeafPtr = getLeafPtr(octant);

                if (subLeafPtr->bb().intersects(tightest))
                {
                    // there might be a better fit inside this subLeaf
                    changed |=
                        subLeafPtr->findTightest
                        (
                            shapes,
                            sample,
                            tightest
                        );
                }
            }
        }
    }

    return changed;
}


template <class Type>
bool treeNode<Type>::findNearest
(
    const Type& shapes,
    const point& sample,
    treeBoundBox& tightest,
    label& tightesti,
    scalar& tightestDist
) const
{
    bool changed = false;


    if (debug & 8)
    {
        Pout<< "In findNearest with sample:" << sample << " cube:"
            << this->bb() << " tightest:" << tightest << endl;
    }

    bool onEdge = false;

    label sampleOctant = this->bb().subOctant(mid(), sample, onEdge);
    
    // Go into all suboctants (one containing sample first) and update tightest.
    // Order of visiting is if e.g. sampleOctant = 5:
    //  5 1 2 3 4 0 6 7
    for(label octanti=0; octanti<8; octanti++)
    {
        label octant;
        if (octanti == 0)
        {
            // Use sampleOctant first
            octant = sampleOctant;
        }
        else if (octanti == sampleOctant)
        {
            octant = 0;
        }
        else
        {
            octant = octanti;
        }

        if (subNodes()[octant])
        {
            if (isNode(octant))
            {
                // Node
                const treeNode<Type>* subNodePtr = getNodePtr(octant);

                if (subNodePtr->bb().intersects(tightest))
                {
                    // there might be a better fit inside this subNode
                    changed |=
                        subNodePtr->findNearest
                        (
                            shapes,
                            sample,
                            tightest,
                            tightesti,
                            tightestDist
                        );
                }
            }
            else
            {
                // Leaf: let leaf::find handle this
                const treeLeaf<Type>* subLeafPtr = getLeafPtr(octant);

                if (subLeafPtr->bb().intersects(tightest))
                {
                    // there might be a better fit inside this subNode
                    changed |=
                        subLeafPtr->findNearest
                        (
                            shapes,
                            sample,
                            tightest,
                            tightesti,
                            tightestDist
                        );
                }
            }
        }
    }

    if (debug & 8)
    {
        Pout<< "Exiting findNearest for sample:" << sample << " cube:"
            << this->bb() << " tightesti:" << tightesti << endl;
    }

    return changed;
}


template <class Type>
bool treeNode<Type>::findNearest
(
    const Type& shapes,
    const linePointRef& ln,
    treeBoundBox& tightest,
    label& tightesti,
    point& linePoint,   // nearest point on line
    point& shapePoint   // nearest point on shape
) const
{
    bool changed = false;

    bool onEdge = false;

    // Estimate for where best to start searching
    label sampleOctant = this->bb().subOctant(mid(), ln.centre(), onEdge);
    
    // Go into all suboctants (one containing sample first) and update tightest.
    // Order of visiting is if e.g. sampleOctant = 5:
    //  5 1 2 3 4 0 6 7
    for(label octanti=0; octanti<8; octanti++)
    {
        label octant;
        if (octanti == 0)
        {
            // Use sampleOctant first
            octant = sampleOctant;
        }
        else if (octanti == sampleOctant)
        {
            octant = 0;
        }
        else
        {
            octant = octanti;
        }

        if (subNodes()[octant])
        {
            if (isNode(octant))
            {
                // Node
                const treeNode<Type>* subNodePtr = getNodePtr(octant);

                if (subNodePtr->bb().intersects(tightest))
                {
                    // there might be a better fit inside this subNode
                    changed |=
                        subNodePtr->findNearest
                        (
                            shapes,
                            ln,
                            tightest,
                            tightesti,
                            linePoint,
                            shapePoint
                        );
                }
            }
            else
            {
                // Leaf: let leaf::find handle this
                const treeLeaf<Type>* subLeafPtr = getLeafPtr(octant);

                if (subLeafPtr->bb().intersects(tightest))
                {
                    // there might be a better fit inside this subNode
                    changed |=
                        subLeafPtr->findNearest
                        (
                            shapes,
                            ln,
                            tightest,
                            tightesti,
                            linePoint,
                            shapePoint
                        );
                }
            }
        }
    }

    return changed;
}


template <class Type>
bool treeNode<Type>::findBox
(
    const Type& shapes,
    const boundBox& box,
    labelHashSet& elements
) const
{
    bool changed = false;

    bool onEdge = false;

    // Estimate for where best to start searching
    point boxMid(0.5*(box.min() + box.max()));
    label sampleOctant = this->bb().subOctant(mid(), boxMid, onEdge);
    
    // Go into all suboctants (one containing sample first) and update tightest.
    // Order of visiting is if e.g. sampleOctant = 5:
    //  5 1 2 3 4 0 6 7
    for(label octanti=0; octanti<8; octanti++)
    {
        label octant;
        if (octanti == 0)
        {
            // Use sampleOctant first
            octant = sampleOctant;
        }
        else if (octanti == sampleOctant)
        {
            octant = 0;
        }
        else
        {
            octant = octanti;
        }

        if (subNodes()[octant])
        {
            if (isNode(octant))
            {
                // Node
                const treeNode<Type>* subNodePtr = getNodePtr(octant);

                if (subNodePtr->bb().intersects(box))
                {
                    // Visit sub node.
                    changed |= subNodePtr->findBox(shapes, box, elements);
                }
            }
            else
            {
                // Leaf: let leaf::find handle this
                const treeLeaf<Type>* subLeafPtr = getLeafPtr(octant);

                if (subLeafPtr->bb().intersects(box))
                {
                    // Visit sub leaf.
                    changed |= subLeafPtr->findBox(shapes, box, elements);
                }
            }
        }
    }

    return changed;
}


// look from <start> in current cube (given by this->bb()).
template <class Type>
const treeLeaf<Type>* treeNode<Type>::findLeafLine
(
    const int level,
    const Type& shapes,
    point& start,
    const point& end
) const
{
    if (debug & 2)
    {
        space(Pout, 2*level);
        Pout<< "findLeafLine : bb:" << this->bb() << "  mid:" << mid()
            << "  start:" << start << endl;
    }

    scalar typDim = this->bb().avgDim();

    const vector direction = end - start;

    // Loop on current level until start has been updated to be outside
    // of this->bb(). Note that max only four subcubes can be crossed so this is
    // check on whether there are any truncation error problems.

    label iter = 0;

    while(true)
    {
        if (!this->bb().contains(direction, start))
        {
            if (debug & 2)
            {
                space(Pout, 2*level);
                Pout<< "findLeafLine : Start not inside bb " << this->bb()
                    << ". Returning with start:" << start << "  subLeaf:"
                    << 0 << endl;
            }
            return NULL;
        }

        // Check if start and <end> equal
        if ((mag(start - end)/typDim) < SMALL)
        {
            if (debug & 2)
            {
                space(Pout, 2*level);
                Pout<< "findLeafLine : start equals end"
                    << ". Returning with start:" << start << "  subLeaf:"
                    << 0 << endl;
            }
            return NULL;
        }

        if (iter >= 4)
        {
            // Too many iterations. Is hanging. Handle outside of loop.
            break;
        }

        bool onEdge = false;
        label octant = this->bb().subOctant(mid(), direction, start, onEdge);

        // Try finding non-empty treeleaf in octant
        const treeLeaf<Type>* leafPtr =
            findLeafLineOctant
            (
                level,
                shapes,
                octant,
                direction,
                start,
                end
            );

        if (leafPtr)
        {
            // Found treeLeaf -> return
            if (debug & 2)
            {
                space(Pout, 2*level);
                Pout<< "findLeafLine : Found treeLeaf"
                        << ". Returning with start:" << start << "  subLeaf:"
                        << long(leafPtr) << endl;
            }

            return leafPtr;
        }

        iter++;
    }

    // Check if is hanging. Max 4 octants can be crossed by a straight line
    FatalErrorIn
    (
        "treeNode<Type>::findLeafLine"
        "(const label, octree<Type>&, point&,"
        " const point&)"
    )   << "Did not leave bb " << this->bb()
        << " after " << iter
        << " iterations of updating starting point."
        << "start:" << start << "  end:" << end
        << abort(FatalError);

    return NULL;
}


template <class Type>
void treeNode<Type>::findLeaves
(
    List<treeLeaf<Type>*>& leafArray,
    label& leafIndex
) const
{
    // Go into all sub boxes
    for(label octant=0; octant<8; octant++)
    {
        if (subNodes()[octant])
        {
            if (isNode(octant))
            {
                // Node: recurse into subnodes
                const treeNode<Type>* subNodePtr = getNodePtr(octant);
                subNodePtr->findLeaves(leafArray, leafIndex);
            }
            else
            {
                // Leaf: store
                treeLeaf<Type>* subLeafPtr = getLeafPtr(octant);
                leafArray[leafIndex++] = subLeafPtr;
            }
        }
    }
}


template <class Type>
void treeNode<Type>::findLeaves
(
    List<const treeLeaf<Type>*>& leafArray,
    label& leafIndex
) const
{
    // Go into all sub boxes
    for(label octant=0; octant<8; octant++)
    {
        if (subNodes()[octant])
        {
            if (isNode(octant))
            {
                // Node: recurse into subnodes
                const treeNode<Type>* subNodePtr = getNodePtr(octant);
                subNodePtr->findLeaves(leafArray, leafIndex);
            }
            else
            {
                // Leaf: store
                treeLeaf<Type>* subLeafPtr = getLeafPtr(octant);
                leafArray[leafIndex++] = subLeafPtr;
            }
        }
    }
}


template <class Type>
void treeNode<Type>::printNode
(
    Ostream& os,
    const label level
) const
{
    space(os, 2*level);

    os << "node:" << this->bb() << endl;

    for(label octant=0; octant<8; octant++)
    {
        label type = getVolType(octant);

        string typeString = octree<Type>::volType(type);

        if (!subNodes_[octant])
        {
            space(os, level);
            os << octant << ":" << typeString << " : null" << endl;
        }
        else if (isNode(octant))
        {            
            space(os, level);
            os << octant << ":" << typeString << " : node" << endl;
            getNodePtr(octant)->printNode(os, level+1);
        }
        else
        {
            space(os, level);
            os << octant << ":" << typeString << " : leaf" << endl;

            treeLeaf<Type>* leafPtr = getLeafPtr(octant);
            leafPtr->printLeaf(os, level+1);
        }
    }
}


template <class Type>
void treeNode<Type>::writeOBJ
(
    Ostream& os,
    const label level,
    label& vertNo
) const
{
    point midPoint(this->bb().mid());

    label midVertNo = vertNo;
    os << "v " << midPoint.x() << " " << midPoint.y() << " "
       << midPoint.z() << endl;
    vertNo++;

    for(label octant=0; octant<8; octant++)
    {
        if (subNodes_[octant])
        {
            if (isNode(octant))
            {
                treeNode<Type>* nodePtr = getNodePtr(octant);

                point subMidPoint(nodePtr->bb().mid());
                os << "v " << subMidPoint.x() << " " << subMidPoint.y() << " "
                   << subMidPoint.z() << endl;
                os << "l " << midVertNo + 1<< " " << vertNo + 1 << endl;
                vertNo++;

                nodePtr->writeOBJ(os, level+1, vertNo);
            }
            else
            {
                treeLeaf<Type>* leafPtr = getLeafPtr(octant);

                point subMidPoint(leafPtr->bb().mid());
                os << "v " << subMidPoint.x() << " " << subMidPoint.y() << " "
                   << subMidPoint.z() << endl;
                os << "l " << midVertNo + 1<< " " << vertNo + 1 << endl;
                vertNo++;

                //leafPtr->writeOBJ(os, level+1, vertNo);
            }
        }
    }
}

 
// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class Type>
Istream& operator>>(Istream& is, treeNode<Type>& oc)
{
    for(label octant = 0; octant < 8; octant++)
    {
        oc.subNodes_[octant] = NULL;
    }

    is >> oc.bb();

    label nPtrs;

    // Read number of entries folllowing
    is >> nPtrs;

    is.readBegin("treeNode");
    for (label octant = 0; octant < nPtrs; octant++)
    {
        label index;
        is >> index;

        if (index >= treeNode<Type>::leafOffset)
        {
            // Leaf recognized by out of range index
            treeLeaf<Type>* leafPtr = new treeLeaf<Type>(is);
            oc.setLeafPtr(index - treeNode<Type>::leafOffset, leafPtr);
        }
        else
        {
            oc.setNodePtr(index, new treeNode<Type>(is));
        }
    }

    // Read end of treeNode list
    is.readEnd("treeNode");

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, treeNode&)");

    return is;
}


template <class Type>
Ostream& operator<<(Ostream& os, const treeNode<Type>& tn)
{
    // Count valid subnodes:
    //   - treeNode
    //   - treeLeafs with non-zero cell list.
    label nPtrs = 0;
    for (label octant = 0; octant < 8; octant++)
    {
        if (tn.subNodes_[octant])
        {
            if
            (
                tn.isNode(octant)
             || (tn.getLeafPtr(octant)->indices().size() != 0)
            )
            {
                nPtrs++;
            }
        }
    }


    // output subnodes as list of length nPtrs
    os << token::SPACE << tn.bb() << token::SPACE << nPtrs
       << token::SPACE << token::BEGIN_LIST;

    for (label octant = 0; octant < 8; octant++)
    {
        if (tn.subNodes_[octant])
        {
            if (tn.isNode(octant))
            {
                const treeNode<Type>* subNodePtr = tn.getNodePtr(octant);

                // Node: output index, value
                os  << token::SPACE << octant << token::SPACE << *subNodePtr
                    << token::NL;
            }
            else if (tn.getLeafPtr(octant)->indices().size() != 0)
            {
                // treeLeaf: mark by putting index invalid
                const treeLeaf<Type>* subLeafPtr = tn.getLeafPtr(octant);

                os  << token::SPACE << octant + treeNode<Type>::leafOffset
                    << token::SPACE << *subLeafPtr
                    << token::NL;
            }
        }
    }

    os << token::SPACE << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
