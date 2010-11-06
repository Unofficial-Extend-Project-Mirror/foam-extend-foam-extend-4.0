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

Class
    ePatch

SourceFiles
    ePatch.C
    newEPatch.C

\*---------------------------------------------------------------------------*/

#ifndef ePatch_H
#define ePatch_H

#include "patchIdentifier.H"
#include "labelList.H"
#include "pointField.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class ePatch Declaration
\*---------------------------------------------------------------------------*/

class ePatch
:
    public patchIdentifier
{
private:

    // Private data

        //- Reference to boundary mesh
        const eBoundaryMesh& boundaryMesh_;

        //- Start index of the patch
        label start_;

        //- Size of the patch
        label size_;


protected:

        // The ePatch geometry initialisation is called by eBoundaryMesh
        friend class eBoundaryMesh;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry()
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry()
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&)
        {}

        //- Correct patch after moving points
        virtual void movePoints(const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh()
        {}

        //- Update of the patch topology
        virtual void updateMesh()
        {}


public:

    typedef eBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("patch");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            ePatch,
            word,
            (
                const word& name,
                const label size,
                const label start,
                const label index,
                const eBoundaryMesh& bm
            ),
            (name, size, start, index, bm)
        );


        declareRunTimeSelectionTable
        (
            autoPtr,
            ePatch,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const eBoundaryMesh& bm
            ),
            (name, dict, index, bm)
        );


    // Constructors

        //- Construct from components
        ePatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const eBoundaryMesh& bm
        );


        //- Construct from dictionary
        ePatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const eBoundaryMesh& bm
        );

        //- Construct as copy
        ePatch(const ePatch&);

        //- Construct as copy, resetting the boundary mesh
        ePatch(const ePatch&, const eBoundaryMesh&);


    // Selectors

        //- Return a pointer to a new patch created on freestore from
        //  components
        static autoPtr<ePatch> New
        (
            const word& patchType,
            const word& name,
            const label size,
            const label start,
            const label index,
            const eBoundaryMesh& bm
        );


        //- Return a pointer to a new patch created
        // on freestore from dictionary
        static autoPtr<ePatch> New
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const eBoundaryMesh& bm
        );


    // Destructor

        virtual ~ePatch();


    // Member Functions

        //- Return the index of this patch in the boundaryMesh
        label index() const
        {
            return patchIdentifier::index();
        }

        //- Return boundaryMesh reference
        const eBoundaryMesh& boundaryMesh() const;

        //- Return true if this patch is coupled
        virtual bool coupled() const
        {
            return false;
        }

        //- Patch start in edge list
        label start() const
        {
            return start_;
        }

        //- Patch size
        virtual label size() const
        {
            return size_;
        }

        //- Slice list to patch
        template<class T>
        typename List<T>::subList patchSlice(const List<T>& l) const
        {
            return typename List<T>::subList(l, size(), start());
        }

        //- Write
        virtual void write(Ostream&) const;

        //- Assignment operator
        void operator=(const ePatch&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const ePatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
