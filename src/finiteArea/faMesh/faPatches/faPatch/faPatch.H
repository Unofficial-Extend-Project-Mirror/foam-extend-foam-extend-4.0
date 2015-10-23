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

Class
    faPatch

Description
    Finite area patch class.  Used for 2-D non-Euclidian finite area method.

SourceFiles
    faPatch.C
    newFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef faPatch_H
#define faPatch_H

#include "patchIdentifier.H"
#include "labelList.H"
#include "pointField.H"
#include "typeInfo.H"
#include "faPatchFieldsFwd.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class faPatch Declaration
\*---------------------------------------------------------------------------*/

class faPatch
:
    public labelList,
    public patchIdentifier
{
private:

    // Private data

        //- Neighbour polyPatch index
        const label ngbPolyPatchIndex_;

        //- Reference to boundary mesh
        const faBoundaryMesh& boundaryMesh_;


        // Demand-driven private data

            //- Edge-face addressing
            mutable labelList::subList* edgeFacesPtr_;

            //- Local points labels
            mutable labelList* pointLabelsPtr_;

            //- Point-edge addressing
            mutable labelListList* pointEdgesPtr_;


    // Private Member Functions

        //- Disallow construct as copy
        faPatch(const faPatch&);

        //- Disallow assignment
        void operator=(const faPatch&);


       //- Clear out topological patch data
       void clearOut();


protected:

        // The faPatch geometry initialisation is called by faBoundaryMesh
        friend class faBoundaryMesh;

        //- Calculate patch point labels
        void calcPointLabels() const;

        //- Calculate patch point-edge addressing
        void calcPointEdges() const;

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
        virtual void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh()
        {}

        //- Update of the patch topology
        virtual void updateMesh()
        {}


public:

    typedef faBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("patch");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            faPatch,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const faBoundaryMesh& bm
            ),
            (name, dict, index, bm)
        );


    // Constructors

        //- Construct from components
        faPatch
        (
            const word& name,
            const labelList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label ngbPolyPatchIndex
        );

        //- Construct from dictionary
        faPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        );

        //- Construct as copy, resetting the boundary mesh
        faPatch(const faPatch&, const faBoundaryMesh&);

        //- Construct and return a clone, resetting the edge list
        //  and boundary mesh
        virtual autoPtr<faPatch> clone
        (
            const faBoundaryMesh& bm,
            const labelList& edgeLabels,
            const label index,
            const label ngbPolyPatchIndex
        ) const
        {
            return autoPtr<faPatch>
            (
                new faPatch(name(), edgeLabels, index, bm, ngbPolyPatchIndex)
            );
        }


    // Selectors

        //- Return a pointer to a new patch created
        // on freestore from dictionary
        static autoPtr<faPatch> New
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        );


    // Destructor

        virtual ~faPatch();


    // Member Functions

        //- Return number of patch points
        label nPoints() const
        {
            return pointLabels().size();
        }

        //- Return neighbour polyPatch index
        label ngbPolyPatchIndex() const;

        //- Return boundaryMesh reference
        const faBoundaryMesh& boundaryMesh() const;

        //- Return true if this patch is coupled
        virtual bool coupled() const
        {
            return false;
        }

        //- Patch start in edge list
        label start() const;

        //- Patch size
        virtual label size() const
        {
            return labelList::size();
        }

        //- Return label of edge in patch from global edge label
        inline label whichEdge(const label l) const
        {
            return l - start();
        }

        //- Slice list to patch
        template<class T>
        typename List<T>::subList patchSlice(const List<T>& l) const
        {
            return typename List<T>::subList(l, size(), start());
        }


        //- Write
        virtual void write(Ostream&) const;


        // Acces functions for geometrical data

            //- Return patch point labels
            const labelList& pointLabels() const;

            //- Return patch point-edge addressing
            const labelListList& pointEdges() const;

            //- Return edge neighbour polyPatch faces
            labelList ngbPolyPatchFaces() const;

            //- Return normals of neighbour polyPatch faces
            tmp<vectorField> ngbPolyPatchFaceNormals() const;

            //- Return normals of neighbour polyPatch joined points
            tmp<vectorField> ngbPolyPatchPointNormals() const;

            //- Return edge-face addressing
            const unallocLabelList& edgeFaces() const;

            //- Return edge centres
            const vectorField& edgeCentres() const;

            //- Return edge length vectors
            const vectorField& edgeLengths() const;

            //- Return edge length magnitudes
            const scalarField& magEdgeLengths() const;

            //- Return edge normals
            tmp<vectorField> edgeNormals() const;

            //- Return neighbour face centres
            tmp<vectorField> edgeFaceCentres() const;

            //- Return cell-centre to face-centre vector
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre vector is returned
            virtual tmp<vectorField> delta() const;


        // Access functions for demand driven data

            //- Make patch weighting factors
            virtual void makeWeights(scalarField&) const;

            //- Return patch weighting factors
            const scalarField& weights() const;

            //- Make patch edge - neighbour face distances
            virtual void makeDeltaCoeffs(scalarField&) const;

            //- Return patch edge - neighbour face distances
            const scalarField& deltaCoeffs() const;


        // Topological changes

            //- Reset edge list
            void resetEdges(const labelList&);


        // Evaluation functions

            //- Return given internal field next to patch as patch field
            template<class Type>
            tmp<Field<Type> > patchInternalField(const UList<Type>&) const;

            //- Return the corresponding patchField of the named field
            template<class GeometricField, class Type>
            const typename GeometricField::PatchFieldType& patchField
            (
                const GeometricField&
            ) const;

            //- Lookup and return the patchField of the named field from the
            //  local objectRegistry.
            //  N.B.  The dummy pointer arguments are used if this function is
            //  instantiated within a templated function to avoid a bug in gcc.
            //  See inletOutletFvPatchField.C and outletInletFvPatchField.C
            template<class GeometricField, class Type>
            const typename GeometricField::PatchFieldType& lookupPatchField
            (
                const word& name,
                const GeometricField* = NULL,
                const Type* = NULL
            ) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const faPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
