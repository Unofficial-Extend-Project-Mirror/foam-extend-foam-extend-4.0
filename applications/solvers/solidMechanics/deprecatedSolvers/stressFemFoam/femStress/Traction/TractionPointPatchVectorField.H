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
    TractionPointPatchVectorField. Traction boundary for stress analysis

Description

SourceFiles
    TractionPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef TractionPointPatchVectorField_H
#define TractionPointPatchVectorField_H

#include "tetPolyPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class TractionPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
class TractionPointPatchVectorField
:
    public PatchField<vector>
{
    // Private data

        //- Traction force
        vectorField traction_;

        //- Pressure force
        scalarField pressure_;


    // Private member functions

        void checkFieldSize() const;


public:

    //- Runtime type information
    TypeName("traction");


    // Constructors

        //- Construct from patch and internal field
        TractionPointPatchVectorField
        (
            const PointPatch&,
            const DimensionedField<vector, tetPointMesh>&
        );

        //- Construct from patch and internal field, traction and pressure
        TractionPointPatchVectorField
        (
            const PointPatch&,
            const DimensionedField<vector, tetPointMesh>&,
            const vectorField&,
            const scalarField&
        );

        //- Construct from patch, internal field and dictionary
        TractionPointPatchVectorField
        (
            const PointPatch&,
            const DimensionedField<vector, tetPointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchVectorField onto a new patch
        TractionPointPatchVectorField
        (
            const TractionPointPatchVectorField
                <PatchField, Mesh, PointPatch, MatrixType>&,
            const PointPatch&,
            const DimensionedField<vector, tetPointMesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<vector> > clone() const
        {
            return autoPtr<PatchField<vector> >
            (
                new TractionPointPatchVectorField
                <PatchField, Mesh, PointPatch, MatrixType>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        TractionPointPatchVectorField
        (
            const TractionPointPatchVectorField
                <PatchField, Mesh, PointPatch, MatrixType>&,
            const DimensionedField<vector, tetPointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<vector> > clone
        (
            const DimensionedField<vector, tetPointMesh>& iF
        ) const
        {
            return autoPtr<PatchField<vector> >
            (
                new TractionPointPatchVectorField
                <PatchField, Mesh, PointPatch, MatrixType>(*this, iF)
            );
        }


    // Member functions

        // Return defining fields

            virtual vectorField& traction()
            {
                return traction_;
            }

            virtual const vectorField& traction() const
            {
                return traction_;
            }


            virtual scalarField& pressure()
            {
                return pressure_;
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }


        // Mapping functions

            //- Reverse map the given PointPatchVectorField onto
            //  this PointPatchVectorField
            virtual void rmap
            (
                const PointPatchField
                    <PatchField, Mesh, PointPatch, MatrixType, vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Add boundary source for gradient-type conditions
            virtual void addBoundarySourceDiag
            (
                MatrixType<vector>&
            ) const;

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "TractionPointPatchVectorField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
