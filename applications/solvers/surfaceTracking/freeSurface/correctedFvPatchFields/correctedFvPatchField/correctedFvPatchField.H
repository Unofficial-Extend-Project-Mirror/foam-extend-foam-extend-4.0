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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    correctedFvPatchField

Description

SourceFiles
    correctedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef correctedFvPatchField_H
#define correctedFvPatchField_H

#include "fvPatchField.H"
#include "volMesh.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class correctedFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class correctedFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        //- Patch correction vectors
        vectorField correctionVectors_;

        //- Field differential along correction vectors
        Field<Type> corrVecGrad_;

        //- Field internal normal gradient
        Field<Type> nGradInternal_;

        //- Patch fvSubMesh
        mutable fvMeshSubset* patchSubMeshPtr_;

        //- Index of sub-mesh patch which correspondes to this patch
        mutable label subMeshPatchID_;

    // Private member functions

        void makePatchSubMesh() const;

        void updateCorrectionVectors();

        void updateCorrVecGrad();

public:

    //- Runtime type information
    TypeName("corrected");


    // Constructors

        //- Construct from patch and internal field
        correctedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        correctedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        correctedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fvPatchField onto a new patch
        correctedFvPatchField
        (
            const correctedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new correctedFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        correctedFvPatchField
        (
            const correctedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new correctedFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

        const GeometricField<Type, fvPatchField, volMesh>& volField() const;

        const Field<Type>& corrVecGrad() const
        {
            return corrVecGrad_;
        }

        const Field<Type>& nGradInternal() const
        {
            return nGradInternal_;
        }

        Field<Type>& nGradInternal()
        {
            return nGradInternal_;
        }

        const fvMeshSubset& patchSubMesh() const;

        label subMeshPatchID() const;

        void movePatchSubMesh();

        const vectorField& correctionVectors() const
        {
            return correctionVectors_;
        }

        //- Return patch-normal gradient
        virtual tmp<Field<Type> > snGrad() const;

        //- Return patch-normal derivative correction
        tmp<Field<Type> > snGradCorrection() const;

        //- Return internal field next to patch as patch field
//         virtual tmp<Field<Type> > patchInternalField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "correctedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
