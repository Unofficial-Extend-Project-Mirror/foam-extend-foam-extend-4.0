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

Description

\*---------------------------------------------------------------------------*/

#include "correctedFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "GeometricField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "gradScheme.H"
#include "gaussGrad.H"
#include "surfaceInterpolationScheme.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
correctedFvPatchField<Type>::correctedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    correctionVectors_(p.size(), vector::zero),
    corrVecGrad_(p.size(), pTraits<Type>::zero),
    nGradInternal_(p.size(), pTraits<Type>::zero),
    patchSubMeshPtr_(NULL),
    subMeshPatchID_(-1)
{
    updateCorrectionVectors();
}


template<class Type>
correctedFvPatchField<Type>::correctedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    fvPatchField<Type>(p, iF, f),
    correctionVectors_(p.size(), vector::zero),
    corrVecGrad_(p.size(), pTraits<Type>::zero),
    nGradInternal_(p.size(), pTraits<Type>::zero),
    patchSubMeshPtr_(NULL),
    subMeshPatchID_(-1)

{
    updateCorrectionVectors();
}


template<class Type>
correctedFvPatchField<Type>::correctedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF),
    correctionVectors_(p.size(), vector::zero),
    corrVecGrad_(p.size(), pTraits<Type>::zero),
    nGradInternal_(p.size(), pTraits<Type>::zero),
    patchSubMeshPtr_(NULL),
    subMeshPatchID_(-1)
{
    updateCorrectionVectors();

    if (dict.found("corrVecGrad"))
    {
        corrVecGrad_ =
            Field<Type>("corrVecGrad", dict, p.size());
    }
}


template<class Type>
correctedFvPatchField<Type>::correctedFvPatchField
(
    const correctedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    correctionVectors_(ptf.correctionVectors_, mapper),
    corrVecGrad_(ptf.corrVecGrad_, mapper),
    nGradInternal_(ptf.nGradInternal_, mapper),
    patchSubMeshPtr_(NULL),
    subMeshPatchID_(-1)
{
    updateCorrectionVectors();
}


template<class Type>
correctedFvPatchField<Type>::correctedFvPatchField
(
    const correctedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    correctionVectors_(ptf.correctionVectors_),
    corrVecGrad_(ptf.corrVecGrad_),
    nGradInternal_(ptf.nGradInternal_),
    patchSubMeshPtr_(NULL),
    subMeshPatchID_(-1)
{
    updateCorrectionVectors();
}


template<class Type>
void correctedFvPatchField<Type>::updateCoeffs()
{
    if(this->patch().boundaryMesh().mesh().moving())
    {
        updateCorrectionVectors();
    }

    updateCorrVecGrad();

    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
const GeometricField<Type, fvPatchField, volMesh>& 
correctedFvPatchField<Type>::volField() const
{
    return this->db().objectRegistry::lookupObject
        <
            GeometricField<Type, fvPatchField, volMesh> 
        >
        (this->dimensionedInternalField().name());
}


// Return gradient at boundary
template<class Type>
tmp<Field<Type> > correctedFvPatchField<Type>::snGrad() const
{
    return fvPatchField<Type>::snGrad()
      - this->patch().deltaCoeffs()*corrVecGrad_;

    // Second order correction
//     return 2.0*fvPatchField<Type>::snGrad()
//       - nGradInternal_
//       - 2.0*this->patch().deltaCoeffs()*corrVecGrad_;
}


// Return internal field next to patch as patch field
// template<class Type>
// tmp<Field<Type> > correctedFvPatchField<Type>::patchInternalField() const
// {
//     return fvPatchField<Type>::patchInternalField() + corrVecGrad_;
// }

// Return gradient at boundary
template<class Type>
tmp<Field<Type> > correctedFvPatchField<Type>::snGradCorrection() const
{
    //Second order correction
//     return fvPatchField<Type>::snGrad()
//       - nGradInternal_
//       - 2*this->patch().deltaCoeffs()*corrVecGrad_;

    return -this->patch().deltaCoeffs()*corrVecGrad_;
}


template<class Type>
void correctedFvPatchField<Type>::updateCorrectionVectors()
{
    const fvPatch& patch = this->patch();
    vectorField delta = patch.delta();

    correctionVectors_ =
        delta - patch.nf()*(patch.nf()&delta);
}


template<class Type>
void correctedFvPatchField<Type>::makePatchSubMesh() const
{
    if (patchSubMeshPtr_)
    {
        FatalErrorIn("correctedFvPatchField<Type>::makePatchSubMesh()")
            << "patch sub-mesh already exists"
            << abort(FatalError);
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    patchSubMeshPtr_ = new fvMeshSubset
    (
        IOobject
        (
            this->patch().name(),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    labelHashSet cellSet;

    const unallocLabelList& faceCells = this->patch().faceCells();

    for(label faceI=0; faceI<faceCells.size(); faceI++)
    {
        cellSet.insert(faceCells[faceI]);
    }

    const labelListList& cellCells = 
        this->patch().boundaryMesh().mesh().cellCells();

    for(label faceI=0; faceI<faceCells.size(); faceI++)
    {
        const labelList& curCells = cellCells[faceCells[faceI]];

        forAll(curCells, cellI)
        {
            if(!cellSet.found(curCells[cellI]))
            {
                cellSet.insert(curCells[cellI]);
            }
        }
    }

    patchSubMeshPtr_->setCellSubset(cellSet);
}


template<class Type>
const fvMeshSubset& correctedFvPatchField<Type>::patchSubMesh() const
{
    bool foundPatchSubMesh =
        this->db().objectRegistry::foundObject<fvMeshSubset>
        (
            this->patch().name()
        );
    
    if(foundPatchSubMesh)
    {
        return this->db().objectRegistry::lookupObject<fvMeshSubset>
        (
            this->patch().name()
        );
    }

    if(!patchSubMeshPtr_)
    {
        makePatchSubMesh();
    }

    return *patchSubMeshPtr_;
}


template<class Type>
label correctedFvPatchField<Type>::subMeshPatchID() const
{
    if(subMeshPatchID_ == -1)
    {
        forAll(patchSubMesh().patchMap(), patchI)
        {
            if(patchSubMesh().patchMap()[patchI] == this->patch().index())
            {
                subMeshPatchID_ = patchI;
                break;
            }
        }

        if(subMeshPatchID_ == -1)
        {
            if(this->patch().size()>0)
            {
                FatalErrorIn("correctedFvPatchField<Type>::subMeshPatchID()")
                    << "Can't determine the subMeshPatchID_"
                        << abort(FatalError);
            }
            else
            {
                subMeshPatchID_ = 0;
            }
        }
    }

    return subMeshPatchID_;
}


template<class Type>
void correctedFvPatchField<Type>::movePatchSubMesh()
{
    if(patchSubMeshPtr_)
    {
        const fvMesh& mesh = this->patch().boundaryMesh().mesh();

        vectorField newPoints(mesh.points(), patchSubMeshPtr_->pointMap());
        
        patchSubMeshPtr_->subMesh().movePoints(newPoints);
    }
}


template<class Type>
void correctedFvPatchField<Type>::updateCorrVecGrad()
{
    GeometricField<Type, fvPatchField, volMesh> subVolField =
        patchSubMesh().interpolate(volField());

    typedef typename 
        outerProduct<vector, typename pTraits<Type>::cmptType>::type 
        GradCmptType;

    vectorField n = this->patch().nf();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        GeometricField<GradCmptType, fvPatchField, volMesh> gradCmpt =
            fvc::grad
            (
                subVolField.component(cmpt),
                "snGradCorr(" + volField().name() + ')'
            );

        if(this->patch().size()>0)
        {
            corrVecGrad_.replace
            (
                cmpt,
                correctionVectors_
               &gradCmpt.boundaryField()[subMeshPatchID()].patchInternalField()
            );

            nGradInternal_.replace
            (
                cmpt,
                n&gradCmpt.boundaryField()[subMeshPatchID()].patchInternalField()
            );
        }
    }
}


// Write
template<class Type>
void correctedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    corrVecGrad_.writeEntry("corrVecGrad", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
