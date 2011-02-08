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
    Finite-Area matrix

\*---------------------------------------------------------------------------*/

#include "areaFields.H"
#include "edgeFields.H"
#include "calculatedFaPatchFields.H"
#include "zeroGradientFaPatchFields.H"
#include "coupledFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void faMatrix<Type>::addToInternalField
(
    const unallocLabelList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorIn
        (
            "faMatrix<Type>::addToInternalField(const unallocLabelList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] += pf[faceI];
    }
}


template<class Type>
template<class Type2>
void faMatrix<Type>::addToInternalField
(
    const unallocLabelList& addr,
    const tmp<Field<Type2> >& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
template<class Type2>
void faMatrix<Type>::subtractFromInternalField
(
    const unallocLabelList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorIn
        (
            "faMatrix<Type>::addToInternalField(const unallocLabelList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] -= pf[faceI];
    }
}


template<class Type>
template<class Type2>
void faMatrix<Type>::subtractFromInternalField
(
    const unallocLabelList& addr,
    const tmp<Field<Type2> >& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
void faMatrix<Type>::addBoundaryDiag
(
    scalarField& diag,
    const direction solveCmpt
) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            internalCoeffs_[patchI].component(solveCmpt),
            diag
        );
    }
}


template<class Type>
void faMatrix<Type>::addCmptAvBoundaryDiag(scalarField& diag) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            cmptAv(internalCoeffs_[patchI]),
            diag
        );
    }
}


template<class Type>
void faMatrix<Type>::addBoundarySource
(
    Field<Type>& source,
    const bool couples
) const
{
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];
        const Field<Type>& pbc = boundaryCoeffs_[patchI];

        if (!ptf.coupled())
        {
            addToInternalField(lduAddr().patchAddr(patchI), pbc, source);
        }
        else if (couples)
        {
            tmp<Field<Type> > tpnf = ptf.patchNeighbourField();
            const Field<Type>& pnf = tpnf();

            const unallocLabelList& addr = lduAddr().patchAddr(patchI);

            forAll(addr, facei)
            {
                source[addr[facei]] += cmptMultiply(pbc[facei], pnf[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
faMatrix<Type>::faMatrix
(
    GeometricField<Type, faPatchField, areaMesh>& psi,
    const dimensionSet& ds
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(ds),
    source_(psi.size(), pTraits<Type>::zero),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(NULL),
    solvingComponent(0)
{
    if (debug)
    {
        Info<< "faMatrix<Type>(GeometricField<Type, faPatchField, areaMesh>&,"
               " const dimensionSet&) : "
               "constructing faMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll (psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

    psi_.boundaryField().updateCoeffs();
}


template<class Type>
faMatrix<Type>::faMatrix(const faMatrix<Type>& fam)
:
    refCount(),
    lduMatrix(fam),
    psi_(fam.psi_),
    dimensions_(fam.dimensions_),
    source_(fam.source_),
    internalCoeffs_(fam.internalCoeffs_),
    boundaryCoeffs_(fam.boundaryCoeffs_),
    faceFluxCorrectionPtr_(NULL),
    solvingComponent(fam.solvingComponent)
{
    if (debug)
    {
        Info<< "faMatrix<Type>::faMatrix(const faMatrix<Type>&) : "
            << "copying faMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (fam.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, faePatchField, edgeMesh>
        (
            *(fam.faceFluxCorrectionPtr_)
        );
    }
}


template<class Type>
faMatrix<Type>::faMatrix
(
    GeometricField<Type, faPatchField, areaMesh>& psi,
    Istream& is
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(is),
    source_(is),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(NULL),
    solvingComponent(0)
{
    if (debug)
    {
        Info<< "faMatrix<Type>(GeometricField<Type, faPatchField, areaMesh>&,"
               " Istream&) : "
               "constructing faMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll (psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

}


template<class Type>
faMatrix<Type>::~faMatrix()
{
    if (debug)
    {
        Info<< "faMatrix<Type>::~faMatrix<Type>() : "
            << "destroying faMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (faceFluxCorrectionPtr_)
    {
        delete faceFluxCorrectionPtr_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set solution in given faces and eliminate corresponding
// equations from the matrix
template<class Type>
void faMatrix<Type>::setValues
(
    const labelList& faceLabels,
    const Field<Type>& values
)
{
    const faMesh& mesh = psi_.mesh();

    //    const faceList& faces = mesh.faces();
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();
    const labelListList& edges=mesh.patch().faceEdges();

    scalarField& Diag = diag();

    forAll(faceLabels, i)
    {
        label facei = faceLabels[i];

        psi_[facei] = values[i];
        source_[facei] = values[i]*Diag[facei];

        if (symmetric() || asymmetric())
        {
	  //            const face& c = faces[facei];
	    const labelList &c=edges[facei];

            forAll(c, j)
            {
                label edgei = c[j];

                if (mesh.isInternalEdge(edgei))
                {
                    if (symmetric())
                    {
                        if (facei == own[edgei])
                        {
                            source_[nei[edgei]] -= upper()[edgei]*values[i];
                        }
                        else
                        {
                            source_[own[edgei]] -= upper()[edgei]*values[i];
                        }

                        upper()[edgei] = 0.0;
                    }
                    else
                    {
                        if (facei == own[edgei])
                        {
                            source_[nei[edgei]] -= lower()[edgei]*values[i];
                        }
                        else
                        {
                            source_[own[edgei]] -= upper()[edgei]*values[i];
                        }

                        upper()[edgei] = 0.0;
                        lower()[edgei] = 0.0;
                    }
                }
                else
                {
                    label patchi = mesh.boundary().whichPatch(edgei);

                    if (internalCoeffs_[patchi].size())
                    {
                        label patchEdgei = 
                            mesh.boundary()[patchi].whichEdge(edgei);

                        internalCoeffs_[patchi][patchEdgei] = 
                            pTraits<Type>::zero;

                        boundaryCoeffs_[patchi][patchEdgei] =
                            pTraits<Type>::zero;
                    }
                }
            }
        }
    }
}


// Set reference level for solution
template<class Type>
void faMatrix<Type>::setReference
(
    const label facei,
    const Type& value
)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            source()[facei] += diag()[facei]*value;
            diag()[facei] += diag()[facei];
        }
    }
}


template<class Type>
void faMatrix<Type>::relax(const scalar alpha)
{
    if (alpha <= 0)
    {
        return;
    }

    Field<Type>& S = source();
    scalarField& D = diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalarField D0(D);

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(D.size(), 0.0);
    sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const unallocLabelList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                const Field<Type>& pCoeffs = boundaryCoeffs_[patchI];

                // For coupled boundaries add the diagonal and
                // off-diagonal contributions
                forAll(pa, face)
                {
                    D[pa[face]] += component(iCoeffs[face], 0);
                    sumOff[pa[face]] += mag(component(pCoeffs[face], 0));
                }
            }
            else
            {
                // For non-coupled boundaries subtract the diagonal
                // contribution off-diagonal sum which avoids having to remove
                // it from the diagonal later.
                // Also add the source contribution from the relaxation
                forAll(pa, face)
                {
                    Type iCoeff0 = iCoeffs[face];
                    iCoeffs[face] = cmptMag(iCoeffs[face]);
                    sumOff[pa[face]] -= cmptMin(iCoeffs[face]);
                    iCoeffs[face] /= alpha;
                    S[pa[face]] +=
                        cmptMultiply(iCoeffs[face] - iCoeff0, psi_[pa[face]]);
                }
            }
        }
    }

    // Ensure the matrix is diagonally dominant...
    max(D, D, sumOff);

    // ... then relax
    D /= alpha;

    // Now remove the diagonal contribution from coupled boundaries
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const unallocLabelList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                forAll(pa, face)
                {
                    D[pa[face]] -= component(iCoeffs[face], 0);
                }
            }
        }
    }

    // Finally add the relaxation contribution to the source.
    S += (D - D0)*psi_.internalField();
}


template<class Type>
void faMatrix<Type>::relax()
{
    scalar alpha = 0;

    if (psi_.mesh().relax(psi_.name()))
    {
        alpha = psi_.mesh().relaxationFactor(psi_.name());
    }

    if (alpha > 0)
    {
        relax(alpha);
    }
}


template<class Type>
tmp<scalarField> faMatrix<Type>::D() const
{
    tmp<scalarField> tdiag(new scalarField(diag()));
    addCmptAvBoundaryDiag(tdiag());
    return tdiag;
}


template<class Type>
tmp<areaScalarField> faMatrix<Type>::A() const
{
    tmp<areaScalarField> tAphi
    (
        new areaScalarField
        (
            IOobject
            (
                "A("+psi_.name()+')',
                psi_.instance(),
                psi_.db()
            ),
            psi_.mesh(),
            dimensions_/psi_.dimensions()/dimArea,
            zeroGradientFaPatchScalarField::typeName
        )
    );

    tAphi().internalField() = D()/psi_.mesh().S();
    tAphi().correctBoundaryConditions();

    return tAphi;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> > faMatrix<Type>::H() const
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > tHphi
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.db()
            ),
            psi_.mesh(),
            dimensions_/dimArea,
            zeroGradientFaPatchScalarField::typeName
        )
    );

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        scalarField psiCmpt = psi_.internalField().component(cmpt);

        scalarField boundaryDiagCmpt(psi_.size(), 0.0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);
        boundaryDiagCmpt.negate();
        addCmptAvBoundaryDiag(boundaryDiagCmpt);

        tHphi().internalField().replace(cmpt, boundaryDiagCmpt*psiCmpt);
    }

    tHphi().internalField() += lduMatrix::H(psi_.internalField()) + source_;
    addBoundarySource(tHphi().internalField());

    tHphi().internalField() /= psi_.mesh().S();
    tHphi().correctBoundaryConditions();

    return tHphi;
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> > faMatrix<Type>::
flux() const
{
    if (!psi_.mesh().fluxRequired(psi_.name()))
    {
        FatalErrorIn("faMatrix<Type>::flux()")
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary of faSchemes."
            << abort(FatalError);
    }

    // construct GeometricField<Type, faePatchField, edgeMesh>
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tfieldFlux
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "flux("+psi_.name()+')',
                psi_.instance(),
                psi_.db()
            ),
            psi_.mesh(),
            dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& fieldFlux = tfieldFlux();

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        fieldFlux.internalField().replace
        (
            cmpt,
            lduMatrix::faceH(psi_.internalField().component(cmpt))
        );
    }

    FieldField<Field, Type> InternalContrib = internalCoeffs_;

    forAll(InternalContrib, patchI)
    {
        InternalContrib[patchI] =
            cmptMultiply
            (
                InternalContrib[patchI],
                psi_.boundaryField()[patchI].patchInternalField()
            );
    }

    FieldField<Field, Type> NeighbourContrib = boundaryCoeffs_;

    forAll(NeighbourContrib, patchI)
    {
        if (psi_.boundaryField()[patchI].coupled())
        {
            NeighbourContrib[patchI] =
                cmptMultiply
                (
                    NeighbourContrib[patchI],
                    psi_.boundaryField()[patchI].patchNeighbourField()
                );
        }
    }

    forAll(fieldFlux.boundaryField(), patchI)
    {
        fieldFlux.boundaryField()[patchI] = 
            InternalContrib[patchI] - NeighbourContrib[patchI];
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }

    return tfieldFlux;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void faMatrix<Type>::operator=(const faMatrix<Type>& famv)
{
    if (this == &famv)
    {
        FatalErrorIn("faMatrix<Type>::operator=(const faMatrix<Type>&)")
            << "attempted to assignment to self"
            << abort(FatalError);
    }

    if (&psi_ != &(famv.psi_))
    {
        FatalErrorIn("faMatrix<Type>::operator=(const faMatrix<Type>&)")
            << "different fields"
            << abort(FatalError);
    }

    lduMatrix::operator=(famv);
    source_ = famv.source_;
    internalCoeffs_ = famv.internalCoeffs_;
    boundaryCoeffs_ = famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, faePatchField, edgeMesh>
        (*famv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void faMatrix<Type>::operator=(const tmp<faMatrix<Type> >& tfamv)
{
    operator=(tfamv());
    tfamv.clear();
}


template<class Type>
void faMatrix<Type>::negate()
{
    lduMatrix::negate();
    source_.negate();
    internalCoeffs_.negate();
    boundaryCoeffs_.negate();

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


template<class Type>
void faMatrix<Type>::operator+=(const faMatrix<Type>& famv)
{
    checkMethod(*this, famv, "+=");

    dimensions_ += famv.dimensions_;
    lduMatrix::operator+=(famv);
    source_ += famv.source_;
    internalCoeffs_ += famv.internalCoeffs_;
    boundaryCoeffs_ += famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, faePatchField, edgeMesh>
        (
            *famv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void faMatrix<Type>::operator+=(const tmp<faMatrix<Type> >& tfamv)
{
    operator+=(tfamv());
    tfamv.clear();
}


template<class Type>
void faMatrix<Type>::operator-=(const faMatrix<Type>& famv)
{
    checkMethod(*this, famv, "+=");

    dimensions_ -= famv.dimensions_;
    lduMatrix::operator-=(famv);
    source_ -= famv.source_;
    internalCoeffs_ -= famv.internalCoeffs_;
    boundaryCoeffs_ -= famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, faePatchField, edgeMesh>
        (-*famv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void faMatrix<Type>::operator-=(const tmp<faMatrix<Type> >& tfamv)
{
    operator-=(tfamv());
    tfamv.clear();
}


template<class Type>
void faMatrix<Type>::operator+=
(
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= su.mesh().S()*su.internalField();
}


template<class Type>
void faMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void faMatrix<Type>::operator-=
(
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(*this, su, "-=");
    source() += su.mesh().S()*su.internalField();
}


template<class Type>
void faMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void faMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    source() -= su.mesh().S()*su;
}


template<class Type>
void faMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    source() += su.mesh().S()*su;
}


template<class Type>
void faMatrix<Type>::operator*=
(
    const areaScalarField& vsf
)
{
    dimensions_ *= vsf.dimensions();
    lduMatrix::operator*=(vsf.internalField());
    source_ *= vsf.internalField();

    forAll(boundaryCoeffs_, patchI)
    {
        scalarField psf = vsf.boundaryField()[patchI].patchInternalField();
        internalCoeffs_[patchI] *= psf;
        boundaryCoeffs_[patchI] *= psf;
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorIn("faMatrix<Type>::operator*=(const tmp<areaScalarField>&)")
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void faMatrix<Type>::operator*=
(
    const tmp<areaScalarField>& tvsf
)
{
    operator*=(tvsf());
    tvsf.clear();
}


template<class Type>
void faMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    lduMatrix::operator*=(ds.value());
    source_ *= ds.value();
    internalCoeffs_ *= ds.value();
    boundaryCoeffs_ *= ds.value();

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ *= ds.value();
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const faMatrix<Type>& fam1,
    const faMatrix<Type>& fam2,
    const char* op
)
{
    if (&fam1.psi() != &fam2.psi())
    {
        FatalErrorIn
        (
            "checkMethod(const faMatrix<Type>&, const faMatrix<Type>&)"
        )   << "incompatible fields for operation "
            << endl << "    "
            << "[" << fam1.psi().name() << "] "
            << op
            << " [" << fam2.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && fam1.dimensions() != fam2.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const faMatrix<Type>&, const faMatrix<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fam1.psi().name() << fam1.dimensions()/dimArea << " ] "
            << op
            << " [" << fam2.psi().name() << fam2.dimensions()/dimArea << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const faMatrix<Type>& fam,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const char* op
)
{
    if (dimensionSet::debug && fam.dimensions()/dimArea != vf.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const faMatrix<Type>&, const GeometricField<Type, "
            "faPatchField, areaMesh>&)"
        )   <<  "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fam.psi().name() << fam.dimensions()/dimArea << " ] "
            << op
            << " [" << vf.name() << vf.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const faMatrix<Type>& fam,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if (dimensionSet::debug && fam.dimensions()/dimArea != dt.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const faMatrix<Type>&, const dimensioned<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fam.psi().name() << fam.dimensions()/dimArea << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
lduSolverPerformance solve
(
    faMatrix<Type>& fam,
    Istream& solverControls
)
{
    return fam.solve(solverControls);
}

template<class Type>
lduSolverPerformance solve
(
    const tmp<faMatrix<Type> >& tfam,
    Istream& solverControls
)
{
    lduSolverPerformance solverPerf = 
        const_cast<faMatrix<Type>&>(tfam()).solve(solverControls);

    tfam.clear();
    return solverPerf;
}


template<class Type>
lduSolverPerformance solve(faMatrix<Type>& fam)
{
    return fam.solve();
}

template<class Type>
lduSolverPerformance solve(const tmp<faMatrix<Type> >& tfam)
{
    lduSolverPerformance solverPerf =
        const_cast<faMatrix<Type>&>(tfam()).solve();

    tfam.clear();
    return solverPerf;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC() += B;
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC() += B;
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<faMatrix<Type> > tC(tB.ptr());
    tC() += A;
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >& tA,
    const tmp<faMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC() += tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().negate();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >& tA
)
{
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().negate();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC() -= B;
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC() -= B;
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<faMatrix<Type> > tC(tB.ptr());
    tC() -= A;
    tC().negate();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >& tA,
    const tmp<faMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}


template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}


template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >& tA,
    const tmp<faMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}


template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>& A,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() -= su.mesh().S()*su.internalField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >& tA,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().S()*su.internalField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() -= su.mesh().S()*su.internalField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const tmp<faMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().S()*su.internalField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu,
    const tmp<faMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>& A,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() += su.mesh().S()*su.internalField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >& tA,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().S()*su.internalField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() += tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const tmp<faMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.mesh().S()*su.internalField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().negate();
    tC().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu,
    const tmp<faMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const tmp<faMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() += su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() += su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const tmp<faMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>& A,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() += su.mesh().S()*su.internalField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >& tA,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().S()*su.internalField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() += tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC().source() += A.psi().mesh().S()*su.value();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC().source() += tC().psi().mesh().S()*su.value();
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator*
(
    const areaScalarField& vsf,
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC() *= vsf;
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator*
(
    const tmp<areaScalarField>& tvsf,
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC() *= tvsf;
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator*
(
    const areaScalarField& vsf,
    const tmp<faMatrix<Type> >& tA
)
{
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC() *= vsf;
    return tC;
}

template<class Type>
tmp<faMatrix<Type> > operator*
(
    const tmp<areaScalarField>& tvsf,
    const tmp<faMatrix<Type> >& tA
)
{
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC() *= tvsf;
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type> > tC(new faMatrix<Type>(A));
    tC() *= ds;
    return tC;
}


template<class Type>
tmp<faMatrix<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const tmp<faMatrix<Type> >& tA
)
{
    tmp<faMatrix<Type> > tC(tA.ptr());
    tC() *= ds;
    return tC;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const faMatrix<Type>& fam)
{
    os  << static_cast<const lduMatrix&>(fam) << nl
        << fam.dimensions_ << nl
        << fam.source_ << nl
        << fam.internalCoeffs_ << nl
        << fam.boundaryCoeffs_ << endl;

    os.check("Ostream& operator<<(Ostream&, faMatrix<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

#include "faMatrixSolve.C"

// ************************************************************************* //
