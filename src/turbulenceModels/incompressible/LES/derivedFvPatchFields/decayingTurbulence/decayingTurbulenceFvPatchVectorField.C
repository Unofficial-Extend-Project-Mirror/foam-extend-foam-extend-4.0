/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "decayingTurbulenceFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "scalarList.H"
#include "vectorList.H"
#include "ListListOps.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decayingTurbulenceFvPatchVectorField::updateVortons()
{
    vectorField& patchField = *this;

    vectorField tloc = patch().Cf();

    List<scalarList> l(Pstream::nProcs());
    List<vectorList> cf(Pstream::nProcs());
    List<vectorList> rf(Pstream::nProcs());

    l[Pstream::myProcNo()].setSize(LField_.size());
    cf[Pstream::myProcNo()].setSize(this->size());
    rf[Pstream::myProcNo()].setSize(refField_.size());

    label i = 0;

    forAll (tloc, faceI)
    {
        l[Pstream::myProcNo()][i]  = LField_[faceI];
        cf[Pstream::myProcNo()][i] = tloc[faceI];
        rf[Pstream::myProcNo()][i] = refField_[faceI];

        ++i;
    }

    Pstream::gatherList(l);
    Pstream::gatherList(cf);
    Pstream::gatherList(rf);

    scalarList L
    (
        ListListOps::combine<scalarList>(l, accessOp<scalarList>())
    );

    vectorList CF
    (
        ListListOps::combine<vectorList>(cf, accessOp<vectorList>())
    );

    vectorList RF
    (
        ListListOps::combine<vectorList>(rf, accessOp<vectorList>())
    );

    List<List<decayingVorton> > v(Pstream::nProcs());

    if (Pstream::master())
    {
        forAll (CF, I)
        {
            if (L[I] > SMALL)
            {
                scalar x = CF[I].x();
                scalar ls = lSpot(L[I]);
                scalar ymin = CF[I].y() - L[I];
                scalar zmin = CF[I].z() - L[I];

                for (label k = 0; k < nVortexes_; k++)
                {
                    vector v
                    (
                        (direction_ > 0) ? x - ls : x + ls,
                        ymin + rndGen_.scalar01()*2*L[I],
                        zmin + rndGen_.scalar01()*2*L[I]
                    );

                    bool allowed = true;

                    for
                    (
                        SLList<decayingVorton>::const_iterator it =
                            vortons_.begin();
                        it != vortons_.end();
                        ++it
                    )
                    {
                        if (mag(v - it().location()) < vortexOverlap_*ls)
                        {
                            allowed = false;
                            break;
                        }
                    }

                    if (!allowed)
                    {
                        continue;
                    }
                    else
                    {
                        vortons_.insert
                        (
                            decayingVorton
                            (
                                L[I],
                                v,
                                RF[I],
                                (direction_ > 0) ? x + ls : x - ls
                            )
                        );
                    }
                }
            }
        }

        v[Pstream::myProcNo()].setSize(vortons_.size());

        i = 0;

        for
        (
            SLList<decayingVorton>::iterator it = vortons_.begin();
            it != vortons_.end();
            ++it
        )
        {
            v[Pstream::myProcNo()][i++] = it();
        }
    }

    Pstream::scatterList(v);

    List<decayingVorton> V
    (
        ListListOps::combine<List<decayingVorton> >
        (
            v,
            accessOp<List<decayingVorton> >()
        )
    );

    if (!patchField.empty())
    {
        vectorField turbulent(refField_.size(), pTraits<vector>::zero);

        forAll (patchField, I)
        {
            vector pt = tloc[I];

            forAll (V, vI)
            {
                const decayingVorton& vorton = V[vI];

                if (mag(vorton.location() - pt) < lSpot(vorton.length()))
                {
                    turbulent[I] += vorton.velocityAt(pt);
                }
            }
        }

        R_ = ((index_ - 1.0)/index_)*R_ + (1.0/index_)*sqr(turbulent);

        tensorField C_(R_.size(), sphericalTensor::I);

        forAll (R_, I)
        {
            scalar D1 = R_[I].xx();

            if (D1 > 0)
            {
                C_[I].xx() = 1.0/sqrt(D1);
            }

            scalar D2 = R_[I].xx()*R_[I].yy() - R_[I].xy()*R_[I].xy();

            if (D1 > 0 && D2 > 0)
            {
                C_[I].yx() = -R_[I].xy()/sqrt(D1*D2);
                C_[I].yy() = sqrt(D1/D2);
            }

            scalar D3 = det(R_[I]);

            if (D2 > 0 && D3 > 0)
            {
                C_[I].zx() = (R_[I].xy()*R_[I].yz() - R_[I].yy()*R_[I].xz())/
                    sqrt(D2*D3);

                C_[I].zy() = -(R_[I].xx()*R_[I].yz()-R_[I].xz()*R_[I].xy())/
                    sqrt(D2*D3);

                C_[I].zz() = sqrt(D2/D3);
            }
        }

        index_++;

        turbulent = C_ & turbulent;
        turbulent = Lund_ & turbulent;

        fixedValueFvPatchVectorField::operator==(refField_+ turbulent);
    }

    if (Pstream::master())
    {
        for
        (
            SLList<decayingVorton>::iterator it = vortons_.begin();
            it != vortons_.end();
            ++it
        )
        {
            it().move(this->db().time().deltaT().value());
        }

        bool modified;

        do
        {
            modified = false;

            for
            (
                SLList<decayingVorton>::iterator it = vortons_.begin();
                it!=vortons_.end();
                ++it
            )
            {
                if (direction_ > 0)
                {
                    if (it().location().x() > it().xmax())
                    {
                        vortons_.remove(it);
                        modified = true;
                        break;
                    }
                }
                else
                {
                    if (it().location().x() < it().xmax())
                    {
                        vortons_.remove(it);
                        modified = true;
                        break;
                    }
                }
            }
        } while (modified);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decayingTurbulenceFvPatchVectorField::
decayingTurbulenceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    LField_(p.size()),
    refField_(p.size()),
    RField_(p.size()),
    Lund_(p.size(), pTraits<tensor>::zero),//HJ???
    R_(p.size()),
    direction_(1),
    vortexOverlap_(0.3),
    nVortexes_(10),
    curTimeIndex_(-1),
    rndGen_(label(0)),
    vortons_(),
    index_(2)
{}


Foam::decayingTurbulenceFvPatchVectorField::
decayingTurbulenceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    LField_("LField", dict, p.size()),
    refField_("refField", dict, p.size()),
    RField_("RField", dict, p.size()),
    Lund_(p.size(), pTraits<tensor>::zero),
    R_(RField_),
    direction_(readScalar(dict.lookup("direction"))),
    vortexOverlap_(dict.lookupOrDefault<scalar>("vortexOverlap", 0.3)),
    nVortexes_(dict.lookupOrDefault<label>("nVortexes", 10)),
    curTimeIndex_(-1),
    rndGen_(label(0)),
    vortons_(),
    index_(dict.lookupOrDefault<label>("index", 2))
{
    // Read optional items
    if (dict.found("vortons"))
    {
        vortons_ = SLList<decayingVorton>(dict.lookup("vortons"));
    }

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(refField_);
    }

    // Check diagonal of RField_: must be positive
    symmTensor minRField = min(RField_);

    if ((minRField.xx() < 0) || (minRField.yy() < 0) || (minRField.zz() < 0))
    {
        FatalErrorIn
        (
            "decayingTurbulenceFvPatchVectorField::"
            "decayingTurbulenceFvPatchVectorField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Negative RMS found for patch " << patch().name()
            << " and field " << dimensionedInternalField().name()
            << ": (" << minRField.xx() << " " << minRField.xx() << " "
            << minRField.xx() << " "
            << abort(FatalError);
    }

    // Set up Lund_
    Lund_.replace
    (
        tensor::XX,
        sqrt(Foam::max(RField_.component(symmTensor::XX), scalar(0)))
    );

    Lund_.replace
    (
        tensor::YX,
        RField_.component(symmTensor::XY)/Lund_.component(tensor::XX)
    );

    Lund_.replace
    (
        tensor::ZX,
        RField_.component(symmTensor::XZ)/Lund_.component(tensor::XX)
    );

    Lund_.replace
    (
        tensor::YY,
        sqrt(Foam::max(RField_.component(symmTensor::YY), scalar(0)))
      - sqr(Lund_.component(tensor::YX))
    );

    Lund_.replace
    (
        tensor::ZY,
        (
            RField_.component(symmTensor::YZ)
          - Lund_.component(tensor::YX)*Lund_.component(tensor::ZX)
        )/Lund_.component(tensor::YY)
    );

    Lund_.replace
    (
        tensor::ZZ,
        sqrt(Foam::max(RField_.component(symmTensor::ZZ), scalar(0))
      - sqr(Lund_.component(tensor::ZX))
      - sqr(Lund_.component(tensor::ZY)))
    );

    if (dict.found("R"))
    {
        R_ = symmTensorField("R", dict, p.size());
    }

    const scalarField& sf = patch().magSf();

    label numsf = 0;

    forAll (sf, faceI)
    {
        if (sf[faceI] > LField_[faceI]*LField_[faceI])
        {
            numsf++;
        }
    }

    // Check for integral length vs face size
    if (numsf > 0)
    {
        WarningIn
        (
            "decayingTurbulenceFvPatchVectorField::"
            "decayingTurbulenceFvPatchVectorField(dict)"
        )   << "For field " << dimensionedInternalField().name()
            << " and patch " << patch().index()
            << " found " << numsf <<" inlet faces (out of " << patch().size()
            << ") where the integral length is smaller than face size."
            << endl;
    }
}


Foam::decayingTurbulenceFvPatchVectorField::
decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    LField_(ptf.LField_, mapper),
    refField_(ptf.refField_, mapper),
    RField_(ptf.RField_, mapper),
    Lund_(ptf.Lund_, mapper),
    R_(ptf.R_, mapper),
    direction_(ptf.direction_),
    vortexOverlap_(ptf.vortexOverlap_),
    nVortexes_(ptf.nVortexes_),
    curTimeIndex_(-1),
    rndGen_(ptf.rndGen_),
    vortons_(ptf.vortons_),
    index_(ptf.index_)
{}


Foam::decayingTurbulenceFvPatchVectorField::
decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    LField_(ptf.LField_),
    refField_(ptf.refField_),
    RField_(ptf.RField_),
    Lund_(ptf.Lund_),
    R_(ptf.R_),
    direction_(ptf.direction_),
    vortexOverlap_(ptf.vortexOverlap_),
    nVortexes_(ptf.nVortexes_),
    curTimeIndex_(-1),
    rndGen_(ptf.rndGen_),
    vortons_(ptf.vortons_),
    index_(ptf.index_)
{}


Foam::decayingTurbulenceFvPatchVectorField::
decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    LField_(ptf.LField_),
    refField_(ptf.refField_),
    RField_(ptf.RField_),
    Lund_(ptf.Lund_),
    R_(ptf.R_),
    direction_(ptf.direction_),
    vortexOverlap_(ptf.vortexOverlap_),
    nVortexes_(ptf.nVortexes_),
    curTimeIndex_(-1),
    rndGen_(ptf.rndGen_),
    vortons_(ptf.vortons_),
    index_(ptf.index_)
{}


void Foam::decayingTurbulenceFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
    LField_.autoMap(m);
    refField_.autoMap(m);
    RField_.autoMap(m);
    R_.autoMap(m);
}


void Foam::decayingTurbulenceFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const decayingTurbulenceFvPatchVectorField& tiptf =
        refCast<const decayingTurbulenceFvPatchVectorField>(ptf);

    LField_.rmap(tiptf.LField_, addr);
    refField_.rmap(tiptf.refField_, addr);
    RField_.rmap(tiptf.RField_, addr);
    R_.rmap(tiptf.R_, addr);
}


void Foam::decayingTurbulenceFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        updateVortons();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::decayingTurbulenceFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    LField_.writeEntry("LField", os);
    refField_.writeEntry("refField", os);
    RField_.writeEntry("RField", os);
    os.writeKeyword("direction")<<direction_<<token::END_STATEMENT<<nl;
    writeEntry("value", os);
    R_.writeEntry("R", os);
    os.writeKeyword("index")<<index_<<token::END_STATEMENT<<nl;

    if (Pstream::master())
    {
        os.writeKeyword("vortons") << vortons_ << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        decayingTurbulenceFvPatchVectorField
    );
}


// ************************************************************************* //
