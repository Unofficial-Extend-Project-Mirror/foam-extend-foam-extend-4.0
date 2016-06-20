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

\*---------------------------------------------------------------------------*/

#include "multiMaterialCohesiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "crackerFvMesh.H"
#include "multiMaterial.H"
#include "constitutiveModel.H"
#include "cohesiveFvPatch.H"
#include "cohesiveLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMaterialCohesiveLaw, 0);
    addToRunTimeSelectionTable
    (
        cohesiveLaw,
        multiMaterialCohesiveLaw,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::multiMaterialCohesiveLaw::indicator
(
    const label i
) const
{
    const scalarField& mat = materials_.internalField();

    tmp<scalarField> tresult(new scalarField(mat.size(), 0.0));
    scalarField& result = tresult();

    forAll (mat, matI)
    {
        if (mat[matI] > i - SMALL && mat[matI] < i + 1 - SMALL)
        {
            result[matI] = 1.0;
        }
    }

    return tresult;
}


Foam::scalar Foam::multiMaterialCohesiveLaw::indicator
(
    const label index,
    const label cellID
) const
{
    const scalar mat = materials_.internalField()[cellID];
    scalar result = 0.0;

    if (mat > index - SMALL && mat < index + 1 - SMALL)
      {
    result = 1.0;
      }

    return result;
}

Foam::scalar Foam::multiMaterialCohesiveLaw::interfaceID
(
    const label mat1,
    const label mat2
) const
{
  word mat1name = (*this)[mat1].name();
  word mat2name = (*this)[mat2].name();

  word interfaceName("interface_"+mat1name+"_"+mat2name);
  label interfaceLawID = -1;
  forAll(interfaceCohesiveLaws_, lawI)
  {
      if (interfaceCohesiveLaws_[lawI].name() == interfaceName)
      {
          interfaceLawID = lawI;
          break;
      }
  }
  if (interfaceLawID == -1)
  {
      // flip name
      interfaceName = word("interface_"+mat2name+"_"+mat1name);
      forAll(interfaceCohesiveLaws_, lawI)
      {
          if (interfaceCohesiveLaws_[lawI].name() == interfaceName)
          {
              interfaceLawID = lawI;
              break;
          }
      }
      if (interfaceLawID == -1)
      {
          FatalError
              << "Cannot find cohesive interfaceLaw "
              << word("interface_"+mat1name+"_"+mat2name) << " or "
              << interfaceName << nl
              << "One of these should be defined!"
              << exit(FatalError);
      }
  }

  return interfaceLawID;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::multiMaterialCohesiveLaw::multiMaterialCohesiveLaw
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    cohesiveLaw(name, sigma, dict),
    PtrList<cohesiveLaw>(),
    materials_
    (
        sigma.mesh().objectRegistry::lookupObject<volScalarField>("materials")
    ),
    // (
    //  IOobject
    //  (
    //   "materials",
    //   mesh().time().timeName(),
    //   mesh(),
    //   IOobject::MUST_READ,
    //   IOobject::AUTO_WRITE
    //   ),
    //  mesh()
    //  ),
    interfaceCohesiveLaws_() //,
    // materialsSurf_
    // (
    //  IOobject
    //  (
    //   "materialsSurf",
    //   mesh().time().timeName(),
    //   mesh(),
    //   IOobject::NO_READ,
    //   IOobject::NO_WRITE
    //   ),
    //  mesh(),
    //  dimensionedScalar("zero", dimless, 0.0)
    //  )
{
    PtrList<cohesiveLaw>& laws = *this;

    PtrList<entry> lawEntries(dict.lookup("laws"));
    laws.setSize(lawEntries.size());

    forAll (laws, lawI)
    {
        laws.set
        (
            lawI,
            cohesiveLaw::New
            (
                lawEntries[lawI].keyword(),
                sigma,
                lawEntries[lawI].dict()
            )
        );
    }

    if
    (
        min(materials_).value() < 0
     || max(materials_).value() > laws.size() + SMALL
    )
    {
        FatalErrorIn
        (
            "multiMaterialCohesiveLaw::multiMaterialCohesiveLaw\n"
            "(\n"
            "    const word& name,\n"
            "    const volSymmTensorField& sigma,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Invalid definition of material indicator field.  "
            << "Number of materials: " << laws.size()
            << " max index: " << max(materials_)
            << ".  Should be " << laws.size() - 1
            << abort(FatalError);
    }

    // cohesive laws must be the same size as rheology laws
    const constitutiveModel& conModel =
        mesh().objectRegistry::lookupObject<constitutiveModel>
        ("rheologyProperties");
    if (conModel.law().type() == "multiMaterial")
      {
    const multiMaterial& mulMatLaw =
        refCast<const multiMaterial>(conModel.law());

    if (laws.size() != mulMatLaw.size())
      {
        FatalError
            << "There should be the same number of cohesive laws "
            << "as rheology laws" << nl
            << "Currently there are " << mulMatLaw.size()
            << " rheology laws "
            << "and " << laws.size() << " cohesive laws!"
            << exit(FatalError);
      }
      }

    // set interfaceCohesiveLaws_
    PtrList<entry> interfaceLawEntries(dict.lookup("interfaceLaws"));
    // if (interfaceLawEntries.size() != int(laws.size()*(laws.size()-1)/2))
    if ( mag(interfaceLawEntries.size() - (laws.size()*(laws.size()-1)/2))
        > SMALL )
      {
    // number of interfaces is a trianular number of number of materials
    // ((n)*(n-1)/2)
          FatalError
              << "There are " << interfaceLawEntries.size()
              << " interface cohesive"
              << " laws defined, but there should be "
              << (laws.size()*(laws.size()-1)/2)
              << "," << nl
              << "as there are " << laws.size()
              << " materials in cohesive laws" << nl
              << exit(FatalError);
      }
    interfaceCohesiveLaws_.setSize(interfaceLawEntries.size());
    forAll (interfaceCohesiveLaws_, lawI)
      {
        interfaceCohesiveLaws_.set
      (
       lawI,
       cohesiveLaw::New
       (
        interfaceLawEntries[lawI].keyword(),
        sigma,
        interfaceLawEntries[lawI].dict()
            )
       );
      }


    // Set materialsSurf
    // materialsSurf_ = fvc::interpolate(materials_);
    // forAll(mesh().boundary(), patchi)
    //   {
    //      materialsSurf_.boundaryField()[patchi] =
    //        materials_.boundaryField()[patchi].patchInternalField();
    //   }
    // // Fix interface values
    // const labelList owner = mesh().owner();
    // const labelList neighbour = mesh().neighbour();
    // forAll (materialsSurf_.internalField(), faceI)
    //   {
    //      // round value to integer and check difference
    //      // if it is small then the face is not on a multi-material
    //      // interface
    //      scalar matIDscalar = materialsSurf_.internalField()[faceI];
    //      label matID = int(matIDscalar);
    //      if (mag(matIDscalar - matID) > SMALL)
    //        {
    //          // find which interface it is on
    //          const label own = owner[faceI];
    //          const label nei = neighbour[faceI];

    //          materialsSurf_.internalField()[faceI] =
    //            interfaceID(materials_[own], materials_[nei]) + laws.size();
    //        }
    //   }

    // philipc
    // processor boundaries are meant to hold the patchNeighbourField
    // but the proc boundaries are interpolated by decomposePar
    // so we must correct them. Now the proc boundaries hold the
    // patchNiehgbourField
    //materials_.correctBoundaryConditions(); // done by rheologyLaw
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiMaterialCohesiveLaw::~multiMaterialCohesiveLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::multiMaterialCohesiveLaw::materials() const
{
    return materials_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::multiMaterialCohesiveLaw::sigmaMax() const
{
    tmp<surfaceScalarField> tresult
    (
        new surfaceScalarField
        (
            IOobject
            (
            "sigmaMax",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
            mesh(),
        dimensionedScalar("zero", dimForce/dimArea, 0)
    )
    );
    surfaceScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<cohesiveLaw>& laws = *this;
    volScalarField indic
      (
       IOobject
       (
    "indic",
    mesh().time().timeName(),
    mesh(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
       mesh(),
       dimensionedScalar("zero", dimless, 0),
       zeroGradientFvPatchScalarField::typeName
       );
    forAll (laws, lawI)
      {
    // interface fields will not be correct but we fix them after
    indic.internalField() = indicator(lawI)();
    surfaceScalarField indicatorSurf = fvc::interpolate(indic);
    // fix boundary fields
    forAll(mesh().boundary(), patchi)
      {
        indicatorSurf.boundaryField()[patchi] =
          indic.boundaryField()[patchi].patchInternalField();
      }
        result += indicatorSurf*laws[lawI].sigmaMax()();

      }

    // fix interfaces
    // forAll surfaces check if surface is a material interface
    // material indicator should read non integer
    // Get the two materials it is an interface of
    // Look up value of sigmaMaxf in dictionary
    // Overwrite existing value on surface with new value

    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();

    forAll(result.internalField(), faceI)
    {
      label ownerMat = label(materials_[owner[faceI]]);
      label neighbourMat = label(materials_[neighbour[faceI]]);

        if (ownerMat != neighbourMat)
        {
      result.internalField()[faceI] = interfaceSigmaMax(ownerMat, neighbourMat);
    }
    }

    forAll(mesh().boundary(), patchI)
    {

        if (mesh().boundary()[patchI].type() == cohesiveFvPatch::typeName)
        {
      const crackerFvMesh& crackerMesh = refCast<const crackerFvMesh>(mesh());
      // we must use owner values
      scalarField localPatchMaterials =
          materials_.boundaryField()[patchI].patchInternalField();
      scalarField globalPatchMaterials =
          crackerMesh.globalCrackField(localPatchMaterials);
      // crackerMesh.globalCrackField(materials_.boundaryField()[patchI]);
      const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
      label globalIndex = crackerMesh.localCrackStart();

      forAll(mesh().boundaryMesh()[patchI], facei)
        {
          label ownerMat = label(globalPatchMaterials[globalIndex]);
          label neighbourMat = label(globalPatchMaterials[gcfa[globalIndex]]);

          if (ownerMat != neighbourMat)
                {
                    result.boundaryField()[patchI][facei] =
                        interfaceSigmaMax(ownerMat, neighbourMat);
                }
          globalIndex++;
            }
        }

    // philipc
    else if (mesh().boundary()[patchI].coupled())
      {
        const scalarField ownerMatField =
            materials_.boundaryField()[patchI].patchInternalField();
        const scalarField neighbourMatField =
            materials_.boundaryField()[patchI].patchNeighbourField();
        //const labelList& faceCells = mesh().boundary()[patchI].faceCells();

        forAll(mesh().boundaryMesh()[patchI], facei)
          {
        label ownerMat = label(ownerMatField[facei]);
        label neighbourMat = label(neighbourMatField[facei]);

        if (ownerMat != neighbourMat)
          {
              // result.boundaryField()[patchI][facei] = sigmaMax();
              result.boundaryField()[patchI][facei] =
                  interfaceSigmaMax(ownerMat, neighbourMat);
          }
          }
      }
    }

    return tresult;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiMaterialCohesiveLaw::tauMax() const
{
    tmp<surfaceScalarField> tresult
    (
        new surfaceScalarField
        (
            IOobject
            (
            "tauMax",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
            mesh(),
        dimensionedScalar("zero", dimForce/dimArea, 0)
    )
    );
    surfaceScalarField& result = tresult();

    const PtrList<cohesiveLaw>& laws = *this;
    volScalarField indic
      (
       IOobject
       (
    "indic",
    mesh().time().timeName(),
    mesh(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
       mesh(),
       dimensionedScalar("zero", dimless, 0),
       zeroGradientFvPatchScalarField::typeName
       );
    forAll (laws, lawI)
      {
    // interface fields will not be correct but we fix them after
    indic.internalField() = indicator(lawI)();
    surfaceScalarField indicatorSurf = fvc::interpolate(indic);
    // fix boundary fields
    forAll(mesh().boundary(), patchi)
      {
        indicatorSurf.boundaryField()[patchi] =
          indic.boundaryField()[patchi].patchInternalField();
      }
        result += indicatorSurf*laws[lawI].tauMax()();
      }

    // forAll surfaces check if surface is a material interface
    // material indicator should read non integer
    // Get the two materials it is an interface of
    // Look up value of tauMaxf in dictionary
    // Overwrite existing value on surface with new value

    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();

    forAll(result.internalField(), faceI)
    {
      label ownerMat = label(materials_[owner[faceI]]);
      label neighbourMat = label(materials_[neighbour[faceI]]);

        if (ownerMat != neighbourMat)
        {
      result.internalField()[faceI] = interfaceTauMax(ownerMat, neighbourMat);
    }
    }

    forAll(mesh().boundary(), patchI)
    {

        if (mesh().boundary()[patchI].type() == cohesiveFvPatch::typeName)
        {
      //label size = (mesh().boundary()[patchI].size())/2;
      const crackerFvMesh& crackerMesh = refCast<const crackerFvMesh>(mesh());
      // we must use owner values
      scalarField localPatchMaterials =
          materials_.boundaryField()[patchI].patchInternalField();
      scalarField globalPatchMaterials =
          crackerMesh.globalCrackField(localPatchMaterials);
      // crackerMesh.globalCrackField(materials_.boundaryField()[patchI]);
      const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
      label globalIndex = crackerMesh.localCrackStart();

      forAll(mesh().boundaryMesh()[patchI], facei)
        {
          label ownerMat = label(globalPatchMaterials[globalIndex]);
          label neighbourMat = label(globalPatchMaterials[gcfa[globalIndex]]);

          if (ownerMat != neighbourMat)
        {
          result.boundaryField()[patchI][facei] =
              interfaceTauMax(ownerMat, neighbourMat);
                }
          globalIndex++;
            }
        }

    // philipc
    else if (mesh().boundary()[patchI].coupled())
      {
        const scalarField ownerMatField =
            materials_.boundaryField()[patchI].patchInternalField();
        const scalarField neighbourMatField =
            materials_.boundaryField()[patchI].patchNeighbourField();
        //const labelList& faceCells = mesh().boundary()[patchI].faceCells();

        forAll(mesh().boundaryMesh()[patchI], facei)
          {
        label ownerMat = label(ownerMatField[facei]);
        label neighbourMat = label(neighbourMatField[facei]);

        if (ownerMat != neighbourMat)
          {
            result.boundaryField()[patchI][facei] =
                interfaceTauMax(ownerMat, neighbourMat);
          }
          }
      }
    }

    return tresult;
}

Foam::tmp<Foam::surfaceScalarField> Foam::multiMaterialCohesiveLaw::GIc() const
{
    tmp<surfaceScalarField> tresult
    (
        new surfaceScalarField
        (
            IOobject
            (
                "GIc",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
        dimensionedScalar("zero", dimForce*dimLength/dimArea, 0)
    )
    );
    surfaceScalarField& result = tresult();

    const PtrList<cohesiveLaw>& laws = *this;
    volScalarField indic
      (
       IOobject
       (
    "indic",
    mesh().time().timeName(),
    mesh(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
       mesh(),
       dimensionedScalar("zero", dimless, 0),
       zeroGradientFvPatchScalarField::typeName
       );
    forAll (laws, lawI)
      {
    // interface fields will not be correct but we fix them after
    indic.internalField() = indicator(lawI)();
    surfaceScalarField indicatorSurf = fvc::interpolate(indic);
    // fix boundary fields
    forAll(mesh().boundary(), patchi)
      {
        indicatorSurf.boundaryField()[patchi] =
          indic.boundaryField()[patchi].patchInternalField();
      }
        result += indicatorSurf*laws[lawI].GIc()();
      }

    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();

    forAll(result.internalField(), faceI)
    {
      label ownerMat = label(materials_[owner[faceI]]);
      label neighbourMat = label(materials_[neighbour[faceI]]);

        if (ownerMat != neighbourMat)
        {
        result.internalField()[faceI] = interfaceGIc(ownerMat, neighbourMat);
    }
    }

    forAll(mesh().boundary(), patchI)
    {

        if (mesh().boundary()[patchI].type() == cohesiveFvPatch::typeName)
        {
      //label size = (mesh().boundary()[patchI].size())/2;
      // const labelList& fCells = mesh().boundary()[patchI].faceCells();
      scalarField localPatchMaterials =
          materials_.boundaryField()[patchI].patchInternalField();
      const crackerFvMesh& crackerMesh = refCast<const crackerFvMesh>(mesh());
      scalarField globalPatchMaterials =
          crackerMesh.globalCrackField(localPatchMaterials);
      const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
      label globalIndex = crackerMesh.localCrackStart();

      forAll(mesh().boundaryMesh()[patchI], facei)
        {
          label ownerMat = label(globalPatchMaterials[globalIndex]);
          label neighbourMat = label(globalPatchMaterials[gcfa[globalIndex]]);

                if (ownerMat != neighbourMat)
                {
                    result.boundaryField()[patchI][facei] =
                        interfaceGIc(ownerMat, neighbourMat);
                }
        globalIndex++;
            }
        }
    else if (mesh().boundary()[patchI].coupled())
      {
        const scalarField ownerMatField =
            materials_.boundaryField()[patchI].internalField();
        const scalarField neighbourMatField =
            materials_.boundaryField()[patchI].patchNeighbourField();

        forAll(mesh().boundaryMesh()[patchI], facei)
        {
            label ownerMat = label(ownerMatField[facei]);
            label neighbourMat = label(neighbourMatField[facei]);

            if (ownerMat != neighbourMat)
            {
                //result.boundaryField()[patchI][facei] = iterGIc();
                result.boundaryField()[patchI][facei] =
                    interfaceGIc(ownerMat, neighbourMat);
            }
        }
      }
    }

    return tresult;
}

Foam::tmp<Foam::surfaceScalarField> Foam::multiMaterialCohesiveLaw::GIIc() const
{
    tmp<surfaceScalarField> tresult
    (
        new surfaceScalarField
        (
            IOobject
            (
                "GIIc",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
        dimensionedScalar("zero", dimForce*dimLength/dimArea, 0)
    )
    );
    surfaceScalarField& result = tresult();

    const PtrList<cohesiveLaw>& laws = *this;
    volScalarField indic
      (
       IOobject
       (
    "indic",
    mesh().time().timeName(),
    mesh(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
       mesh(),
       dimensionedScalar("zero", dimless, 0),
       zeroGradientFvPatchScalarField::typeName
       );
    forAll (laws, lawI)
      {
    // interface fields will not be correct but we fix them after
    indic.internalField() = indicator(lawI)();
    surfaceScalarField indicatorSurf = fvc::interpolate(indic);
    // fix boundary fields
    forAll(mesh().boundary(), patchi)
      {
        indicatorSurf.boundaryField()[patchi] =
          indic.boundaryField()[patchi].patchInternalField();
      }
        result += indicatorSurf*laws[lawI].GIIc()();
      }

    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();

    forAll(result.internalField(), faceI)
    {
      label ownerMat = label(materials_[owner[faceI]]);
      label neighbourMat = label(materials_[neighbour[faceI]]);

        if (ownerMat != neighbourMat)
        {
      result.internalField()[faceI] = interfaceGIIc(ownerMat, neighbourMat);
    }
    }

    forAll(mesh().boundary(), patchI)
    {

        if (mesh().boundary()[patchI].type() == cohesiveFvPatch::typeName)
        {
      //label size = (mesh().boundary()[patchI].size())/2;
      // const labelList& fCells = mesh().boundary()[patchI].faceCells();
      scalarField localPatchMaterials =
          materials_.boundaryField()[patchI].patchInternalField();
      const crackerFvMesh& crackerMesh = refCast<const crackerFvMesh>(mesh());
      scalarField globalPatchMaterials =
          crackerMesh.globalCrackField(localPatchMaterials);
      const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
      label globalIndex = crackerMesh.localCrackStart();

      forAll(mesh().boundaryMesh()[patchI], facei)
        {
          label ownerMat = label(globalPatchMaterials[globalIndex]);
          label neighbourMat = label(globalPatchMaterials[gcfa[globalIndex]]);

          if (ownerMat != neighbourMat)
          {
              result.boundaryField()[patchI][facei] =
                  interfaceGIIc(ownerMat, neighbourMat);
          }
          globalIndex++;
        }
        }

    else if (mesh().boundary()[patchI].coupled())
      {
        const scalarField ownerMatField =
            materials_.boundaryField()[patchI].internalField();
        const scalarField neighbourMatField =
            materials_.boundaryField()[patchI].patchNeighbourField();

        forAll(mesh().boundaryMesh()[patchI], facei)
          {
                label ownerMat = label(ownerMatField[facei]);
                label neighbourMat = label(neighbourMatField[facei]);

                if (ownerMat != neighbourMat)
          {
              result.boundaryField()[patchI][facei] =
                  interfaceGIIc(ownerMat, neighbourMat);
          }
          }
      }

    }

    return tresult;
}

Foam::scalar Foam::multiMaterialCohesiveLaw::interfaceSigmaMax
(
    double mat1,
    double mat2
) const
{
  label interfaceLawID = interfaceID(mat1, mat2);

  return interfaceCohesiveLaws_[interfaceLawID].sigmaMaxValue();
}

Foam::scalar Foam::multiMaterialCohesiveLaw::
interfaceTauMax
(
    double mat1,
    double mat2
) const
{
  label interfaceLawID = interfaceID(mat1, mat2);

  return interfaceCohesiveLaws_[interfaceLawID].tauMaxValue();
}

Foam::scalar Foam::multiMaterialCohesiveLaw::interfaceGIc
(
    double mat1,
    double mat2
) const
{
  label interfaceLawID = interfaceID(mat1, mat2);

  return interfaceCohesiveLaws_[interfaceLawID].GIcValue();
}

Foam::scalar Foam::multiMaterialCohesiveLaw::interfaceGIIc
(
    double mat1,
    double mat2
) const
{
  label interfaceLawID = interfaceID(mat1, mat2);

  return interfaceCohesiveLaws_[interfaceLawID].GIIcValue();
}


void Foam::multiMaterialCohesiveLaw::damageTractions
(
 scalar& tN,
 scalar& tS,
 const scalar deltaN,
 const scalar deltaS,
 const scalar GI,
 const scalar GII,
 const label faceID,
 const scalarField& globalPatchMaterials
 ) const
{
  // Find out which cohesive law does the face belong to
  const crackerFvMesh& crackerMesh = refCast<const crackerFvMesh>(mesh());
  const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
  label ownerMat = label(globalPatchMaterials[faceID]);
  label neighbourMat = label(globalPatchMaterials[gcfa[faceID]]);

  if (ownerMat != neighbourMat)
    {
      label matID = interfaceID(ownerMat, neighbourMat);

      // face is on multi-material interface
      interfaceCohesiveLaws_[matID].damageTractions
          (tN, tS, deltaN, deltaS, GI, GII, faceID, globalPatchMaterials);
    }
  else
    {
      // face is within one material
      // call material law function
      label matID = ownerMat;
      (*this)[matID].damageTractions
          (tN, tS, deltaN, deltaS, GI, GII, faceID, globalPatchMaterials);
    }
}


Foam::tmp<Foam::surfaceVectorField>
Foam::multiMaterialCohesiveLaw::interfaceTraction
(
 surfaceVectorField n,
 volVectorField U,
 volTensorField gradU,
 volScalarField mu,
 volScalarField lambda
 ) const
{
  tmp<surfaceVectorField> tresult
    (
        new surfaceVectorField
        (
            IOobject
            (
            "traction",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
            mesh(),
        dimensionedVector("zero", dimForce/dimArea, vector::zero)
    )
    );
    surfaceVectorField& result = tresult();

    surfaceScalarField mu1 = fvc::interpolate(mu);
    surfaceScalarField lambda1 = fvc::interpolate(lambda);

    surfaceTensorField sGradU =
        ((I - n*n)&fvc::interpolate(gradU));

    vectorField UI = U.internalField();

    // result =
    //   (2*mu1 + lambda1)*fvc::snGrad(U)
    //   - (mu1 + lambda1)*(n&sGradU)
    //   + mu1*(sGradU&n)
    //   + lambda1*tr(sGradU)*n;
    result =
        2*mu1
        *(n&fvc::interpolate(symm(gradU)))
        + lambda1*n*fvc::interpolate(tr(gradU));

    //    const labelList& owner = mesh().owner();
    //const labelList& neighbour = mesh().neighbour();

    /*        forAll(result.internalField(), faceI)
    {
            label ownerMat = materials_[owner[faceI]];
            label neighbourMat = materials_[neighbour[faceI]];

            if (ownerMat != neighbourMat)
            {
        vector ownN = n.internalField()[faceI];
        vector ngbN = -n.internalField()[faceI];

        tensor ownSGradU =
        ((I-ownN*ownN)&gradU.internalField()[owner[faceI]]);
        tensor ngbSGradU =
        ((I-ngbN*ngbN)&gradU.internalField()[neighbour[faceI]]);

        scalar ownTrSGradUt = tr(ownSGradU&(I-ownN*ownN));
        scalar ngbTrSGradUt = tr(ngbSGradU&(I-ngbN*ngbN));

        vector ownSGradUn = (ownSGradU&ownN);
        vector ngbSGradUn = (ngbSGradU&ngbN);

        scalar ownMu = mu.internalField()[owner[faceI]];
        scalar ngbMu = mu.internalField()[neighbour[faceI]];

        scalar ownLambda = lambda.internalField()[owner[faceI]];
        scalar ngbLambda = lambda.internalField()[neighbour[faceI]];

        vector ownUt = ((I-ownN*ownN)&UI[owner[faceI]]);
        vector ngbUt = ((I-ngbN*ngbN)&UI[neighbour[faceI]]);

        vector ownUn = ownN*(ownN&U.internalField()[owner[faceI]]);
        vector ngbUn = ngbN*(ngbN&U.internalField()[neighbour[faceI]]);

        scalar ownDn = mesh().weights().internalField()[faceI]
                 * (1.0/mesh().deltaCoeffs().internalField()[faceI]);
        scalar ngbDn =
        (1.0/mesh().deltaCoeffs().internalField()[faceI]) - ownDn;

        scalar own2ML = (2*ownMu + ownLambda);
        scalar ngb2ML = (2*ngbMu + ngbLambda);

        scalar face2ML = 1.0/((1.0-mesh().weights().internalField()[faceI])
                   / own2ML + mesh().weights().internalField()[faceI]/ngb2ML);
        scalar faceMu = 1.0/((1.0-mesh().weights().internalField()[faceI])
                  / ownMu + mesh().weights().internalField()[faceI]/ngbMu);

        vector tractionN =
                    face2ML*(ownUn - ngbUn)/(ownDn+ngbDn)
                  + ((own2ML*ngbDn*ngbLambda*ngbN*ngbTrSGradUt)
                  + (ngb2ML*ownDn*ownLambda*ngbN*ownTrSGradUt))
                  / ((own2ML*ngbDn) + (ngb2ML*ownDn));

        vector tractionT =
                faceMu*(ownUt - ngbUt)/(ownDn+ngbDn)
                  + ((ownMu*ngbMu*ngbDn*ngbSGradUn)
                  + (ownMu*ngbMu*ownDn*ownSGradUn))
                  / ((ownMu*ngbDn) + (ngbMu*ownDn));

        //      Info << tractionN << " *** " << tractionT <<endl;

        result[faceI] = tractionN + tractionT;
        //vector traction = tractionN + tractionT;

        //vector ownDelta =
        mesh().C()[neighbour[faceI]] - mesh().C()[owner[faceI]];
        }
        }
    */
    // philipc
    // cracks and bi-material interfaces may be on processor boundaries
    // so we must correct tractions there too
    forAll(mesh().boundary(), patchI)
      {
        /*
 //if (mesh().boundaryMesh()[patchI].type() == processorPolyPatch::typeName)
        if (mesh().boundary()[patchI].coupled())
          {
        const scalarField ownerMatField =
        materials_.boundaryField()[patchI].patchInternalField();
        const scalarField neighbourMatField =
        materials_.boundaryField()[patchI].patchNeighbourField();
        const labelList& faceCells = mesh().boundary()[patchI].faceCells();
        const tensorField ngbGradUField =
        gradU.boundaryField()[patchI].patchNeighbourField();
        const scalarField ngbMuField =
        mu.boundaryField()[patchI].patchNeighbourField();
        const scalarField ngbLambdaField =
        lambda.boundaryField()[patchI].patchNeighbourField();
        const vectorField ngbUField =
        U.boundaryField()[patchI].patchNeighbourField();
        const scalarField weights = mesh().weights().boundaryField()[patchI];

        forAll(mesh().boundaryMesh()[patchI], faceI)
          {
            label faceCelli = faceCells[faceI];
            label ownerMat = ownerMatField[faceI];
            label neighbourMat = neighbourMatField[faceI];

            if (ownerMat != neighbourMat)
              {
            const vector& ownN = n.boundaryField()[patchI][faceI];
            const vector& ngbN = -n.boundaryField()[patchI][faceI];
            const tensor& ngbGradU = ngbGradUField[faceI];
            const vector& ngbU = ngbUField[faceI];

            //tensor ownSGradU =
//            ((I-ownN*ownN)&gradU.boundaryField()[patchI][owner[faceI]]);
            tensor ownSGradU = ((I-ownN*ownN)&gradU.internalField()[faceCelli]);
            //tensor ngbSGradU =
            // ((I-ngbN*ngbN)&gradU.boundaryField()[patchI][neighbour[faceI]]);
            tensor ngbSGradU = ((I-ngbN*ngbN)&ngbGradU);

            scalar ownTrSGradUt = tr(ownSGradU&(I-ownN*ownN));
            scalar ngbTrSGradUt = tr(ngbSGradU&(I-ngbN*ngbN));

            vector ownSGradUn = (ownSGradU&ownN);
            vector ngbSGradUn = (ngbSGradU&ngbN);

            scalar ownMu = mu.internalField()[faceCelli];
            //scalar ngbMu = mu.boundaryField()[patchI][neighbour[faceI]];
            scalar ngbMu = ngbMuField[faceI];

            scalar ownLambda = lambda.internalField()[faceCelli];
            //scalar ngbLambda =
            //lambda.boundaryField()[patchI][neighbour[faceI]];
            scalar ngbLambda = ngbLambdaField[faceI];

            vector ownUt = ((I-ownN*ownN)&UI[owner[faceI]]);
            //vector ngbUt = ((I-ngbN*ngbN)&UI[neighbour[faceI]]);
            vector ngbUt = ((I-ngbN*ngbN)&ngbU);

            vector ownUn = ownN*(ownN&U.internalField()[faceCelli]);
            //vector ngbUn =
            //ngbN*(ngbN&U.boundaryField()[patchI][neighbour[faceI]]);
            vector ngbUn = ngbN*(ngbN&ngbU);

            //scalar ownDn = mesh().weights().boundaryField()[patchI][faceI]
            scalar ownDn = weights[faceI]
              * (1.0/mesh().deltaCoeffs().boundaryField()[patchI][faceI]);
            scalar ngbDn =
            (1.0/mesh().deltaCoeffs().boundaryField()[patchI][faceI]) - ownDn;

            scalar own2ML = (2*ownMu + ownLambda);
            scalar ngb2ML = (2*ngbMu + ngbLambda);

            // scalar face2ML =
            // 1.0/((1.0-mesh().weights().boundaryField()[patchI][faceI])
            // / own2ML
            //+ mesh().weights().boundaryField()[patchI][faceI]/ngb2ML);
            scalar face2ML = 1.0/((1.0-weights[faceI])
                          / own2ML + weights[faceI]/ngb2ML);
            // scalar faceMu =
            // 1.0/((1.0-mesh().weights().boundaryField()[patchI][faceI])
            // / ownMu + mesh().weights().boundaryField()[patchI][faceI]/ngbMu);
            scalar faceMu = 1.0/((1.0-weights[faceI])
                         / ownMu + weights[faceI]/ngbMu);

            vector tractionN =
              face2ML*(ownUn - ngbUn)/(ownDn+ngbDn)
              + ((own2ML*ngbDn*ngbLambda*ngbN*ngbTrSGradUt)
              + (ngb2ML*ownDn*ownLambda*ngbN*ownTrSGradUt))
              / ((own2ML*ngbDn) + (ngb2ML*ownDn));

            vector tractionT =
              faceMu*(ownUt - ngbUt)/(ownDn+ngbDn)
              + ((ownMu*ngbMu*ngbDn*ngbSGradUn)
              + (ownMu*ngbMu*ownDn*ownSGradUn))
              / ((ownMu*ngbDn) + (ngbMu*ownDn));

            result.boundaryField()[patchI][faceI] = tractionN + tractionT;
              }
          }
          }

          else*/
          if (mesh().boundary()[patchI].type() == cohesiveFvPatch::typeName)
          {
              // I think mu and lambda are wrong on new crack faces
              // so use internal cell values
              const scalarField muPatch =
                  mu.boundaryField()[patchI].patchInternalField();
              const scalarField lambdaPatch =
                  lambda.boundaryField()[patchI].patchInternalField();
              result.boundaryField()[patchI] =
                  (2*muPatch*n.boundaryField()[patchI]
                   &symm(gradU.boundaryField()[patchI]))
                  + (lambdaPatch*tr(gradU.boundaryField()[patchI])
                     *n.boundaryField()[patchI]);
          }

      }

     return tresult;
}

void Foam::multiMaterialCohesiveLaw::correct()
{
    PtrList<cohesiveLaw>& laws = *this;

    forAll (laws, lawI)
    {
        laws[lawI].correct();
    }
}


// ************************************************************************* //
