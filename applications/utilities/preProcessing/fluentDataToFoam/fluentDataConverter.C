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

#include "fluentDataConverter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluentDataConverter::fluentDataConverter
(
    const fvMesh& mesh,
    const SLList<label>& fieldID,
    const SLList<label>& zoneID,
    const SLList<label>& firstID,
    const SLList<label>& lastID,
    const SLPtrList<FieldField<Field, scalar> >& zoneData
)
:
    mesh_(mesh),
    zoneToPatchName_
    (
        IOobject
        (
            "zoneToPatchName",
            mesh.time().constant(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    fieldID_(fieldID),
    zoneID_(zoneID),
    firstID_(firstID),
    lastID_(lastID),
    zoneData_(zoneData)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluentDataConverter::convertField
(
    const word& fieldName,
    const label unitNumber,
    const dimensionedScalar& defaultValue
) const
{
    // Create field
    Info << "Creating field " << fieldName << " for unit "<< unitNumber << endl;

    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            defaultValue
        )
    );
    volScalarField& result = tresult();

    SLList<label>::const_iterator fieldIDIter = fieldID_.begin();
    SLList<label>::const_iterator zoneIDIter = zoneID_.begin();
    SLList<label>::const_iterator firstIDIter = firstID_.begin();
    SLList<label>::const_iterator lastIDIter = lastID_.begin();
    SLPtrList<FieldField<Field, scalar> >::const_iterator zoneDataIter =
        zoneData_.begin();

    for
    (
        ;
        fieldIDIter != fieldID_.end();

        ++fieldIDIter,
        ++zoneIDIter,
        ++firstIDIter,
        ++lastIDIter,
        ++zoneDataIter
    )
    {
        // Look for field index
        if (fieldIDIter() == unitNumber)
        {
            Info<< "Found field ID for zone " << zoneIDIter();

            word patchName = zoneToPatchName_[zoneIDIter()];

            // Internal Field
            if
            (
                patchName == "unknown"
             && zoneDataIter()[0].size() == mesh_.nCells()
            )
            {
                Info<< " internal cell zone.  Size = "
                    << mesh_.nCells() << endl;
                result.internalField() = zoneDataIter()[0];
            }
            else
            {
                label patchID = mesh_.boundaryMesh().findPatchID(patchName);

                if (patchID > -1)
                {
                    Info<< " and patch " << patchName
                        << " with id " << patchID << endl;

                    result.boundaryField()[patchID] == zoneDataIter()[0];
                }
                else
                {
                    Info<< " and patch not found" << endl;
                }
            }
        }
    }

    return tresult;
}


// ************************************************************************* //
