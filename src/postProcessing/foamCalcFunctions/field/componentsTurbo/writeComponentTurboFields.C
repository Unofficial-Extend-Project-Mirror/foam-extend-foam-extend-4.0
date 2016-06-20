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

template <class Type>
void Foam::calcTypes::componentsTurbo::writeComponentFields
(
    const IOobject& header,
    const fvMesh& mesh,
    bool& processed
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (header.headerClassName() == fieldType::typeName)
    {
        const char* turboComponentNames[] = {"radial", "tangential", "z"};

        Info<< "    Reading " << header.name() << endl;
        fieldType field(header, mesh);
        fieldType turboField(field);

        vectorField cellCentres(mesh.cellCentres());
        cellCentres.replace(vector::Z, scalar(0));

        const scalarField r = mag(cellCentres);

        forAll(turboField, fI)
        {
            turboField[fI].replace(vector::X,
                (
                    cellCentres[fI].component(vector::X)*field[fI].component(vector::X) +
                    cellCentres[fI].component(vector::Y)*field[fI].component(vector::Y)
                ) / r[fI]
            );

            turboField[fI].replace(vector::Y,
                (
                    -cellCentres[fI].component(vector::Y)*field[fI].component(vector::X) +
                     cellCentres[fI].component(vector::X)*field[fI].component(vector::Y)
                ) / r[fI]
            );
        }

        // Convert boundary fields
        forAll(turboField.boundaryField(), patchI)
        {
            fvPatchField<Type>& bf  = field.boundaryField()[patchI];
            fvPatchField<Type>& tbf = turboField.boundaryField()[patchI];

            vectorField faceCentres(tbf.patch().Cf());

            faceCentres.replace(vector::Z, scalar(0));

            const scalarField r = mag(faceCentres);

            forAll(tbf, tbfI)
            {
                tbf[tbfI].replace(vector::X,
                    (
                        faceCentres[tbfI].component(vector::X)*bf[tbfI].component(vector::X) +
                        faceCentres[tbfI].component(vector::Y)*bf[tbfI].component(vector::Y)
                    ) / r[tbfI]
                );

                tbf[tbfI].replace(vector::Y,
                    (
                        -faceCentres[tbfI].component(vector::Y)*bf[tbfI].component(vector::X) +
                         faceCentres[tbfI].component(vector::X)*bf[tbfI].component(vector::Y)
                    ) / r[tbfI]
                );
            }
        }

        for (direction i=0; i<Type::nComponents; i++)
        {
            Info<< "    Calculating " << header.name()
                << turboComponentNames[i] << endl;

            volScalarField componentField
            (
                IOobject
                (
                    header.name() + word(turboComponentNames[i]),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                turboField.component(i)
            );
            componentField.write();
        }

        processed = true;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
