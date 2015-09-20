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

Description
    Selects a cell set through a dictionary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "volFields.H"

using namespace Foam;

template<class GeoField>
void setFieldType
(
    const fvMesh& mesh,
    const labelList& selectedCells,
    Istream& fieldValueStream
)
{
    // Read field and value together; otherwise there will be an input error
    // when a field is not found.  HJ, 3/Aug/2011
    word fieldName(fieldValueStream);

    typename GeoField::value_type value
    (
        static_cast<const typename GeoField::value_type&>
        (
            pTraits<typename GeoField::value_type>(fieldValueStream)
        )
    );

    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check field exists
    if (fieldHeader.headerOk())
    {
        Info<< "    Setting " << fieldHeader.headerClassName()
            << " " << fieldName << endl;

        GeoField field(fieldHeader, mesh);

        if (selectedCells.size() == field.size())
        {
            field.internalField() = value;
        }
        else
        {
            forAll (selectedCells, celli)
            {
                field[selectedCells[celli]] = value;
            }
        }

        forAll (field.boundaryField(), patchi)
        {
            // Forced patch assignment.  HJ, 1/Aug/2010
            field.boundaryField()[patchi] ==
                field.boundaryField()[patchi].patchInternalField();
        }

        field.write();
    }
    else
    {
        WarningIn
        (
            "void setFieldType"
            "(const fvMesh& mesh, const labelList& selectedCells,"
            "Istream& fieldValueStream)"
        ) << "Field " << fieldName << " not found" << endl;
    }
}


class setField
{

public:

    setField()
    {}

    autoPtr<setField> clone() const
    {
        return autoPtr<setField>(new setField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedCells_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedCells)
        :
            mesh_(mesh),
            selectedCells_(selectedCells)
        {}

        autoPtr<setField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if (fieldType == "volScalarFieldValue")
            {
                setFieldType<volScalarField>
                    (mesh_, selectedCells_, fieldValues);
            }
            else if (fieldType == "volVectorFieldValue")
            {
                setFieldType<volVectorField>
                    (mesh_, selectedCells_, fieldValues);
            }
            else if (fieldType == "volSphericalTensorFieldValue")
            {
                setFieldType<volSphericalTensorField>
                    (mesh_, selectedCells_, fieldValues);
            }
            else if (fieldType == "volSymmTensorFieldValue")
            {
                setFieldType<volSymmTensorField>
                    (mesh_, selectedCells_, fieldValues);
            }
            else if (fieldType == "volTensorFieldValue")
            {
                setFieldType<volTensorField>
                    (mesh_, selectedCells_, fieldValues);
            }
            else
            {
                WarningIn("setField::iNew::operator()(Istream& is)")
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setField>(new setField());
        }
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Time = " << runTime.timeName() << endl;

#   include "createMesh.H"

    Info<< "Reading setFieldsDict\n" << endl;

    IOdictionary setFieldsDict
    (
        IOobject
        (
            "setFieldsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (setFieldsDict.found("defaultFieldValues"))
    {
        Info<< "Setting field default values" << endl;
        PtrList<setField> defaultFieldValues
        (
            setFieldsDict.lookup("defaultFieldValues"),
            setField::iNew(mesh, labelList(mesh.nCells()))
        );
        Info<< endl;
    }


    Info<< "Setting field region values" << endl;

    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    forAll (regions, regionI)
    {
        const entry& region = regions[regionI];

        autoPtr<topoSetSource> cellSelector =
            topoSetSource::New(region.keyword(), mesh, region.dict());

        cellSet selectedCellSet
        (
            mesh,
            "cellSet",
            mesh.nCells()/10+1  // Reasonable size estimate.
        );

        cellSelector->applyToSet
        (
            topoSetSource::NEW,
            selectedCellSet
        );

        PtrList<setField> fieldValues
        (
            region.dict().lookup("fieldValues"),
            setField::iNew(mesh, selectedCellSet.toc())
        );
    }

    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //
