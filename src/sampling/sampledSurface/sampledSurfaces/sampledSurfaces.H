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
    Foam::sampledSurfaces

Description
    Set of surfaces to sample.

    The write() method is used to sample and write files.

SourceFiles
    sampledSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSurfaces_H
#define sampledSurfaces_H

#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "volFieldsFwd.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                      Class sampledSurfaces Declaration
\*---------------------------------------------------------------------------*/

class sampledSurfaces
:
    public PtrList<sampledSurface>
{
    // Private classes

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public wordList
        {
        public:

            //- Surface formatter
            autoPtr<surfaceWriter<Type> > formatter;

            //- Construct null
            fieldGroup()
            :
                wordList(0),
                formatter(NULL)
            {}

            //- Construct for a particular surface format
            fieldGroup(const word& writeFormat)
            :
                wordList(0),
                formatter(surfaceWriter<Type>::New(writeFormat))
            {}

            //- Construct for a particular surface format and a list of field names
            fieldGroup
            (
                const word& writeFormat,
                const wordList& fieldNames
            )
            :
                wordList(fieldNames),
                formatter(surfaceWriter<Type>::New(writeFormat))
            {}

            void operator=(const word& writeFormat)
            {
                formatter = surfaceWriter<Type>::New(writeFormat);
            }

            void operator=(const wordList& fieldNames)
            {
                wordList::operator=(fieldNames);
            }
        };


        //- Class used for surface merging information
        class mergeInfo
        {
        public:
            pointField points;
            faceList   faces;
            labelList  pointsMap;
        };


    // Static data members

        //- output verbosity
        static bool verbose_;


    // Private data

        //- Name of this set of surfaces,
        //  Also used as the name of the sampledSurfaces directory.
        word name_;

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        //- output path
        fileName outputPath_;

        //- pointMesh for interpolation
        autoPtr<pointMesh> pMeshPtr_;

        //- volPointInterpolation for interpolation
        autoPtr<volPointInterpolation> pInterpPtr_;


        // Read from dictonary

            //- Names of fields to sample
            wordList fieldNames_;

            //- Interpolation scheme to use
            word interpolationScheme_;

            //- Output format to use
            word writeFormat_;


        // surfaces

            //- Information for merging surfaces
            List<mergeInfo> mergeList_;


        // Calculated

            //- Categorized scalar/vector/tensor fields
            fieldGroup<scalar> scalarFields_;
            fieldGroup<vector> vectorFields_;
            fieldGroup<sphericalTensor> sphericalTensorFields_;
            fieldGroup<symmTensor> symmTensorFields_;
            fieldGroup<tensor> tensorFields_;


    // Private Member Functions

        //- Classify field types, return true if nFields > 0
        bool checkFieldTypes();

        //- Merge points on surfaces
        void mergeSurfaces();

        //- Find the fields in the list of the given type, return count
        template<class Type>
        label grep
        (
            fieldGroup<Type>& fieldList,
            const wordList& fieldTypes
        ) const;

        //- Set interpolator for the field
        template<class Type>
        autoPtr<interpolation<Type> > setInterpolator
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const surfaceWriter<Type>& formatter
        );

        //- Sample and write all the fields of the given type
        template <class Type>
        void sampleAndWrite(fieldGroup<Type>&);

        //- Disallow default bitwise copy construct and assignment
        sampledSurfaces(const sampledSurfaces&);
        void operator=(const sampledSurfaces&);


public:

    //- Runtime type information
    TypeName("surfaces");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        sampledSurfaces
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~sampledSurfaces();


    // Member Functions

        //- Return name of the set of probes
        virtual const word& name() const
        {
            return name_;
        }

        //- set verbosity level
        void verbose(const bool verbosity = true);

        //- Sample and write
        virtual void write();

        //- Read the sampledSurfaces
        virtual void read(const dictionary&);

        //- Correct for mesh changes
        void correct();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const pointField&);

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sampledSurfacesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
