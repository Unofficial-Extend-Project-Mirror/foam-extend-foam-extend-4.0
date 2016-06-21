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

Class
    Foam::fvSchemes

Description
    Selector class for finite volume differencing schemes.
    fvMesh is derived from fvSchemes so that all fields have access to the
    fvSchemes from the mesh reference they hold.

SourceFiles
    fvSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef fvSchemes_H
#define fvSchemes_H

#include "IOdictionary.H"
#include "debugSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvSchemes Declaration
\*---------------------------------------------------------------------------*/

class fvSchemes
:
    public IOdictionary
{
    // Private data

        dictionary ddtSchemes_;
        ITstream defaultDdtScheme_;

        dictionary d2dt2Schemes_;
        ITstream defaultD2dt2Scheme_;

        dictionary interpolationSchemes_;
        ITstream defaultInterpolationScheme_;

        dictionary divSchemes_;
        ITstream defaultDivScheme_;

        dictionary gradSchemes_;
        ITstream defaultGradScheme_;

        dictionary snGradSchemes_;
        ITstream defaultSnGradScheme_;

        dictionary laplacianSchemes_;
        ITstream defaultLaplacianScheme_;

        mutable dictionary fluxRequired_;
        bool defaultFluxRequired_;


    // Private Member Functions

        //- Clear the dictionaries and streams before reading
        void clear();

        //- Disallow default bitwise copy construct
        fvSchemes(const fvSchemes&);

        //- Disallow default bitwise assignment
        void operator=(const fvSchemes&);


public:

    //- Debug switch
    static debug::debugSwitch debug;


    // Constructors

        //- Construct for objectRegistry
        fvSchemes(const objectRegistry& obr);


    // Member Functions

        // Access

            const dictionary& schemesDict() const;

            ITstream& ddtScheme(const word& name) const;

            ITstream& d2dt2Scheme(const word& name) const;

            ITstream& interpolationScheme(const word& name) const;

            ITstream& divScheme(const word& name) const;

            ITstream& gradScheme(const word& name) const;

            ITstream& snGradScheme(const word& name) const;

            ITstream& laplacianScheme(const word& name) const;

            void setFluxRequired(const word& name) const;

            bool fluxRequired(const word& name) const;


        // Edit

            //- Return access to ddt schemes
            dictionary& ddtSchemes()
            {
                return ddtSchemes_;
            }

            //- Return access to d2dt2 schemes
            dictionary& d2dt2Schemes()
            {
                return d2dt2Schemes_;
            }

            //- Return access to interpolation schemes
            dictionary& interpolationSchemes()
            {
                return interpolationSchemes_;
            }

            //- Return access to div schemes
            dictionary& divSchemes()
            {
                return divSchemes_;
            }

            //- Return access to grad schemes
            dictionary& gradSchemes()
            {
                return gradSchemes_;
            }

            //- Return access to snGrad schemes
            dictionary& snGradSchemes()
            {
                return snGradSchemes_;
            }

            //- Return access to laplacian schemes
            dictionary& laplacianSchemes()
            {
                return laplacianSchemes_;
            }

            //- Return access to flux required
            dictionary& fluxRequired()
            {
                return fluxRequired_;
            }


        // Read

            //- Read the fvSchemes
            bool read();


        // Write

            //- WriteData function required for regIOobject write operation
            virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
