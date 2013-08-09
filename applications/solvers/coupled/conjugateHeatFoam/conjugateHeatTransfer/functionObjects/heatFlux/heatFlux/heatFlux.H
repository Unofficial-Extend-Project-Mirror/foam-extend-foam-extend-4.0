/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::heatFlux

Description

SourceFiles
    heatFlux.C
    heatFluxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef heatFlux_H
#define heatFlux_H

#include "volFieldsFwd.H"
#include "pointFieldFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class heatFluxItem;
class OFstream;
template<class Type>
class List;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class heatFlux Declaration
\*---------------------------------------------------------------------------*/

class heatFlux
{
protected:

    // Private data

        //- Name of this set of heat flux check
        word name_;

        //- name of fluid thermal diffusivity
        word Kfluid_;

        //- Database this class is registered to
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

    // Private Member Functions

        // Calculation functions

            //- Write
            virtual void calcAndPrint();

        // Functions to be over-ridden from IOoutputFilter class

            //- Update mesh
            virtual void updateMesh(const mapPolyMesh&);

            //- Move points
            virtual void movePoints(const Field<point>&);


        //- Disallow default bitwise copy construct
        heatFlux(const heatFlux&);

        //- Disallow default bitwise assignment
        void operator=(const heatFlux&);


public:

    //- Runtime type information
    TypeName("heatFlux");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        heatFlux
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~heatFlux();


    // Member Functions

        //- Return name of the set of field averages
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field average data
        virtual void read(const dictionary&);

        //- Execute the averaging
        virtual void execute();

        //- Execute the averaging at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the field average data and write
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
