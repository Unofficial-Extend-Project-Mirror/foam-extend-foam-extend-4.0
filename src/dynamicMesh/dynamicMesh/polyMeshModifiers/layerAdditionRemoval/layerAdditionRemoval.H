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
    Foam::layerAdditionRemoval

Description
    Cell layer addition mesh modifier

SourceFiles
    layerAdditionRemoval.C
    addCellLayer.C
    removeCellLayer.C

\*---------------------------------------------------------------------------*/

#ifndef layerAdditionRemoval_H
#define layerAdditionRemoval_H

#include "polyMeshModifier.H"
#include "primitiveFacePatch.H"
#include "ZoneIDs.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class layerAdditionRemoval Declaration
\*---------------------------------------------------------------------------*/

class layerAdditionRemoval
:
    public polyMeshModifier
{
    // Private data

        //- Master face zone ID
        faceZoneID faceZoneID_;

        //- Min thickness of extrusion layer.  Triggers layer removal
        mutable scalar minLayerThickness_;

        //- Max thickness of extrusion layer.  Triggers layer addition
        mutable scalar maxLayerThickness_;

        //- Layer thickness from previous step
        //  Used to decide whether to add or remove layers
        mutable scalar oldLayerThickness_;

        //- Point pairing
        mutable labelList* pointsPairingPtr_;

        //- Face pairing
        mutable labelList* facesPairingPtr_;

        //- Layer removal trigger time index
        mutable label triggerRemoval_;

        //- Layer addition trigger  time index
        mutable label triggerAddition_;

        //- cell zone for new cells
        label cellZone_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        layerAdditionRemoval(const layerAdditionRemoval&);

        //- Disallow default bitwise assignment
        void operator=(const layerAdditionRemoval&);

        //- Check validity of construction data
        void checkDefinition();


        // Topological changes

            //- Check for valid layer
            bool validCollapse() const;

            //- Set layer pairing.  Return true if a valid layer exists
            bool setLayerPairing() const;

            //- Return points pairing in a layer (not automatic!)
            const labelList& pointsPairing() const;

            //- Return faces pairing in a layer (not automatic!)
            const labelList& facesPairing() const;

            //- Add a layer of cells
            void addCellLayer(polyTopoChange&) const;

            //- Remove a layer of cells
            void removeCellLayer(polyTopoChange&) const;

            //- Clear addressing
            void clearAddressing() const;

        // Helpers

            //- Optionally read old thickness
            static scalar readOldThickness(const dictionary&);


    // Static data members

        //- Motion detection fraction: if the motion is smaller than
        // motionDelta*minLayerThickness_, it is assumed that the mesh is
        // not layering
        static const debug::tolerancesSwitch motionDelta_;

        //- Thickness insertion fraction for the pre-motion
        static const debug::tolerancesSwitch addDelta_;


public:

    //- Runtime type information
    TypeName("layerAdditionRemoval");


    // Constructors

        //- Construct from components
        layerAdditionRemoval
        (
            const word& name,
            const label index,
            const polyTopoChanger& mme,
            const word& zoneName,
            const scalar minThickness,
            const scalar maxThickness,
            const label cellZone = -1
        );

        //- Construct from dictionary
        layerAdditionRemoval
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyTopoChanger& mme
        );


    //- Destructor
    virtual ~layerAdditionRemoval();


    // Member Functions

        //- Set layer removal
        void setRemoval();

        //- Set addition
        void setAddition();

        //- Check for topology change
        virtual bool changeTopology() const;

        //- Insert the layer addition/removal instructions
        //  into the topological change
        virtual void setRefinement(polyTopoChange&) const;

        //- Modify motion points to comply with the topological change
        virtual void modifyMotionPoints(pointField& motionPoints) const;

        //- Force recalculation of locally stored data on topological change
        virtual void updateMesh(const mapPolyMesh&);


        // Edit

            //- Return min layer thickness which triggers removal
            scalar minLayerThickness() const
            {
                return minLayerThickness_;
            }

            //- Set min layer thickness which triggers removal
           void setMinLayerThickness(const scalar t) const;

            //- Return max layer thickness which triggers removal
            scalar maxLayerThickness() const
            {
                return maxLayerThickness_;
            }

            //- Set max layer thickness which triggers removal
           void setMaxLayerThickness(const scalar t) const;


        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
