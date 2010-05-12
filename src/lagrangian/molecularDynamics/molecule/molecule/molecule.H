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
    molecule

Description

SourceFiles
    moleculeI.H
    molecule.C
    moleculeIO.C

\*---------------------------------------------------------------------------*/

#ifndef molecule_H
#define molecule_H

#include "Particle.H"
#include "IOstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class moleculeCloud;

/*---------------------------------------------------------------------------*\
                           Class molecule Declaration
\*---------------------------------------------------------------------------*/

class molecule
:
    public Particle<molecule>
{
    // Private data

        //- Be careful with the ordering of data.  It has an impact on binary
        //  transfer:
        //  1) Put the largest data members 1st
        //  2) Pair up labels,
        //  3) Don't go scalar-label, scalar-label, becasue in 64bit mode,
        //     the labels will be padded by 4bytes.

        // - mass of molecule
        scalar mass_;

        // - Velocity of molecule
        vector U_;

        // - Acceleration of molecule
        vector A_;

        // - Tether position
        vector tetherPosition_;

        // - Potential energy that this molecules posseses
        scalar potentialEnergy_;

        // - r_ij f_ij, stress dyad
        tensor rf_;

        // - Is the molecule tethered?
        label tethered_;

        // - id (type) of molecule
        label id_;


public:

    friend class Cloud<molecule>;

    //- Class used to pass tracking data to the trackToFace function
    class trackData
    :
        public Particle<molecule>::trackData
    {
        //- Reference to the cloud containing this particle
        moleculeCloud& molCloud_;

        // label specifying which part of the integration algorithm is taking
        // place (i.e. leapfrog 1 or leapfrog 2. Predictor or Corrector)
        label part_;


    public:

        // Constructors

            inline trackData
            (
                moleculeCloud& molCloud,
                label part
            );

        // Member functions

            inline moleculeCloud& molCloud();

            inline label part() const;
    };

    // Constructors

        //- Construct from components
        inline molecule
        (
            const Cloud<molecule>& c,
            const vector& position,
            const label celli,
            const scalar mass,
            const vector& U,
            const vector& A,
            const vector& tetherPosition,
            const label tethered,
            const label id
        );

        //- Construct from Istream
        molecule
        (
            const Cloud<molecule>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<molecule> clone() const
        {
            return autoPtr<molecule>(new molecule(*this));
        }


    // Member Functions

        void transformProperties(const tensor& T);

        void transformProperties(const vector& separation);

        // Access

            //- Return id
            inline label id() const;

            //- Return mass
            inline scalar mass() const;

            //- Return velocity
            inline const vector& U() const;
            inline vector& U();

            //- Return acceleration
            inline const vector& A() const;
            inline vector& A();

            //- Return potential energy
            inline scalar potentialEnergy() const;
            inline scalar& potentialEnergy();

            //- Return stress contribution
            inline const tensor& rf() const;
            inline tensor& rf();

            //- Return tethered
            inline label tethered() const;

            //- Return tetherPosition
            inline const vector& tetherPosition() const;
            inline vector& tetherPosition();


        //- Tracking
        bool move(trackData&);


    // Member Operators

        //- Overridable function to handle the particle hitting a
        //  processorPatch
        void hitProcessorPatch
        (
            const processorPolyPatch&,
            molecule::trackData& td
        );

        //- Overridable function to handle the particle hitting a
        //  processorPatch without trackData
        void hitProcessorPatch
        (
            const processorPolyPatch&,
            int&
        );

        //- Overridable function to handle the particle hitting a wallPatch
        void hitWallPatch
        (
            const wallPolyPatch&,
            molecule::trackData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        //- without trackData
        void hitWallPatch
        (
            const wallPolyPatch&,
            int&
        );

        //- Overridable function to handle the particle hitting a polyPatch
        void hitPatch
        (
            const polyPatch&,
            molecule::trackData& td
        );

        //- Overridable function to handle the particle hitting a polyPatch
        //- without trackData
        void hitPatch
        (
            const polyPatch&,
            int&
        );

        static void readFields(moleculeCloud& mC);

        static void writeFields(const moleculeCloud& mC);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const molecule&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "moleculeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
