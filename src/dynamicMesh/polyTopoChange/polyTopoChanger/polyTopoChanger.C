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

#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyTopoChanger, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::polyTopoChanger::readModifiers()
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        PtrList<polyMeshModifier>& modifiers = *this;

        // Read modifiers
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        modifiers.setSize(patchEntries.size());

        forAll(modifiers, modifierI)
        {
            modifiers.set
            (
                modifierI,
                polyMeshModifier::New
                (
                    patchEntries[modifierI].keyword(),
                    patchEntries[modifierI].dict(),
                    modifierI,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "polyTopoChanger::polyTopoChanger"
            "(const IOobject&, const polyMesh&)"
        );

        close();
    }
}


// Read constructor given IOobject and a polyMesh reference
Foam::polyTopoChanger::polyTopoChanger
(
    const IOobject& io,
    polyMesh& mesh
)
:
    PtrList<polyMeshModifier>(),
    regIOobject(io),
    mesh_(mesh)
{
    readModifiers();
}


// Read constructor given IOobject and a polyMesh reference
Foam::polyTopoChanger::polyTopoChanger(polyMesh& mesh)
:
    PtrList<polyMeshModifier>(),
    regIOobject
    (
        IOobject
        (
            "meshModifiers",
            mesh.time().findInstance
            (
                mesh.meshDir(),
                "meshModifiers",
                IOobject::READ_IF_PRESENT
            ),
            mesh.meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    morphIndex_(mesh.time().timeIndex()),
    mesh_(mesh)
{
    readModifiers();
}


// Return a list of modifier types
Foam::wordList Foam::polyTopoChanger::types() const
{
    const PtrList<polyMeshModifier>& modifiers = *this;

    wordList t(modifiers.size());

    forAll (modifiers, modifierI)
    {
        t[modifierI] = modifiers[modifierI].type();
    }

    return t;
}


// Return a list of modifier names
Foam::wordList Foam::polyTopoChanger::names() const
{
    const PtrList<polyMeshModifier>& modifiers = *this;

    wordList t(modifiers.size());

    forAll (modifiers, modifierI)
    {
        t[modifierI] = modifiers[modifierI].name();
    }

    return t;
}


// Is topology change required
bool Foam::polyTopoChanger::changeTopology() const
{
    // Go through all mesh modifiers and accumulate the morphing information
    const PtrList<polyMeshModifier>& topoChanges = *this;

    bool triggerChange = false;

    forAll (topoChanges, morphI)
    {
        if (topoChanges[morphI].active())
        {
            bool curTriggerChange = topoChanges[morphI].changeTopology();

            if (debug)
            {
                Pout<< "Modifier " << morphI << " named "
                    << topoChanges[morphI].name();

                if (curTriggerChange)
                {
                    Pout << " morphing" << endl;
                }
                else
                {
                    Pout << " unchanged" << endl;
                }
            }

            triggerChange = triggerChange || curTriggerChange;
        }
        else
        {
            if (debug)
            {
                Pout<< "Modifier " << morphI  << " named "
                    << topoChanges[morphI].name() << " inactive" << endl;
            }
        }
    }

    return triggerChange;
}


// Return topology change request
Foam::autoPtr<Foam::polyTopoChange>
Foam::polyTopoChanger::topoChangeRequest() const
{
    // Collect changes from all modifiers
    const PtrList<polyMeshModifier>& topoChanges = *this;

    polyTopoChange* refPtr(new polyTopoChange(mesh()));
    polyTopoChange& ref = *refPtr;

    forAll (topoChanges, morphI)
    {
        if (topoChanges[morphI].active())
        {
            topoChanges[morphI].setRefinement(ref);
        }
    }

    return autoPtr<polyTopoChange>(refPtr);
}


// Correct polyTopoChanger after moving points
void Foam::polyTopoChanger::modifyMotionPoints(pointField& p) const
{
    if (debug)
    {
        Pout<< "void polyTopoChanger::modifyMotionPoints("
            << "pointField& p) const : "
            << "Updating motion points by mesh modifiers" << endl;
    }

    const PtrList<polyMeshModifier>& topoChanges = *this;

    forAll (topoChanges, morphI)
    {
        if (topoChanges[morphI].active())
        {
            topoChanges[morphI].modifyMotionPoints(p);
        }
    }
}


// Force recalculation of locally stored data on topological change
void Foam::polyTopoChanger::update(const mapPolyMesh& m)
{
    if (debug)
    {
        Pout<< "void polyTopoChanger::update(const mapPolyMesh& m) : "
            << "Updating mesh modifiers" << endl;
    }

    // Go through all mesh modifiers and accumulate the morphing information
    PtrList<polyMeshModifier>& topoChanges = *this;

    forAll (topoChanges, morphI)
    {
        topoChanges[morphI].updateMesh(m);
    }

    // Force the mesh modifiers to auto-write.  This allows us to
    // preserve the current state of modifiers corresponding with
    // the mesh.
    writeOpt() = IOobject::AUTO_WRITE;
    instance() = mesh_.time().timeName();
}


// Add mesh modifiers to the morph engine
void Foam::polyTopoChanger::addTopologyModifiers
(
    const List<polyMeshModifier*>& tm
)
{
    setSize(tm.size());

    // Copy the patch pointers
    forAll (tm, tmI)
    {
        if (tm[tmI]->topoChanger() != *this)
        {
            FatalErrorIn
            (
                "void polyTopoChanger::addTopologyModifiers("
                "const List<polyMeshModifier*>& tm)"
            )   << "Mesh modifier created with different mesh reference."
                << abort(FatalError);
        }
        set(tmI, tm[tmI]);
    }

    writeOpt() = IOobject::AUTO_WRITE;
}


Foam::label Foam::polyTopoChanger::findModifierID
(
    const word& modName
) const
{
    const PtrList<polyMeshModifier>& topoChanges = *this;

    forAll (topoChanges, morphI)
    {
        if (topoChanges[morphI].name() == modName)
        {
            return morphI;
        }
    }

    // Modifier not found
    if (debug)
    {
        Pout<< "label polyTopoChanger::::findModifierID(const word& "
            << "modName) const"
            << "Modifier named " << modName << " not found.  "
            << "List of available modifier names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


Foam::label Foam::polyTopoChanger::morphIndex() const
{
    return morphIndex_;
}


// writeData member function required by regIOobject
bool Foam::polyTopoChanger::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::polyTopoChanger::operator!=(const polyTopoChanger& me) const
{
    return &me != this;
}


bool Foam::polyTopoChanger::operator==(const polyTopoChanger& me) const
{
    return &me == this;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyTopoChanger& mme)
{
    os  << mme.size() << nl << token::BEGIN_LIST;

    forAll(mme, mmeI)
    {
        mme[mmeI].writeDict(os);
    }

    os  << token::END_LIST;

    return os;
}


// ************************************************************************* //
