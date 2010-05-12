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

#include "dictionary.H"
#include "primitiveEntry.H"
#include "dictionaryEntry.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

defineTypeNameAndDebug(Foam::dictionary, 0);

const Foam::dictionary Foam::dictionary::null;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionary::dictionary()
:
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary
(
    const dictionary& parentDict,
    const dictionary& dict
)
:
    IDLList<entry>(dict, *this),
    name_(dict.name()),
    parent_(parentDict)
{
    for
    (
        IDLList<entry>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        hashedEntries_.insert(iter().keyword(), &iter());
    }
}


Foam::dictionary::dictionary
(
    const dictionary& dict
)
:
    IDLList<entry>(dict, *this),
    name_(dict.name()),
    parent_(dictionary::null)
{
    for
    (
        IDLList<entry>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        hashedEntries_.insert(iter().keyword(), &iter());
    }
}


Foam::autoPtr<Foam::dictionary> Foam::dictionary::clone() const
{
    return autoPtr<dictionary>(new dictionary(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dictionary::~dictionary()
{
    // cerr<< "~dictionary() " << name() << " " << long(this) << std::endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::dictionary::startLineNumber() const
{
    if (size())
    {
        return first()->startLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::label Foam::dictionary::endLineNumber() const
{
    if (size())
    {
        return last()->endLineNumber();
    }
    else
    {
        return -1;
    }
}


bool Foam::dictionary::found(const word& keyword, bool recursive) const
{
    if (hashedEntries_.found(keyword))
    {
        return true;
    }
    else if (recursive && &parent_ != &dictionary::null)
    {
        return parent_.found(keyword, recursive);
    }
    else
    {
        return false;
    }
}


const Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive
) const
{
    HashTable<entry*>::const_iterator iter = hashedEntries_.find(keyword);

    if (iter == hashedEntries_.end())
    {
        if (recursive && &parent_ != &dictionary::null)
        {
            return parent_.lookupEntryPtr(keyword, recursive);
        }
        else
        {
            return NULL;
        }
    }

    return iter();
}


Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive
)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(keyword);

    if (iter == hashedEntries_.end())
    {
        if (recursive && &parent_ != &dictionary::null)
        {
            return const_cast<dictionary&>(parent_).lookupEntryPtr
            (
                keyword,
                recursive
            );
        }
        else
        {
            return NULL;
        }
    }

    return iter();
}


const Foam::entry& Foam::dictionary::lookupEntry
(
    const word& keyword,
    bool recursive
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive);

    if (entryPtr == NULL)
    {
        FatalIOErrorIn
        (
            "dictionary::lookupEntry(const word& keyword) const",
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return *entryPtr;
}


Foam::ITstream& Foam::dictionary::lookup
(
    const word& keyword,
    bool recursive
) const
{
    return lookupEntry(keyword, recursive).stream();
}


bool Foam::dictionary::isDict(const word& keyword) const
{
    const entry* entryPtr = lookupEntryPtr(keyword);

    if (entryPtr)
    {
        return entryPtr->isDict();
    }
    else
    {
        return false;
    }
}


const Foam::dictionary* Foam::dictionary::subDictPtr(const word& keyword) const
{
    const entry* entryPtr = lookupEntryPtr(keyword);

    if (entryPtr)
    {
        return &entryPtr->dict();
    }
    else
    {
        return NULL;
    }
}


const Foam::dictionary& Foam::dictionary::subDict(const word& keyword) const
{
    const entry* entryPtr = lookupEntryPtr(keyword);
    if (entryPtr == NULL)
    {
        FatalIOErrorIn
        (
            "dictionary::subDict(const word& keyword) const",
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }
    return entryPtr->dict();
}


Foam::dictionary& Foam::dictionary::subDict(const word& keyword)
{
    entry* entryPtr = lookupEntryPtr(keyword);
    if (entryPtr == NULL)
    {
        FatalIOErrorIn
        (
            "dictionary::subDict(const word& keyword)",
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }
    return entryPtr->dict();
}


Foam::wordList Foam::dictionary::toc() const
{
    wordList keys(size());

    label i = 0;
    for
    (
        IDLList<entry>::const_iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        keys[i++] = iter().keyword();
    }

    return keys;
}


bool Foam::dictionary::add(entry* entryPtr, bool mergeEntry)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(entryPtr->keyword());

    if (mergeEntry && iter != hashedEntries_.end())
    {
        // merge dictionary with dictionary
        if (iter()->isDict() && entryPtr->isDict())
        {
            iter()->dict().merge(entryPtr->dict());
            delete entryPtr;

            return true;
        }
        else
        {
            // replace existing dictionary with entry or vice versa
            IDLList<entry>::replace(iter(), entryPtr);
            delete iter();
            hashedEntries_.erase(iter);

            if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
            {
                entryPtr->name() = name_ + "::" + entryPtr->keyword();
                return true;
            }
            else
            {
                IOWarningIn("dictionary::add(entry*)", (*this))
                    << "problem replacing entry "<< entryPtr->keyword()
                    << " in dictionary " << name() << endl;

                IDLList<entry>::remove(entryPtr);
                delete entryPtr;
                return false;
            }
        }
    }

    if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
    {
        entryPtr->name() = name_ + "::" + entryPtr->keyword();
        IDLList<entry>::append(entryPtr);

        return true;
    }
    else
    {
        IOWarningIn("dictionary::add(entry* entryPtr)", (*this))
            << "attempt to add entry "<< entryPtr->keyword()
            << " which already exists in dictionary " << name()
            << endl;

        delete entryPtr;
        return false;
    }
}


void Foam::dictionary::add(const entry& e, bool mergeEntry)
{
    add(e.clone(*this).ptr(), mergeEntry);
}

void Foam::dictionary::add(const word& k, const word& w, bool overwrite)
{
    add(new primitiveEntry(k, token(w)), overwrite);
}

void Foam::dictionary::add(const word& k, const Foam::string& s, bool overwrite)
{
    add(new primitiveEntry(k, token(s)), overwrite);
}

void Foam::dictionary::add(const word& k, const label l, bool overwrite)
{
    add(new primitiveEntry(k, token(l)), overwrite);
}

void Foam::dictionary::add(const word& k, const scalar s, bool overwrite)
{
    add(new primitiveEntry(k, token(s)), overwrite);
}

void Foam::dictionary::add(const word& k, const dictionary& d, bool mergeEntry)
{
    add(new dictionaryEntry(k, *this, d), mergeEntry);
}


void Foam::dictionary::set(entry* entryPtr)
{
    entry* existingPtr = lookupEntryPtr(entryPtr->keyword());

    // clear dictionary so merge acts like overwrite
    if (existingPtr && existingPtr->isDict())
    {
        existingPtr->dict().clear();
    }
    add(entryPtr, true);
}


void Foam::dictionary::set(const entry& e)
{
    set(e.clone(*this).ptr());
}

void Foam::dictionary::set(const word& k, const dictionary& d)
{
    set(new dictionaryEntry(k, *this, d));
}


bool Foam::dictionary::remove(const word& Keyword)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(Keyword);

    if (iter != hashedEntries_.end())
    {
        IDLList<entry>::remove(iter());
        delete iter();
        hashedEntries_.erase(iter);

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::dictionary::changeKeyword
(
    const word& oldKeyword,
    const word& newKeyword,
    bool forceOverwrite
)
{
    // no change
    if (oldKeyword == newKeyword)
    {
        return false;
    }

    HashTable<entry*>::iterator iter = hashedEntries_.find(oldKeyword);

    // oldKeyword not found - do nothing
    if (iter == hashedEntries_.end())
    {
        return false;
    }

    HashTable<entry*>::iterator iter2 = hashedEntries_.find(newKeyword);

    // newKeyword already exists
    if (iter2 != hashedEntries_.end())
    {
        if (forceOverwrite)
        {
            IDLList<entry>::replace(iter2(), iter());
            delete iter2();
            hashedEntries_.erase(iter2);
        }
        else
        {
            WarningIn("dictionary::changeKeyword(const word&, const word&)")
                << "cannot rename keyword "<< oldKeyword
                << " to existing keyword " << newKeyword
                << " in dictionary " << name() << endl;
            return false;
        }
    }

    // change name and HashTable, but leave DL-List untouched
    iter()->keyword() = newKeyword;
    iter()->name() = name_ + "::" + newKeyword;
    hashedEntries_.erase(oldKeyword);
    hashedEntries_.insert(newKeyword, iter());

    return true;
}


bool Foam::dictionary::merge(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::merge(const dictionary&)")
            << "attempted merge to self for dictionary " << name()
            << abort(FatalError);
    }

    bool changed = false;

    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        const word& keyword = iter().keyword();

        HashTable<entry*>::iterator iter2 = hashedEntries_.find(keyword);

        if (iter2 != hashedEntries_.end())
        {
            // Recursively merge sub-dictionaries
            // TODO: merge without copying
            if (iter2()->isDict() && iter().isDict())
            {
                if (iter2()->dict().merge(iter().dict()))
                {
                    changed = true;
                }
            }
            else
            {
                add(iter().clone(*this).ptr(), true);
                changed = true;
            }
        }
        else
        {
            // not found - just add
            add(iter().clone(*this).ptr());
            changed = true;
        }
    }

    return changed;
}


void Foam::dictionary::clear()
{
    IDLList<entry>::clear();
    hashedEntries_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::ITstream& Foam::dictionary::operator[](const word& keyword) const
{
    return lookup(keyword);
}


void Foam::dictionary::operator=(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::operator=(const dictionary&)")
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalError);
    }

    name_ = dict.name();
    clear();

    // Create clones of the entries in the given dictionary
    // resetting the parentDict to this dictionary
    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        IDLList<entry>::append(iter().clone(*this).ptr());
    }

    for
    (
        IDLList<entry>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        hashedEntries_.insert(iter().keyword(), &iter());
    }
}


void Foam::dictionary::operator+=(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::operator+=(const dictionary&)")
            << "attempted addition assignment to self for dictionary " << name()
            << abort(FatalError);
    }

    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        add(iter().clone(*this).ptr());
    }
}


void Foam::dictionary::operator|=(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::operator|=(const dictionary&)")
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalError);
    }

    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        if (!found(iter().keyword()))
        {
            add(iter().clone(*this).ptr());
        }
    }
}


void Foam::dictionary::operator<<=(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::operator<<=(const dictionary&)")
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalError);
    }

    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        set(iter().clone(*this).ptr());
    }
}


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

Foam::dictionary Foam::operator+
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary sum(dict1);
    sum += dict2;
    return sum;
}


Foam::dictionary Foam::operator|
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary sum(dict1);
    sum |= dict2;
    return sum;
}


// ************************************************************************* //
