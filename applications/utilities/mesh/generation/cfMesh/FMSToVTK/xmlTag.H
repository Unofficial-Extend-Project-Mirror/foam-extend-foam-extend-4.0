/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::xmlTag

Description
    Simple XML tag class allowing child tags and attributes. Specialized
    output stream operators are provided to display or write the XML
    structure.

Author
    Ivor Clifford <ivor.clifford@psi.ch>

\*---------------------------------------------------------------------------*/

#ifndef XMLtag_H
#define XMLtag_H

#include "OStringStream.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class xmlTag Declaration
\*---------------------------------------------------------------------------*/

class xmlTag
:
    public OStringStream
{
    // Private data
    
        //- Tag name
        word name_;
        
        //- Attributes
        HashTable<string> attributes_;
        
        //- Child tags
        DynamicList<xmlTag> children_;
        
public:

    // Constructors

        //- Null constructor
        xmlTag()
        :
            OStringStream(),
            name_("unknown"),
            attributes_(),
            children_()
        {}
        
        //- Construct given tag name
        xmlTag(const word& name)
        :
            OStringStream(),
            name_(name),
            attributes_(),
            children_()
        {}

        //- Construct as copy
        xmlTag(const xmlTag& tag)
        :
            OStringStream(tag),
            name_(tag.name_),
            attributes_(tag.attributes_),
            children_(tag.children_)
        {}

    //- Destructor
    ~xmlTag()
    {}


    // Member Functions
    
        //- Add an attribute
        template<class T>
        void addAttribute(const keyType& key, const T& value)
        {
            OStringStream os;
            os << value;
            attributes_.insert(key, os.str());
        };
        
        //- Add a fileName attribute
        void addAttribute(const keyType& key, const fileName& value)
        {
            attributes_.insert(key, value);
        };
        
        //- Add a string attribute
        void addAttribute(const keyType& key, const string& value)
        {
            attributes_.insert(key, value);
        };

        //- Add a word attribute
        void addAttribute(const keyType& key, const word& value)
        {
            attributes_.insert(key, value);
        };
        
        //- Add a child
        xmlTag& addChild(const xmlTag& tag)
        {
            children_.append(tag);
            
            return children_[children_.size()-1];
        };
        
        //- Create and add a new child
        xmlTag& addChild(const word& name)
        {
            return addChild(xmlTag(name));
        }

    // Member Operators

        void operator=(const xmlTag& tag)
        {
            name_ = tag.name_;
            attributes_ = tag.attributes_;
            children_ = tag.children_;            
            OStringStream::rewind();
            Foam::operator<<(*this, tag.str().c_str());
        };
        
    // Friend IOstream Operators

        friend Ostream& operator<<(Ostream&, const xmlTag&);
        
        template<class Form, class Cmpt, int nCmpt>
        friend xmlTag& operator<<(xmlTag&, const VectorSpace<Form, Cmpt, nCmpt>&);
        
        friend xmlTag& operator<<(xmlTag&, const labelledTri&);
        
        template<class T, unsigned Size>
        friend xmlTag& operator<<(xmlTag&, const FixedList<T, Size>&);

        template<class T>
        friend xmlTag& operator<<(xmlTag&, const LongList<T>&);
        
        template<class T>
        friend xmlTag& operator<<(xmlTag&, const UList<T>&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write the tag in XML format to the supplied output stream
Ostream& operator<<(Ostream& os, const xmlTag& tag)
{
    // Tag name
    os << indent << '<' << tag.name_;
    
    // Attributes and text
    for
    (
        HashTable<string>::const_iterator iter = tag.attributes_.cbegin();
        iter != tag.attributes_.cend();
        ++iter
    )
    {
        os << token::SPACE << iter.key() << '=' << iter();
    }
    
    if (tag.str().size() || tag.children_.size())
    {
        os << '>' << nl;
                
        // Children
        os.incrIndent();
        
        forAll(tag.children_, i)
        {
            os << tag.children_[i];
        }
        
        os.decrIndent();
        
        // Tag text
        os << tag.str().c_str();
        
        // Close tag
        os << indent << "</" << tag.name_ << '>' << endl;
    }
    else
    {
        // Empty element tag
        os << "/>" << endl;
    }
    return os;
}


//- Append the supplied data to the tag text
template<class T>
xmlTag& operator<<(xmlTag& tag, const UList<T>& data)
{
    forAll(data, i)
    {
        tag << data[i] << token::SPACE;
    }
    
    tag << nl;
    
    return tag;
}

//- Append the supplied data to the tag text
template<class T>
xmlTag& operator<<(xmlTag& tag, const LongList<T>& data)
{
    forAll(data, i)
    {
        tag << data[i] << token::SPACE;
    }
    
    tag << nl;
    
    return tag;
}


//- Append the supplied data to the tag text
template<class Form, class Cmpt, int nCmpt>
xmlTag& operator<<(xmlTag& tag, const VectorSpace<Form, Cmpt, nCmpt>& data)
{
    forAll(data, i)
    {
        tag << data[i] << token::SPACE;
    }
    
    tag << nl;
    
    return tag;
}

//- Append the supplied data to the tag text
template<class T, unsigned Size>
xmlTag& operator<<(xmlTag& tag, const FixedList<T, Size>& data)
{
    forAll(data, i)
    {
        tag << data[i] << token::SPACE;
    }
    
    tag << nl;
    
    return tag;
}

//- Append the supplied data to the tag text
xmlTag& operator<<(xmlTag& tag, const labelledTri& data)
{
    const triFace& tFace = data;
    
    return tag << tFace;
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
