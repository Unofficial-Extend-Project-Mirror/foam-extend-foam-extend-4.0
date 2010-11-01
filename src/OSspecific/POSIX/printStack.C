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

#include "error.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "readHexLabel.H"

#include <cxxabi.h>
#ifndef darwin
#include <execinfo.h>
#endif
#include <dlfcn.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

string pOpen(const string &cmd, label line=0)
{
    const int MAX = 1000;

    FILE *cmdPipe = popen(cmd.c_str(), "r");

    if (cmdPipe)
    {
        // Read line number of lines
        for (label cnt = 0; cnt <= line; cnt++)
        {
            char buffer[MAX];

            char* s = fgets(buffer, MAX-1, cmdPipe);

            if (s == NULL)
            {
#ifdef darwin
                // workaround for the Python-Script
                for(int i=0;i<MAX;i++) {
                    if(buffer[i]=='\n') {
                        buffer[i]='\0';
                    }
                }
                return buffer;
#else
                return "";
#endif
            }

            if (cnt == line)
            {
                string str(buffer);
                return str.substr(0, str.size()-1);
            }
        }
        pclose(cmdPipe);
    }

    return "";
}


// use popen to call addr2line (using bfd.h directly would have
// meant relinking everything)

void printSourceFileAndLine
(
    Ostream& os,
    const HashTable<label, fileName>& addressMap,
    const fileName& filename,
    const word& address
)
{
    word myAddress = address;

#ifndef darwin
    if (filename.ext() == "so")
#else
    if (filename.ext() == "dylib")
#endif
    {
        // Convert offset into .so into offset into executable.

        void *addr;
        sscanf(myAddress.c_str(), "%p",&addr);

        Dl_info info;

        dladdr(addr, &info);

        unsigned long offset = reinterpret_cast<unsigned long>(info.dli_fbase);

        IStringStream addressStr(address.substr(2));
        label addressValue = readHexLabel(addressStr);
        label relativeAddress = addressValue-offset;

        // Reconstruct hex word from address
        OStringStream nStream;
        nStream << "0x" << hex << relativeAddress;
        myAddress = nStream.str();
    }

    if (filename[0] == '/')
    {
        string line = pOpen
        (
#ifndef darwin
            "addr2line -f --demangle=auto --exe "
#else
            //            "gaddr2line -f --inline --demangle=auto --exe "
            "addr2line4Mac.py "
#endif
          + filename
          + " "
          + myAddress,
            1
        );

        if (line == "")
        {
            os << " addr2line failed";
        }
        else if (line == "??:0")
        {
            os << " in " << filename;
        }
        else
        {
            string cwdLine(line.replaceAll(cwd() + '/', ""));

            string homeLine(cwdLine.replaceAll(home(), '~'));

            os << " at " << homeLine.c_str();
        }
    }
}

#ifdef darwin

// Trying to emulate the original backtrace and backtrace_symbol from the glibc
// After an idea published by Rush Manbert at http://lists.apple.com/archives/xcode-users/2006/Apr/msg00528.html

template<int level>
void *getStackAddress() 
{
    const unsigned int stackLevel=level;
    return (
        __builtin_frame_address(level) 
        ? __builtin_return_address(stackLevel)
        : (void *)0
    );
};

#define GET_STACK_ADDRESS(lvl)              \
    case lvl: {return getStackAddress<lvl>(); break; }

// please don't laugh. For some reason this is necessary (the compiler won't accept it otherwise)
void *getStackAddress(int level)
{
    switch(level) {
        GET_STACK_ADDRESS(0);
        GET_STACK_ADDRESS(1);
        GET_STACK_ADDRESS(2);
        GET_STACK_ADDRESS(3);
        GET_STACK_ADDRESS(4);
        GET_STACK_ADDRESS(5);
        GET_STACK_ADDRESS(6);
        GET_STACK_ADDRESS(7);
        GET_STACK_ADDRESS(8);
        GET_STACK_ADDRESS(9);
        GET_STACK_ADDRESS(10);
        GET_STACK_ADDRESS(11);
        GET_STACK_ADDRESS(12);
        GET_STACK_ADDRESS(13);
        GET_STACK_ADDRESS(14);
        GET_STACK_ADDRESS(15);
        GET_STACK_ADDRESS(16);
        GET_STACK_ADDRESS(17);
        GET_STACK_ADDRESS(18);
        GET_STACK_ADDRESS(19);
        GET_STACK_ADDRESS(20);
        GET_STACK_ADDRESS(21);
        default:
            return (void *)0;
            break;
    }
}

unsigned backtrace(void **bt, unsigned maxAddrs) 
{
    unsigned valid=0;
    bool ok=true;

    for(int level=0;level<maxAddrs;level++) {
        if(ok) {
            bt[level]=getStackAddress(level);
            
            if(bt[level]!=(void *)0) {
                valid=level;
            } else {
                ok=false;
            }
        } else {
            bt[level]=(void *)0;
        }
    }

    return valid;
}

// This function is a potential memory leak. But I don't care because the program is terminating anyway
char **backtrace_symbols(void **bt,unsigned nr) 
{
    char **strings=(char **)malloc(sizeof(char *)*nr);

    for(unsigned i=0;i<nr;i++) {
        Dl_info info;
        int result=dladdr(bt[i],&info);

        char tmp[1000];
#ifdef darwinIntel64
        sprintf(tmp,"%s(%s+%p) [%p]",info.dli_fname,info.dli_sname,(void *)((unsigned long)bt[i]-(unsigned long)info.dli_saddr),bt[i]);
#else
        sprintf(tmp,"%s(%s+%p) [%p]",info.dli_fname,info.dli_sname,(void *)((unsigned int)bt[i]-(unsigned int)info.dli_saddr),bt[i]);
#endif
        strings[i]=(char *)malloc(strlen(tmp)+1);
        strcpy(strings[i],tmp);
    }

    return strings;
}

#endif

void getSymbolForRaw
(
    Ostream& os,
    const string& raw,
    const fileName& filename,
    const word& address
)
{
    if (filename.size() && filename[0] == '/')
    {
        string fcnt = pOpen
        (
#ifndef darwin
            "addr2line -f --demangle=auto --exe "
#else
            //            "gaddr2line -f --inline --demangle=auto --exe "
            "addr2line4Mac.py "
#endif
          + filename
          + " "
          + address
        );

        if (fcnt != "")
        {
            os << fcnt.c_str();
            return;
        }
    }
    os << "Uninterpreted: " << raw.c_str();
}

void error::printStack(Ostream& os)
{
    // Do not print anything if FOAM_ABORT is not set
    if (!env("FOAM_ABORT"))
    {
        return;
    }

    // Reads the starting addresses for the dynamically linked libraries
    // from the /proc/pid/maps-file
    // I'm afraid this works only for Linux 2.6-Kernels (may work on 2.4)
    // Note2: the filenames in here will have softlinks resolved so will
    // go wrong when having e.g. OpenFOAM installed under a softlink.

    HashTable<label, fileName> addressMap;
    {
        IFstream is("/proc/" + name(pid()) + "/maps");

        while(is.good())
        {
            string line;
            is.getLine(line);

            string::size_type space = line.rfind(' ') + 1;
            fileName libPath = line.substr(space, line.size()-space);

            if (libPath.size() && libPath[0] == '/')
            {
                string offsetString(line.substr(0, line.find('-')));
                IStringStream offsetStr(offsetString);
                addressMap.insert(libPath, readHexLabel(offsetStr));
            }
        }
    }

    // Get raw stack symbols
    void *array[100];
    size_t size = backtrace(array, 100);
    char **strings = backtrace_symbols(array, size);

    // See if they contain function between () e.g. "(__libc_start_main+0xd0)"
    // and see if cplus_demangle can make sense of part before +
    // HJ, formatting of stack backtrace.  17/Dec/2008
    os << nl;
    for (size_t i = 0; i < size; i++)
    {
        string msg(strings[i]);
        fileName programFile;
        word address;

        os << '#' << label(i) << "  ";
        //os << "Raw   : " << msg << "\n\t";
        {
            string::size_type lPos = msg.find('[');
            string::size_type rPos = msg.find(']');

            if (lPos != string::npos && rPos != string::npos && lPos<rPos)
            {
                address = msg.substr(lPos+1, rPos-lPos-1);
                msg = msg.substr(0, lPos);
            }

            string::size_type bracketPos = msg.find('(');
            string::size_type spacePos = msg.find(' ');
            if (bracketPos != string::npos || spacePos != string::npos)
            {
                programFile = msg.substr(0, min(spacePos, bracketPos));

                // not an absolute path
                if (programFile[0] != '/')
                {
                    string tmp = pOpen("which " + programFile);
                    if (tmp[0] == '/' || tmp[0] == '~')
                    {
                        programFile = tmp;
                    }
                }
            }
        }

        string::size_type bracketPos = msg.find('(');

        if (bracketPos != string::npos)
        {
            string::size_type start = bracketPos+1;

            string::size_type plusPos = msg.find('+', start);

            if (plusPos != string::npos)
            {
                string cName(msg.substr(start, plusPos-start));

                int status;
                char* cplusNamePtr = abi::__cxa_demangle
                (
                    cName.c_str(),
                    NULL,                   // have it malloc itself
                    0,
                    &status
                );

                if (status == 0 && cplusNamePtr)
                {
                    os << cplusNamePtr;
                    free(cplusNamePtr);
                }
                else
                {
                    os << cName.c_str();
                }
            }
            else
            {
                string::size_type endBracketPos = msg.find(')', start);

                if (endBracketPos != string::npos)
                {
                    string fullName(msg.substr(start, endBracketPos-start));

                    os << fullName.c_str() << nl;
                }
                else
                {
                    // Print raw message
                    getSymbolForRaw(os, msg, programFile, address);
                }
            }
        }
        else
        {
            // Print raw message
            getSymbolForRaw(os, msg, programFile, address);
        }

        printSourceFileAndLine(os, addressMap, programFile, address);

        os << nl;
    }

    free(strings);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
