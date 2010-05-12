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

InNamespace
    Foam

Description
    Functions used by OpenFOAM which are specific to the UNIX operating system
    and need to be replaced or emulated on other systems.

SourceFiles
    Unix.C

\*---------------------------------------------------------------------------*/

#ifndef OSspecific_H
#define OSspecific_H

#include "fileNameList.H"
#include "long.H"

#include <sys/types.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return the PID of this process
pid_t pid();

//- Return the parent PID of this process
pid_t ppid();

//- Return the group PID of this process
pid_t pgid();

//- Return true if environment variable of given name is defined
bool env(const word&);

//- Return environment variable of given name
//  Return string::null if the environment is undefined
string getEnv(const word& name);

//- Set an environment variable
bool setEnv(const word& name, const string& value, const bool overwrite);

//- Return the system's host name
word hostName();

//- Return the user's login name
word userName();

//- Return home directory path name for the current user
fileName home();

//- Return home directory path name for a particular user
fileName home(const word& userName);

//- Return current working directory path name
fileName cwd();

//- Change the current directory to the one given and return true,
//  else return false
bool chDir(const fileName& dir);

//- Search for @em name in the following hierarchy:
//  -# personal settings:
//    - ~/.OpenFOAM/\<VERSION\>/
//      <em>for version-specific files</em>
//    - ~/.OpenFOAM/
//      <em>for version-independent files</em>
//  -# site-wide settings:
//    - $WM_PROJECT_INST_DIR/site/\<VERSION\>
//      <em>for version-specific files</em>
//    - $WM_PROJECT_INST_DIR/site/
//      <em>for version-independent files</em>
//  -# shipped settings:
//    - $WM_PROJECT_DIR/etc/
//
//  @return the full path name or fileName::null if the name cannot be found
fileName dotFoam(const fileName& name);

//- Make a directory and return an error if it could not be created
//  and does not already exist
bool mkDir(const fileName&, mode_t=0777);

//- Set the file mode
bool chmod(const fileName&, const mode_t);

//- Return the file mode
mode_t mode(const fileName&);

//- Return the file type: FILE or DIRECTORY
fileName::Type type(const fileName&);

//- Does the name exist in the filing system?
bool exists(const fileName& name);

//- Does the file exist?
bool file(const fileName&);

//- Does the directory exist?
bool dir(const fileName&);

//- Return size of file
off_t size(const fileName&);

//- Return time of last file modification
time_t lastModified(const fileName&);

//- Read a directory and return the entries as a string list
fileNameList readDir
(
    const fileName&,
    const fileName::Type=fileName::FILE,
    const bool filtergz=true
);

//- Copy, recursively if necessary, the source to the destination
bool cp(const fileName& srcFile, const fileName& destFile);

//- Create a softlink. destFile should not exist. Returns true if successful.
bool ln(const fileName& srcFile, const fileName& destFile);

//- Rename srcFile destFile
bool mv(const fileName& srcFile, const fileName& destFile);

//- Remove a file returning true if successful otherwise false
bool rm(const fileName&);

//- Remove a dirctory and its contents
bool rmDir(const fileName&);

//- Sleep for the specified number of seconds
unsigned int sleep(const unsigned int);

//- Close file descriptor
void fdClose(const int);

//- Check if machine is up by pinging given port
bool ping(const word&, const label port, const label timeOut);

//- Check if machine is up by ping port 22 = ssh and 222 = rsh
bool ping(const word&, const label timeOut=10);

//- Executes a command specified in command
int system(const string& command);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
