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

Application
    fluentDataToFoam

Description
    Converts Fluent data to OpenFOAM format

\*---------------------------------------------------------------------------*/

%{

#undef yyFlexLexer

/* ------------------------------------------------------------------------- *\
   ------ local definitions
\* ------------------------------------------------------------------------- */

#include "fvCFD.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


label dimensionOfGrid = 0;

// Machine config
label mc1 = 0;
label mc2 = 0;
label mc3 = 0;
label mc4 = 0;
label mc5 = 0;
label mc6 = 0;
label mc7 = 0;
label mc8 = 0;
label mc9 = 0;
label mc10 = 0;
label mc11 = 0;

label nPoints = 0;
label nFaces = 0;
label nCells = 0;

SLList<label> fieldID;
SLList<label> zoneID;
SLList<label> nEntriesPerObject;
SLList<label> firstID;
SLList<label> lastID;

SLPtrList<FieldField<Field, scalar> > zoneData;

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
#if YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}

%}

one_space                  [ \t\f]
space                      {one_space}*
some_space                 {one_space}+
cspace                     ","{space}
spaceNl                    ({space}|\n|\r)*

alpha                      [_[:alpha:]]
digit                      [[:digit:]]
decDigit                   [[:digit:]]
octalDigit                 [0-7]
hexDigit                   [[:xdigit:]]

lbrac                      "("
rbrac                      ")"
quote                      \"
dash                       "-"
comma                      ","
dotColonDash               [.:-]

schemeSpecialInitial       [!$%&*/:<=>?~_^#.]
schemeSpecialSubsequent    [.+-]
schemeSymbol               (({some_space}|{alpha}|{quote}|{schemeSpecialInitial})({alpha}|{quote}|{digit}|{schemeSpecialInitial}|{schemeSpecialSubsequent})*)


identifier                 {alpha}({alpha}|{digit})*
integer                    {decDigit}+
label                      [1-9]{decDigit}*
hexLabel                   {hexDigit}+
zeroLabel                  {digit}*

word                       ({alpha}|{digit}|{dotColonDash})*

exponent_part              [eE][-+]?{digit}+
fractional_constant        [-+]?(({digit}*"."{digit}+)|({digit}+".")|({digit}))

double                     ((({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))|0)

x                          {double}
y                          {double}
z                          {double}
scalar                     {double}
labelListElement           {space}{zeroLabel}
hexLabelListElement        {space}{hexLabel}
scalarListElement          {space}{double}
schemeSymbolListElement    {space}{schemeSymbol}
labelList                  ({labelListElement}+{space})
hexLabelList               ({hexLabelListElement}+{space})
scalarList                 ({scalarListElement}+{space})
schemeSymbolList           ({schemeSymbolListElement}+{space})

starStar                   ("**")
text                       ({space}({comma}*{word}*{space})*)

dateDDMMYYYY               ({digit}{digit}"/"{digit}{digit}"/"{digit}{digit}{digit}{digit})
dateDDMonYYYY              ((({digit}{digit}{space})|({digit}{space})){alpha}*{space}{digit}{digit}{digit}{digit})
time                       ({digit}{digit}":"{digit}{digit}":"{digit}{digit})

versionNumber              ({digit}|".")*

comment                    {spaceNl}"(0"{space}
header                     {spaceNl}"(1"{space}
dimension                  {spaceNl}"(2"{space}
machineConfig              {spaceNl}"(4"{space}
gridSize                   {spaceNl}"(33"{space}
variables                  {spaceNl}"(37"{space}
fieldData                  {spaceNl}"(300"{space}

endOfSection               {space}")"{space}



 /* ------------------------------------------------------------------------- *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------- */

%option stack

%x readComment
%x embeddedCommentState
%x readHeader
%x readDimension
%x readMachineConfigHeader
%x readMachineConfig
%x readGridSizeHeader
%x readGridSize

%x readFieldData
%x readFieldDataHeader
%x readFieldSize
%x readField

%x unknownBlock
%x embeddedUnknownBlock
%%

%{
    // Data index
    label curNEntriesPerObject = 0;

    // Dummy, not used
    label nTimeLevels;
    label nPhases;

    label curFirstID = 0;
    label curLastID = 0;
    label objI = 0;

    FieldField<Field, scalar>* curZdPtr = NULL;
%}


 /* ------------------------------------------------------------------------- *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------- */

 /*                      ------ Reading control header ------                 */

{comment} {
        yy_push_state(readComment);
    }


<readComment>{quote}{text}{quote} {
    }


<readComment>{spaceNl}{endOfSection} {
        yy_pop_state();
    }

{header} {
        BEGIN(readHeader);
    }

<readHeader>{quote}{text}{quote} {
        Info<< "Reading header: " << YYText() << endl;
    }


{dimension} {
        BEGIN(readDimension);
    }

<readDimension>{space}{label}{space} {
        IStringStream dimOfGridStream(YYText());

        dimensionOfGrid = readLabel(dimOfGridStream);

        Info<< "Dimension of grid: " << dimensionOfGrid << endl;
    }


{machineConfig} {
        yy_push_state(readMachineConfigHeader);
    }

<readMachineConfigHeader>{spaceNl}{lbrac} {
        yy_push_state(readMachineConfig);
    }

<readMachineConfig>{space}{labelList} {
        IStringStream machineConfigStream(YYText());

        mc1 = readLabel(machineConfigStream);
        mc2 = readLabel(machineConfigStream);
        mc3 = readLabel(machineConfigStream);
        mc4 = readLabel(machineConfigStream);
        mc5 = readLabel(machineConfigStream);
        mc6 = readLabel(machineConfigStream);
        mc7 = readLabel(machineConfigStream);
        mc8 = readLabel(machineConfigStream);
        mc9 = readLabel(machineConfigStream);
        mc10 = readLabel(machineConfigStream);
        mc11 = readLabel(machineConfigStream);

        Info<< "Machine config: "
            << mc1 << mc2 << mc3 << mc4 << mc5 << mc6
            << mc7 << mc8 << mc9 << mc10 << mc11
            << endl;
    }


{gridSize} {
        yy_push_state(readGridSizeHeader);
    }


<readGridSizeHeader>{spaceNl}{lbrac} {
        yy_push_state(readGridSize);
    }


<readGridSize>{space}{labelList} {
        IStringStream gridSizeStream(YYText());

        nCells = readLabel(gridSizeStream);
        nFaces = readLabel(gridSizeStream);
        nPoints = readLabel(gridSizeStream);

        Info<< "Grid size: nCells = " << nCells << " nFaces = " << nFaces
            << " nPoints = " << nPoints << endl;
    }


{variables} {
        yy_push_state(readComment);
    }


{fieldData} {
        yy_push_state(readFieldDataHeader);
    }


<readFieldDataHeader>{spaceNl}{lbrac} {
        BEGIN(readFieldSize);
    }


<readFieldSize>{space}{labelList} {
        IStringStream fieldDataStream(YYText());

        // Field and zone ID
        fieldID.append(readLabel(fieldDataStream));
        zoneID.append(readLabel(fieldDataStream));

        // Number of entries per object (cell/face)
        curNEntriesPerObject = readLabel(fieldDataStream);
        nEntriesPerObject.append(curNEntriesPerObject);

        // Dummy, not used
        nTimeLevels = readLabel(fieldDataStream);
        nPhases = readLabel(fieldDataStream);

        // Start and end of list
        curFirstID = readLabel(fieldDataStream);
        firstID.append(curFirstID);
        curLastID = readLabel(fieldDataStream);
        lastID.append(curLastID);

        // Create field for data
//         Info << "Create : " << curNEntriesPerObject << endl;
        curZdPtr = new FieldField<Field, scalar>(curNEntriesPerObject);

//         Info<< "Setting field with " << curNEntriesPerObject
//             << " components of length " << curLastID - curFirstID +1 << endl;

        for (label cmptI = 0; cmptI < curNEntriesPerObject; cmptI++)
        {
            curZdPtr->set(cmptI, new scalarField(curLastID - curFirstID + 1));
        }

        zoneData.append(curZdPtr);
    }


<readFieldSize>{endOfSection} {
        BEGIN(readFieldData);
    }


<readFieldData>{spaceNl}{lbrac} {
//         Info<< "Reading field. nEntries = " << curNEntriesPerObject
//             << " firstID = " << curFirstID << " lastID = " << curLastID << endl;

        // Reset the counter
        objI = 0;

        yy_push_state(readField);
    }


<readField>{spaceNl}{scalarList} {

        IStringStream fieldStream(YYText());

        FieldField<Field, scalar>& zd = *curZdPtr;

        // Reading field
        for (label entryI = 0; entryI < curNEntriesPerObject; entryI++)
        {
            zd[entryI][objI] = readScalar(fieldStream);
        }
        objI++;
    }


<readMachineConfig,readMachineConfigHeader,readGridSize,readGridSizeHeader,readFieldData,readField>{spaceNl}{endOfSection} {

        yy_pop_state();
    }

 /*             ------ Reading end of section and others ------               */

<readHeader,readDimension>{spaceNl}{endOfSection} {
        BEGIN(INITIAL);
    }

 /*    ------ Reading unknown type or non-standard comment ------             */


{lbrac}{label} {
//         Info<< "Found unknown block:" << YYText() << endl;
        yy_push_state(unknownBlock);
    }

<readComment,unknownBlock,embeddedUnknownBlock>{spaceNl}{schemeSymbol} {
    }

<readComment,unknownBlock,embeddedUnknownBlock>{spaceNl}{lbrac} {
//         Info<< "Embedded blocks in comment or unknown:" << YYText() << endl;
        yy_push_state(embeddedUnknownBlock);

    }

<readComment,unknownBlock,embeddedUnknownBlock>{spaceNl}{endOfSection} {
//         Info<< "Found end of section in unknown:" << YYText() << endl;
        yy_pop_state();
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{labelList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{hexLabelList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{scalarList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{schemeSymbolList} {
    }

<unknownBlock,embeddedUnknownBlock>{spaceNl}{text} {
    }


 /* ------ Ignore remaining space and \n s.  Any other characters are errors. */

 /*  ------ On EOF return to previous file, if none exists terminate. ------  */

<<EOF>> {
            yyterminate();
    }
%%


#include "fileName.H"
#include "fluentDataConverter.H"

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("Fluent data file");

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

#   include "createTime.H"

    fileName fluentFile(args.additionalArgs()[0]);
    std::ifstream fluentStream(fluentFile.c_str());

    if (!fluentStream)
    {
        FatalErrorIn("fluentToFoam::main(int argc, char *argv[])")
            << args.executable()
            << ": file " << fluentFile << " not found"
            << exit(FatalError);
    }

    yyFlexLexer lexer(&fluentStream);
    while(lexer.yylex() != 0)
    {}

    Info<< "\n\nFINISHED LEXING\n\n\n";

    // Re-package the fields into OpenFOAM data

    // Warning: syncronous iterators, check size
    // HJ, 25/May/2009
    const label iterSize = fieldID.size();
    Info << "Number of entries read: " << iterSize << nl << endl;

    if
    (
        zoneID.size() != iterSize
     || nEntriesPerObject.size() != iterSize
     || firstID.size() != iterSize
     || lastID.size() != iterSize
     || zoneData.size() != iterSize
    )
    {
        FatalErrorIn(args.executable())
            << "Problem in reading: incorrect iterator size: " << iterSize
            << abort(FatalError);
    }

    labelHashSet fields;

    for
    (
        SLList<label>::const_iterator fieldIDIter = fieldID.begin();
        fieldIDIter != fieldID.end();
        ++fieldIDIter
    )
    {
        fields.insert(fieldIDIter());
    }

    // Info: Available units
    {
        labelList u = fields.toc();
        sort(u);
        Info<< "Available units: " << u << endl;
    }

    // Create a mesh
#   include "createMesh.H"

    // Create a converter
    fluentDataConverter fdc
    (
        mesh,
        fieldID,
        zoneID,
        firstID,
        lastID,
        zoneData
    );

    if (fields.found(1))
    {
        Info << "Converting pressure field" << endl;
        fdc.convertField
        (
            "p",   // field name
            1,     // fluent Unit number
            dimensionedScalar("zero", dimPressure, 0)
        )().write();
    }

    if (fields.found(2))
    {
        Info << "Converting momentum field" << endl;
        fdc.convertField
        (
            "momentum",   // field name
            2,            // fluent Unit number
            dimensionedScalar("zero", dimDensity*dimVelocity, 0)
        )().write();
    }

    if (fields.found(3))
    {
        Info << "Converting temperature field" << endl;
        fdc.convertField
        (
            "T",   // field name
            3,     // fluent Unit number
            dimensionedScalar("zero", dimTemperature, 273.15)
        )().write();
    }

    if (fields.found(4))
    {
        Info << "Converting enthalpy field" << endl;
        fdc.convertField
        (
            "h",   // field name
            4,     // fluent Unit number
            dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0)
        )().write();
    }

    if (fields.found(5))
    {
        Info << "Converting turbulence kinetic energy field" << endl;
        fdc.convertField
        (
            "k",   // field name
            5,     // fluent Unit number
            dimensionedScalar("zero", sqr(dimVelocity), 0)
        )().write();
    }

    // Species not done yet.  HJ, 29/May/2009

    if (fields.found(8))
    {
        Info << "Converting G field" << endl;
        fdc.convertField
        (
            "G",   // field name
            8,           // fluent Unit number
            dimensionedScalar("zero", dimless, 0)
        )().write();
    }

    if (fields.found(15))
    {
        Info << "Converting body force field" << endl;
        fdc.convertField
        (
            "bodyForce",   // field name
            15,            // fluent Unit number
            dimensionedScalar("zero", dimAcceleration, 0)
        )().write();
    }

    if (fields.found(101))
    {
        Info << "Converting density field" << endl;
        fdc.convertField
        (
            "rho",   // field name
            101,            // fluent Unit number
            dimensionedScalar("zero", dimDensity, 0)
        )().write();
    }

    if (fields.found(102))
    {
        Info << "Converting laminar viscosity field" << endl;
        fdc.convertField
        (
            "muLam",   // field name
            102,       // fluent Unit number
            dimensionedScalar("zero", dimPressure*dimTime, 0)
        )().write();
    }

    if (fields.found(103))
    {
        Info << "Converting turbulent field" << endl;
        fdc.convertField
        (
            "muT",   // field name
            103,     // fluent Unit number
            dimensionedScalar("zero", dimPressure*dimTime, 0)
        )().write();
    }

    if (fields.found(104))
    {
        Info << "Converting  field" << endl;
        fdc.convertField
        (
            "Cp",   // field name
            104,    // fluent Unit number
            dimensionedScalar("zero", dimSpecificHeatCapacity, 0)
        )().write();
    }

    // Velocity field
    if (fields.found(111) || fields.found(112) || fields.found(113))
    {
        Info << "Converting velocity field" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimVelocity, vector::zero)
        );

        if (fields.found(111))
        {
            U.replace
            (
                vector::X, 
                fdc.convertField
                (
                    "Ux",   // field name
                    111,           // fluent Unit number
                    dimensionedScalar("zero", dimVelocity, 0)
                )
            );
        }

        if (fields.found(112))
        {
            U.replace
            (
                vector::Y, 
                fdc.convertField
                (
                    "Uy",   // field name
                    112,    // fluent Unit number
                    dimensionedScalar("zero", dimVelocity, 0)
                )
            );
        }

        if (fields.found(113))
        {
            U.replace
            (
                vector::Z, 
                fdc.convertField
                (
                    "Uz",   // field name
                    113,    // fluent Unit number
                    dimensionedScalar("zero", dimVelocity, 0)
                )
            );
        }

        U.write();
    }


    Info<< nl << "End" << endl;
    return 0;
}


 /* ------------------------------------------------------------------------- *\
    ------ End of fluentDataToFoam.L
 \* ------------------------------------------------------------------------- */
