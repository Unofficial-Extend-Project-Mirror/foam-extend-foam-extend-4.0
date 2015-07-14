import os, sys, string, subprocess, re
from os.path import join, getsize

origString = "#include \"OF_dlLibraryTable.H\""
newString = ""

def CollectSourceDirectories(rootDirectory) :
    sourceDirectories = []
    for root,dirs,files in os.walk(rootDirectory) :
        for dirname in dirs :
            if dirname.find("Make") != -1 :
                sourceDirectories.append(root)

    return sourceDirectories
                
def ReplaceStrings(sourceDir) :
    global newString, origString
    try :
        # read the content of the file
        options = []
        for root,dirs,files in os.walk(sourceDir) :
            for file in files :
                if file == "options" or file == "files":
                    fName = join(root, file)
                    file = open(fName, 'raU')
                    lines = file.readlines()
                    file.close()

                    hasString = False
                    for line in lines :
                        if origString in line :
                            hasString = True
                            break
                            
                    if hasString == False :
                        continue
                        
                    newlines = []
                    for line in lines :
                        newlines.append(string.replace(line, origString, newString))

                    # write processed file
                    newFile = open(fName, 'wa')
                    for line in newlines :
                        newFile.write(line)
                    newFile.close()
    except :
        print "Cannot open options file",options

# start of the workflow
# find the root of the installation
if len(sys.argv) == 4 :
    startDir = sys.argv[1]
    origString = sys.argv[2]
    newString = sys.argv[3]
else :
    startDir = sys.argv[1]

print "Installation directory is",startDir
print "Orig string",origString
print "New string",newString

# find directories containing source files
print "Finding directories containing source files"
sourceDirectories = CollectSourceDirectories(startDir)
print "Number of directories with source files is", len(sourceDirectories)

# check which souce directories contain declareRunTimeSelectionTable
for dir in sourceDirectories :
    ReplaceStrings(dir)
