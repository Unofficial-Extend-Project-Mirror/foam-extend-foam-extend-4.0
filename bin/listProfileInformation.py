#! /usr/bin/env python

# Lists the profiling information in time directories (uniform/profilingInfo)
# in a human readable form

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import sys

pf=ParsedParameterFile(sys.argv[1],
                       treatBinaryAsASCII=True)

data={}
children={}
root=None

for p in pf["profilingInfo"]:
    if p["id"] in data:
        print "Duplicate definition of",p["id"]
        sys.exit(-1)
    if p["description"][0]=='"':
        p["description"]=p["description"][1:]
    if p["description"][-1]=='"':
        p["description"]=p["description"][:-1]

    data[p["id"]]=p
    if "parentId" in p:
        if p["parentId"] in children:
            children[p["parentId"]].append(p["id"])
        else:
            children[p["parentId"]]=[p["id"]]
    else:
        if root!=None:
            print "Two root elements"
            sys-exit(-1)
        else:
            root=p["id"]

def depth(i):
    if i in children:
        return max([depth(j) for j in children[i]])+1
    else:
        return 0

#make sure that children are printed in the correct order
for i in children:
    children[i].sort()

maxdepth=depth(root)

depths={}

def nameLen(i,d=0):
    depths[i]=d
    maxi=len(data[i]["description"])
    if i in children:
        maxi=max(maxi,max([nameLen(j,d+1) for j in children[i]]))
    return maxi+3

maxLen=nameLen(root)

format=" %5.1f%% (%5.1f%%) - %5.1f%% | %8d %9.4gs %9.4gs"
totalTime=data[root]["totalTime"]

header=" "*(maxLen)+" | parent  (total ) -  self  |    calls      total      self "
print header
print "-"*len(header)

def printItem(i):
    result=""
    if depths[i]>1:
        result+="  "*(depths[i]-1)
    if depths[i]>0:
        result+="|- "
    result+=data[i]["description"]
    result+=" "*(maxLen-len(result)+1)+"| "

    parentTime=data[i]["totalTime"]
    if "parentId" in data[i]:
        parentTime=data[data[i]["parentId"]]["totalTime"]

    tt=data[i]["totalTime"]
    ct=data[i]["childTime"]

    result+=format % (100*tt/parentTime,
                      100*(tt-ct)/totalTime,
                      100*(tt-ct)/tt,
                      data[i]["calls"],
                      tt,
                      tt-ct)
    print result
    if i in children:
        for c in children[i]:
            printItem(c)

printItem(root)
