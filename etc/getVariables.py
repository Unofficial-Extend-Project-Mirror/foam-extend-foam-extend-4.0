#! /usr/bin/env python

# This file sources the etc bashrc, finds the differences to the old
# configuration and makes them available in a format that the calling
# shell can understand
#

import sys

# This part is lifted from six.py (https://pythonhosted.org/six/) to
# make sure that this script runs with Python 2 and Python 3

# True if we are running on Python 3.
PY3 = sys.version_info[0] == 3

if PY3:
    import builtins
    print_ = getattr(builtins, "print")
    del builtins
else:
    def print_(*args, **kwargs):
        """The new-style print function."""
        fp = kwargs.pop("file", sys.stdout)
        if fp is None:
            return
        def write(data):
            if not isinstance(data, basestring):
                data = str(data)
            fp.write(data)
        want_unicode = False
        sep = kwargs.pop("sep", None)
        if sep is not None:
            if isinstance(sep, unicode):
                want_unicode = True
            elif not isinstance(sep, str):
                raise TypeError("sep must be None or a string")
        end = kwargs.pop("end", None)
        if end is not None:
            if isinstance(end, unicode):
                want_unicode = True
            elif not isinstance(end, str):
                raise TypeError("end must be None or a string")
        if kwargs:
            raise TypeError("invalid keyword arguments to print()")
        if not want_unicode:
            for arg in args:
                if isinstance(arg, unicode):
                    want_unicode = True
                    break
        if want_unicode:
            newline = unicode("\n")
            space = unicode(" ")
        else:
            newline = "\n"
            space = " "
        if sep is None:
            sep = space
        if end is None:
            end = newline
        for i, arg in enumerate(args):
            if i:
                write(sep)
            write(arg)
        write(end)

def printDebug(*args):
    print_(*args,file=sys.stderr)

# the actual work starts here
from os import environ,path

if sys.version_info<(2,6):
    from popen2 import popen4
else:
    from subprocess import Popen,PIPE,STDOUT

# only for development
verbose="FOAM_NEW_STARTUP_DEBUG" in environ

# Location of this script
here=path.dirname(path.abspath(sys.argv[0]))
if verbose:
    printDebug("Using scripts in",here)

# For which shell
destShell=sys.argv[1]

# Variables like WM_COMPILER=Gcc46
additional=sys.argv[2:]

if verbose:
    printDebug("Target shell",destShell)
    printDebug("Additional settings:",additional)

# Certain bashrc-s fail if these are set
for v in ["FOAM_INST_DIR",
          "WM_THIRD_PARTY_DIR",
          "WM_PROJECT_USER_DIR",
          "OPAL_PREFIX"]:
    try:
        del environ[v]
    except KeyError:
        pass

# To be executed in bash
cmd ='echo "=== Export pre";export;echo "=== Alias pre";alias;'
cmd+='echo "=== Script";. '+path.join(here,"bashrc")+' '+' '.join(additional)+';'
cmd+='echo "=== Export post";export;echo "=== Alias post";alias;'
cmd+='echo "=== End"'

cmd="bash -c '"+cmd+"'"

if verbose:
    printDebug("Cmd:",cmd)

# Execute the shell commands
if sys.version_info<(2,6):
    # for old machines (RHEL 5)
    raus,rein = popen4(cmd)
    ret=0 # standin for a proper implementation
else:
    p = Popen(cmd, shell=True,
              stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    (rein,raus)=(p.stdin,p.stdout)
    ret=p.wait()

lines=[l.strip().decode() for l in raus.readlines()]
rein.close()
raus.close()

if verbose or ret!=0:
    printDebug("========= Script output start")
    for l in lines:
        printDebug(l)
    printDebug("========= Script output end")

if ret!=0:
    printDebug("Command return code",ret,". There seems to be a problem")

def extractVariables(lines):
    vars={}
    # different forms if bash is called as sh
    prefixes=["export ","declare -x "]
    for l in lines:
        pref=None
        for p in prefixes:
            if l.find(p)==0:
                pref=p
                break
        pos=l.find("=")
        name=l[len(pref):pos]
        if pos>0:
            val=l[pos+1:]
            if val[0]=='"' and val[-1]=='"':
                val=val[1:-1]
        else:
            val=""
        vars[name]=val

    return vars

def extractAliases(lines):
    aliases={}
    for l in lines:
        pref=""  # if called as sh
        if l.find("alias ")==0:
            pref="alias "
        pos=l.find("=")
        aliases[l[len(pref):pos]]=l[pos+2:-1]

    return aliases

def changedVars(old,new):
    changed={}
    for k in new:
        if k not in old:
            changed[k]=new[k]
        elif old[k]!=new[k]:
            changed[k]=new[k]
    return changed

def splitPaths(orig):
    new={}
    for k in orig.keys():
        if k.find("PATH")>=0 and orig[k].find(":")>=0:
            new[k]=orig[k].split(":")
        else:
            new[k]=orig[k]
    return new

try:
    vars=splitPaths(
        changedVars(extractVariables(lines[lines.index("=== Export pre")+1:
                                           lines.index("=== Alias pre")]),
                    extractVariables(lines[lines.index("=== Export post")+1:
                                           lines.index("=== Alias post")])))
    aliases=changedVars(extractAliases(lines[lines.index("=== Alias post")+1:
                                             lines.index("=== Script")]),
                        extractAliases(lines[lines.index("=== Alias post")+1:
                                             lines.index("=== End")]))

    scriptOutput=lines[lines.index("=== Script")+1:
                       lines.index("=== Export post")]
except Exception:
    e = sys.exc_info()[1] # Needed because python 2.5 does not support 'as e'
    err, detail, tb = sys.exc_info()
    printDebug("Output of",cmd,"not successfully parsed")
    import traceback
    traceback.print_exception(err,detail,tb,file=sys.stderr)
    sys.exit(-1)

# Pass output of the script to stderr for debugging
if len(scriptOutput)>0:
    for l in scriptOutput:
        print_(l,file=sys.stderr)

class ShellConvert(object):
    def __call__(self,vars,aliases):
        result=""
        for v in sorted(vars.keys()):
            result+=self.toVar(v,vars[v])+"\n"
        for a in sorted(aliases.keys()):
            result+=self.toAlias(a,aliases[a])+"\n"
        return result

class BashConvert(ShellConvert):
    def toVar(self,n,v):
        if type(v)==list:
            val=":".join(v)
        else:
            val=v
        if val.find(" "):
            return 'export %s="%s"' % (n,val)
        else:
            return 'export %s=%s' % (n,val)

    def toAlias(self,n,v):
        return "alias %s='%s'" % (n,v)

class CshConvert(ShellConvert):
    def __init__(self,sName="csh"):
        self.shellName=sName

    def toVar(self,n,v):
        if type(v)==list:
            val=":".join(v)
        else:
            val=v
        result='setenv %s "%s"' % (n,val)
        if n=="PATH":
            result+="\nset path=(%s)" % " ".join(v)
        return result

    def toAlias(self,n,v):
        val=v.replace(" . ","source ").replace("bash",self.shellName)
        if val.find("unset ")>=0:
            val=val.replace("unset ","unsetenv ")
        # Make sure that more than one export is possible and no wrong = is replaced
        while val.find("export ")>=0:
            pos=val.find("export ")
            val=val.replace("export ","setenv ",1)
            val=val[:pos]+val[pos:].replace("="," ",1)

        # Prefix with a variable that tells the script that it is sourced from an alias
        if val.find("source")>=0:
            val=("setenv FOAM_SOURCED_FROM_ALIAS %s/tcshrc ; " % here)+val

        return "alias %s '%s'" % (n,val)

class TcshConvert(CshConvert):
    def __init__(self):
        CshConvert.__init__(self,"tcsh")

class FishConvert(ShellConvert):
    def toVar(self,n,v):
        if type(v)==list:
            val=":".join(v)
        else:
            val=v

        if n=="PATH":
            result="\nset -x PATH %s" % " ".join(v)
        else:
            result='set -x %s "%s"' % (n,val)

        return result

    def toAlias(self,n,v):
        val=v.replace(" . ","source ").replace("bash","fish")
        if val.find("unset ")>=0:
            val=val.replace("unset ","set -e ")
        # Make sure that more than one export is possible and no wrong = is replaced
        while val.find("export ")>=0:
            pos=val.find("export ")
            val=val.replace("export ","set -x ",1)
            val=val[:pos]+val[pos:].replace("="," ",1)
        return "alias %s '%s'" % (n,val)

class ZshConvert(BashConvert):
    def toAlias(self,n,v):
        return BashConvert.toAlias(self,n,v).replace("bash","zsh")

shells={"bash" : BashConvert,
        "zsh"  : ZshConvert,
        "fish" : FishConvert,
        "csh"  : CshConvert,
        "tcsh" : TcshConvert}

try:
    converter=shells[destShell]()
except KeyError:
    printDebug("No converter for shell",destShell,
               "Available:"," ".join(shells.keys()))
    sys.exit(-1)

try:
    result=converter(vars,aliases)
except Exception:
    e = sys.exc_info()[1] # Needed because python 2.5 does not support 'as e'
    err, detail, tb = sys.exc_info()
    printDebug("Problem while converting output for",destShell)
    import traceback
    traceback.print_exception(err,detail,tb,file=sys.stderr)
    sys.exit(-1)

if PY3:
    osException=OSError
else:
    osException=IOError

try:
    open(path.abspath(sys.argv[0])+"."+destShell,"w").write(result)
except osException:
    # We are not allowed to write here. Try tmp
    try:
        open(path.join("/tmp",
                       path.basename(sys.argv[0])+"."+destShell),"w").write(result)
    except osException:
        # Nobody wants this
        pass

# This is the only part that goes to stdout to be sourced by the calling script
print_(result)
