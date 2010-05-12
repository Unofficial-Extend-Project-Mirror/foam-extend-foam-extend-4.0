# Clears a release directory before untarring a Release-Tarball into it
# Leaves the directory structure and all .svn-directories intact

find $1 \( \( -name .svn -prune \) -or \( -type f -print -exec rm -f {} \; \) \)
