
FREESTEAM
=========

Note: We ripped out what we needed, here is the full description!

This is a package for calculating the properties of water and steam
using the IAPWS-IF97 industry-standard steam properties correlations.

freesteam includes a thorough set of test cases; you are encouraged to inspect
these for yourself to assure yourself that freesteam is accurate for the cases
you required.

freesteam is not a stand-alone program in its own right. Rather, it is a 'shared
library' that implements a range of functions that you can use to work out steam
properties in your own program that you write yourself.

This version 2.x release of freesteam is a complete re-write of freesteam in 
pure C. We have removed all the complicated C++ template code (excessively
complicated) as well as the units-of-measurement code (cute, but also
rather hard to maintain). This has the immediate advantage that the resulting
DLLs/libraries can be linked to from all C/C++ compilers, something which was 
not possible with the previous C++ version.

For more information, see
http://freesteam.sourceforge.net/


