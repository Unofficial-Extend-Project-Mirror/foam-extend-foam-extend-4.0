/*
freesteam - IAPWS-IF97 steam tables library
Copyright (C) 2004-2009  John Pye

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#ifndef FREESTEAM_ZEROIN_H
#define FREESTEAM_ZEROIN_H

#include "common.h"

/**
	Any function that you want to solve using the 'zeroin_solve' function
	needs to conform to this function prototype. Any parameters required
	internally by the function can be passed in to the function as a pointer
	to a struct, array, etc, via the void 'user_data' pointer. Your subject
	function can then cast the pointer to its correct type, and access the
	necesssary parameter data, for example:

	-- in your main code
	typedef struct{double param1, param2} MyData;
	MyData D = {100., 3.141};
	double myfunc(double x, void *user_data){
		MyData *D1 = (MyData *)user_data;
		return x * D1->param1 + D1->param2;
	}
	zeroin_solve(&myfunc, &D, ...);
*/
typedef double ZeroInSubjectFunction(double, void *user_data);

/**
	Attempt to solve the function y = f(x) = 0 by varying x between
	a lower and upper bound, using the Brent algorithm.

	Originally based on brent solver from netlib, then converted to C++ for
	used in earlier freesteam versions, and now converted back to pure C again.
	@see brent.shar at http://www.netlib.org/c/

	@param func the function being solved, must be a ZeroInSubjectFunction.
	@param lowerbound the lower bound of the range in which a root is sought
	@param upperbound the upper bound of the range in which a root is sought
	@param tol maximum permissible magnitude of the function at the solved root location
	@param user_data additional data that will be passed to the subject function func.
	@param solution (returned) the value of 'x' at the solution
	@param error (returned) the value of 'y' at the solution.
	@return 0 on success
*/

FREESTEAM_DLL char zeroin_solve(ZeroInSubjectFunction *func, void *user_data, double lowerbound, double upperbound, double tol, double *solution, double *error);

#endif

