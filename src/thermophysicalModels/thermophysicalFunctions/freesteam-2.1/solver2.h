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

#ifndef FREESTEAM_SOLVER2_H
#define FREESTEAM_SOLVER2_H

#include "common.h"
#include "steam.h"

/*
	two-way solver interface, for unusual property pairs that are not
	implemented in official IAPWS releases.

	this solver will use provided 'first guess' estimates, meaning that
	if correlations, even approximate, can be provided, then this will
	decrease unneeded iteration.

	these solvers first require that you know which region you're in.
	This is because the correlation variables are different for the different
	regions.

	'freesteam_region_XX functions will calculate the region for you; but
	these need to be hand-coded at this stage.
*/

/**
	Two-variable Newton solver for Region 3, using finite difference
	approximations for derivatives. The names of the provided variables are
	given as X and Y (eg 'T' and 's') then the require state values are given
	as x and y (eg 300 and 5500). The solver will attempt to find a state
	that satisfies the provided input, and return it as output.

	On return, the variable status will be set to 0 on success, or an error code
	if something has gone wrong, eg failed to converge.

	A first guess must be provided via the 'guess' parameter.

	FIXME provide 'default' guess functions.
*/
FREESTEAM_DLL SteamState freesteam_solver2_region1(FREESTEAM_CHAR X, FREESTEAM_CHAR Y, double x, double y, SteamState guess, int *status);
FREESTEAM_DLL SteamState freesteam_solver2_region2(FREESTEAM_CHAR X, FREESTEAM_CHAR Y, double x, double y, SteamState guess, int *status);
FREESTEAM_DLL SteamState freesteam_solver2_region3(FREESTEAM_CHAR X, FREESTEAM_CHAR Y, double x, double y, SteamState guess, int *status);
FREESTEAM_DLL SteamState freesteam_solver2_region4(FREESTEAM_CHAR X, FREESTEAM_CHAR Y, double x, double y, SteamState guess, int *status);

#endif
