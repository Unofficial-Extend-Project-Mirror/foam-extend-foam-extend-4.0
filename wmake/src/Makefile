#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright held by original author
#    \\/     M anipulation  |
#-------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation; either version 2 of the License, or (at your
#     option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM; if not, write to the Free Software Foundation,
#     Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
#
# Script
#     Makefile
#
# Description
#     Makefile used to compile wmake utilities.
#
#------------------------------------------------------------------------------
# Set shell the Makefile uses to the Bourne shell
#------------------------------------------------------------------------------

SHELL      = /bin/sh

#------------------------------------------------------------------------------
# Initialise suffices list
#------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .o

RULES_DIR = $(WM_DIR)/rules/$(WM_ARCH)$(WM_COMPILER)

all:  $(RULES_DIR)/dirToString $(RULES_DIR)/wmkdep

$(RULES_DIR)/dirToString: dirToString.c
	$(cc) $(cFLAGS) dirToString.c -o $(RULES_DIR)/dirToString

$(RULES_DIR)/wmkdep: wmkdep.l
	flex wmkdep.l ; $(cc) $(cFLAGS) lex.yy.c -o $(RULES_DIR)/wmkdep

#------------------------------------------------------------------------------
# set compilation and dependency building rules
#------------------------------------------------------------------------------

GENERAL_RULES = $(WM_DIR)/rules/General
RULES         = $(RULES_DIR)

include $(RULES)/general
include $(RULES)/$(WM_LINK_LANGUAGE)

#------------------------------------------------------------------------------
