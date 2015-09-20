#!/usr/bin/perl
#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     | Version:     3.2
#   \\  /    A nd           | Web:         http://www.foam-extend.org
#    \\/     M anipulation  | For copyright notice see file Copyright
#------------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     parmetis_includes_hack.pl
#
# Description
#     Adds OpenMPI lib/includes dirs to CMake-generated GCC options.
#
# Author:
#     Cesare Guardino, Alstom Power Ltd., (2015)
#
#------------------------------------------------------------------------------

use strict;
use warnings;

use File::Find;

my $MPI_ROOTDIR = $ENV{"MPI_ROOTDIR"};
$MPI_ROOTDIR =~ s/^\/\w//;
my $drive_letter = $&;
$drive_letter =~ s/\///;
$MPI_ROOTDIR = uc($drive_letter) . ":" . $MPI_ROOTDIR;

my @dirs = (".");
find(\&wanted, @dirs);

sub wanted
{
    my $file = $_;
    my $path = $File::Find::name;

    if ($file eq "linklibs.rsp" or $file eq "includes_C.rsp")
    {
        open (FILE, '<', $file) or die ("ERROR: Can't open '$path' [$!]");
        my @contents = <FILE>;
        close (FILE);

        my $string = ($file eq "linklibs.rsp") ? "-L$MPI_ROOTDIR/lib -lmpi" : "-I$MPI_ROOTDIR/include";
        open (FILE, '>', $file) or die ("ERROR: Can't open '$path' [$!]");
        foreach my $line (@contents)
        {
            chomp($line);
            print FILE $line . $string;
        }
        close (FILE);
    }
}
