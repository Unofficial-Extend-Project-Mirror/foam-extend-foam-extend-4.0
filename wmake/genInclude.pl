#!/usr/bin/perl
#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     | Version:     4.0
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
#     genInclude.pl
#
# Description
#     Generates include files into lnInclude, as an alternative to creating
#     symbolic links (useful for Windows)
#
# Usage:
#     genInclude.pl <PATH_TO_SOURCE_FILE> .
#
# Author:
#     Cesare Guardino, Alstom Power Ltd., (2015)
#
#------------------------------------------------------------------------------

use strict;
use warnings;

use Cwd;
use Cwd 'abs_path';
use File::Basename;

my $cwd = cwd();

my $source = abs_path($ARGV[0]);
my $fileName = basename($source);
$cwd =~ s/lnInclude//;
$source =~ s/$cwd//;

my $dir = $cwd;
if ($dir =~ /^.*\/src\//) {
    $dir =~ s/^.*\/src\///;
}
elsif ($dir =~ /^.*\/applications\//) {
    $dir =~ s/^.*\/applications\///;
}

my $link = $dir . $source;
open (FILE, '>', $fileName) or die ("ERROR: Can't open '$fileName' [$!]");
print FILE "#include \"$link\"\n";
close (FILE);

