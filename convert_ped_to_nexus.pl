#!/usr/bin/perl
use strict;
use warnings;

die "This program converts PED to NEXUS.\nAuther: Woody
Usage: $0 <in.ped> <in.map> <out.nex>\n" if @ARGV < 3;

open NTAX, "-|", "wc -l $ARGV[0]";
my $wc = <NTAX>;
my $ntax = (split / /, $wc)[0];
close NTAX;

open NCHAR, "-|", "wc -l $ARGV[1]";
$wc = <NCHAR>;
my $nchar = (split / /, $wc)[0];
close NCHAR;

open PED, "<", $ARGV[0];
open OUT, ">", $ARGV[2];

my %base = ("A A","A", "C C","C", "G G","G", "T T","T", "A G","R", "C T","Y", "A C","M", "G T","K", "G C","S", "A T","W", "0 0","N");

print OUT "#NEXUS
BEGIN DATA;
Dimensions ntax=$ntax nchar=$nchar;
Format datatype=DNA gap=- missing=?;
Matrix\n";

while (<PED>) {
	chomp;
	my @a = split /\t/;
	shift @a;
	my $id = shift @a;
	splice @a, 0, 4;
	die "$id error: inconsistent number of SNPs between $ARGV[0] and $ARGV[1].\n" if @a != $nchar;
	print OUT "$id\t", join("", @base{@a}), "\n";
}

print OUT ";\nEND;\n";
close PED;
close OUT;
