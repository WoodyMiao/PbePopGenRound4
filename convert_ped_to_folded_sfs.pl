#!/usr/bin/perl
use strict;
use warnings;

die "
This program converts PED to folded SFS of DaDi (2 outgroups; 2 or 3 ingroup populations).\n
Only biallelic loci without missing allele are counted.\n
Auther: Woody\n
Usage: $0 <a one-line list file of PopID|\"Excluded\" of each line in the ped file> <in.ped> <out.sfs>\n\n" if @ARGV < 3;

open POP, "<", $ARGV[0];
open PED, "<", $ARGV[1];
open SFS, ">", $ARGV[2];

my @pop;
while (<POP>) {
	chomp;
	@pop = split /\s+/;
}
close POP;

my @ped;
my @exc;
for (0 .. @pop - 1) {
	my $a = <PED>;
	if ($pop[$_] eq "Excluded") {
		unshift @exc, $_;
		next;
	}
	chomp $a;
	my @b = split /\t/, $a;
	my @c = splice @b, 0, 6;
	push @ped, \@b;
#	warn "Read $c[1] complete.\n";
}
close PED;

if (@exc) {
	splice @pop, $exc[$_], 1 for (0 .. @exc - 1);
}

my %pop;
++$pop{$_} for @pop;
my $nip = keys %pop; # number of ingroup pouplations
my @ipl = sort keys %pop; # ingroup population list
my @dim; # dimensions of SFS
push @dim, 2*$pop{$_} for @ipl;
my @sfs;

if ($nip == 2) {
	for my $x (0 .. $dim[0]) {
		for my $y (0 .. $dim[1]) {
			$sfs[$x][$y] = 0;
		}
	}
} elsif ($nip == 3) {
	for my $x (0 .. $dim[0]) {
		for my $y (0 .. $dim[1]) {
			for my $z (0 .. $dim[2]) {
				$sfs[$x][$y][$z] = 0;
			}
		}
	}
} else {
	die "Error: there are $nip populations.\n";
}

my $nfs = 0; # number of filtered SNPs
my $ncs = 0; # number of counted SNPs
for my $j (0 .. @{$ped[0]} - 1) { # SNP (column) index
	my %alle; # overall allele count
	for my $i (0 .. @ped - 1) { # Sample (row) index
		my @a = split / /, $ped[$i][$j];
		++$alle{$a[0]};
		++$alle{$a[1]};
	}
	if ($alle{"0"} or keys %alle != 2) {
		++$nfs;
		next;
	}
	++$ncs;
	my $minor; # minor allele
	for (keys %alle) {
		if ($alle{$_} <= @pop) {
			$minor = $_;
			next;
		}
	}
	my %mac; # minor allele count
	$mac{$_} = 0 for @ipl;
	for my $i (0 .. @ped - 1) {
		my @a = split / /, $ped[$i][$j];
		for (@a) {
			++$mac{$pop[$i]} if $_ eq $minor;
		}
	}
	if ($nip == 2) {
		++$sfs[$mac{$ipl[0]}][$mac{$ipl[1]}];
	} else {
		++$sfs[$mac{$ipl[0]}][$mac{$ipl[1]}][$mac{$ipl[2]}];
	}
#	warn "Counting: $nfs SNPs were filtered, $ncs SNPs were counted.\n" if $ncs =~ /00000$/;
}

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
print SFS "# This SFS was generated from $ARGV[1], by $0, saved at $now.\n";
print SFS 2*$pop{$_}+1, " " for @ipl;
print SFS "unfolded";
print SFS " \"$_\"" for @ipl;
print SFS "\n";
my $num_ele;
if ($nip == 2) {
	$num_ele = ($dim[0]+1) * ($dim[1]+1) - 2;
	for my $x (0 .. $dim[0]) {
		for my $y (0 .. $dim[1]) {
			print SFS $sfs[$x][$y], " ";
		}
	}
} else {
	$num_ele = ($dim[0]+1) * ($dim[1]+1) * ($dim[2]+1) - 2;
	for my $x (0 .. $dim[0]) {
		for my $y (0 .. $dim[1]) {
			for my $z (0 .. $dim[2]) {
				print SFS $sfs[$x][$y][$z], " ";
			}
		}
	}
}
print SFS "\n1 ", "0 " x $num_ele, "1\n";
close SFS;
warn "$ARGV[2] was done: $nfs SNPs were filtered, $ncs SNPs were counted.\n";
