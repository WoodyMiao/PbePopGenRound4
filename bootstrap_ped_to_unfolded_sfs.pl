#!/usr/bin/perl
use strict;
use warnings;

die "
This program make unfolded SFS bootstraps (100 times, output from 00.sfs to 99.sfs in current folder) from PED file (2 outgroups; 2|3|4 ingroup populations).\n
Only biallelic loci without missing allele are counted. If both outgroup individuals are monomorphic for the same allele, this allele is considered ancestral (Nater 2015).\n
Auther: Woody\n
Usage: $0 <control> <in.ped> <in.map> 

The 1st line of the control file is a list of PopID|\"Outgroup\"|\"Excluded\" for each line in the ped file.
The 2nd line of the control file is a list of PopIDs, which defines the order in the output SFS.\n\n" if @ARGV < 2;

my $nbr = 100; # number of bootstrap replicates

open CTR, "<", $ARGV[0];
open PED, "<", $ARGV[1];
open MAP, "<", $ARGV[2];

my $control = <CTR>;
chomp $control;
my @pop = split /\s+/, $control;
$control = <CTR>;
chomp $control;
my @ipl = split /\s+/, $control;
close CTR;

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
}
close PED;
warn "Reading $ARGV[1] completed.\n";

my @map;
while (<MAP>) {
	chomp;
	my @a = split /\t/;
	push @map, [$a[0], $a[3]];
}
close MAP;
warn "Reading $ARGV[2] completed.\n";

if (@exc) {
	splice @pop, $exc[$_], 1 for (0 .. @exc - 1);
}

my %pop;
++$pop{$_} for @pop;
delete $pop{Outgroup};
my $nip = keys %pop; # number of ingroup pouplations
die "Error: inconsistent number of ingroup populations!" if @ipl != $nip;
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
} elsif ($nip == 4) {
	for my $x (0 .. $dim[0]) {
		for my $y (0 .. $dim[1]) {
			for my $z (0 .. $dim[2]) {
				$sfs[$x][$y][$z][$_] = 0 for (0 .. $dim[3]);
			}
		}
	}
} else {
	die "Error: there are $nip populations.\n";
}

my @bs; # array of 1 Mb windows of @ped;
my @bs_ance;
my $chrnum = 1;
my $coordi = 1000000;
my $bsi = 0; # @bs array index;
for my $j (0 .. @{$ped[0]} - 1) { # SNP (column) index
	my %alle; # overall allele count
	my %oual; # outgroup allele count
	for my $i (0 .. @ped - 1) { # Sample (row) index
		my @a = split / /, $ped[$i][$j];
		++$alle{$a[0]};
		++$alle{$a[1]};
		if ($pop[$i] eq "Outgroup") {
			++$oual{$a[0]};
			++$oual{$a[1]};
		}
	}
	if ($alle{"0"} or keys %alle != 2 or keys %oual != 1) {
		next;
	}
	my ($ance) = keys %oual; # ancestral allele
	if ($map[$j][0] != $chrnum) {
		$chrnum += 1;
		$coordi = 1000000;
		++$bsi;
	}
	if ($map[$j][1] > $coordi) {
		$coordi += 1000000;
		++$bsi;
	}
	push @{$bs_ance[$bsi]}, $ance;
	for my $i (0 .. @ped-1) {
		my @a = split / /, $ped[$i][$j];
		push @{$bs[$bsi][$i]}, \@a;
	}
}
++$bsi;
warn "Picking SNP completed. Number of bootstrap units: $bsi\n";
undef @ped;
undef @map;

for (0 .. $nbr-1) {
	my $file = sprintf("%02d.fs", $_);
	open my $out, ">", $file;
	my @s = @sfs;
	for (1 .. $bsi) {
		my $r = int(rand($bsi));
		if (!$bs[$r][0][0][0]) {
			warn "There is no SNP in the No. $r bootstrap units.\n";
			next;
		}
		for my $j (0 .. @{$bs[$r][0]}-1) {
			my %dac; # derived allele count
			$dac{$_} = 0 for @ipl;
			for my $i (0 .. @{$bs[0]}-1) {
				next if $pop[$i] eq "Outgroup";
				for (@{$bs[$r][$i][$j]}) {
					++$dac{$pop[$i]} if $_ ne $bs_ance[$r][$j];
				}
			}
			if ($nip == 2) {
				++$s[$dac{$ipl[0]}][$dac{$ipl[1]}];
			} elsif ($nip == 3) {
				++$s[$dac{$ipl[0]}][$dac{$ipl[1]}][$dac{$ipl[2]}];
			} else {
				++$s[$dac{$ipl[0]}][$dac{$ipl[1]}][$dac{$ipl[2]}][$dac{$ipl[3]}];
			}
		}
	}
	print $out 2*$pop{$_}+1, " " for @ipl;
	print $out "unfolded";
	print $out " \"$_\"" for @ipl;
	print $out "\n";
	my $num_ele;
	if ($nip == 2) {
		$num_ele = ($dim[0]+1) * ($dim[1]+1) - 2;
		for my $x (0 .. $dim[0]) {
			for my $y (0 .. $dim[1]) {
				print $out $s[$x][$y], " ";
			}
		}
	} elsif ($nip == 3) {
		$num_ele = ($dim[0]+1) * ($dim[1]+1) * ($dim[2]+1) - 2;
		for my $x (0 .. $dim[0]) {
			for my $y (0 .. $dim[1]) {
				for my $z (0 .. $dim[2]) {
					print $out $s[$x][$y][$z], " ";
				}
			}
		}
	} else {
		$num_ele = ($dim[0]+1) * ($dim[1]+1) * ($dim[2]+1) * ($dim[3]+1) - 2;
		for my $x (0 .. $dim[0]) {
			for my $y (0 .. $dim[1]) {
				for my $z (0 .. $dim[2]) {
					print $out $s[$x][$y][$z][$_], " " for (0 .. $dim[3]);
				}
			}
		}
	}
	print $out "\n1 ", "0 " x $num_ele, "1\n";
	close $out;
}
