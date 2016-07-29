#!/usr/bin/perl
use strict;
use warnings;

die "
This program makes nexus alignment from PED file using only SNPs for counting unfolded SFS.\n
Loci satisfying all the following conditions simultaneously are outputed: 1. no missing allele, 2. overall biallelic, 3. outgroup monoallelic, 4. ingroup biallelic\n
Auther: Woody\n
Usage: $0 <control> <in.ped> <out.nex>

The one line control file is a list of PopID|\"Outgroup\"|\"Excluded\" for each line in the ped file.\n\n" if @ARGV < 3;

open POP, "<", $ARGV[0];
open PED, "<", $ARGV[1];

my $control = <POP>;
chomp $control;
my @pop = split /\s+/, $control;
close POP;

my @ped;
my @exc;
my @id;
for (0 .. @pop - 1) {
	my $a = <PED>;
	if ($pop[$_] eq "Excluded") {
		unshift @exc, $_;
		next;
	}
	chomp $a;
	my @b = split /\t/, $a;
	shift @b;
	my $c = shift @b;
	push @id, $c;
	splice @b, 0, 4;
	push @ped, \@b;
}
close PED;

if (@exc) {
	splice @pop, $exc[$_], 1 for (0 .. @exc - 1);
}

my %base = ("A A","A", "C C","C", "G G","G", "T T","T", "A G","R", "C T","Y", "A C","M", "G T","K", "G C","S", "A T","W", "0 0","N");
my $ntax = @id;
my $nchar;
my @nucleotide;

for my $j (0 .. @{$ped[0]} - 1) { # SNP (column) index
	my %alle; # overall allele count
	my %oual; # outgroup allele count
	my %inal; # ingroup allele count
	for my $i (0 .. @ped-1) { # Sample (row) index
		my @a = split / /, $ped[$i][$j];
		++$alle{$a[0]};
		++$alle{$a[1]};
		if ($pop[$i] eq "Outgroup") {
			++$oual{$a[0]};
			++$oual{$a[1]};
		} else {
			++$inal{$a[0]};
			++$inal{$a[1]};
		}
	}
	if ($alle{"0"} or keys %alle != 2 or keys %oual != 1 or keys %inal == 1) {
		next;
	}
	++$nchar;
	for my $i (0 .. @ped-1) {
		$nucleotide[$i] .= $base{$ped[$i][$j]};
	}
}

open OUT, ">", $ARGV[2];
print OUT "#NEXUS
BEGIN DATA;
Dimensions ntax=$ntax nchar=$nchar;
Format datatype=DNA gap=- missing=?;
Matrix\n";
print OUT $id[$_], "\t", $nucleotide[$_], "\n" for (0 .. @id-1);
print OUT ";\nEND;\n";
close OUT;
