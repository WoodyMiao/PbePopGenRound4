#!/usr/bin/perl

use strict;
use warnings;

# This program generate pseudo PVI reference genome using FCA reference and BSNP output.

open REF, "<", "/bak/archive/projects/LeopardCat/reference/felCat8_gene_masked_auto.fa";
open SNP, "-|", "zcat /bak/archive/projects/LeopardCat/2.unBQSR_BSNP/PVIP0012.SNP.gz";
open NSN, "<", "/bak/archive/projects/LeopardCat/2.unBQSR_BSNP/PVIP0012.nSNP";
open OUT, "|-", "gzip -9c >pseudo_pvi_ref.fa.gz";

my %ref;
my %len;
while (<REF>) {
	chomp;
	s/>//;
	my $seq = <REF>;
	chomp $seq;
	$ref{$_} = $seq;
	$len{$_} = length $seq;
}
close REF;

my %out;
foreach (keys %ref) {
       $out{$_} = "N" x $len{$_};
}       

my %base_allele = ("R", ["A","G"], "Y", ["C","T"], "M", ["A","C"], "K", ["G","T"], "S", ["G","C"], "W", ["A","T"]);
<SNP>;
while (<SNP>) {
	next unless /^chr\w\d/;
	chomp;
	my @a = split /\t/;
	next if $a[1] >= $len{$a[0]};
	next if substr($ref{$a[0]}, $a[1], 1) eq "N";
	if ($a[4] =~ /[ACGT]/) {
		substr $out{$a[0]}, $a[1], 1, $a[4];
	} else {
		my @b = split //, $a[31];
		my %base;
		++$base{$_} for @b;
		if ($base{$base_allele{$a[4]}[0]} > $base{$base_allele{$a[4]}[1]}) {
			substr $out{$a[0]}, $a[1], 1, $base_allele{$a[4]}[0];
		} elsif ($base{$base_allele{$a[4]}[1]} > $base{$base_allele{$a[4]}[0]}) {
			substr $out{$a[0]}, $a[1], 1, $base_allele{$a[4]}[1];
		} else {
			my $r = rand;
			if ($r < 0.5) {
				substr $out{$a[0]}, $a[1], 1, $base_allele{$a[4]}[0];
			} else {
				substr $out{$a[0]}, $a[1], 1, $base_allele{$a[4]}[1];
			}
		}
	}
}
close SNP;

<NSN>;
while (<NSN>) {
	next unless /^chr\w\d/;
	my @a = split / +/;
	my @b = split //, $a[4];
	my @c = split //, $a[3];
	foreach (0 .. $a[2]-1) {
		my $coor = $a[1] + $_; # 0 based coordinate
		next if $coor >= $len{$a[0]};
		substr $out{$a[0]}, $coor, 1, substr($ref{$a[0]}, $coor, 1);
	}
}
close NSN;

print OUT ">$_\n$out{$_}\n" foreach sort keys %out;
close OUT;
