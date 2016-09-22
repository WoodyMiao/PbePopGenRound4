#!/usr/bin/perl
use strict;
use warnings;

die "This program generate fasta using reference and BSNP output.\nAuther: Woody
Usage: $0 <reference.fa> <SNP.gz> <nSNP> <output.fa.gz> <minimum Phred quality> [min coverage] [max coverage]\n" if @ARGV < 5;

open REF, "<", $ARGV[0];
open SNP, "-|", "zcat $ARGV[1]";
open NSN, "<", $ARGV[2];
open OUT, "|-", "gzip -9c >$ARGV[3]";
my $pSNP = 1 - 10**(-$ARGV[4]/10); # Min SNP posterior probability 
my $pNSN = $ARGV[4]/2 + 33; # Min ASCII value of ProbSNP in nSNP file 
my $minc = $ARGV[5] if $ARGV[5]; # Min coverage
my $maxc = $ARGV[6] if $ARGV[6]; # Max coverage

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

<SNP>;
while (<SNP>) {
	my @a = split /\t/;
	if ($a[1] >= $len{$a[0]}) {
		warn "$ARGV[1] $a[0] $a[1] >= $len{$a[0]}\n";
		next;
	}
	next if substr($ref{$a[0]}, $a[1], 1) eq "N";
	if ($minc and $maxc) {
		my $cover = $a[29]+$a[30];
		next if ($cover < $minc) or ($cover > $maxc);
	}
	foreach (9 .. 18) {
		if ($a[$_] > $pSNP) {
       			substr $out{$a[0]}, $a[1], 1, $a[4];
			last;
		}
	}
}
close SNP;

<NSN>;
while (<NSN>) {
	my @a = split / +/;
	my @b = split //, $a[4];
	my @c = split //, $a[3];
	foreach (0 .. $a[2]-1) {
		if ($minc and $maxc) {
			my $cover = ord($c[$_]) - 33;
			next if ($cover < $minc) or ($cover > $maxc);
		}
		if (ord($b[$_]) > $pNSN) {
			my $coor = $a[1] + $_; # 0 based coordinate
			if ($coor >= $len{$a[0]}) {
				warn "$ARGV[2] $a[0] $coor >= $len{$a[0]}\n";
				next;
			}
			substr $out{$a[0]}, $coor, 1, substr($ref{$a[0]}, $coor, 1);
		}
	}
}
close NSN;

print OUT ">$_\n$out{$_}\n" foreach sort keys %out;
close OUT;
