#!/usr/bin/perl
use strict;
use warnings;

# This program makes pseudo PBE reference genome using gzipped fasta inputs. The reference allele is the major allele among all inputs.

open LIST, "-|", "ls -1 /share/users/miaolin/6.Pbe_genomics/a1.fasta_CpG_gene_masked/PBEP*.auto.fa.gz";
open O, "|-", "gzip -9c >pseudo_pbe_ref.fa.gz";

my @file;
while (<LIST>) {
	open my $i, "-|", "zcat $_";
	push @file, $i;
}
close LIST;

my %base_allele = ("A", ["A","A"], "C", ["C","C"], "G", ["G","G"], "T", ["T","T"], "R", ["A","G"], "Y", ["C","T"], "M", ["A","C"], "K", ["G","T"], "S", ["G","C"], "W", ["A","T"], "N", ["N","N"]);
while (readline $file[0]) {
	my $chr0 = $_;
	my $seq0 = readline $file[0];
	chomp $seq0;
	my $len = length $seq0;
	my @seq;
	push @seq, $seq0;
	for (1 .. @file-1) {
		my $chr1 = readline $file[$_];
		die "Error: inconsistent chromosome order among inputs.\n" if $chr1 ne $chr0;
		my $seq1 = readline $file[$_];
		push @seq, $seq1;
	}
	my $ref_seq;
	for my $coo (0 .. $len-1) {
		my %base;
		for (@seq) {
			++$base{$base_allele{substr($_, $coo, 1)}[0]};
			++$base{$base_allele{substr($_, $coo, 1)}[1]};
		}
		delete $base{N} if $base{N};
		my $ref_alle = "N";
		my $num_alle = 0;
		for (keys %base) {
			if ($base{$_} > $num_alle) {
				$ref_alle = $_;
				$num_alle = $base{$_};
			}
		}
		my @max_alle;
		for (keys %base) { 
			if ($base{$_} == $num_alle) {
				push @max_alle, $_;
			}
		}
		if (@max_alle and @max_alle > 1) {
			my $r = rand;
			if ($r < 1/@max_alle) {
				$ref_alle = $max_alle[0];
			} elsif ($r >= 1/@max_alle and $r < 2/@max_alle) {
				$ref_alle = $max_alle[1];
			} elsif ($r >= 2/@max_alle and $r < 3/@max_alle) {
				$ref_alle = $max_alle[2];
			} elsif ($r >= 3/@max_alle and $r < 4/@max_alle) {
				$ref_alle = $max_alle[3];
			} else {
				warn;
			}
		}
		$ref_seq .= $ref_alle;
	}
	print O $chr0, $ref_seq, "\n";
}
close O;
