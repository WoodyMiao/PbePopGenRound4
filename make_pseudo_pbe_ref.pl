#!/usr/bin/perl
use strict;
use warnings;

die "\nThis program makes pseudo reference genome using gzipped fasta inputs. The reference allele is the major allele among all inputs.
\nAuther: Woody\n
Usage: $0 <in.fa.gz.list> <out.fa.gz>\n\n" if @ARGV < 2;

open LIST, "<", $ARGV[0];
open O, "|-", "gzip -9c >$ARGV[1]";

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
		$ref_seq .= $ref_alle;
	}
	print O $chr0, $ref_seq, "\n";
}
close O;
