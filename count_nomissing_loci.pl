#!/usr/bin/perl
use strict;
use warnings;

# This program counts number of nucleotide loci without missing data using gzipped fasta files.
# A list of fa.gz should be passed to STDIN.

my @file;
while (<>) {
	chomp;
	open my $i, "-|", "zcat $_";
	push @file, $i;
}

my %miss;
while (readline $file[0]) {
	chomp;
	s/^>chr//;
	my $chr0 = $_;
	my $seq0 = readline $file[0];
	chomp $seq0;
	my $len = length $seq0;
	my @seq;
	push @seq, $seq0;
	for (1 .. @file-1) {
		my $chr1 = readline $file[$_];
		chomp $chr1;
		die "Error: different chromosome order!\n" if $chr1 ne ">chr$chr0";
		my $seq1 = readline $file[$_];
		chomp $seq1;
		push @seq, $seq1;
	}
	for my $coo (0 .. $len-1) {
		my %base;
		for (@seq) {
			my $b = substr($_, $coo, 1);
			++$base{$b};
		}
		if ($base{N}) {
			++$miss{$base{N}};
		} else {
			++$miss{"0"};
		}
	}
}

print "#missing\t#loci\n";
for (sort {$a <=> $b} keys %miss) {
	print "$_\t$miss{$_}\n";
}
