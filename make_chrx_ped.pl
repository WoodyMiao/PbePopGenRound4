#!/usr/bin/perl
use strict;
use warnings;

die "This program makes chrX PED file using gzipped fasta files.\nAuther: Woody
Usage: $0 <list (column1: sample, column2: file, column3: 1=male|2=female, column4: population)> <output prefix>\n" if @ARGV < 2;

open LIST, "<", "$ARGV[0]";
open GT, ">", "$ARGV[1].ped";
open SNP, ">", "$ARGV[1].pedsnp";
open IND, ">", "$ARGV[1].pedind";

my $num_miss = 10; # max number of missing individual for each locus
my @id;
my %sex;
my %pop;
my @seq;
while (<LIST>) {
	chomp;
	my @a = split / /;
	push @id, $a[0];
	$sex{$a[0]} = $a[2];
	$pop{$a[0]} = $a[3];
	open my $in, "-|", "zcat $a[1]";
	<$in>;
	my $s = <$in>;
	push @seq, $s;
	close $in;
}
close LIST;

my %chr_gt;
my %chr_coo;
my $chrXlen = length $seq[0];
for my $coo (0 .. $chrXlen-1) {
	my @base;
	my %base;
	for (@seq) {
		my $b = substr($_, $coo, 1);
		push @base, $b;
		++$base{$b};
	}
	my $nk = keys %base;
	if ($base{N}) {
		next if $base{N} > $num_miss;
		next if $nk == 2;
	} else {
		next if $nk == 1;
	}
	for (0 .. @id-1) {
		$chr_gt{X}[$_] .= $base[$_];
	}
	push @{$chr_coo{X}}, $coo;
}

my %base_gt = ("A","A A", "C","C C", "G","G G", "T","T T", "R","A G", "Y","C T", "M","A C", "K","G T", "S","G C", "W","A T", "N","0 0");
for my $i (0 .. @id-1) {
	print IND "0\t$id[$i]\t0\t0\t$sex{$id[$i]}\t$pop{$id[$i]}\n";
	print GT "0\t$id[$i]\t0\t0\t$sex{$id[$i]}\t$pop{$id[$i]}";
	for my $c (sort keys %chr_gt) {
		my $len = length $chr_gt{$c}[$i];
		for (0 .. $len-1) {
			my $b = substr $chr_gt{$c}[$i], $_, 1;
			print GT "\t$base_gt{$b}";
		}
	}
	print GT "\n";
}
close IND;
close GT;
for my $c (sort keys %chr_coo) {
	for (0 .. @{$chr_coo{$c}}-1) {
		print SNP "$c\tchr${c}snp";
		printf SNP "%08d", 1+$_;
		print SNP "\t0\t$chr_coo{$c}[$_]\n";
	}
}
close SNP;
