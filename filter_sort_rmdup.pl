#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

die "\nUsage: $0 <in.geneious.bam>\n\n" if @ARGV != 1;

my $bn = basename $ARGV[0];
$bn =~ s/.geneious.bam$//;

open I, "-|", "samtools-0.1.7 view -h $ARGV[0]";
open O, "|-", "samtools-0.1.7 view -Shb - >$bn.filter.bam";

my ($num_in, $num_out);

while (<I>) {
	if (/^\@/) {
		if (/^\@HD/) {
			print O "\@HD\tVN:1.3\n";
		} elsif (/^\@SQ/) {
			print O;
		} elsif (/^\@RG/) {
			next;
		} 
	} else {
		++$num_in;
		my @a = split /\t/;
		if (($a[1] == 83 or $a[1] == 99 or $a[1] == 147 or $a[1] == 163) and ($a[6] eq "=")) {
				++$num_out;
				my @b = split / /, $a[0];
				shift @a;
				pop @a;
				print O "$b[0]\t", join("\t", @a), "\n";
		}
	}
}
warn "$bn\tInput: $num_in reads\tOutput: $num_out reads\n";
close I;
close O;

warn "sorting $bn ...\n";
system "samtools-0.1.7 sort $bn.filter.bam $bn.sort";

warn "rmduping $bn ...\n";
system "samtools-0.1.7 rmdup $bn.sort.bam $bn.rmdup.bam";

warn "$bn was done!\n";
