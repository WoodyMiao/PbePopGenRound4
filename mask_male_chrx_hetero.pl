#!/usr/bin/perl
use strict;
use warnings;

#my @input = qw/
#PVIP0012.chrx.fa.gz
#PBEP0005.chrx.fa.gz
#PBEP0009.chrx.fa.gz
#PBEP0010.chrx.fa.gz
#PBEP0023.chrx.fa.gz
#PBEP0025.chrx.fa.gz
#PBEP0028.chrx.fa.gz
#PBEP0036.chrx.fa.gz
#PBEP0039.chrx.fa.gz
#PBEP0065.chrx.fa.gz
#PBEP0067.chrx.fa.gz
#PBEP0068.chrx.fa.gz
#PBEP0069.chrx.fa.gz
#/;
my @input = qw/
PBEP0009.chrx.fa.gz
PBEP0065.chrx.fa.gz
/;

for (@input) {
	open my $in, "-|", "zcat $_";
	my $chr = <$in>;
	my $seq = <$in>;
	close $in;
	$seq =~ tr/RYMKSW/N/;
	s/.fa.gz$//;
	open my $out, "|-", "gzip -9c >$_.hm.fa.gz";
	print $out $chr, $seq;
	close $out;
}
