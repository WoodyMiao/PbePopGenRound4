#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

die "This program converts PED to ARLEQUIN project file.\nAuther: Woody
Usage: $0 <in.ped> <out.arp> <auto|chrx> [N (pick 1 snp per N snps)]\n" if @ARGV < 3;

open PED, "<", $ARGV[0];
open ARP, ">", $ARGV[1];

my %ped;
my %sex;
while (<PED>) {
	chomp;
	my @a = split /\t/;
	shift @a;
	my $id = shift @a;
	shift @a;
	shift @a;
	my $sex = shift @a;
	my $pop = shift @a;
	$sex{$id} = $sex;
	$ped{$pop}{$id} = \@a;
}	
close PED;

my $title = basename $ARGV[1];
$title =~ s/.arp$//;
my $nbs = keys %ped;
print ARP "[Profile]
	Title = \"$title\"
	DataType = DNA
	NbSamples = $nbs
	GameticPhase = 0
	GenotypicData = 1
	LocusSeparator = NONE

[Data]
[[Samples]]\n";

#for my $i (sort keys %ped) {
for my $i ("Borneo", "Malay", "Cambodia", "Amur") { # Order for my special case.
	my $ss = keys %{$ped{$i}};
	print ARP "
	SampleName = \"$i\"
	SampleSize = $ss
	SampleData = {\n";
	for my $j (sort keys %{$ped{$i}}) {
		my $h0;
		my $h1;
		if ($ARGV[3]) {
			for (0 .. @{$ped{$i}{$j}}-1) {
				next if $_ % $ARGV[3];
				my @a = split / /, $ped{$i}{$j}[$_];
				$h0 .= $a[0];
				$h1 .= $a[1];
			}
		} else {
			for (@{$ped{$i}{$j}}) {
				my @a = split / /;
				$h0 .= $a[0];
				$h1 .= $a[1];
			}
		}
		if ($ARGV[2] eq "chrx" and $sex{$j} == 1) {
			my $l = length $h0;
			$h1 = "?" x $l;
		}
		print ARP "\t\t$j\t1\t$h0\n\t\t\t\t$h1\n";
	}
	print ARP "\t}\n";
}
close ARP;
