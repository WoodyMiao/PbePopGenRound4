#!/usr/bin/perl
use strict;
use warnings;

die "This program picks no-missing SNPs from PED format.\nAuther: Woody
Usage: $0 <input prefix> <output prefix>\n" if @ARGV < 2;

open PEDI, "<", "$ARGV[0].ped";
open MAPI, "<", "$ARGV[0].map";

open PEDO, ">", "$ARGV[1].ped";
open MAPO, ">", "$ARGV[1].map";

my @gt;
my @id;
while (<PEDI>) {
	chomp;
	my @a = split /\t/;
	my @b = splice @a, 0, 6;
	push @gt, \@a;
	push @id, \@b;
}
close PEDI;

my @snp;
while (<MAPI>) {
	push @snp, $_;
}
close MAPI;

my $n = @{$gt[0]};
die "Error: inconsistent number of SNPs!\n" if $n != @snp;

my @index;
for my $i (0 .. $n-1) {
	my %allele;
	++$allele{"${$_}[$i]"} for @gt;
	next if $allele{"0 0"};
	push @index, $i;
	print MAPO $snp[$i];
}
close MAPO;

for my $i (0 .. @id-1) {
	print PEDO join("\t", @{$id[$i]});
	print PEDO "\t$gt[$i][$_]" for @index;
	print PEDO "\n";
}
close PEDO;

my $n2 = @index;
warn "Pick $n2 no-missing SNPs out of $n.\n";
