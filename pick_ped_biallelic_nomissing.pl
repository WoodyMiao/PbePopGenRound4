#!/usr/bin/perl
use strict;
use warnings;

die "This program picks no-missing biallelic SNPs from PED format.\nAuther: Woody
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
	for (@gt) {
		my @a = split / /, ${$_}[$i];
		++$allele{"$a[0]"};
		++$allele{"$a[1]"};
	}
	if ($allele{"0"}) {
		next;
	} else {
		if (keys %allele == 2) {
			push @index, $i;
			print MAPO $snp[$i];
		}
	}
}
close MAPO;

for my $i (0 .. @id-1) {
	print PEDO join("\t", @{$id[$i]});
	for (@index) {
		print PEDO "\t$gt[$i][$_]";
	}
	print PEDO "\n";
}
close PEDO;

my $n2 = @index;
warn "Pick $n2 no-missing biallelic SNPs out of $n.\n";
#system "cp $ARGV[0].pedind $ARGV[1].pedind";
#system "ln $ARGV[1].map $ARGV[1].pedsnp";
#system "p-link --file $ARGV[1] --out $ARGV[1] --make-bed";
