#!/usr/bin/perl
use strict;
use warnings;

die "\nThis program makes DaDi SNP data file using gzipped fasta files.\n\nAuther: Woody\n
Usage: $0 <list (column1: sampleID, column2: file, column3: \"Reference\"|\"Outgroup\"|PopID)> <output>\n\n" if @ARGV < 2;

open LIST, "<", $ARGV[0];
open O, ">", $ARGV[1];

my $max_miss = 7; # max number of missing individuals
my $C1;
my $C2;
my @id;
my @po;
my %pop;
my %id_pop;
my %file;
while (<LIST>) {
	chomp;
	my @a = split /\s+/;
	open $file{$a[0]}, "-|", "zcat $a[1]";
	if ($a[2] eq "Reference") {
		$C1 = $a[0];
	} elsif ($a[2] eq "Outgroup") {
		$C2 = $a[0];
	} else {
		push @id, $a[0];
		push @po, $a[2];
		$pop{$a[2]} = 1;
	}
}
close LIST;
my %base_allele = ("A", ["A","A"], "C", ["C","C"], "G", ["G","G"], "T", ["T","T"], "R", ["A","G"], "Y", ["C","T"], "M", ["A","C"], "K", ["G","T"], "S", ["G","C"], "W", ["A","T"], "N", ["N","N"]);
my %chr_gt;
print O "$C1\t$C2\tAllele1\t", join("\t", sort keys %pop), "\tAllele2\t", join("\t", sort keys %pop), "\tChrID\tOneBasedCoordinate\n";
while (readline $file{$C1}) {
	chomp;
	s/^>//;
	last unless /^chr\w\d$/;
	my $chr0 = $_;
	my $seq0 = readline $file{$C1};
	chomp $seq0;
	my $len = length $seq0;
	my @seq;
	push @seq, $seq0;
	readline $file{$C2};
	$seq0 = readline $file{$C2};
	chomp $seq0;
	push @seq, $seq0;
	for (@id) {
		my $chr1 = readline $file{$_};
		chomp $chr1;
		warn if $chr1 ne ">$chr0";
		my $seq1 = readline $file{$_};
		push @seq, $seq1;
	}
	for my $coo (0 .. $len-1) {
		my @base;
		my %base;
		for (2 .. @seq-1) {
			my @b = @{$base_allele{substr($seq[$_], $coo, 1)}};
			push @base, \@b;
			++$base{$b[0]};
			++$base{$b[1]};
		}
		if ($base{N}) {
			if ($base{N} > 2*$max_miss or keys %base != 3) {
				next;
			} else {
				delete $base{N};
			}
		} else {
			next if keys %base != 2;
		}
		my ($allele0, $allele1) = keys %base;
		my %pop_allele;
		for (0 .. @base-1) {
			++$pop_allele{$po[$_]}{$base[$_][0]};
			++$pop_allele{$po[$_]}{$base[$_][1]};
		}
		my @num_allele0;
		my @num_allele1;
		for (sort keys %pop) {
			$pop_allele{$_}{$allele0} = 0 unless $pop_allele{$_}{$allele0};
			$pop_allele{$_}{$allele1} = 0 unless $pop_allele{$_}{$allele1};
			push @num_allele0, $pop_allele{$_}{$allele0};
			push @num_allele1, $pop_allele{$_}{$allele1};
		}
		my $column1 = substr($seq[0], $coo-1, 3);
		$column1 =~ s/N/-/g;
		my $column2 = substr($seq[1], $coo-1, 3);
		$column2 =~ s/N/-/g;
		print O "$column1\t$column2\t$allele0\t";
		print O join("\t", @num_allele0);
		print O "\t$allele1\t";
		print O join("\t", @num_allele1);
		print O "\t$chr0\t", $coo+1, "\n";
	}
}
close O;
