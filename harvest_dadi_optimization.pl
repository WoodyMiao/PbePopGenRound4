#!/usr/bin/perl
use strict;
use warnings;

for (@ARGV) {
	my $i = $_;
	open I, "<", $i;
	my $finish;
	while (<I>) {
		if (/^Finshed/) {
			$finish = 1;
			last;
		}
	}
	next if !$finish;
	$/ = "]";
	my $a = <I>;
	chomp $a;
	$/ = "\n";
	<I>;
	$a =~ s/^Best-fit parameters: \[ +//;
	my @p = split /\s+/, $a;
	$a = <I>;
	chomp $a;
	$a =~ s/^Maximum log composite likelihood: //;
	push @p, $a;
	$a = <I>;
	chomp $a;
	$a =~ s/^Optimal value of theta: //;
	push @p, $a;
	print $i, "\t", join("\t", @p), "\n";
	close I;
}
