#!/usr/bin/perl
use strict;
use warnings;

while (<>) {
	next if !/^\d+ \d+ \d+ \d+ \d+/;
	my @a = split / /;
	shift @a;
	pop @a;
	my $t;
	$t += $_ for @a;
	print "sum = $t\n";
	last;
}
