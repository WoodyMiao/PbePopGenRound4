#!/usr/bin/perl
use strict;
use warnings;
use threads;

die "\nThis program calculates number of segregating sites matrix from nexus input. Auther: Woody\n
Usage: $0 <input.nex> <output_prefix>\n
This program reads all sequences into RAM, then counts all pairs simultaneously, so make sure there is enough RAM before running.\n\n" if @ARGV < 2;

open I, "<", $ARGV[0];
open O1, ">", "$ARGV[1].txt";
open O2, ">", "$ARGV[1].nex";

my @id;
my %seq;
while (<I>) {
	last if /^matrix/i;
}
while (<I>) {
	last if /^;/;
	chomp;
	my @a = split /\t/;
	push @id, $a[0];
	$seq{$a[0]} = $a[1];
}
close I;
warn "Read sequences complete!\n";

my %threads;
for my $i1 (0 .. @id-2) {
	for my $i2 ($i1+1 .. @id-1) {
		my $pair = "$id[$i1]_$id[$i2]";
		warn "Comparing $id[$i1] and $id[$i2] ...\n";
		$threads{$pair} = threads->new(\&segregating_sites, $seq{$id[$i1]}, $seq{$id[$i2]});
	}
}
my %distance;
$distance{$_} = $threads{$_}->join for keys %threads;
warn "Compare sequences complete!\n";

print O1 "#SamplePair\t#0identical\t#1identical\t#2identical\t#BothNotN\t#SegregatingSites\n"; 
foreach (sort keys %distance) {
	print O1 $_;
	foreach (@{$distance{$_}}) {
		print O1 "\t$_";
	}
	print O1 "\n";
}
close O1;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
print O2 "#NEXUS\n
[Number of segregating sites matrix calculated by segregating_sites.pl, saved at $now]\n
Begin taxa;
\tDimensions ntax=", @id+0, ";
\tTaxLabels\n\t\t", join("\n\t\t", @id), "\n\t\t;\nEnd;\n\n";

print O2 "Begin distances;
\tFormat triangle=lower labels nodiagonal;
\tMatrix\n";
for my $i1 (0 .. @id-1) {
	print O2 $id[$i1];
	for my $i2 (0 .. @id-1) {
		my $pair = "$id[$i2]_$id[$i1]";
		print O2 "\t", $distance{$pair}[4] if $distance{$pair};
	}
	print O2 "\n";
}
print O2 "\t;\nEnd;\n\n";

print O2 "Begin paup;
\tSet autoclose=yes increase=auto warntree=no warnreset=no;
\tSet criterion=distance;
\tNj brlens=yes treefile=$ARGV[1].tre;
\tQuit;
End;\n";
close O2;
system "paup $ARGV[1].nex";

sub segregating_sites {
	my ($seq1, $seq2) = @_;
	my $len1 = length $seq1;
	my $len2 = length $seq2;
	die "Different sequence length!" if $len1 != $len2;
	my $two = 0; # both alleles are same.
	my $one = 0; # only one allele is same.
	my $zero = 0; # both alleles are different.
	my $total; # total informative bp;
	foreach my $i (0 .. $len1 - 1) {
		my $s1 = substr $seq1, $i, 1;
		my $s2 = substr $seq2, $i, 1;
		if (($s1 eq "N") or ($s2 eq "N")) {
			next;
		} else {
			++$total;
			if ($s1 eq $s2) {
				++$two;
			} elsif ((($s1 eq "A") and ($s2 =~ /[WMR]/))
			or (($s1 eq "C") and ($s2 =~ /[SMY]/))
			or (($s1 eq "G") and ($s2 =~ /[SKR]/))
			or (($s1 eq "T") and ($s2 =~ /[WKY]/))
			or (($s1 eq "W") and ($s2 =~ /[ATMKRY]/))
			or (($s1 eq "S") and ($s2 =~ /[CGMKRY]/))
			or (($s1 eq "M") and ($s2 =~ /[ACWSRY]/))
			or (($s1 eq "K") and ($s2 =~ /[GTWSRY]/))
			or (($s1 eq "R") and ($s2 =~ /[AGWSMK]/))
			or (($s1 eq "Y") and ($s2 =~ /[CTWSMK]/))) {
				++$one;
			} else {
				++$zero;
			}
		}
	}
	my $segregating = $zero + $one/2;
	return [$zero, $one, $two, $total, $segregating];
}
