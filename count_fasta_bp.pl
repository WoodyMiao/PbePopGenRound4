#!/usr/bin/perl
use strict;
use warnings;
use threads;
use File::Basename;

my @list;
open $list[0], "-|", "ls -1 ../a1.fasta_CpG_gene_masked/*.fa.gz";
open $list[1], "-|", "ls -1 ../a2.fasta_onlyCDS_unmasked/*.fa.gz";

open O, ">", "count_fasta_bp.txt";

my %length = (
	chrA1	=>	240380223,
	chrA2	=>	168638799,
	chrA3	=>	140925898,
	chrB1	=>	206538554,
	chrB2	=>	152998503,
	chrB3	=>	148068395,
	chrB4	=>	142431058,
	chrC1	=>	222198629,
	chrC2	=>	159252932,
	chrD1	=>	115468741,
	chrD2	=>	88096124,
	chrD3	=>	94101111,
	chrD4	=>	94492513,
	chrE1	=>	61081816,
	chrE2	=>	61960243,
	chrE3	=>	41224383,
	chrF1	=>	70119229,
	chrF2	=>	83953389,
	chrX	=>	127282370,
);

my %sex = (
	FCAP0072	=>	"2",
	PBEP0005	=>	"1",
	PBEP0008	=>	"2",
	PBEP0009	=>	"1",
	PBEP0010	=>	"1",
	PBEP0013	=>	"2",
	PBEP0014	=>	"2",
	PBEP0023	=>	"1",
	PBEP0025	=>	"1",
	PBEP0026	=>	"2",
	PBEP0028	=>	"1",
	PBEP0036	=>	"1",
	PBEP0038	=>	"2",
	PBEP0039	=>	"1",
	PBEP0064	=>	"2",
	PBEP0065	=>	"1",
	PBEP0066	=>	"2",
	PBEP0067	=>	"1",
	PBEP0068	=>	"1",
	PBEP0069	=>	"1",
	PVIP0012	=>	"1",
);

print O "#File\tSex\t%NonN\t%GC\t%Ti\t%Tv\tTi/Tv\t%SNP\t#A\t#M\t#R\t#W\t#C\t#S\t#Y\t#G\t#K\t#T\t#N\t#ACGT\t#Ti\t#Tv\t#NonN\t#Total\n";
foreach my $l (@list) {
	my %file;
	while (<$l>) {
		chomp;
		$file{$_} = openfile($_);
	}
	close $l;
	print STDERR "Read file list complete!\n";

	my %thread;
	foreach (sort keys %file) {
		$thread{$_} = threads->new(\&countfasta, $file{$_});
	}

	my %count;
	foreach (sort keys %thread) {
		$count{$_} = $thread{$_}->join;
		print STDERR "Count $_ complete!\n";
	}

	foreach my $f (sort keys %count) {
		my $g = basename $f;
		$g =~ s/.fa.gz$//;
		$g =~ /^(\w+)/;
		print O "$g\t$sex{$1}\t", join("\t", @{$count{$f}}[16 .. 21], @{$count{$f}}[0 .. 15]), "\n";
	}
}
close O;

sub openfile {
	my $filename = shift;
	my $infile;
	open $infile, "-|", "zcat $filename" or die "Error opening $filename: $!\n";
	return $infile;
}

sub countfasta {
	my $in = shift;
	my %chrcount;
	while (<$in>) {
		chomp;
		s/>//;
		my $seq = <$in>;
		chomp $seq;
		my $A = ($seq =~ tr/A//); #AA
		my $M = ($seq =~ tr/M//); #AC
		my $R = ($seq =~ tr/R//); #AG
		my $W = ($seq =~ tr/W//); #AT
		my $C = ($seq =~ tr/C//); #CC
		my $S = ($seq =~ tr/S//); #CG
		my $Y = ($seq =~ tr/Y//); #CT
		my $G = ($seq =~ tr/G//); #GG
		my $K = ($seq =~ tr/K//); #GT
		my $T = ($seq =~ tr/T//); #TT
		my $N = ($seq =~ tr/N//); #N
		my $ho = $A + $C + $G + $T; #homozygote
		my $ti = $R + $Y; #transition
		my $tv = $M + $W + $S + $K; #transversion
		my $nn = $ho + $ti + $tv; # not N
		my $total = $N + $nn;
		warn "$_\t$total" if $total != $length{$_};
#		my $pnn = sprintf("%.2f", $nn / $total);
#		my $gcc = sprintf("%.2f", (($M + $R + $Y + $K) / 2 + $C + $G + $S) / $nn); #GC content
#		my $pti = sprintf("%.4f", $ti / $nn); # %transition
#		my $ptv = sprintf("%.4f", $tv / $nn); # %transversion
#		my $titv = sprintf("%.2f", $ti / $tv);
#		my $psnp = sprintf("%.4f", ($ti + $tv) / $nn); # %SNP
		$chrcount{$_} = [$A, $M, $R, $W, $C, $S, $Y, $G, $K, $T, $N, $ho, $ti, $tv, $nn, $total];
	}
	close $in;
	my @sum;
	foreach my $c (keys %chrcount) {
		foreach my $n (0 .. 15) {
			$sum[$n] += $chrcount{$c}[$n];
		}
	}
	$sum[16] = sprintf("%.6f", $sum[14] / $sum[15]);
	$sum[17] = sprintf("%.6f", (($sum[1] + $sum[2] + $sum[6] + $sum[8]) / 2 + $sum[4] + $sum[5] + $sum[7]) / $sum[14]);
	$sum[18] = sprintf("%.6f", $sum[12] / $sum[14]);
	$sum[19] = sprintf("%.6f", $sum[13] / $sum[14]);
	if ($sum[13]) {
		$sum[20] = sprintf("%.6f", $sum[12] / $sum[13]);
	} else {
		$sum[20] = "N/A";
	}
	$sum[21] = sprintf("%.6f", ($sum[12] + $sum[13]) / $sum[14]);;
	return \@sum;
}
