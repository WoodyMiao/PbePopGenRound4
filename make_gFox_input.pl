#!/usr/bin/perl
use strict;
use warnings;

my @list = (
"../../a1.fasta_CpG_gene_masked/FCAP0072.auto.fa.gz",
"../../a1.fasta_CpG_gene_masked/PVIP0002.auto.fa.gz",
"../../a1.fasta_CpG_gene_masked/PVIP0012.auto.fa.gz",
"../../a1.fasta_CpG_gene_masked/PBEP0025.auto.fa.gz",
"../../a1.fasta_CpG_gene_masked/PBEP0064.auto.fa.gz",
"../../a1.fasta_CpG_gene_masked/PBEP0065.auto.fa.gz",
"../../a1.fasta_CpG_gene_masked/PBEP0066.auto.fa.gz",
"../../a1.fasta_CpG_gene_masked/PBEP0067.auto.fa.gz",
"../../a1.fasta_CpG_gene_masked/PBEP0068.auto.fa.gz",
"../../a1.fasta_CpG_gene_masked/PBEP0069.auto.fa.gz");

my %chr = (
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
);

open OUT, ">", "gFox.input";

my %seq;
for (@list) {
	/(\w{8}).auto.fa.gz$/;
	my $id = $1;
	open my $in, "-|", "zcat $_";
	while (<$in>) {
		chomp;
		s/^>//;
		my $s = <$in>;
		chomp $s;
		$seq{$id}{$_} = $s;
	}
	close $in;
}
my @id = sort keys %seq;
my %loci;
my $num_loci;
for my $u (sort keys %chr) {
	for (my $i=0; $i < $chr{$u}; $i +=10) {
		my $w = substr($seq{FCAP0072}{$u}, $i, 1000);
		next if (length $w) < 1000;
		next if $w =~ /^N{10}/;
		my $n = ($w =~ tr/N//);
		next if $n > 200;
		$w = substr($seq{PVIP0002}{$u}, $i, 1000);
		next if $w =~ /^N{10}/;
		$n = ($w =~ tr/N//);
		next if $n > 200;
		$w = substr($seq{PVIP0012}{$u}, $i, 1000);
		next if $w =~ /^N{10}/;
		$n = ($w =~ tr/N//);
		next if $n > 200;
		$w = substr($seq{PBEP0025}{$u}, $i, 1000);
		next if $w =~ /^N{10}/;
		$n = ($w =~ tr/N//);
		next if $n > 200;
		$w = substr($seq{PBEP0068}{$u}, $i, 1000);
		next if $w =~ /^N{10}/;
		$n = ($w =~ tr/N//);
		next if $n > 200;
		my @subseq;
		for my $v (@id) {
			push @subseq, substr($seq{$v}{$u}, $i, 1000);
		}
		push @{$loci{$u}}, \@subseq;
		$i += 50990;
	}
	my $w = @{$loci{$u}};
	warn "$u\t$w loci\n";
	$num_loci += $w;
}
warn "Pick loci complete!\n";

my $num_sample = @id;
print OUT "$num_loci\n\n";
for my $u (sort keys %loci) {
	my $v;
	for my $w (@{$loci{$u}}) {
		++$v;	
		print OUT "${u}_$v\t$num_sample\t1000\n";
		for (0 .. $num_sample-1) {
			print OUT "$id[$_]\t${$w}[$_]\n";
		}
		print OUT "\n";
	}
}
close OUT;
warn "Done!\n";
