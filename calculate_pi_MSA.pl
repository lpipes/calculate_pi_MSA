#!/usr/bin/perl
use strict;
use warnings;

my %sequences = ();
my %msa = ();
my $length = 0;
my $number_of_sequences = 0;
if ( scalar(@ARGV) < 1){
	print "Usage: calculate_pi_MSA.pl MSA.fasta\n";
	print "Supply a multiple sequence alignment fasta file as the first argument and the script outputs pi\n";
	print "Script also outputs a file pi_per_site.txt which contains the pi per site\n";
	exit(1);
}
my %per_site = ();
open(MSA,$ARGV[0]) || die("cannot open file!");
while(<MSA>){
	my $line = $_;
	chomp($line);
	if ( $line !~ /^\>/ ){
		$sequences{$line}++;
		my @spl = split(//,$line);
		for(my $i=0; $i<scalar(@spl); $i++){
			push(@{$msa{$line}},$spl[$i]);
			if ( $spl[$i] ne '-' and $spl[$i] ne 'n' and $spl[$i] ne 'N' ){
				$per_site{$i}{$spl[$i]}++;
			}
		}
		$length = scalar(@spl);
	}else{
		$number_of_sequences++;
	}
}
close(MSA);

my @key_order;
foreach my $key (sort keys %msa){
	push(@key_order,$key);
}
my $sum=0;
my @pi_per_site;
for(my $i=0; $i<$length; $i++){
	$pi_per_site[$i]=0;
}
for(my $i=0; $i<scalar(@key_order); $i++){
	for(my $j=$i+1; $j<scalar(@key_order); $j++){
		my $differences=0;
		for(my $k=0; $k<$length; $k++){
			my $pi_per_site_total = 0;
			foreach my $key (keys %{$per_site{$k}} ){
				$pi_per_site_total += $per_site{$k}{$key};
			}
			if ( @{$msa{$key_order[$i]}}[$k] ne @{$msa{$key_order[$j]}}[$k] and @{$msa{$key_order[$i]}}[$k] ne '-' and @{$msa{$key_order[$j]}}[$k] ne '-' and @{$msa{$key_order[$i]}}[$k] ne 'N' and @{$msa{$key_order[$j]}}[$k] ne 'N'){
				$differences++;
				$pi_per_site[$k] += ( $per_site{$k}{ @{$msa{$key_order[$i]}}[$k] } / $pi_per_site_total ) * ( $per_site{$k}{ @{$msa{$key_order[$j]}}[$k] } / $pi_per_site_total);
			}
		}
		$sum += ($sequences{$key_order[$i]}/$number_of_sequences) * ($sequences{$key_order[$j]}/$number_of_sequences) * ($differences/$length);
	}
}
my @pi_per_site_sum;
open(PI_PER_SITE,">","pi_per_site.txt") || die("cannot open file!");
for(my $i=0; $i<$length; $i++){
	my $pi_per_site_total = 0;
	foreach my $key (keys %{$per_site{$i}} ){
		$pi_per_site_total += $per_site{$i}{$key};
	}
	$pi_per_site_sum[$i] = ( $pi_per_site_total / ($pi_per_site_total-1) )*$pi_per_site[$i];
	print PI_PER_SITE "$i\t$pi_per_site_sum[$i]\n";
}
close(PI_PER_SITE);
#print "SUM is $sum\n";
print "Number of sequences: $number_of_sequences\n";
my $pi = ($number_of_sequences/($number_of_sequences-1))*$sum;
print "PI is $pi\n";

