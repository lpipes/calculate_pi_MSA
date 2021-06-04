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
	exit(1);
}
open(MSA,$ARGV[0]) || die("cannot open file!");
while(<MSA>){
	my $line = $_;
	chomp($line);
	if ( $line !~ /^\>/ ){
		$sequences{$line}++;
		my @spl = split(//,$line);
		for(my $i=0; $i<scalar(@spl); $i++){
			push(@{$msa{$line}},$spl[$i]);
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
for(my $i=0; $i<scalar(@key_order); $i++){
	for(my $j=$i+1; $j<scalar(@key_order); $j++){
		my $differences=0;
		for(my $k=0; $k<$length; $k++){
			if ( @{$msa{$key_order[$i]}}[$k] ne @{$msa{$key_order[$j]}}[$k] and @{$msa{$key_order[$i]}}[$k] ne '-' and @{$msa{$key_order[$j]}}[$k] ne '-' and @{$msa{$key_order[$i]}}[$k] ne 'N' and @{$msa{$key_order[$j]}}[$k] ne 'N'){
				#print "difference at $k between @{$msa{$key_order[$i]}}[$k] and @{$msa{$key_order[$j]}}[$k]\n";
				$differences++;
			}
		}
		$sum += ($sequences{$key_order[$i]}/$number_of_sequences) * ($sequences{$key_order[$j]}/$number_of_sequences) * ($differences/$length);
	}
}
#print "SUM is $sum\n";
print "Number of sequences: $number_of_sequences\n";
my $pi = ($number_of_sequences/($number_of_sequences-1))*$sum;
print "PI is $pi\n";
