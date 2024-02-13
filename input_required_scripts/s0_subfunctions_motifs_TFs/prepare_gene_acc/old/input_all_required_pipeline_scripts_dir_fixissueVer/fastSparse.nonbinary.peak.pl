#!/usr/bin/perl
use strict;
use warnings;
use Sort::Naturally;

my %hash;
my $its = 0;

open F, $ARGV[0] or die;
while(<F>){
	$its++;
	if(($its % 1000000)==0){
		print STDERR " - read $its reads ... \n";
	}
	chomp;
	my @col = split("\t",$_);
	$col[0] = 'chr' . $col[0];
	if($col[0] =~ /B73V/){
		$col[0] =~ s/chrB73V4_ctg/chrB73V4ctg/g;
	}
	my $id = @col[8];
	if(exists $hash{$col[3]}{$id}){
		$hash{$col[3]}{$id}++;
	}else{
		$hash{$col[3]}{$id} = 1;
	}
}
close F;

my @cells = nsort keys %hash;
for (my $i = 0; $i < @cells; $i++){
	my @sites = keys %{$hash{$cells[$i]}};
	for (my $j = 0; $j < @sites; $j++){
		print "$sites[$j]\t$cells[$i]\t$hash{$cells[$i]}{$sites[$j]}\n";
	}
}
