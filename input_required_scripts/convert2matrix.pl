#!/usr/bin/perl
use strict;
use warnings;

my %tss;
my %tss_save;

open F, $ARGV[0] or die;
while(<F>){
	chomp;
	my @col = split("\t",$_);
	my $id = join("\t",@col[0..2]);
	$tss_save{$col[3]} = $id;
	my $str = $col[1];
	my $end = $col[2];
	for (my $i = $str; $i <= $end; $i++){
		if($i >= $col[5] && $i <= $col[6]){
			$tss{$col[3]}{$i}++;
		}
	}
}
close F;

my @gene = keys %tss;
for (my $i = 0; $i < @gene; $i++){
	my @sites = keys %{$tss{$gene[$i]}};
	my @coord = split("\t",$tss_save{$gene[$i]});
	if(exists $tss{$gene[$i]}{$coord[1]}){
		print "$tss{$gene[$i]}{$coord[1]}";
	}else{
		print "0";
	}
	for (my $j = $coord[1]+1; $j <= $coord[2]; $j++){
		if(exists $tss{$gene[$i]}{$j}){
			print "\t$tss{$gene[$i]}{$j}";
		}else{
			print "\t0";
		}
	}
	print "\n";
}
