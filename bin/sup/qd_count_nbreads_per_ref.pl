#!/usr/bin/perl

use strict;
use warnings;

my $sam = shift @ARGV;

my %h_ref; # nb reads per ref

open(SAM, $sam) or die "cannot open file >$sam<\n";
while(my $line=<SAM> )
{
	if ($line !~ /^@/ )
	{
		undef my @a_line;
		@a_line=split(/\t/, $line);
		$h_ref{$a_line[2]}++;
	}
}
close(SAM);

my $out = "NbReadsPerTargets.csv";
open( OUT, " > $out " );
foreach my $ref (sort keys %h_ref)
{
	print OUT $ref."\t".$h_ref{$ref}."\n";
}
close(OUT);
