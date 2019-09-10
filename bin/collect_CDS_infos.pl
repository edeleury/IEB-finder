#!/usr/bin/perl


use strict;
use warnings; 
#use POSIX qw(strftime);
use Getopt::Long;


## ==================================================
## Parametres 
## ==================================================
my ($fasta, $gff3, $output, $help);
my $ini = GetOptions("f=s" => \$fasta,
										 "g=s" => \$gff3,
	                   "o=s" => \$output,
                     "help"=> \$help );

if ( $help || !$fasta || !$output )
{
	print "USAGE : perl collect_CDS_infos.pl -f <CDS.fasta> [-g <gff3>] [-h] -o <output.csv> \n\n" ; 
	print "Collect CDS information  : CDS lengths and, if intron-exon structure is known (GFF available), exon positions on CDS sequence. Probably must be adapted for each GFF file (different nomage according to the species). \n\n";

	print "arguments:\n";
	print "  -f\t  fasta file of CDS sequences\n";
	print "  -o\t  output file \n";
	print "\noptional arguments:\n";
	print "  -h\t  show the help message\n";
	print "  -g\t  GFF file with intron-exon gene annotation\n";
	
	
	exit 0;
}
# default values 

## ==================================================
## MAIN
## ==================================================

# -- recup Accession et Longueur des CDS pour lesquels on veut connaitre les IEB 
my $rh_ref = &recup_info_CDS_seq($fasta);

# -- recup positions des exons si evaluation de la prediction des IEB 
if ( $gff3 && ($gff3 ne "") ) 
{
	$rh_ref = &recup_exon_positions_from_gff3($rh_ref, $gff3);
}

# -- write output 
open(OUT, " > $output ");
# entete
print OUT "ACCESSION\tLGTH";
if ( $gff3 && ($gff3 ne "") ) { print OUT "\tEXON_POSITIONS"; }
print OUT "\n";
# one output line per CDS
foreach my $acc ( sort keys %{ $rh_ref } ) 
{
	print OUT $acc."\t".$rh_ref->{$acc}->{"LGTH"};
	# si IE connue: position des exons sur CDS 
	if ( $gff3 && ($gff3 ne "") ) 
	{
		# test si longueur de la seq CDS dans le fasta et la meme que celle de la somme des longueurs des exons
		if (  $rh_ref->{$acc}->{"EXONS_LGTH"} == $rh_ref->{$acc}->{"LGTH"} ) 
		{
			print OUT "\t".$rh_ref->{$acc}->{"EXONS"};
		}
		else
		{
			print OUT "\tUNK";
			print "WARNING : ".$acc." CDS_LGTH=".$rh_ref->{$acc}->{"LGTH"}." DIFF SUM_EXON_LGTH=".$rh_ref->{$acc}->{"EXONS_LGTH"}."!!\n";
		}
	}
	print OUT "\n";

}
close(OUT);



## #########################################################################
## ###################### SUBROUTINES ######################################
## #########################################################################

# ================================================================
# sub recup_info_CDS_seq
#
# This function parses a fasta file and calculates the length of the CDS sequences.
# ================================================================
sub recup_info_CDS_seq
{
 	my ($fasta) = @_ ;
	my %h_ref; 

	my %h_ref_seq; 
	my $acc = ""; 
	open(FA, $fasta) or die "cannot open file >$fasta<!\n" ; 
	while(my $line=<FA>)
	{
		chomp $line; 
		if ($line =~ /^>([^\s]+)/ )
		{
			$acc = $1 ;
			#$acc =~ s/\.\d+$//; # delete version in accession number
			$h_ref_seq{$acc} = "" ;  
		}	
		else
		{
			$h_ref_seq{$acc} .= $line ; 
		}
	}
	close(FA);


	foreach my $acc (keys %h_ref_seq) 
	{
		$h_ref_seq{$acc} =~ s/ //g;
		$h_ref_seq{$acc} =~ s/\s//g;
		my $ref_lgth = length($h_ref_seq{$acc});
		$h_ref{$acc} = {"LGTH"=>$ref_lgth, "EXONS"=>"", "EXONS_LGTH"=>0, "CDSfile"=>1, "GFF3file"=>0};
	}

	return(\%h_ref);
}


# ================================================================
# sub recup_exon_positions_from_gff3
#
# This function parses a gff3 file to retrieve the length of the exons on the genome and transposes this length to the CDS (taking into account the direction) to position the exons on the CDS. 
# ================================================================
sub recup_exon_positions_from_gff3
{
	my ($rh_ref, $gff3) = @_ ;
	my %h_ref_listexon;

	# -- read gff3 to identify exon positions
	open(GFF, $gff3) or die "cannot open GFF3 file >$gff3< !\n"; 
	while(my $line=<GFF>)
	{
		chomp $line; 
		if ($line !~ /^#/ )
		{
			undef my @a_line ; 
			@a_line = split(/\t/, $line); 
			if ($a_line[2] eq "CDS") 
			{
				my $lg = $a_line[4] - $a_line[3] +1; # exon length
				my $sens = $a_line[6]; 
				$a_line[8] =~ /Parent=transcript:(\w+);/; #### WARNING: NOT GENERIC !!!!!!
				my $ref_acc = $1;
				if ( ! exists($h_ref_listexon{$ref_acc}) ) 
				{
					$h_ref_listexon{$ref_acc}="";# ini
				}

				if ( $sens eq "+" ) 
				{ $h_ref_listexon{$ref_acc} .= $lg.";;"; }
				elsif ( $sens eq "-") 
				{ $h_ref_listexon{$ref_acc}  = $lg.";;".$h_ref_listexon{$ref_acc}; }
			}
		}	
	}
	close(GFF);

	# -- position of exons on each CDS sequence
	foreach my $ref_acc ( keys %h_ref_listexon) 
	{
		$h_ref_listexon{$ref_acc} =~ s/;;$// ;
		undef my @a_listexons ; 
		@a_listexons = split(/;;/, $h_ref_listexon{$ref_acc}) ;
		my $exon_lgth_sum=0;
		
		# -- liste des exons codants d'un CDS
		my $h_exon_positions = "";
		my $pos = 0 ; 
		foreach my $exon_lg ( @a_listexons )
		{
			$exon_lgth_sum += $exon_lg;
			# -- exon begin
			$pos+=1; 		
			$h_exon_positions .= $pos."-" ; 
			# -- exon end
			$pos += $exon_lg -1 ;
			$h_exon_positions .= $pos.";" ; 
		}

		if ( exists($rh_ref->{$ref_acc}) )
		{
			$rh_ref->{$ref_acc}->{"EXONS"} = $h_exon_positions;
			$rh_ref->{$ref_acc}->{"GFF3file"} = 1;
			$rh_ref->{$ref_acc}->{"EXONS_LGTH"} = $exon_lgth_sum;
			
		}
		else
		{
			$rh_ref->{$ref_acc} = {"LGTH"=>"", "EXONS"=>$h_exon_positions, "EXONS_LGTH"=>$exon_lgth_sum, "CDSfile"=>0, "GFF3file"=>1};
		}
	}


	return($rh_ref);
}


