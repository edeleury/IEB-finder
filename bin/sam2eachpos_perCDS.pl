#!/usr/bin/perl 


use strict;
use warnings; 
#use POSIX qw(strftime);
use Statistics::Descriptive;
use Getopt::Long;

## ==================================================
## Parametres 
## ==================================================
my ($refLg_exonPos, $sam_local_file, $output_dir, $help);
my $ini = GetOptions("i=s" => \$refLg_exonPos,
	                   "s=s" => \$sam_local_file,
										 "o=s" => \$output_dir,
                     "help"=> \$help );

if ( $help || !$refLg_exonPos || !$sam_local_file  || !$output_dir )
{
	print "USAGE : perl sam2eachpos_perCDS.pl [-h] -i <CDSinfos> -s <SAM> -o <dir_output>\n\n" ;

	print "arguments:\n";
	print "  -i\t  CDS information file (step1 output)\n";
	print "  -s\t  SAM file, output of mapping genomic reads on CDS sequences (local alignments; step2)\n";
	print "  -o\t  output directory\n";
	print "\noptional arguments:\n";
	print "  -h\t  show the help message\n";

	exit 0;
}
# default values 

                    
## ==================================================
## MAIN
## ==================================================

# -- recov. of CDS information 
my $rh_ref = &recup_trLgth_exonPos($refLg_exonPos);

# -- split the SAM output file into one file per CDS in a dedicated directory 
my $rh_file = &split_SAMfile_perCDS($rh_ref,$sam_local_file,$output_dir);


my $stats_output = $output_dir."/../CDS_STATS_.csv";
open(STATS," > $stats_output ");
print STATS "#ACC\tLGTH\tCOV_MEAN\tCOV_MED\tnbBE_MEAN\tnbBE_MED\tLGTH_COVERED\tCOV_MEAN2\tCOV_MED2\tnbBE_MEAN2\tnbBE_MED2\n";

# -- for each CDS, convert SAM file to EACHPOS file
foreach my $acc (keys %{ $rh_file } ) 
{
	my $sam = $rh_file->{$acc};
	my $acc_lgth = $rh_ref->{$acc}->{"LGTH"};
	my $acc_exon_pos = $rh_ref->{$acc}->{"EXONS"};
	my $str_ontr = &sam2eachpos_oneCDS($acc, $acc_lgth, $acc_exon_pos, $sam);
	print STATS $str_ontr."\n";
}
close(STATS);




## #########################################################################
## ###################### SUBROUTINES ######################################
## #########################################################################

# ==========================================================================
# sub recup_trLgth_exonPos
#
# This function recovers the length of the CDS sequences, and the exon positions on each CDS if known.
# ==========================================================================
sub recup_trLgth_exonPos
{
 	my ($ref_infos) = @_ ;
	my %h_ref; 
	
	open(CSV, $ref_infos) or die "cannot open file >$ref_infos<!\n" ; 
	while(my $line=<CSV>)
	{
		chomp $line;
		if ( $line !~ /ACCESSION/ )
		{
			undef my @a_line;
			@a_line = split(/\t/, $line);
			my $acc = $a_line[0];
			my $acc_lg = $a_line[1];
			my $exon_pos = "UNK";
			if ( $a_line[2] ) {  $exon_pos = $a_line[2]; }
	
			$h_ref{$acc} = {"LGTH"=>$acc_lg, "EXONS"=>$exon_pos};
		}
	}
	close(CSV);

	return (\%h_ref);
}

# ==========================================================================
# sub split_SAMfile_perCDS
#
# This function cuts the SAM file by CDS sequences. CAUTION, the input SAM file must be sorted. 
# ==========================================================================
sub split_SAMfile_perCDS
{
 	my ($rh_ref, $sorted_sam, $output_dir) = @_ ; 
	my %h_file;
	
	my $acc_previous = ""; my $i=0;
	my $out="";
	open(SORTED_SAM, $sorted_sam) or die "cannot open file >>$sorted_sam<< !\n" ;
	while(my $line=<SORTED_SAM>)
	{
		if ( $line !~/^\@/)
		{
			undef my @a_line;
			@a_line = split(/\t/, $line);
			my $ref_acc = $a_line[2]; # CDS on which the read is mapped
			if ( ($ref_acc eq "*") || (! exists($rh_ref->{$ref_acc}) ) ) { next; }

			if ($ref_acc ne $acc_previous )
			{
				if ( $acc_previous ne "") 
				{
					close(OUT);
				}
				my $out_dir = $output_dir."/".$ref_acc;
				if ( ! (-d $out_dir) ) 
				{
					mkdir($out_dir,0755) or die "cannot create dir >>".$out_dir."<<\n";
				}
				$out = $out_dir."/".$ref_acc.".sam";
				$h_file{$ref_acc} = $out;
				open(OUT, " > $out ");
			}
			print OUT $line; $i++;
			$acc_previous = $ref_acc;
		}
	}
	close(OUT);
	close(SORTED_SAM);


	return (\%h_file);
}


# ================================================================
# sub sam2eachpos_oneCDS
#
# This function converts a SAM file into an EACHPOS file. 
# ================================================================
sub sam2eachpos_oneCDS
{
	my ($acc, $acc_lgth, $acc_exon_pos, $sam) = @_;
	
	# -- ini du hash avec toutes les positions du CDS
	my %h_pos; my $rh_pos=\%h_pos;
	my $c=0;
	while($c < $acc_lgth)
	{
		$c++;
		$h_pos{$c} = {'EE'=>0, 'cov'=>0, 'nbBE'=>0};
	}
	# -- si positions des exons connue
	if ( ($acc_exon_pos ne "") && ($acc_exon_pos ne "UNK"))
	{
		$acc_exon_pos =~s/;$//;
		undef my @a_exons;
		@a_exons = split(/;/, $acc_exon_pos);
		foreach my $exon (@a_exons)
		{
			$exon =~/^(\d+)-(\d+)$/;
			$rh_pos->{$1}->{'EE'}=1;
			$rh_pos->{$2}->{'EE'}=1;
		}
	}

	# -- lecture du SAM
	open(SAM, $sam) or die "cannot open SAM file >$sam< !\n" ; 
	while ( my $line=<SAM>) 
	{
		if ( $line !~ /^@/ ) 
		{
			chomp $line;
			undef my @a_line ; 
			@a_line = split(/\t/, $line);
	
			my $pos_beg_mapping = $a_line[3] ; 
			# -- calculation of the end position of the read mapping on the CDS using the CIGAR code 
			my $cigar = $a_line[5] ; 
			my @a_occur_M = ( $cigar =~ /(\d+)M/ ) ; # match or mismatch
			my @a_occur_D = ( $cigar =~ /(\d+)D/ ) ; # deletion
			my @a_occur_M_D = ( @a_occur_M, @a_occur_D) ; 
			my $length_align = 0 ; 
			foreach my $pb ( @a_occur_M_D ) 
			{
				$length_align += $pb ; 
			}
			my $pos_end_mapping = $pos_beg_mapping + $length_align -1 ;

			# -- start position
			$rh_pos->{$pos_beg_mapping}->{'nbBE'}++;
			# -- end position
			$rh_pos->{$pos_end_mapping}->{'nbBE'}++;

			# -- couverage at each position for the read # TODO: must be improved considering the CIGAR format 
			# a read is considered as covering the CDS sequence from its start position to its end position 
			my $pos = $pos_beg_mapping; # de 1 à n de l'alignement 
 			while ( $pos <= $pos_end_mapping ) 
			{	
				$rh_pos->{$pos}->{'cov'}++;
				$pos++; 
			}
			
		}
	}
	close(SAM); 

	# -- output : EACHPOS file 
	if ( $acc_exon_pos eq "")
	{
		$acc_exon_pos = "UNK";
	}
	my $out = $sam;
	$out =~ s/([^\/]+)$/eachpos/;
	open(OUT, " > $out ");
	print OUT "REF\tPOS\tEE\tCOV\tnbBE\n";
	print OUT "#".$acc." LGTH=".$acc_lgth." EXONS=".$acc_exon_pos."\n";
	my @a_cov; my @a_nbBE; 
	my @a_cov_onlybasecovered; my @a_nbBE_onlybasecovered;
	my $nbBases_covered=0;
	foreach my $pos (sort {$a<=>$b} keys %h_pos)
	{
		print OUT $acc."\t".$pos."\t".$rh_pos->{$pos}->{'EE'}."\t".$rh_pos->{$pos}->{'cov'}."\t".$rh_pos->{$pos}->{'nbBE'}."\n" ;
		
		# -- pour calcul mean et median de cov et nbBE pour le CDS
		push(@a_cov, $rh_pos->{$pos}->{'cov'});
		push(@a_nbBE, $rh_pos->{$pos}->{'nbBE'});
		if ( $rh_pos->{$pos}->{'cov'} > 0 )
		{
			$nbBases_covered++;
			push(@a_cov_onlybasecovered, $rh_pos->{$pos}->{'cov'});
			push(@a_nbBE_onlybasecovered, $rh_pos->{$pos}->{'nbBE'});
		}
		
	}
	my $cov_mean = &calcul_mean_from_array(@a_cov);
	my $cov_med = &calcul_median_from_array(@a_cov);
	my $nbBE_mean =  &calcul_mean_from_array(@a_nbBE);
	my $nbBE_med = &calcul_median_from_array(@a_nbBE);
	my $cov_mean_2 = &calcul_mean_from_array(@a_cov_onlybasecovered);
	my $cov_med_2 = &calcul_median_from_array(@a_cov_onlybasecovered);
	my $nbBE_mean_2 =  &calcul_mean_from_array(@a_nbBE_onlybasecovered);
	my $nbBE_med_2 = &calcul_median_from_array(@a_nbBE_onlybasecovered);

	print OUT "#".$acc." LGTH=".$acc_lgth." COV_MEAN=".$cov_mean." COV_MED=".$cov_med." nbBE_MEAN=".$nbBE_mean." nbBE_MED=".$nbBE_med." LGTH_COVERED=".$nbBases_covered." COV_MEAN2=".$cov_mean." COV_MED2=".$cov_med." nbBE_MEAN2=".$nbBE_mean." nbBE_MED2=".$nbBE_med ;
	my $str_onetr= $acc."\t".$acc_lgth."\t".$cov_mean."\t".$cov_med."\t".$nbBE_mean."\t".$nbBE_med."\t".$nbBases_covered."\t".$cov_mean_2."\t".$cov_med_2."\t".$nbBE_mean_2."\t".$nbBE_med_2 ;
	close(OUT);

	# -- on supprime le fichier SAM 
	my $cmd = "rm ".$sam ;
	system($cmd);
	
	return $str_onetr;
}

# ================================================================
# function calcul_median_from_array
# ================================================================
sub calcul_median_from_array
{
 	my (@a_values) = @_ ;
	my $median;

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@a_values);
	$median = sprintf("%.2f",$stat->median() );

	return $median;
}


# ================================================================
# calcul_mean
# Cette fonction calcule la moyenne
# ================================================================
sub calcul_mean_from_array
{
 	my (@a_values) = @_ ;
	my $mean;

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@a_values);
	$mean = sprintf("%.2f", $stat->mean() );

	return $mean;
}





# by emeline.deleury@inra.fr
# date : 2016-06-30 
# date last modified : 2019-04-29

# ===================================================================================
# script : qd_02_sam2eachpos_perCDS.pl

# Ce script travaille sur les CDS présents dans un fichier fasta. Il découpe le fichier 
# SAM SORTED passé en entrée, en un fichier SAM par CDS. Puis il parse pour chaque CDS,
# le fichier SAM associé en fichier EACHPOS (pos, cov, nbBE). Il parse le fichier SAM
# pour calculer à chaque position du CDS, la couverture (cov) et le nombre de reads qui
# commencent ou finissent de mapper à cette position (nbBE). 
# Si les positions des exons sont connues sur le CDS, alors on les renseigne en mettant 
# un "1" dans la colonne "EE" du fichier EACHPOS ("1" quand extrémité d'exon). 

### 
# DETAIL CALCUL DES POSITIONS DE DEBUT ET DE FIN DE MAPPING D UN READ SUR LA REFERENCE 
# La position de debut de mapping du read sur la ref/transcript est donnée dans la colonne 2 du SAM.
# La position de fin du mapping du read sur la ref est calculée avec le code CIGAR qui decrit l'alignement du read sur la reference. On compte le nb de match ou mismatch (M) et le nb de deletion (D) et on les ajoute à la position de début de mapping du read sur la reference.
#
# Exemple M : 	REF :  XXXXCCCGTCTATAGCATCGGTAAXXXX
#													 ||||||||||| ||||||||
#						  	READ:	     CCCGTCTATAGAATCGGTAA         CIGAR : 20M
# 
# Exemple D:	REF :  XXXXCCCGTCTAGCATACGCATGAXXXX
#								         ||||||||-|||||||||||
#				    	READ:	     CCCGTCTA*CATACGCATGA         CIGAR : 8M 1D 11M 
# 
# Exemple I:	REF :  XXXXCCCGTCTAG**CATACGCATGAXXXX
#								         |||||||||  |||||||||||
#				    	READ:	     CCCGTCTAGATCATACGCATGA		 CIGAR : 9M 2I 11M
# 
# Exemple N:	REF :  AGCTAGCATCGTGTCGCCCGTCTAGCATACGCATGATCGACTGTCAGCTAGTCAGACTAGTCGATC
#										           |||||||||...................................|||||||||
#					    READ:	           GTGTAACCC...................................GACTAGTCG
# 																		 CIGAR : 9M 35N 9M 
#				=> pas pris en compte dans le script
# 				TODO : checker si on a des formats CIGAR avec des N 
# 				cut -f 6 <_clean.sam> | sort | uniq | grep -e N 
#
# Exemple S:	REF :  XXXXXXATCGTGTCGCCCGTCTAGCATACXXXXXXXXXX
#									    ||||||||||||||||
#					    READ:	       gggGTGTCGCCC-TCTAGCgggg	 CIGAR : 3S 9M 1D 6M 4S
# Note : La colonne 4 du sam (position de début du mapping) tient déjà compte des 'Soft Clipped' (S). La position 
# de debut sur la ref commence bien au 'G' dans l'exemple ci-dessous, donc on n'a pas à se préoccuper de cette info. 

# TODO : 
	# format SAM : 
	# 		- checker les flags de la col 2 du SAM 
	# 		- meilleure prise en compte du format CIGAR (verif si cas avec des N ) et checker sur plusieurs ex
# ===================================================================================





