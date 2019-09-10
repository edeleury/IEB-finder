#!/usr/bin/perl

use strict;
use warnings; 
#use POSIX qw(strftime);
use Statistics::Descriptive;
use Getopt::Long;


## ==================================================
## Parametres 
## ==================================================
my ($dir_input, $cov_min, $lg_min, $wind_nbBE_median, $wind_cov_mean, $cutoff, $signaldeport, $nbBases_betw2signal, $evaluation, $help);
my $ini = GetOptions("d=s" => \$dir_input,
	                   "c=i" => \$cov_min,
										 "l=i" => \$lg_min,
										 "w1=i"=> \$wind_nbBE_median,
										 "x=i" => \$cutoff,
							       "n=i" => \$signaldeport,
										 "e=i" => \$evaluation,
                     "help"=> \$help );

if ( $help || ( !$dir_input && !$cutoff) )
{
	print "USAGE : IEB_predict.pl [-h] [...] -d <DIR_ALL_COVERED_CDS> -x <cutoff_toDetectSignal>\n\n" ; 

	print "arguments:\n";
	print "  -d\t  Directory with covered CDS to annote (defined in step 3, ie output directory for all CDS)\n";
	print "  -c\t  Minimum coverage to predict IEB (default: 10)\n";
	print "  -l\t  Minimum region length with coverage < c to consider a region not covered (default: 10)\n";
  print "    \t  A region is considered not covered if it has more than x consecutive bases with a coverage below c.\n"; 
	print "  -w1\t  Number of bases to be taken from each side of a position to calculate the median of the nbBE\n"; 
  print "    \t  in this window (default: 5 for a window of 11 bases centered on the position)\n";
	print "  -x\t  Cutoff to detect a signal. A signal is detected at a position if nbBE at that position is x times larger\n";
	print "    \t  than the local median of nbBE.\n";
	print "  -n\t  Number of bases on which the signal is extended on either side of the measured signal (default: 3)\n";
	print "  -e\t  Evaluation of the predicter's performance (-e 1), IE structure of CDS need to be known (default: 0)\n";

	print "\noptional arguments:\n";
	print "  -h\t  show the help message\n";

	exit 0;
}
# default values 
if ( !$cov_min ) 							{ $cov_min = 10; }
if ( !$lg_min ) 							{ $lg_min = 10; }
if ( !$wind_nbBE_median )			{ $wind_nbBE_median = 5; }
if ( !$signaldeport ) 				{ $signaldeport = 3; }
if ( !$evaluation ) 					{ $evaluation = 0; }

my $str_param = "";
$str_param .= "input_dir   : $dir_input\n";
$str_param .= "cov_min     : $cov_min\n";
$str_param .= "lg_min      : $lg_min\n";
$str_param .= "lgW1_median : $wind_nbBE_median\n";
$str_param .= "cutoff      : $cutoff\n";
$str_param .= "near        : $signaldeport\n";
$str_param .= "evaluation  : $evaluation\n";
my $outputfile_suffixeparam = "_c".$cov_min."_x".$cutoff."_n".$signaldeport ;
                    
## ==================================================
## MAIN
## ==================================================

## -- initialization / header of output files

my %h_perf = ("TP"=>0, "FN"=>0, "FP"=>0, "TN"=>0, "SN"=>"", "SP"=>"");
my $rh_perf = \%h_perf;

# output file with the predictions of regions (region +=>IEB, region -=>EXON, no=>no_prediction)
my $output_allCDS_regions = $dir_input."/../regions_allCDS".$outputfile_suffixeparam.".csv";
open(REGION, " > $output_allCDS_regions ");
print REGION "ACC\tNUM\tBEGIN-END\tLGTH\tMEAN_COV\tSIGNAL_NB\tPREDICTION";

# output file with evaluation if known IE structure and evaluation requested 
if ( $evaluation == 1)
{
	my $output_allCDS_eval = $dir_input."/../evaluation_allCDS".$outputfile_suffixeparam.".csv";
	open(EVAL, " > $output_allCDS_eval " );
	print EVAL "#CDS\tLGTH\tLGTH_covered\tCOV_MEAN\tCOV_MEAN_covered\tmed_nbBE\tNbRegCov0\tnbIEB\tnbIEBcov\tTP\tFN\tFP\tTN\tSN\tSP\n";

	print REGION "\tPRED\tEVAL";
}
print REGION "\n";


## -- IEB prediction for each covered CDS found in directory 
my $cmd = "ls -l $dir_input | grep ^d ";
open(LIST, " $cmd | ");
while(my $res =<LIST>)
{
	$res =~ /([^\s]+)$/;
	my $eachpos_onetr = $dir_input."/".$1."/eachpos" ;

	# -- data recovery in 2 hashs (%h_tr, %h_pos) 
	undef my $rh_pos; undef my $rh_tr;
	($rh_pos, $rh_tr ) = &eachpos_file2hash( $eachpos_onetr );
	my $tr_lgth=$rh_tr->{'lgth'};

	# -- condition for predicting IEBs:
	if ( $rh_tr->{"tr_cov_mean_covered"} > $cov_min )
	{
		# -- search for regions of more than *l* consecutive bases with no or very low coverage
		$rh_pos = &zone_lowcov($tr_lgth, $rh_pos, $cov_min, $lg_min);

		# -- calculation of the median of nbBE around each position and signal detection
		$rh_pos = &signal_detection( $tr_lgth, $rh_pos, $wind_nbBE_median, $cutoff );

		# -- extension of the measured signal to the surrounding positions
		$rh_pos = &signal_offset( $rh_pos, $signaldeport, $tr_lgth );

		# -- division of the CDS sequence into regions
		undef my $rh_region;
		$rh_region = &division_into_regions( $rh_pos, $tr_lgth );

		# -- prediction of IE structure of the CDS
		$rh_region = &prediction_IE_struct($rh_region);
	
		# -- Possible evaluation of the IEB predictions
		#    only if the CDS has known intron-exon structure
		my $str_evalIEB="";
		if ( ($evaluation == 1) && ($rh_tr->{"exons_pos"} ne "") && ($rh_tr->{"exons_pos"} ne "UNK") )
		{
			# -- region assessment (TP,FP,FN,TN)
			my ($rh_region, $rh_eval_regions) = &regions_evaluation($rh_region);
			# -- sensitivity (SN) and specificity (SP) calculation
			($str_evalIEB, $rh_perf) = &calcul_SN_SP($rh_eval_regions, $rh_perf ); 
			# output in 'EVAL' file for one CDS (one line per CDS)
			print EVAL $rh_tr->{'acc'}."\t".$tr_lgth."\t".$rh_tr->{'lgth_covered'}."\t".$rh_tr->{'tr_cov_mean'}."\t".$rh_tr->{'tr_cov_mean_covered'}."\t".$rh_tr->{'tr_nbBE_median'}."\t".$str_evalIEB."\n"; 		
		}

		# -- output in 'REGION' file for one CDS (several lines per CDS)
		print REGION "#".$rh_tr->{'acc'}." LGTH=".$tr_lgth." EXONS=".$rh_tr->{'exons_pos'}." COV_MEAN=".$rh_tr->{'tr_cov_mean'}." MEDIAN_nbBE=".$rh_tr->{'tr_nbBE_median'}."\n";
		print REGION &regions_hash2tab_output($rh_tr,$rh_region,$evaluation, 0);

		# -- detailed output file in the directory of each CDS 
		# 	region predicted and possible evaluation were reported at the end of the file
		my $outputfile = $dir_input."/".$rh_tr->{'acc'}."/".$rh_tr->{'acc'}."_eachpos__IEBpred".$outputfile_suffixeparam.".csv";
		my $res = &make_outputfile_onetr( $outputfile, $rh_tr, $rh_pos, $rh_region, $evaluation, $str_evalIEB);

	}
}
close(REGION);


## -- Global assessment of IEB predictions on all covered CDS
if ( $evaluation == 1) 
{
	if ( ( $h_perf{'TP'} + $h_perf{'FN'} ) == 0 ) 
	{ $h_perf{'SN'} = "NA"; }
	else
	{ $h_perf{'SN'} = sprintf("%.3f",$h_perf{'TP'} / ( $h_perf{'TP'} + $h_perf{'FN'} ) ); }

	if ( ( $h_perf{'TP'} + $h_perf{'FP'} )== 0 ) 
	{ $h_perf{'SP'} = "NA"; }
	else
	{ $h_perf{'SP'} = sprintf("%.3f", $h_perf{'TP'} / ( $h_perf{'TP'} + $h_perf{'FP'} ) ); }

	print EVAL "\n\t\t\t\t\t\t\t\t".$h_perf{'TP'}."\t".$h_perf{'FN'}."\t".$h_perf{'FP'}."\t".$h_perf{'TN'}."\t".$h_perf{'SN'}."\t".$h_perf{'SP'}."\n";
	print EVAL "\n\n\n\n## PARAMETERS SCRIPT :\n".$str_param."\n";
	close(EVAL);
}




## #########################################################################
## ###################### SUBROUTINES ######################################
## #########################################################################

# ==========================================================================
# sub eachpos_file2hash
#
# This function transforms an EACHPOS file into 2 hashes, %h_pos with the 
# information at each position, and %h_tr with the information on the transcript.
# It calculates some stats for the transcript (cov_mean, median_nbBE), and 
# verifies the length of the transcript passed in header. 
# FR: Cette fonction transforme un fichier EACHPOS en 2 hashs, %h_pos avec les 
# infos à chaque position, et %h_tr avec les infos sur le transcrit. Elle
# calcule quelques stats pour le transcrit (cov_mean, median_nbBE), et 
# verifie la longueur du transcrit passée en entete. 
# ==========================================================================
sub eachpos_file2hash
{
 	my ($eachpos_onetr) = @_ ;
	my %h_pos; my $rh_pos=\%h_pos;
	my %h_tr; my $rh_tr=\%h_tr;
	my $nb_line = 0;# nb of position lines in EACHPOS file
	my $str_nbBE ="";

	# -- read eachpos file 
	open(EACHPOS, $eachpos_onetr) or die "cannot open file >$eachpos_onetr< ! \n";
	while( my $line=<EACHPOS>)
	{
		chomp $line;
		undef my @a_line; 
		if ( $line =~ /^#([^ ]+) LGTH=(\d+) EXONS=([^ ]+)/) # eg: #ENST00000007735 LGTH=1215 EXONS=UNK
		{
			%h_tr = ( "acc"=>$1, "lgth"=>$2, "lgth_covered"=>0, "exons_pos"=>$3, "line"=>$line, "nbEE"=>0, "sum_cov"=>0, "nbBaseWithCov0"=>0, "nbEECov0"=>0, "tr_cov_mean"=>0, "tr_cov_mean_covered"=>0, "tr_nbBE_median"=>"", "warning"=>"" );	
		}
		elsif ( ($line !~ /^REF/) && ($line !~ /^#/) )
		{
			$nb_line++;
			@a_line = split(/\t/, $line);
			my $acc = $a_line[0]; # transcript accession
			my $pos = $a_line[1]; # position on transcrit 
			my $EE = $a_line[2]; # exon end (1) or not (0)
			my $cov = $a_line[3]; # nb of mapped reads at that position (ie coverage)
			my $nbBE = $a_line[4]; # nb of reads that begin or end their alignment at that position 
		
			$h_pos{$pos} = {"EE"=>$EE , "cov"=>$cov, "nbBE"=>$nbBE, "nbBEmedian"=>"", "nbBE_nbBEmedian"=>"", "pred"=>"", "signal"=>"", "covstatut"=>""} ;
	
			# to calculate tr_cov_mean
			$h_tr{"sum_cov"} += $cov;
			# nb of exon ends 
			if ( $EE == 1 )	{ $h_tr{"nbEE"}++; }
			# nb of not covered bases
			if ( $cov == 0) 
			{ 
				$h_tr{"nbBaseWithCov0"}++; 
				if ($EE == 1) { $h_tr{"nbEECov0"}++; }
			}
			# to calculate the median of the nbBE on the transcript for only the covered bases
			if ( $cov > 0 ) 
			{ 
				$str_nbBE .= $nbBE.";" ; 
				$h_tr{"lgth_covered"}++;
			}

		}
	}
	close(EACHPOS);

	# -- test if the CDS length in the header is the same as the number of counted positions
	if ( $h_tr{"lgth"} != $nb_line )
	{
		print "WARNING: TR=".$h_tr{"acc"}." TR_LGTH_ENTETE: ***".$h_tr{"lgth"}."*** != ".$nb_line." (NB LINES EACHPOS FILE)\n";
		$h_tr{"lgth"} = $nb_line;
		$h_tr{"warning"} = "TR_LGTH: ".$h_tr{"lgth"}."!=".$nb_line ;
	}
	
	# -- calculation of the mean coverage of the CDS
	# -- calculation of the median of the nbBE of the CDS
	if ( $h_tr{"sum_cov"} > 0)
	{
		$h_tr{"tr_cov_mean"} = sprintf("%.1f", $h_tr{"sum_cov"} / $h_tr{"lgth"} ); # arrondi
		$h_tr{"tr_cov_mean_covered"} = sprintf("%.1f", $h_tr{"sum_cov"} / $h_tr{"lgth_covered"} );
		$h_tr{"tr_nbBE_median"} = 0;# ini 
		if ( $str_nbBE ne "")
		{
			$h_tr{"tr_nbBE_median"} = &calcul_median_from_str($str_nbBE);
		}
	}

	return ($rh_pos, $rh_tr);
}


# ==========================================================================
# function calcul_median_from_str
#
# This function calculates the median of a list of values passed as a parameter 
# via a string with the separator";".
# FR: Cette fonction calcule la mediane d'une liste de valeurs passée en parametre
# via une chaine de caracteres avec le separateur ";".
# ==========================================================================
sub calcul_median_from_str
{
 	my ($list) = @_ ;
	my @a_values;
	my $median;

	$list =~ s/;$//;
	my @a_list = split(/;/, $list);
	foreach my $value (@a_list)
	{
		push(@a_values, $value);
	}
	$median = &calcul_median_from_array(\@a_values , 0 );

	return $median;
}

# ==========================================================================
# function calcul_median_from_array
#
# This function calculates the median of a list of values passed as parameters
# via an array. CAUTION when calculating the median of the nbBE around a position,
# if the median value is less than 1, then replace the value with 1 (ie $action=1).
# In this way, when calculating at a given position the ratio nbBE(pos)/median_nbBE(window_around_pos),
# dividing by 1 will not find a signal when the median value of the nbBE is too small.
# And avoid the case of dividing by zero. 
# FR: Cette fonction calcule la mediane d'une liste de valeurs passée en parametre 
# via un array. ATTENTION dans le cas du calcul de la mediane du nbBE autour 
# d'une position, si la valeur de la mediane est inferieure à 1, alors on 
# remplace la valeur par 1. De cette facon, quand on calculera à une position
# donnée le ratio nbBE(pos)/median_nbBE(fenetre_autour_pos), diviser par 1 ne 
# permettra pas de trouver un signaux quand la valeur mediane du nbBE est trop 
# petit. Et evite le cas de diviser par zero. 
# ==========================================================================
sub calcul_median_from_array
{
 	my ($ra_values, $action) = @_; # $action is 0 or 1 (if 1 then a local median value <1 is replaced by 1 )
	my $median;
	my @a_values = @$ra_values; # on dereference le array

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@$ra_values);
	$median = sprintf("%.2f", $stat->median() );

	# CHOICE: if mediane < 1, then this value is replaced by 1 
	# divided by 1, will not create signals  
	if ( ($median < 1) && ($action ==1) ) { $median=1; }	

	return $median;
}


# ==========================================================================
# function zone_lowcov
#
# This function searches for region with low or no coverage of more than *l* bases 
# (default 10) on a transcript. For these areas, the tool will not attempt 
# to predict IEBs (and will not evaluate these regions). This function divides CDS
# into covered and uncovered regions. It calls the "delete_shortzone_withnocov" 
# function to delete very small areas that are not covered. Once done, the function
# assigns to each position of the transcript a covstatut: 0 if not or slightly covered
# (cov < *c*) and 1 if not.
# FR: Cette fonction recherche des zones pas ou peu couvertes de plus de x bases 
# (default 10) sur un transcrit. Pour ces zones, l'outil ne tentera pas de 
# predire les IEB (et n'évaluera pas ces zones). Cette fonction decoupe en 
# zones couvertes et non couverte le transcrit. Elle appelle la fonction 
# "delete_shortzone_withnocov" pour supprimer les tres petites zones non 
# couvertes. Une fois fait, la fonction attribue à chaque position du transcrit
# un covstatut: 0 si pas ou peu couvert (cov < cov_min) et 1 sinon. #### NOTE: 
# les zones couvertes et non couvertes ne sont pas identifiées par une position 
# de debut et une position de fin (on ne conserve pas l'info sous cette forme).
# ==========================================================================
sub zone_lowcov
{
	my ($tr_lgth, $rh_pos, $cov_min, $lg_lowmin) = @_;
	my %h_zone_covstatut; 
	my $rh_zone_covstatut=\%h_zone_covstatut;
	my $z=0;# autoincrement ID
	my $cov_previous = 100;
	my $min_lgth=100000;

	# -- we go through the CDS sequence to know the positions where cov < c
  #    and we cut the CDS sequence into areas (covered and uncovered)
	foreach my $pos ( sort {$a<=>$b} keys %{ $rh_pos } ) 
	{
		my $covstatut=""; 
		my $cov = $rh_pos->{$pos}->{'cov'};
		if ( $cov < $cov_min )	{ $covstatut = 0; }
		else 										{ $covstatut = 1; }
	
		if ( $covstatut != $cov_previous ) 
		{
			if ( $pos == 1 ) { goto NEWR; }

			# -- treatment of the previous zone 
			$rh_zone_covstatut->{$z}->{"end"} = $pos -1 ; 
			$rh_zone_covstatut->{$z}->{"lgth"} = $pos -1 - $rh_zone_covstatut->{$z}->{"beg"} +1; 
			if ( ($rh_zone_covstatut->{$z}->{"covstatut"} == 0) && ($rh_zone_covstatut->{$z}->{"lgth"} < $min_lgth) )
			{
				$min_lgth = $rh_zone_covstatut->{$z}->{"lgth"};
			}

			# -- new zone
			NEWR:$z++; 
			$h_zone_covstatut{$z} = {"beg"=>$pos, "end"=>"", "covstatut"=>$covstatut, "lgth"=>""};

		}

		# -- treatment of the last zone 
		if ( $pos == $tr_lgth )
		{
			$rh_zone_covstatut->{$z}->{"end"} = $pos; 
			$rh_zone_covstatut->{$z}->{"lgth"} = $pos - $rh_zone_covstatut->{$z}->{"beg"} +1;
			if ( ($rh_zone_covstatut->{$z}->{"covstatut"} == 0) && ($rh_zone_covstatut->{$z}->{"lgth"} < $min_lgth) )
			{
				$min_lgth = $rh_zone_covstatut->{$z}->{"lgth"};
			}
		}

		$cov_previous = $covstatut;
	}


	# -- possibly, removal of very short areas not covered (default = 10 bases)
 	if ( $min_lgth < $lg_lowmin ) # at least one area to be deleted 
	{
			$rh_zone_covstatut = &delete_shortzone_withnocov($rh_zone_covstatut, $lg_lowmin, $z);
	}

	# -- we go through the covered and uncovered areas and report the covstatut in the hash %h_pos 
	#    covstatut = 0 if position belongs to an area with low or no coverage of more than 10 bases
	foreach my $z ( sort {$a<=>$b} keys %{ $rh_zone_covstatut } ) 
	{
		my $beg = $rh_zone_covstatut->{$z}->{"beg"};
		my $end = $rh_zone_covstatut->{$z}->{"end"};
		my $covstatut = $rh_zone_covstatut->{$z}->{"covstatut"};
		for ( my $pos=$beg; $pos<=$end; $pos++)
		{
			$rh_pos->{$pos}->{"covstatut"} = $covstatut;
		}
	}

	return ($rh_pos);
}


# ==========================================================================
# function delete_shortzone_withnocov
#
# This function rewrites the hash passed at the input by 'deleting' areas 
# with weakly or no coverage of less than x bp (default 10). Short, 
# poorly/not covered areas are'merged' with neighbouring covered areas, 
# as in the example below. 
# FR: Cette fonction reecrit le hash passé en entree en 'supprimant' les zones peu 
# ou pas couvertes de moins de x pb (defaut 10). Les zones courtes peu/pas 
# couvertes sont 'refusionnées' avec les zones couvertes voisines, comme sur
# l'exemple présenté ci-dessous. 

# -- Example : # ENST00000438665
# input :
#		Z:1	1-17			0	17		
#		Z:2	18-463		1	446
#		Z:3	464-473		0	10
#		Z:4	474-850		1	377
#		Z:5	851-852		0	2		  <----- short area with little coverage
#		Z:6	853-2310	1	1458
#		Z:7	2311-2311	0	1	    <----- short area with little coverage
#		Z:8	2312-3135	1	824		
#
# output :
#		Z:1	1-17			0	17
#		Z:2	18-463		1	446
#		Z:3	464-473		0	10
#		Z:4	474-3135	1	2662
# ==========================================================================
sub delete_shortzone_withnocov
{
	my ($rh_zone, $lgth, $z) = @_; 
	my %h_new_zone; my $rh_new_zone =\%h_new_zone; 
	my $n = 0 ; # key of the new hash

	my $zone1_statut=0;
	my $w = 0 ; 
	while ( $w < $z ) # we go through the old hash
	{ 	
		$w++;
		# if short area not covered or little covered, then merge with the 2 surrounding areas
		if ( ($rh_zone->{$w}->{"covstatut"} == 0) && ($rh_zone->{$w}->{"lgth"} < 10) )
		{
			if ( $w == 1) # first area on the CDS is a short uncovered area
			{
				$zone1_statut = "1" ;
			}
			else
			{
				my $after= $w+1;
				if ( exists($rh_zone->{$after}->{"end"}) ) 
				{ $rh_new_zone->{$n}->{"end"} = $rh_zone->{$after}->{"end"}; }
				else
				{ $rh_new_zone->{$n}->{"end"} = $rh_zone->{$w}->{"end"}; }
					
				$w++;
			}
		}
		else
		{
				$n++;
				$h_new_zone{$n} = {"beg"=>$rh_zone->{$w}->{"beg"}, "end"=>$rh_zone->{$w}->{"end"}, "covstatut"=>$rh_zone->{$w}->{"covstatut"} };
				# case if first area on the CDS is a short uncovered area
				if ( $zone1_statut ==1 )
				{
					$rh_new_zone->{$n}->{'beg'} = 1;
				}
		}	
	}

	return ($rh_new_zone);
}


# ==========================================================================
# function signal_detection
#
# This function searches for signals. It goes through the sequence of the CDS,
# at each position, it calculates for a window of given length (default 11, 
# ie 5 bases upstream and 5 bases downstream of the position) the median of
# the nbBE (nbBE_median) and the ratio "nbBE / nbBE_median". If this ratio is 
# greater than or equal to the cutoff value (parameter x chosen by the user),
# then a signal, noted +, of the presence of an IEB at this position is measured.
# ATTENTION: a signal detected at the ends of the CDS is not considered because,
# by construction, the reads stop mapping to the ends of the CDS (direct mapping 
# on CDS and not on genome).
# FR: Cette fonction recherche les signaux. Elle parcourt la sequence du transcrit,
# à chaque position, elle calcule pour une fenetre de longueur donnée (defaut 11,
# ie 5 bases en amont et 5 bases en aval de la position) la mediane du nbBE 
# (nbBE_median) et le ratio "nbBE / nbBE_median". Si ce ratio est superieure ou
# egale à la valeur cutoff (parametre choisi par l'utilisateur), alors un signal, 
# noté +, de la presence d'un IEB à cette position est détecté/mesuré.
# ATTENTION: un signal détecté aux extremites du CDS n'est pas considéré car,
# par construction, les reads s'arretent de mapper aux extremites du transcrit 
# (mapping direct sur CDS et non sur génome).
# ==========================================================================
sub signal_detection
{
	my ($tr_lgth, $rh_pos, $w_lg, $cutoff ) = @_; 
	my @a_values_in_window;

	my $q=0;
	# -- initialization of the array at the 1st positions
	while($q < $w_lg)
	{
		$q++;
		push(@a_values_in_window, $rh_pos->{$q}->{'nbBE'} );
	}

	foreach my $p ( sort {$a<=>$b} keys %{ $rh_pos } ) 
	{
		# we remove an element at the beginning of the array 
		if ( $p > ($w_lg +1) )
		{
			shift(@a_values_in_window);
		}

		# we add an element at the end of the array (q advance)
		if ( $p < $tr_lgth - $w_lg +1 )
		{
			$q++;
			push(@a_values_in_window, $rh_pos->{$q}->{'nbBE'} );
		}

		# -- calculation, for position p, of the median nbBE of the window 
		#   NOTE if the median is <1, then it is replaced by 1
		#   NOTE the value nbBE at the considered position is included in the calculation of the median
		my $median = &calcul_median_from_array(\@a_values_in_window, 1);
		$rh_pos->{$p}->{'nbBEmedian'} = $median;

		# -- calculation of ratio : nbBE(p) / median_nbBE(Window)
		my $ratio = &calcul_ratio_nbBE_nbBEmedian( $rh_pos->{$p}->{'nbBE'} , $median );
		if ( $ratio != 1000000000 ) 
		{
			$rh_pos->{$p}->{'nbBE_nbBEmedian'} = $ratio;
			$rh_pos->{$p}->{'pred'} = "-";
			# -- signal prediction if the ratio is higher than the cutoff value X 
			#    *** ATTENTION: no signal prediction at the ends of the CDS  ***
			if ( ($ratio >= $cutoff) && ($p!=1) && ($p!=$tr_lgth) )
			{	
				$rh_pos->{$p}->{'pred'} = "+";
			}	
		}
		else
		{
			print "WARNING: division by zero !!! "; exit 0; 
		}
	}

 	return ($rh_pos);
}

# ==========================================================================
# function calcul_ratio_nbBE_nbBEmedian
# ==========================================================================
sub calcul_ratio_nbBE_nbBEmedian
{
	my ($val, $val_wind) = @_;
	
	if ( $val_wind !=0 )
	{
		my $ratio =sprintf("%.1f", $val / $val_wind ); # arrondi
		return $ratio;
	}
	else
	{
		return 1000000000;
	}	
}


# ==========================================================================
# function signal_offset
#
# This function extends a measured signal, noted "+" found at a position 
# on n bases upstream and n bases downstream to take into account a possible 
# offset of the signal (related to the alignment of a small piece of intron 
# on the CDS sequence). The bases that "benefit" from a signal offset are 
# noted "+n" (n for near). Care should be taken not to replace a measured 
# signal "+" with a remote signal "+n", in order not to lose the information 
# of the positions of the measured signals.  
# FR: Cette fonction transpose un signal mesuré, noté '+' trouvé à une position sur
# n bases en amont et n bases en aval pour tenir compte d'un déport possible
# du signal (lié à l'alignement d'un petit bout d'intron sur la sequence
# du transcrit). Les bases qui "beneficient" d'un déport de signal sont notées 
# "+n" (n pour near). On fait attention de ne pas remplacer un signal mesuré
# "+" par un signal deporté "+n", pour ne pas perdre l'information des positions
# des signaux mesurés.  
# ==========================================================================
sub signal_offset
{
	my ($rh_pos, $n_lg, $tr_lgth ) = @_;
	my @a_values_in_window;

	foreach my $pos ( sort {$a<=>$b} keys %{ $rh_pos } ) 
	{
		my $pred = $rh_pos->{$pos}->{'pred'}; 
		if ( $pred eq "+" ) # position with a measured signal 
		{
			# -- the signal is extended to n upstream bases 
			for (my $i = 1; $i <= $n_lg ; $i++)
			{
				my $near_pos = $pos - $i ; if ( $near_pos <= 0) { last; }
				if ( $rh_pos->{$near_pos}->{'pred'} ne "+" )
				{
					$rh_pos->{$near_pos}->{'pred'} = "+n";
				}
			}
			# -- the signal is extended to n downstream bases
			for (my $i = 1; $i <= $n_lg ; $i++)
			{
				my $near_pos = $pos + $i ; if ( $near_pos > $tr_lgth) { last; }
				if ( $rh_pos->{$near_pos}->{'pred'} ne "+" )
				{
					$rh_pos->{$near_pos}->{'pred'} = "+n";
				}
			}

		}
	}

 	return ($rh_pos);
}

# ==========================================================================
# function division_into_regions
#
# This function scans the CDS sequence and cuts it into different regions, 
# i.e:
# - uncovered regions for which the tool cannot predict 
# - short regions + (prediction of an IEB) in which a signal (or several) 
#   has been measured 
# - regions - (exon prediction) in which no signal has been detected.
# This breakdown into regions first takes into account the coverage of the
# regions.  For each region, it recovers the number of exon extremities present
# in the region. Please note that exons ends at the ends of the CDS
# are not counted. 
#
# FR: Cette fonction parcourt le transcrit et le découpe en differentes regions, 
# i.e :
# - des regions non couvertes pour lesquelles l'outil ne peut pas predire
# - des regions + (prediction d'un IEB) courtes dans lesquelles un signal 
#   (ou plusieurs) a été mesuré 
# - des regions - (prediction d'un exon) dans lesquelles aucun signal n'a été
#   détecté
# Ce decoupage en regions tient compte en premier de la couverture des regions.
# Pour chaque region, il recupere le nombre d'extremites d'exons présentes 
# dans la region. Attention les EE aux extremites du transcrit ne sont pas 
# comptabilisées. 
# 
# -- Example pour HA_10003, découpage en regions : 
#	Z:1		1-19::no_prediction		0.0
#	Z:2		20-23::+							2273.7
#	Z:3		24-166::-							3555.9
#	Z:4		167-173::+						3111.1
#	Z:5		174-683::-						4458.6
#	Z:6		684-692::+						1901.3
# [...]
# ==========================================================================
sub division_into_regions
{
	my ($rh_pos, $tr_lgth) = @_ ;
	my %h_region; 
	my $rh_region = \%h_region;

	my $pred_previous="";
	my $covstatus_previous=100; # initial value different from 0 and 1
	my $r=0; # autoincremente ID
	
	foreach my $pos ( sort {$a<=>$b} keys %{ $rh_pos } ) 
	{
		my $EE=$rh_pos->{$pos}->{'EE'};
		my $cov=$rh_pos->{$pos}->{'cov'};
		my $pred=$rh_pos->{$pos}->{'pred'};
		if ( $pred=~/^\+/ ) { $pred = "+"; }
		my $covstatut=$rh_pos->{$pos}->{'covstatut'};

		# -- if we "pass" in an UNcovered region ( 1 => 0  ou 100 => 0 )
		# -- ---------------------------------------------
		if ( ($covstatut==0) && ($covstatus_previous !=0) )
		{
			# if first region of the CDS
			if ( $r == 0) { goto PASSE; }

			# -- previous region treatment
			$rh_region->{$r}->{"end"} = $pos-1;

			# -- initialization of a new region
			PASSE:$r++;
			$h_region{$r} = {"beg"=>$pos, "end"=>"", "pred"=>"no_prediction", "nbEE"=>0, "prediction"=>"", "nbsignaux"=>0, "sum_cov"=>0, "lgth"=>0, "cov_mean"=>"", "eval"=>""};
		}

		# -- if we "stay" in an UNcovered region ( 0 => 0 ) 
		# -- -------------------------------------------
		#elsif ( ($covstatut==0) && ($covstatus_previous ==0) ) { }

		# -- if we "pass" in a covered region ( 0 => 1  ou 100 => 1 ) 
		# -- -------------------------------------------
		elsif ( ($covstatut==1) && ($covstatus_previous !=1) )
		{
			# if first region of the CDS
			if ( $r == 0) { goto LALA; }

			# -- previous region treatment
			$rh_region->{$r}->{"end"} = $pos-1;

			# -- initialization of a new region
			LALA:$r++;
			$h_region{$r} = {"beg"=>$pos, "end"=>"", "pred"=>$pred, "nbEE"=>0, "prediction"=>"", "nbsignaux"=>0, "sum_cov"=>0, "lgth"=>0, "cov_mean"=>"", "eval"=>""};
		}


		# -- if we stay in a covered region
		# -- ------------------------------------------
		elsif (  ($covstatut==1) && ($covstatus_previous ==1) )
		{
			# --  if change of region type (  + => - ou - => + )
			if ( (($pred eq "-") && ($pred_previous ne "-")) || (($pred eq "+") && ($pred_previous ne "+"))  )
			{
				# -- previous region treatment
				$rh_region->{$r}->{"end"} = $pos-1;

				# -- initialization of a new region
				$r++;
				$h_region{$r} = {"beg"=>$pos, "end"=>"", "pred"=>$pred, "nbEE"=>0, "prediction"=>"", "nbsignaux"=>0, "sum_cov"=>0, "lgth"=>0, "cov_mean"=>"", "eval"=>""}; #-# &ini_region_info($pos,$pred);
			}
		}

		# -- treatment of the last region of the CDS
		# -- ------------------------------------------
		if ( $pos == $tr_lgth)
		{
			$rh_region->{$r}->{"end"}=$pos;
		}

		# -- per region, calculation of : nb exon ends, cumulative coverage, nb signals 
		# the exon ends at the CDS end are not counted 
		if ( ($pos==1) || ($pos==$tr_lgth) ) { $EE = 0; }
		$rh_region->{$r}->{"nbEE"} += $EE;
		$rh_region->{$r}->{"sum_cov"} += $cov;
		if ( $rh_pos->{$pos}->{'pred'} eq "+" ) { $rh_region->{$r}->{"nbsignaux"}++; } # pas de "+" aux extremités du tr 

		# -- we continue to go through the CDS
		$covstatus_previous = $covstatut;
		if ($covstatut==0) { $pred_previous = 0; }
		else { $pred_previous = $pred; }
		
	}

	# -- calculation of region lgth and region mean coverage
	foreach my $r (sort keys %{ $rh_region } )
	{
		$rh_region->{$r}->{"lgth"}=  $rh_region->{$r}->{"end"} - $rh_region->{$r}->{"beg"} +1 ;
		$rh_region->{$r}->{"cov_mean"}= $rh_region->{$r}->{"sum_cov"} / $rh_region->{$r}->{"lgth"} ;

		$rh_region->{$r}->{"cov_mean"} =sprintf("%.2f",$rh_region->{$r}->{"cov_mean"}); 
	}
	
	if ( $r > 0 ) { return $rh_region; }
	else { return 0; }
}


# ==========================================================================
# function prediction_IE_struct
#
# This function more explicitly names the regions + and regions -. 
# FR: Cette fonction nomme de facon plus explicite les region + et -. 
# ==========================================================================
sub prediction_IE_struct
{
	my ( $rh_region ) = @_; 

	foreach my $z ( sort {$a<=>$b} keys %{ $rh_region } )
	{
		# region non suffisamment couverte
		if ( $rh_region->{$z}->{"pred"} eq "no_prediction" )
		{
			$rh_region->{$z}->{"prediction"} = "no_prediction (region with low coverage)" ;
		}
		
		# region + == IEB
		elsif ($rh_region->{$z}->{"pred"} eq "+")
		{
			$rh_region->{$z}->{"prediction"} = "IEB prediction (".$rh_region->{$z}->{"nbsignaux"}." signals)" ;
		}
		# region - == exon
		elsif ($rh_region->{$z}->{"pred"} eq "-")
		{
			if ( $rh_region->{$z}->{"lgth"} > 20 ) 
			{
				$rh_region->{$z}->{"prediction"} = "EXON" ;
			}
			else 
			{
				$rh_region->{$z}->{"prediction"} = "exon";
			}
		}
	}

	return ( $rh_region );
}

# ==========================================================================
# function regions_evaluation
# 
# This function is only called if the evaluation of predictions is requested, 
# i.e. only if the I/E structure of the transcript is known.
# 
# For all regions + or - of the CDS, this function evaluates whether the region is:
# - a true positive (TP), if region + containing at least 1 EE (true IEB well predicted)
# - a false negative (FN), if region - containing at least 1 EE (true unpredicted IEB)
# - a false positive (FP), if region + while 0 EE (prediction while no true IEB)
# - a true negative (TN), if region - and 0 EE (true exon prediction)
# Note: Regions not covered are not evaluated. 
#
# FR: Cette fonction n'est appelée que si l'évaluation des predictions est demandée,
# donc uniquement si la structure I/E du transcrit est connue.
#
# Pour toutes les regions + ou - du transcrit, cette fonction evalue si la region est :
# - un vrai positif (TP), si region + contenant au moins 1 EE (vrai IEB bien predit)
# - un faux negatif (FN), si region - contenant au moins 1 EE (vrai IEB non predit)
# - un faux positif (FP), si region + alors que 0 EE (prediction alors que pas de vrai IEB)
# - un vrai negatif (TN), si region - et 0 EE (vrai prediction d'un exon)
# Note: les regions non couvertes ne sont pas évaluées. 
# Rappel : on a "annulé" la prediction de signal mesuré aux extrémites des CDS, 
# donc pas non plus de deport de ces signaux mesurés non pris en compte, donc pas 
# de regions + aux extrémités du CDS (lié à un signal à la position 1 ou à la pos
# fin du CDS). Une region + vers l'extremité d'un CDS sera le fait d'un signal 
# mesuré non loin de l'extrémité (mais pas à l'exact extrémité) du CDS.
# Rappel : les EE aux extrémités du CDS ne sont pas non plus comptees.  

# NOTE : une IEB correspond à 2 EE consecutives. De ce fait, la somme des TP et 
# des FN peut etre superieure au nb de vrais IEB. 

# RAPPEL schématique :
# 
# 												   Region
#											Region+				Region-
#				Realite			
# 				vrai EE>=1		TP						FN
# 				vrai EE=0 		FP						TN
# ==========================================================================
sub regions_evaluation
{
	my ($rh_region) = @_ ;
	my %h_eval_regions = ( "TP"=>0, "FN"=>0, "FP"=>0, "TN"=>0, "no_prediction"=>0, "nbIEB"=>0, "nbIEBcovered"=>0 );
	
	foreach my $r (sort {$a<=>$b} keys %{ $rh_region } )
	{
		my $region_eval = "";		
		$h_eval_regions{"nbIEB"}+= $rh_region->{$r}->{"nbEE"};# nbIEB = nbEE/2 
						# Reminder: EE at the ends of the CDS have not been recorded
		$h_eval_regions{"nbIEBcovered"}+= $rh_region->{$r}->{"nbEE"}; # if the region is not covered, we will deduct after

		# --  case of an uncovered region => no prediction
		if ( $rh_region->{$r}->{"pred"} eq "no_prediction" ) 
		{ 
			$region_eval = "no_prediction";
			$h_eval_regions{"nbIEBcovered"}-= $rh_region->{$r}->{"nbEE"};
		}

		# -- cas of a region + 
		elsif ( $rh_region->{$r}->{"pred"} eq "+")
		{
				if ( $rh_region->{$r}->{"nbEE"} >0 ) { $region_eval = "TP"; } # (*)
				else																 { $region_eval = "FP"; }
		}
		## (*) REFLEXION: si 1 signal, soit doute, soit superposition signaux

		# -- cas of a region - 
		elsif ( $rh_region->{$r}->{"pred"} eq "-")
		{
			if ( $rh_region->{$r}->{"nbEE"} ==0 ) { $region_eval = "TN"; }
			else																  { $region_eval = "FN"; } # (**)
		}
		## (**) REFLEXION: si region_lg < 20pb et entouré de 2 regions + alors probable erreur
		##    du à un déport plus important que prévu du signal. 
		
		$rh_region->{$r}->{"eval"} = $region_eval;

		# -- for general statistics for the CDS :
		$h_eval_regions{ $region_eval }++;
	}
	# nbIEB = nbEE/2 (NOTE: EE at the CDS ends have not been recorded)
	$h_eval_regions{"nbIEB"} = $h_eval_regions{"nbIEB"}/2; 
	$h_eval_regions{"nbIEBcovered"} = $h_eval_regions{"nbIEBcovered"}/2; #-# discutable !?
	
	return($rh_region, \%h_eval_regions);
}


# ==========================================================================
# function calcul_SN_SP
# 
# This function is only called if the evaluation of predictions is requested, 
# i.e. only if the I/E structure of the transcript is known.
# 
# Calculation of the SeNsibility (SN) and SPecificity (SP) of the prediction of IEBs
# on a CDS. SN: predict sufficiently (ie all true IEBs) but SP: not too much 
# (do not over predict).
# WARNING: we are trying to predict IEBs, something that does not necessarily
# exist on one CDS. Not predicting IEB when the gene has no intron is correct.
# In this case, the number of TP is 0, and by calculation SP=SP=SP=0 (so problem).
# If CDS without intron, SN is NA, and SP=100 if FP=0 or SP=0 if FP>=1 .
# 
#
# FR: Cette fonction n'est appelée que si l'évaluation des predictions est demandée,
# donc uniquement si la structure I/E du transcrit est connue.
#
# Calcul de la SeNsibilite et SPecificite de la prediction des IEB sur un CDS.
# SN prédire suffisamment (ie tous les IEB) mais SP pas trop (ne pas surpredire)
#
# ATTENTION: on essaye de prédire des IEB, quelque chose qui n'existe pas 
# forcement. Ne pas predire d'IEB quand le gene n'a pas d'intron, est correct.
# Dans ce cas, le nb de TP vaut 0, et par calcul SP=SP=0 (donc probleme).
# 
# Regles : 
# 	Si CDS sans intron, SN vaut NA, et SP=100 si FP=0 ou SP=0 si FP>=1
# ==========================================================================
sub calcul_SN_SP 
{
	my ($rh_eval_regions, $rh_eval_all) = @_ ; # TODO: verif que des hash pas vide
	my $sn = ""; my $sp = "";

	# -- special case: transcribed without IEB and no predicted IEB 
	if ( ($rh_eval_regions->{"TP"} + $rh_eval_regions->{"FN"}) == 0 )  # nb real intron = TP + FN
	{	
		# -- SeNsibility TP/(TP+FN)
		$sn = "NA";
		# -- SPecificity TP/(TP+FP)
		if ( $rh_eval_regions->{"FP"} == 0  ) # no false IEBs are predicted
		{ $sp = 100;}
		else 
		{ $sp = 0;}
	}

	# -- transcrit with introns 
	else 
	{
		# -- SeNsibility
		$sn = $rh_eval_regions->{"TP"} *100 / ($rh_eval_regions->{"TP"} + $rh_eval_regions->{"FN"}); 
		$sn = sprintf("%.2f", $sn);
		
		# -- SPecificity 
		if ( ($rh_eval_regions->{"TP"} + $rh_eval_regions->{"FP"}) == 0 )
		{	$sp = "NA" ;} # TP=FP=0 nothing is predicted (so nothing too much) => SP=0 ou SP="NA"
		else 
		{
			$sp = $rh_eval_regions->{"TP"} *100 / ($rh_eval_regions->{"TP"} + $rh_eval_regions->{"FP"});
			$sp = sprintf("%.2f", $sp);
		} 
	}
	$rh_eval_regions->{'SN'} = $sn; 
	$rh_eval_regions->{'SP'} = $sp; 

	# -- output : evaluation of prediction of one CDS 
	my $str_evalIEB = $rh_eval_regions->{'no_prediction'}."\t".$rh_eval_regions->{'nbIEB'}."\t".$rh_eval_regions->{'nbIEBcovered'}."\t".$rh_eval_regions->{'TP'}."\t".$rh_eval_regions->{'FN'}."\t".$rh_eval_regions->{'FP'}."\t".$rh_eval_regions->{'TN'}."\t".$rh_eval_regions->{'SN'}."\t".$rh_eval_regions->{'SP'} ;
	
	# -- for the calculation of performance for all CDS
	$rh_eval_all->{'TP'} += $rh_eval_regions->{'TP'};
	$rh_eval_all->{'FN'} += $rh_eval_regions->{'FN'};
	$rh_eval_all->{'FP'} += $rh_eval_regions->{'FP'};
	$rh_eval_all->{'TN'} += $rh_eval_regions->{'TN'};

	return ($str_evalIEB, $rh_eval_all);
}


# ==========================================================================
# function regions_hash2tab_output
#
# This function transforms the hash %h_regions (the division into regions of
# a CDS) into a tabular result for output in the REGION file for all transcripts 
# (without header) or the output file of one CDS (with header). 
# FR: Cette fonction transforme le hash %h_regions (le découpage en regions d'un
# transcrit) en resultat tabulaire pour la sortie dans le fichier REGION pour 
# tous les transcrits (sans header) ou le fichier de sortie du transcrit (avec
# header). 
# ==========================================================================
sub regions_hash2tab_output
{
	my ($rh_tr, $rh_region, $evaluation, $header_bool ) = @_ ; 
	my $str_output = "";

	# -- display or not of the table header 
	if ( $header_bool == 1)
	{
		$str_output = "# -- PREDICTIONS\n"; 
		$str_output.= "#ACC\tNUM\tBEGIN-END\tLGTH\tMEAN_COV\tSIGNAL_NB\tPREDICTION";
		if ( $evaluation == 1)
		{
			$str_output.= "\tPRED\tEVAL";
		}
		$str_output.= "\n";
	}

	# -- hash -> tabulated output
	my $num_exon=0; 
	foreach my $i (sort {$a<=>$b} keys %{ $rh_region } )
	{
		$str_output .= $rh_tr->{'acc'} ;
		if ( $rh_region->{$i}->{'prediction'} =~ /exon/i )
		{
			$num_exon++; # numerotation des predicted exons
			$str_output .= "\t".$num_exon;
		}
		else # any other region is not numbered
		{
			$str_output .= "\t-";
		}

		# line for each region of the transcript 
		$str_output .= "\t".$rh_region->{$i}->{'beg'}."-".$rh_region->{$i}->{'end'}."\t".$rh_region->{$i}->{"lgth"}."\t".$rh_region->{$i}->{"cov_mean"}."\t".$rh_region->{$i}->{'nbsignaux'}."\t".$rh_region->{$i}->{'prediction'};

 		if ( $evaluation == 1 ) 
		{
			$str_output .= "\t".$rh_region->{$i}->{'pred'}."\t".$rh_region->{$i}->{'eval'};
		}
		$str_output .= "\n";
		
	}

	return $str_output;
}


# ==========================================================================
# function make_outputfile_onetr
#
# This function writes the detailed output file for one CDS. This includes:
# - header with CDS information
# - the EACHPOS part with median_calculation, ratio, covstatut, pred
# - the predicted regions (IEB, exons, no_prediction)
# - if IE structure is known, the evaluation of the prediction if requested 
#   (option -e 1)
# FR: Cette fonction écrit le fichier detaillé de sortie pour un CDS. Cela
# comprend :
# - entete avec infos sur le CDS 
# - la partie EACHPOS avec calcul_median, ratio, covstatut, pred
# - le decoupage en regions prédites (IEB, exons, no_prediction)
# - si struct IE connue, l'evaluation de la prediction si demandé (option -e 1)
# ==========================================================================
sub make_outputfile_onetr
{
	my ( $outputfile, $rh_tr, $rh_pos, $rh_region, $evaluation, $str_evalIEB) = @_;
	
	# -- output file
	open(OUT, " > $outputfile ");
	# header with CDS information
	print OUT "#".$rh_tr->{'acc'}." LGTH=".$rh_tr->{'lgth'}." EXONS=".$rh_tr->{'exons_pos'}." COV_MEAN=".$rh_tr->{'tr_cov_mean'}." COV_MEAN_COVERED_BASES=".$rh_tr->{'tr_cov_mean_covered'}." LGTH_COVERD=".$rh_tr->{'lgth_covered'}." MEDIAN_nbBE=".$rh_tr->{'tr_nbBE_median'}."\n\n"; 

	# -- the EACHPOS part (features at each CDS position)
	print OUT &eachpos_hash2tab($rh_tr, $rh_pos);
	
	# -- the prediction part (ie predicted regions)
	print OUT "\n\n".&regions_hash2tab_output($rh_tr,$rh_region,$evaluation,1);

	# -- the evaluation part : one summary line
	if ( $evaluation == 1)
	{
		print OUT "\n\n# -- PREDICTION EVALUATION\n";
		print OUT "#CDS\tLGTH\tLGTH_COVERED\tCOV_MEAN\tCOV_MEAN_covered\tmednbBE\tNbRegCov0\tnbIEB\tnbIEBcov\tTP\tFN\tFP\tTN\tSN\tSP\n";#-# NbRegCDSEE\t
		print OUT $rh_tr->{'acc'}."\t".$rh_tr->{'lgth'}."\t".$rh_tr->{'lgth_covered'}."\t".$rh_tr->{'tr_cov_mean'}."\t".$rh_tr->{'tr_cov_mean_covered'}."\t".$rh_tr->{'tr_nbBE_median'}."\t".$str_evalIEB."\n";
	}
	
	return 1;
}


# ==========================================================================
# sub eachpos_hash2tab
#
# This function transforms the hash eachpos into a string of characters for 
# printing in a tabular file (csv).
# FR: Cette fonction transforme le hash eachpos en une chaine de caractere
# pour impression dans un fichier tabulaire (csv).
# ==========================================================================
sub eachpos_hash2tab
{
 	my ($rh_tr, $rh_pos) = @_ ;
	
	my $str_eachpos = "#ACC\tPOS\tEE\tCOV\tnbBE\tCovstatut\tMED\tRATIO\tPRED\n";

	# one line per CDS position
	foreach my $pos ( sort {$a<=>$b} keys %{ $rh_pos } )
	{
		$str_eachpos .= $rh_tr->{'acc'}."\t".$pos."\t".$rh_pos->{$pos}->{'EE'}."\t".$rh_pos->{$pos}->{'cov'}."\t".$rh_pos->{$pos}->{'nbBE'}."\t".$rh_pos->{$pos}->{'covstatut'}."\t".$rh_pos->{$pos}->{'nbBEmedian'}."\t".$rh_pos->{$pos}->{'nbBE_nbBEmedian'}."\t".$rh_pos->{$pos}->{'pred'}."\n";
	}

	return ($str_eachpos);
}


