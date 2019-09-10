#!/usr/bin/perl 

use strict;
use warnings;
use GD::Simple;
use Getopt::Long;

## ==================================================
## Parameters
## ==================================================
my ($eachpos, $tr_acc, $help);
my $ini = GetOptions("i=s" => \$eachpos,
										 "t=s" => \$tr_acc,
                     "help"=> \$help );

if ( $help || !$eachpos || !$tr_acc )
{
	print "USAGE : qd_make_img_1_transcrit.pl [options] -i <EACHPOS> -t <CDS_accession>\n\n" ;

	print "arguments:\n";
	print "  -i\t  Input file in eachpos format (with or without predictions)\n";
	print "  -t\t  CDS accession\n";

	print "\noptional arguments:\n";
	print "  -h\t  Show the help message\n";

	exit 0;
}
# default values (aucune)


## ==================================================
## Main
## ==================================================

# -- recup des info necessaires à la representation par defaut (CDS lgth, IE structure if known)
my ($tr_lgth, $rh_exons_lg) = &eachpos_tr_lgth_IEstruct($eachpos,$tr_acc);
my ($rh_cov, $rh_nbBE, $cov_max, $nbBE_max, $bool_predIEB) = &eachpos_extract_stats($eachpos,$tr_acc);

# -- PNG output file (in the same directory as the input file)
my $png = $eachpos;
if ( $png=~/\.csv$/ ) {	$png =~ s/\.csv/.png/; }
else { $png .=".png"; } 

# -- creer l'image avec la prediction de la structure IE
if ( $bool_predIEB ==1 ) 
{
	my $rh_predieb = &region_tr_predIEB($eachpos,$tr_acc);
	&make_graph_tr_with_cov_and_nbBEreads($tr_acc, $tr_lgth, $rh_exons_lg, $rh_cov, $rh_nbBE, $cov_max, $nbBE_max, $png, $rh_predieb);
}
#  -- creer l'image sans predition de la structure IE
else 
{
	&make_graph_tr_with_cov_and_nbBEreads($tr_acc, $tr_lgth, $rh_exons_lg, $rh_cov, $rh_nbBE, $cov_max, $nbBE_max, $png);
}



## #########################################################################
## ###################### SUBROUTINES ######################################
## #########################################################################

# ================================================================
# sub eachpos_tr_lgth_IEstruct
#
# This function retrieves in an EACHPOS file the CDS length and if known 
# its IE structure (i.e. exon positions on CDS sequence)
# ================================================================
sub eachpos_tr_lgth_IEstruct
{
	my ($eachpos, $tr_acc) = @_ ;
	my $tr_lgth;
	my %h_exons_lg;

	my $cmd = "grep -E '^#".$tr_acc."' ".$eachpos." ";
	my $tr_info = `$cmd`;
	if ( $tr_info =~ /^#$tr_acc LGTH=(\d+)/ )
	{
		$tr_lgth = $1;
		# if the IE structure is known
		if ( $tr_info !~ /EXONS=UNK/ ) 
		{
			$tr_info =~ /EXONS=([^\s]+)/;
			my $list_exons = $1; 
			$list_exons =~ s/;$//;
			my @a_exons = split(/;/, $list_exons) ; 
			my $i=0; 
			foreach my $exon_pos (@a_exons)
			{
				$i++;
				$exon_pos =~ /(\d+)-(\d+)/; 
				$h_exons_lg{$i} = $2 - $1 +1 ;
			}
		}
	}

	FFF:return ($tr_lgth, \%h_exons_lg);
}

# ================================================================
# sub region_tr_predIEB
#
# This function retrieves in the IEB prediction output file the 
# positions of the predicted regions. 
# ================================================================
sub region_tr_predIEB
{
	my ($predieb_file, $tr_acc) = @_ ;
	my $tr_lgth;
	my %h_regions_lg;
	my $i=0;
	my $bool="n";

	open(CSV, $predieb_file) or die "Error: cannot open file >$predieb_file< !!\n" ;
	while(my $line=<CSV>)
	{
		chomp $line;
		
		if ($bool eq "y")
		{
			if ( $line =~ /^$tr_acc\t/ )
			{
				undef my @a_line;
				@a_line = split(/\t/, $line);
				my $region_lg = $a_line[3];
				my $pred="";
				if ( !$a_line[7] )
				{
					if ( $a_line[6] =~/^exon/i ) { $pred="-";}
					elsif ($a_line[6] =~/^IEB/i ) { $pred="+";}
					else  { $pred="no_prediction";}
				}
				else
				{ 
					 $pred=$a_line[7];
				}
				$i++;
				$h_regions_lg{$i} = $region_lg."::".$pred;
			}
			elsif (  $line =~ /^# -- PREDICTION EVALUATION/ ) { $bool ="n"; goto GGG; }
		}

		if ( $line =~/^# -- PREDICTIONS/ ) { $bool ="y";}
	
	}
	GGG:close(CSV);

	return (\%h_regions_lg);
}

# ================================================================
# sub eachpos_extract_stats
#
# Retrieves the features (cov, nbBE) at each CDS position. It identifies
# the maximum coverage and the maximum nbBE to adjust the image size.
# ================================================================
sub eachpos_extract_stats
{
 	my ($eachpos, $tr_acc ) = @_ ;
	my %h_cov; my %h_nbBE;
	my $cov_max=0;
	my $nbBE_max=0;
	my $bool_predIEB = 0;

	my $tr_lg;
	my $bool = "y";

	open(EACHPOS, $eachpos) or die "Error: cannot open file >$eachpos< !!\n" ; 
	while(my $line=<EACHPOS>)
	{
		chomp $line;
		# -- At each CDS position : cov, nbBE
		if ( ($line =~ /^$tr_acc/) && ($bool eq "y") ) 
		{
			undef my @a_line; 
			@a_line = split(/\t/, $line); 
			$h_cov{$a_line[1]} = $a_line[3]; # pos => cov
			$h_nbBE{$a_line[1]} = $a_line[4]; # pos => nbBE
			# -- cov_max, nbBE_max (to define png size) 
			if ( $a_line[3] > $cov_max) { $cov_max = $a_line[3];}
			if ( $a_line[4] > $nbBE_max) { $nbBE_max = $a_line[4];}
			# -- stop condition	
			if ( $a_line[1] == $tr_lg ) { $bool = "n"; }
		
	 	}
		elsif ( ($line=~/LGTH=(\d+)/)  && ($bool eq "y")) { $tr_lg=$1; }
		elsif ( ($bool eq "n") && ($line=~/^# -- PREDICTIONS/) ) { $bool_predIEB = 1; goto END; }
	
	}
	END:close(EACHPOS);

	return(\%h_cov, \%h_nbBE, $cov_max, $nbBE_max, $bool_predIEB );
}

# ================================================================
# sub make_graph_tr_with_cov_and_nbBEreads()
# ================================================================
=pod

=head1 make_graph_tr_with_cov_and_nbBEreads() 

 DESCRIPTION

  This function creates an image of the CDS with the information of read depth and nbBE along the sequence.

 USAGE
   make_graph_tr_with_cov_and_nbBEreads($tr_acc, $tr_lgth, $rh_exons_lg, $rh_cov, $rh_nbBE, $cov_max, $nbBE_max, $png); 
 
 ATTRIBUT  
   --tr_acc						CDS/Transcrit accession
   --tr_lgth 					CDS length
   --rh_exons_info		Length of ordered exons
	 --rh_cov 					Couverage at each CDS position 
	 --rh_nbBE					nbBE at each CDS position 
   --cov_max					Maximum coverage (for calculation of image size)
   --nbBE_max					Maximum nbBE (for calculation of image size)
   --png							PNG output file name
	 --rh_predieb				Length of ordered predicted regions (region+, region-, noprediction) 
 
 RETURN
   returns a PNG image

=cut

sub make_graph_tr_with_cov_and_nbBEreads
{	
	my ($tr_acc, $tr_lgth, $rh_exons_lg, $rh_cov, $rh_nbBE, $cov_max, $nbBE_max, $png, $rh_predieb) = @_;

	# -- max(cov_max, nbBE_max)
	my $max = ($cov_max, $nbBE_max)[$cov_max < $nbBE_max];

	# -- dimensions of the different parts of the image and representation units
	my $x_left = 120 ; # px
	my $x_right = 30; # px
	my $x_base = 2 ; # 2 px / base 

	my $y_top = 50 ;
	my $y_down_predieb = 70;
	if ( !$rh_predieb ) { $y_down_predieb = 0; }
	my $y_down = 20 ;
	my $y_tr = 30 ; # CDS representation as a boxe
	my $y_unite = 2 ; # 2 px / read (read depth or nbBE)
 
	# -- image size : $x et $y  
	my $x_graph = ($tr_lgth * $x_base) ; # graph width  = CDS lgth x nb_px/base
	my $x = $x_left + $x_graph + $x_right ; # image width  
	my $y_graph = ($max +20) * $y_unite ; # graph height = valeur max * nb_px/read
	my $y = $y_top + $y_tr + $y_graph + $y_down_predieb + $y_down ; # image height
	
	# -- create image
	my $img = new GD::Image($x,$y) or die "Error: create image !\n";
	$img->colorAllocate(255, 255, 255); # fond blanc

	# -- color allocation
	my $black = $img->colorAllocate(0, 0, 0);
	my $white = $img->colorAllocate(255, 255, 255);
	my $gray =  $img->colorAllocate(166,166,166) ; 
	my $gray1 = $img->colorAllocate(230,230,230) ; #D8D8D8 light grey
	my $gray2 = $img->colorAllocate(179,179,179) ; #585858 dark grey
	my %h_color = (1=>$gray1, 2=>$gray2, 3=>$white, 4=>$black);


	# -- scale representation to know position along the CDS (every 100pb)
	$img = 	&draw_scale_bp($img, $x_left, $y_top, $tr_lgth, $x_base, $gray2);

	# -- representation of the CDS, ie each exon with alternating color 
	$img = &draw_CDS($img, $x_left, $y_top+5, $tr_acc, $tr_lgth, $x_base, $rh_exons_lg, \%h_color );

	# -- representation of the CDS (ie each exon with alternating color if IE structure known)
	$y = $y_top + 30;
	foreach my $pos ( sort {$a<=>$b} keys %{ $rh_cov } ) 
	{
		my $x1 = $x_left + ( $pos*$x_base);
		my $x2 = $x1 + $x_base;
		my $y1 = $y ; 
		my $y2 = $y + ( $rh_cov->{$pos} * $y_unite ); 
		$img->filledRectangle($x1, $y1 , $x2 , $y2 , $gray2);# cov à une position du transcrit
	}

	# -- scale of read depth (every 50X or every 10X)
	my $i=0; my $j=0; my $pas;
	if ( $max =~ /^(\d)(\d)\d\.?/ )
	{
		$pas=50;
		$j = $1*2; # 100X = 2 * 50X (2 traits d'echelle par 100X )
		if ( $2 >= 3 && $2 < 8 ) { $j++ ; }
		#elsif ( $2 >= 8) { $j+=2 ; }
	}
	elsif ( $max =~ /^(\d)(\d)\.?/ )
	{
		$pas=10;
		$j=$1; # 1 trait d'echelle tous les 10X 
		if ( $2 >= 5) { $j++; }
	}
	
	while ($i < $j)
	{
		$i++; 
		my $ybis = $y + ($i * $pas * $y_unite); 
		# trait pour materialiser i * $pas
		$img->filledRectangle($x_left, $ybis , $x_left + ($tr_lgth * $x_base) , $ybis+1 , $gray2); 
		my $nb = $i * $pas ;
		$img->string(gdMediumBoldFont,$x_left-40,$ybis-10,"$nb X ",$gray2); # scale information
	}

	# -- representation of nbBE value at each CDS position
	$y = $y_top + 30 ;
	foreach my $pos ( sort {$a<=>$b} keys %{ $rh_nbBE } ) 
	{
		my $x1 = $x_left + ( $pos*$x_base);
		my $x2 = $x1 + $x_base;
		my $y1 = $y ; 
		my $y2 = $y + ( $rh_nbBE->{$pos} * $y_unite ); 
		$img->filledRectangle($x1, $y1 , $x2 , $y2, $black); # nbBE
	}

	if ( $rh_predieb )
	{
		my $Y = $y_top+$y_tr+$y_graph;
		# -- scale representation to know position along the CDS (every 100pb)
		$img = 	&draw_scale_bp($img, $x_left, $Y+15, $tr_lgth, $x_base, $gray2);

		# -- representation of the CDS (ie each exon with alternating color if IE structure known)
		$img = &draw_CDS($img, $x_left, $Y+20, $tr_acc, $tr_lgth, $x_base, $rh_exons_lg, \%h_color );

		# -- representation of the predicted regions (IEB prediction, exon and no_prediction region)
		$img = &draw_structIE($img, $x_left, $Y+50, $x_base, $rh_predieb, \%h_color );
	}


	# -- image printing 
	open(IMG, "> $png " ) or die ("Error: cannot open file >$png<");
	binmode IMG;
	print IMG $img->png;
	close(IMG);

	return 1;
}



# ================================================================
# sub draw_scale_bp
#
# Draws the scale in bp along the sequence, every 50 bp with a dash
# and every 100 bp with a dash and value.
# ================================================================
sub draw_scale_bp 
{
	my ($img, $x, $y, $tr_lgth, $x_base, $color) = @_;

	my $pas2 = 100;
	$tr_lgth=~ /^(\d+)\d\d/;
	my $g = $1;
	my $p = 0;
	while ($p < $g)
	{
		$p++;
		my $xbis = $x + ($p * $pas2 * $x_base);
		# trait pour materialiser p * $pas2
		$img->filledRectangle($xbis, $y-5 , $xbis+1, $y , $color); # a dash every 100 bp
		my $nb = $p * $pas2 ;
		$img->string(gdMediumBoldFont, $xbis-6, $y-20, $nb , $color); # scale value
		# trait pour materialiser 50X plus loin
		$xbis = $xbis + (50 * $x_base);
		$img->filledRectangle($xbis, $y-5 , $xbis+1, $y , $color); #  a dash every 50 bp
	
	}

	return $img;
}

# ================================================================
# sub draw_CDS
#
# draws the CDS, each exon with alternating color if the IE structure 
# is known. 
# ================================================================
sub draw_CDS
{
	my ($img, $x, $y, $tr_acc, $tr_lgth, $x_base, $rh_exons_lg, $rh_color) = @_;

	# if IE struct known
	my $nb_exons=0;
	foreach my $exon (keys %{ $rh_exons_lg } ) 
	{
		$nb_exons++;
	}

	my $black = $img->colorAllocate(0, 0, 0);
	my $white = $img->colorAllocate(250, 250, 250);
	my $graph_x1 = 0; # ini coordonnées x1 sur le graphe
	my $graph_x2 = 0; # ini coordonnées x2 sur le graphe
	my $color_cle = 1;

	my $str_tr = $tr_acc." (".$tr_lgth." pb)";
	$img->string(gdMediumBoldFont,$x,($y-5)-40 ,$str_tr,$black); # CDS Accession 

	if ( $nb_exons ==0 ) # case: exon positions unknown
	{
			$graph_x1 = 0 ;
			# calculation of x2 position according to CDS length 
			$graph_x2 = $graph_x2 + ($tr_lgth * $x_base) ;  
			# draw
			$img->filledRectangle($x, $y-2 , $x+$graph_x2-1, $y+15, $rh_color->{"2"}); # CDS representation 
			$img->string(gdMediumBoldFont,$x+$graph_x1+($tr_lgth*$x_base/2)-5,$y,$tr_lgth,$black); # CDS length
	}
	else # case: exon positions known -> successive exons with alternating color 
	{
		foreach my $exon_num ( sort {$a<=>$b} keys %{ $rh_exons_lg } ) 
		{
			$graph_x1 = $graph_x2 +1 ; # start position of the new exon= end position of the previous exon -1
			# calculation of x2 position according to exon length 
			$graph_x2 = $graph_x2 + ($rh_exons_lg->{$exon_num} * $x_base) ;  

			# alterning color
			if ( $color_cle == 1) { $color_cle = 2; }
			else { $color_cle = 1; }
	
			# draw
			$img->filledRectangle($x+$graph_x1, $y-2 , $x+$graph_x2-1, $y+15, $rh_color->{$color_cle});#exon representation
			$img->string(gdMediumBoldFont,$x+$graph_x1+($rh_exons_lg->{$exon_num}*$x_base/2)-5,$y,$rh_exons_lg->{$exon_num},$white);#exon length
		}
	}

	return $img;
}

# ================================================================
# sub draw_structIE
#
# draws the predicted regions (IEB, exon or no_prediction regions)
# ================================================================
sub draw_structIE
{
	my ($img, $x, $y, $x_base, $rh_regions_lg, $rh_color) = @_;

	my $black = $img->colorAllocate(0, 0, 0);
	my $white = $img->colorAllocate(250, 250, 250);
	my $graph_x1 = 0; # ini coordonnées x1 sur le graphe
	my $graph_x2 = 0; # ini coordonnées x2 sur le graphe
	my $color_cle = 1;

  $img->string(gdMediumBoldFont,$x-95, $y ,"Predictions: ",$black);# draw
	foreach my $region_num ( sort {$a<=>$b} keys %{ $rh_regions_lg } ) 
	{
		$graph_x1 = $graph_x2 +1 ; # pos debut nouvel exon = pos fin region precedant +1
		# calculation of x2 position according to region length 
		$rh_regions_lg->{$region_num} =~ /^(\d+)::(.*)/; 
		my $lgth = $1;
		my $pred = $2;
		$graph_x2 = $graph_x2 + ($lgth * $x_base) ;  

		# colour according to the type of region
		my $color_string=$rh_color->{2};
		if ( $pred eq "+") { $color_cle = 4; $color_string=$white; }
		elsif ( $pred eq "-") { $color_cle = 1;}
		else { 	$color_cle = 3; }

		# draw
		$img->filledRectangle($x+$graph_x1, $y-2 , $x+$graph_x2-1, $y+15, $rh_color->{$color_cle});# region representation
		$img->string(gdMediumBoldFont,$x+$graph_x1+($lgth*$x_base/2)-5,$y,$lgth,$color_string);# region length
	
	}

	return $img;
}


