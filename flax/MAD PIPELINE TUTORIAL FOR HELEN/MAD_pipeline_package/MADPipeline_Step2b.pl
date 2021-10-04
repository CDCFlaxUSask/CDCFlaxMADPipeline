#!/usr/bin/perl -w
############################################################################################
# Program name: MAD_analysis_Step2b.pl
# Version:      1.0
# Copyright:    Copyright 2013 by Dr. Frank You. All rights reserved
# Modified:     March, 2013
# Function:     Summarize ANOVA results from SAS for joint analysis over multiple enviroments.
#               The prequisites are the same design with the same plot and subplot controls in
#               all experiments.

#     The program has the following assumptions for effect model of ANOVA:
#     (1)	For joint analysis of multiple years and locations, two effect models are used
#           (see Table 9 in the paper of You et al. 2013):
#           a.	Fixed model: all effects of Year, Location and Genotype are fixed.
#           b.	Mixed model: the Location and Genotype effects are fixed and the Year effect is fixed.  
#     (2)	For joint analysis of multiple locations in one year, the Genotype and Location
#           effects are fixed (see Table 8 in the paper of You et al. 2013).
#     (3)	For joint analysis of multiple years in one year, the Genotype and Year effects a
#           re fixed (see Table 8 in the paper of You et al. 2013, but ‘Year’ is replaced by ‘Location’.

# Input:        Three files will be automatically read into the program. So the following three
#               files must be present:
#               (1) multi_envirom_anova_stat_YL.txt
#               (2) multi_envirom_anova_stat_by_loc.txt
#               (3) multi_envirom_anova_stat_by_year.txt

# Output:
#   An ANOVA summary file will be generated:
#   MAD_multi_envirom_ANOVA_result_summary.txt

############################################################################################
# Author: Dr. Frank M. You
# Affiliation: Agricultural and Agri-Food Canada (AAFC), Winnipeg
# Contact: frank.you@agr.gc.ca
############################################################################################

use strict;

use Statistics::Distributions;       # statistical test
use Getopt::Std;

# Files of MAD data analysis Step 2 from SAS data analysis 
my $multi_anova_file = "multi_envirom_anova_stat_YL.txt";
my $multi_anova_year_file = "multi_envirom_anova_stat_by_loc.txt";
my $multi_anova_location_file = "multi_envirom_anova_stat_by_year.txt";

# output file
my $anova_output_file = "MAD_multi_envirom_ANOVA_result_summary.txt";

# global variables
my ($year, $location, $yl, $genotype, $gy, $gl, $gyl, $error);

if (-e $multi_anova_file || -e $multi_anova_year_file || -e $multi_anova_location_file) {
	open (OUT, ">$anova_output_file") or die ("Can't open the file $anova_output_file!\n");

	# read ANOVA results
	if ( -e $multi_anova_file) {
		($year, $location, $yl, $genotype, $gy, $gl, $gyl, $error) = &read_multi_anova_stats($multi_anova_file);
		# Generate ANOVA tables
		&print_anova_table($anova_output_file);
	} else {
		print "\nThere is no $multi_anova_file available.\n";	
	}

	if ( -e $multi_anova_year_file) {
		($year, $genotype, $gy, $error) = &read_multi_anova_stats_by_loc($multi_anova_year_file);
		# Generate ANOVA tables
		&print_anova_table_by_loc($anova_output_file);
	} else {
		print "\nThere is no $multi_anova_year_file available.\n";	
	}

	if ( -e $multi_anova_location_file) {
		($location, $genotype, $gl, $error) = &read_multi_anova_stats_by_year($multi_anova_location_file);
		# Generate ANOVA tables
		&print_anova_table_by_year($anova_output_file);
	} else {
		print "\nThere is no $multi_anova_location_file available.\n";	
	}
	
	close OUT;
	
	print "The ANOVA summary is saved in the file $anova_output_file.\n";
}




#############################  End of the main program ############################

###################################################################################
# Print ANOVA table for data of multiple years and locations
# Two effect modles: (a) Fixed model
#                    (b) Mixed model: Year is fixed and Location and Genotype are fixed 
###################################################################################
sub print_anova_table {
	my $anova_output_file = shift;

	print "Print MAD ANOVA summary for multiple years and locations...\n";
	
	print OUT "ANOVA summary\n";
	print OUT "\nANOVA in multiple years and locations";
	print OUT "\n";
	print OUT "Trait\tSource\tDF\tSS\tMS\tFixed model\t\t\tMixed model (Y random, and L and G fixed)\t\t\t\n";
	print OUT "\t\t\t\t\tF\tProb>F\tSignificance\tF\tProb>F\tSignificance\n";
	print OUT "\n";

	my @ids = sort(keys(%$error));
	
	for (my $i = 0; $i < @ids; $i++) {
		#my ($year, $location, $trait) = split(/\//, $ids[$i]);
		my $trait = $ids[$i];

		my $significance;
		my $error_ms = $error->{$ids[$i]}{'ERROR'}{'MS'};
		my $dfe = $error->{$ids[$i]}{'ERROR'}{'DF'};
		
		############ Year #################################
		print OUT "$trait\tYear\t";
		print OUT $year->{$ids[$i]}{'Year'}{'DF'} . "\t";
		print OUT $year->{$ids[$i]}{'Year'}{'SS'} . "\t";
		print OUT $year->{$ids[$i]}{'Year'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		# for both fixed model and mixed model
		my $f = $year->{$ids[$i]}{'Year'}{'MS'}/$error_ms;
		my $fprob=Statistics::Distributions::fprob($year->{$ids[$i]}{'Year'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
        # fixed model
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\t";
        # mixed model
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Location #################################
		print OUT "$trait\tLocation\t";
		print OUT $location->{$ids[$i]}{'Location'}{'DF'} . "\t";
		print OUT $location->{$ids[$i]}{'Location'}{'SS'} . "\t";
		print OUT $location->{$ids[$i]}{'Location'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		# fixed model
		$f = $location->{$ids[$i]}{'Location'}{'MS'}/$error_ms;;
		$fprob=Statistics::Distributions::fprob($location->{$ids[$i]}{'Location'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\t";

		# mixed model
		$f = $location->{$ids[$i]}{'Location'}{'MS'}/$yl->{$ids[$i]}{'Year*Location'}{'MS'};
		$fprob=Statistics::Distributions::fprob($location->{$ids[$i]}{'Location'}{'DF'},$yl->{$ids[$i]}{'Year*Location'}{'DF'},$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Year*Location #################################
		print OUT "$trait\tYear*Location\t";
		print OUT $yl->{$ids[$i]}{'Year*Location'}{'DF'} . "\t";
		print OUT $yl->{$ids[$i]}{'Year*Location'}{'SS'} . "\t";
		print OUT $yl->{$ids[$i]}{'Year*Location'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		$f = $yl->{$ids[$i]}{'Year*Location'}{'MS'}/$error_ms;
		$fprob=Statistics::Distributions::fprob($yl->{$ids[$i]}{'Year*Location'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		# fixed model
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\t";

		# mixed model
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Genotype #################################
		print OUT "$trait\tGenotype\t";
		print OUT $genotype->{$ids[$i]}{'Genotype'}{'DF'} . "\t";
		print OUT $genotype->{$ids[$i]}{'Genotype'}{'SS'} . "\t";
		print OUT $genotype->{$ids[$i]}{'Genotype'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		# fixed model
		$f = $genotype->{$ids[$i]}{'Genotype'}{'MS'}/$error_ms;
		$fprob=Statistics::Distributions::fprob($genotype->{$ids[$i]}{'Genotype'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\t";

		# mixed model
		$f = $genotype->{$ids[$i]}{'Genotype'}{'MS'}/$gy->{$ids[$i]}{'Year*Genotype'}{'MS'};
		$fprob=Statistics::Distributions::fprob($genotype->{$ids[$i]}{'Genotype'}{'DF'},$gy->{$ids[$i]}{'Year*Genotype'}{'DF'},$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";


		############ Year*Genotype #################################
		print OUT "$trait\tYear*Genotype\t";
		print OUT $gy->{$ids[$i]}{'Year*Genotype'}{'DF'} . "\t";
		print OUT $gy->{$ids[$i]}{'Year*Genotype'}{'SS'} . "\t";
		print OUT $gy->{$ids[$i]}{'Year*Genotype'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		#fix model
		$f = $gy->{$ids[$i]}{'Year*Genotype'}{'MS'}/$error_ms;
		$fprob=Statistics::Distributions::fprob($gy->{$ids[$i]}{'Year*Genotype'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\t";

		# mixed model
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Location*Genotype #################################
		print OUT "$trait\tLocation*Genotype\t";
		print OUT $gl->{$ids[$i]}{'Location*Genotype'}{'DF'} . "\t";
		print OUT $gl->{$ids[$i]}{'Location*Genotype'}{'SS'} . "\t";
		print OUT $gl->{$ids[$i]}{'Location*Genotype'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		# fixed model
		$f = $gl->{$ids[$i]}{'Location*Genotype'}{'MS'}/$error_ms;
		$fprob=Statistics::Distributions::fprob($gl->{$ids[$i]}{'Location*Genotype'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\t";

		# mixed model
		$f = $gl->{$ids[$i]}{'Location*Genotype'}{'MS'}/$gyl->{$ids[$i]}{'Year*Locatio*Genotyp'}{'MS'};
		$fprob=Statistics::Distributions::fprob($gl->{$ids[$i]}{'Location*Genotype'}{'DF'},$gyl->{$ids[$i]}{'Year*Locatio*Genotyp'}{'DF'},$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";


		############ Year*Location*Genotype #################################
		print OUT "$trait\tYear*Location*Genotype\t";
		print OUT $gyl->{$ids[$i]}{'Year*Locatio*Genotyp'}{'DF'} . "\t";
		print OUT $gyl->{$ids[$i]}{'Year*Locatio*Genotyp'}{'SS'} . "\t";
		print OUT $gyl->{$ids[$i]}{'Year*Locatio*Genotyp'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		#fixed model
		$f = $gyl->{$ids[$i]}{'Year*Locatio*Genotyp'}{'MS'}/$error_ms;
		$fprob=Statistics::Distributions::fprob($gyl->{$ids[$i]}{'Year*Locatio*Genotyp'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\t";
		
		# mixed model
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Error #################################
		print OUT "$trait\tERROR\t";
		print OUT $error->{$ids[$i]}{'ERROR'}{'DF'} . "\t";
		print OUT $error->{$ids[$i]}{'ERROR'}{'SS'} . "\t";
		print OUT $error->{$ids[$i]}{'ERROR'}{'MS'} . "\n";
		print OUT "\n";
	}
		
	print OUT "\n\n";
}

###################################################################################
# Print ANOVA table for data of multiple years (fixed model)
###################################################################################
sub print_anova_table_by_loc {
	my $anova_output_file = shift;

	print "Print MAD ANOVA summary for multiple years ...\n";
	
	print OUT "ANOVA summary\n";
	print OUT "\nANOVA in multiple years (fixed model)\n";
	print OUT "\n";
	print OUT "Location\tTrait\tSource\tDF\tSS\tMS\tF\tProb>F\tSignificance\n";
	print OUT "\n";

	my @ids = sort(keys(%$error));
	
	for (my $i = 0; $i < @ids; $i++) {
		my $location_trait = $ids[$i];
		my @tmps = split(/\//, $location_trait);
		my $location = $tmps[0];
		my $trait = $tmps[1];

		my $significance;
		my $error_ms = $error->{$ids[$i]}{'ERROR'}{'MS'};
		my $dfe = $error->{$ids[$i]}{'ERROR'}{'DF'};
		
		############ Year #################################
		print OUT "$location\t$trait\tYear\t";
		print OUT $year->{$ids[$i]}{'Year'}{'DF'} . "\t";
		print OUT $year->{$ids[$i]}{'Year'}{'SS'} . "\t";
		print OUT $year->{$ids[$i]}{'Year'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		# for both fixed model and mixed model
		my $f = $year->{$ids[$i]}{'Year'}{'MS'}/$error_ms;
		my $fprob=Statistics::Distributions::fprob($year->{$ids[$i]}{'Year'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
        # fixed model
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Genotype #################################
		print OUT "$location\t$trait\tGenotype\t";
		print OUT $genotype->{$ids[$i]}{'Genotype'}{'DF'} . "\t";
		print OUT $genotype->{$ids[$i]}{'Genotype'}{'SS'} . "\t";
		print OUT $genotype->{$ids[$i]}{'Genotype'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		# fixed model
		$f = $genotype->{$ids[$i]}{'Genotype'}{'MS'}/$error_ms;
		$fprob=Statistics::Distributions::fprob($genotype->{$ids[$i]}{'Genotype'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Year*Genotype #################################
		print OUT "$location\t$trait\tYear*Genotype\t";
		print OUT $gy->{$ids[$i]}{'Year*Genotype'}{'DF'} . "\t";
		print OUT $gy->{$ids[$i]}{'Year*Genotype'}{'SS'} . "\t";
		print OUT $gy->{$ids[$i]}{'Year*Genotype'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		#fix model
		$f = $gy->{$ids[$i]}{'Year*Genotype'}{'MS'}/$error_ms;
		$fprob=Statistics::Distributions::fprob($gy->{$ids[$i]}{'Year*Genotype'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Error #################################
		print OUT "$location\t$trait\tERROR\t";
		print OUT $error->{$ids[$i]}{'ERROR'}{'DF'} . "\t";
		print OUT $error->{$ids[$i]}{'ERROR'}{'SS'} . "\t";
		print OUT $error->{$ids[$i]}{'ERROR'}{'MS'} . "\n";
		print OUT "\n";
	}
		
	print OUT "\n\n";
}

###################################################################################
# Print ANOVA table for data of multiple locations (fixed model)
###################################################################################
sub print_anova_table_by_year {
	my $anova_output_file = shift;

	print "Print MAD ANOVA summary for multiple locations...\n";
	
	print OUT "ANOVA summary\n";
	print OUT "\nANOVA in multiple locations";
	print OUT "\n";
	print OUT "Year\tTrait\tSource\tDF\tSS\tMS\tF\tProb>F\tSignificance\n";
	print OUT "\n";

	my @ids = sort(keys(%$error));
	
	for (my $i = 0; $i < @ids; $i++) {
		my $year_trait = $ids[$i];
		my @tmps = split(/\//, $year_trait);
		my $year = $tmps[0];
		my $trait = $tmps[1];

		my $significance;
		my $error_ms = $error->{$ids[$i]}{'ERROR'}{'MS'};
		my $dfe = $error->{$ids[$i]}{'ERROR'}{'DF'};
		
		############ Location #################################
		print OUT "$year\t$trait\tLocation\t";
		print OUT $location->{$ids[$i]}{'Location'}{'DF'} . "\t";
		print OUT $location->{$ids[$i]}{'Location'}{'SS'} . "\t";
		print OUT $location->{$ids[$i]}{'Location'}{'MS'} . "\t";
			
		# need consider model: fixed or random
		# fixed model
		my $f = $location->{$ids[$i]}{'Location'}{'MS'}/$error_ms;;
		my $fprob=Statistics::Distributions::fprob($location->{$ids[$i]}{'Location'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Genotype #################################
		print OUT "$year\t$trait\tGenotype\t";
		print OUT $genotype->{$ids[$i]}{'Genotype'}{'DF'} . "\t";
		print OUT $genotype->{$ids[$i]}{'Genotype'}{'SS'} . "\t";
		print OUT $genotype->{$ids[$i]}{'Genotype'}{'MS'} . "\t";
			
		# fixed model
		$f = $genotype->{$ids[$i]}{'Genotype'}{'MS'}/$error_ms;
		$fprob=Statistics::Distributions::fprob($genotype->{$ids[$i]}{'Genotype'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Location*Genotype #################################
		print OUT "$year\t$trait\tLocation*Genotype\t";
		print OUT $gl->{$ids[$i]}{'Location*Genotype'}{'DF'} . "\t";
		print OUT $gl->{$ids[$i]}{'Location*Genotype'}{'SS'} . "\t";
		print OUT $gl->{$ids[$i]}{'Location*Genotype'}{'MS'} . "\t";
			
		# fixed model
		$f = $gl->{$ids[$i]}{'Location*Genotype'}{'MS'}/$error_ms;
		$fprob=Statistics::Distributions::fprob($gl->{$ids[$i]}{'Location*Genotype'}{'DF'},$dfe,$f);
		$significance = &get_significane($fprob);
		print OUT  "$f\t$fprob\t";
		print OUT "$significance\n";

		############ Error #################################
		print OUT "$year\t$trait\tERROR\t";
		print OUT $error->{$ids[$i]}{'ERROR'}{'DF'} . "\t";
		print OUT $error->{$ids[$i]}{'ERROR'}{'SS'} . "\t";
		print OUT $error->{$ids[$i]}{'ERROR'}{'MS'} . "\n";
		print OUT "\n";
	}
		
	print OUT "\n\n";
}

sub get_significane {
	my $prob = shift;
	my $significance;
	if ($prob <= 0.01) {
		$significance = "**";
	} elsif ($prob <= 0.05) {
		$significance = "*";
	} else {
		$significance = "ns";
	}
	return $significance;
}


sub read_multi_anova_stats {
	my $anova_file = shift;
	
	print "Read ANOVA results for multiple years and locations  ...\n";
	
	my %year = ();
	my %location = ();
	my %yl = ();
	my %genotype = ();
	my %gy = ();
	my %gl = ();
	my %gyl = ();
	my %error = ();
	
	open (IN, "<$anova_file") or die ("Can not open the file $anova_file: $!\n");
	while (my $line = <IN>) {
		chomp($line);
		$line = &trim($line);
		next if ($line =~ /^\s*$/);

		my @cols = split(/\t/, $line);
		my $id = $cols[0];        #trait
		my $source = $cols[1];
		my $df     = $cols[3];
		my $ss     = $cols[4];
		my $f      = $cols[5];
		#my $p      = $cols[6];
		if ($source eq 'ERROR') {
			$error{$id}{$source}{'DF'} = $df;
			$error{$id}{$source}{'SS'} = $ss;
			$error{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^year$/i) {
			$year{$id}{$source}{'DF'} = $df;
			$year{$id}{$source}{'SS'} = $ss;
			$year{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^location$/i) {
			$location{$id}{$source}{'DF'} = $df;
			$location{$id}{$source}{'SS'} = $ss;
			$location{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /year\*location/i) {
			$yl{$id}{$source}{'DF'} = $df;
			$yl{$id}{$source}{'SS'} = $ss;
			$yl{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^genotype$/i) {
			$genotype{$id}{$source}{'DF'} = $df;
			$genotype{$id}{$source}{'SS'} = $ss;
			$genotype{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^year\*genotype$/i) {
			$gy{$id}{$source}{'DF'} = $df;
			$gy{$id}{$source}{'SS'} = $ss;
			$gy{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^location\*genotype$/i) {
			$gl{$id}{$source}{'DF'} = $df;
			$gl{$id}{$source}{'SS'} = $ss;
			$gl{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^year\*locatio\*genotyp$/i) {
			$gyl{$id}{$source}{'DF'} = $df;
			$gyl{$id}{$source}{'SS'} = $ss;
			$gyl{$id}{$source}{'MS'} = $ss/$df;
		}
		else {
			print "Error: $line\n";
		}
	}

	close IN;
	
	return (\%year,\%location,\%yl, \%genotype, \%gy, \%gl, \%gyl, \%error);
}

sub read_multi_anova_stats_by_loc {
	my $anova_file = shift;
	
	print "Read ANOVA results for multiple years  ...\n";
	
	my %year = ();
	my %genotype = ();
	my %gy = ();
	my %error = ();
	
	open (IN, "<$anova_file") or die ("Can not open the file $anova_file: $!\n");
	while (my $line = <IN>) {
		chomp($line);
		$line = &trim($line);
		next if ($line =~ /^\s*$/);

		my @cols = split(/\t/, $line);
		my $id = $cols[0] . "/" . $cols[1];        #Location/Trait
		my $source = $cols[2];
		my $df     = $cols[4];
		my $ss     = $cols[5];
		my $f      = $cols[6];
		#my $p      = $cols[6];
		if ($source eq 'ERROR') {
			$error{$id}{$source}{'DF'} = $df;
			$error{$id}{$source}{'SS'} = $ss;
			$error{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^year$/i) {
			$year{$id}{$source}{'DF'} = $df;
			$year{$id}{$source}{'SS'} = $ss;
			$year{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^genotype$/i) {
			$genotype{$id}{$source}{'DF'} = $df;
			$genotype{$id}{$source}{'SS'} = $ss;
			$genotype{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^year\*genotype$/i) {
			$gy{$id}{$source}{'DF'} = $df;
			$gy{$id}{$source}{'SS'} = $ss;
			$gy{$id}{$source}{'MS'} = $ss/$df;
		}
		else {
			print "Error: $line\n";
		}
	}

	close IN;
	
	return (\%year, \%genotype, \%gy, \%error);
}

sub read_multi_anova_stats_by_year {
	my $anova_file = shift;
	
	print "Read ANOVA results for multiple locations  ...\n";
	
	my %location = ();
	my %genotype = ();
	my %gl = ();
	my %error = ();
	
	open (IN, "<$anova_file") or die ("Can not open the file $anova_file: $!\n");
	while (my $line = <IN>) {
		chomp($line);
		$line = &trim($line);
		next if ($line =~ /^\s*$/);

		my @cols = split(/\t/, $line);
		my $id = $cols[0] . "/" . $cols[1];        #Location/Trait
		my $source = $cols[2];
		my $df     = $cols[4];
		my $ss     = $cols[5];
		my $f      = $cols[6];

		if ($source eq 'ERROR') {
			$error{$id}{$source}{'DF'} = $df;
			$error{$id}{$source}{'SS'} = $ss;
			$error{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^location$/i) {
			$location{$id}{$source}{'DF'} = $df;
			$location{$id}{$source}{'SS'} = $ss;
			$location{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^genotype$/i) {
			$genotype{$id}{$source}{'DF'} = $df;
			$genotype{$id}{$source}{'SS'} = $ss;
			$genotype{$id}{$source}{'MS'} = $ss/$df;
		}
		elsif ($source =~ /^location\*genotype$/i) {
			$gl{$id}{$source}{'DF'} = $df;
			$gl{$id}{$source}{'SS'} = $ss;
			$gl{$id}{$source}{'MS'} = $ss/$df;
		}
		else {
			print "Error: $line\n";
		}
	}

	close IN;
	
	return (\%location, \%genotype, \%gl, \%error);
}

# Trim function to remove leading and trailing whitespaces
sub trim {
    my $string = shift;
    if ($string =~ /^\s+/) {
	$string =~ s/^\s+//;
    }
    if ($string =~ /\s+$/) {
	$string =~ s/\s+$//;
    }
    return $string;
}

