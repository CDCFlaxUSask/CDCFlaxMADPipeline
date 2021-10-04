#!/usr/bin/perl -w
############################################################################################
# Program name: transpose_MAD_data.pl
# Version:      1.0
# Copyright:    Copyright 2013 by Dr. Frank You. All rights reserved
# Modified:     March, 2013
# Function:     Summarize ANOVA results from SAS and perform adjustment based on ANOVA results and
#               relative efficient estimates
# Input:        (1) MAD phenotypic data file which is used for ANOVA
#               (2) unadjusted_subplot_anova_stats.txt from SAS output
#               (3) unadjusted_plot_anova_stats.txt from SAS output

# Output:
#   (1)	MAD_ANOVA_result_summary.txt: Summary of complete MDA ANOVA for each trait in each
#       individual experiment. 
#   (2)	MAD_control_total_means.txt: Total means of plot controls for all traits and experiments.
#   (3)	MAD_adjusted_all_data_suggested.txt: Values of test genotypes plus controls adjusted
#       by the best adjustment method. This file has the same data format as the raw phenotypic
#       data used in Step 2. This file can be also used in the next step for joint ANOVA if the
#       same design with the same plot and subplot controls. 
#   (4)	MAD_adjusted_all_data_Method1.txt: Adjusted values by Method 1.
#   (5)	MAD_adjusted_all_data_Method3.txt: Adjusted values by Method 3.
#   (6)	MAD_adjusted_all_data_Method13.txt: Adjusted values by Method 1+3.
#   (7)	MAD_adjusted_data_by_three_methods.txt: Adjusted values by all three methods, Method 1,
#       Method 3 and Method1+3.
#   (8)	MAD_adjusted_data_for_genotypes_at_multi_locations.txt: Adjusted values by the best
#       method are arranged by experiments in columns.
#   (9)	MAD_adjusted_subplot_controls_data.txt: Adjusted values for subplot controls (using 
#      the best adjustment method).
#   (10)MAD_adjusted_RE.txt: Relative efficiency comparison using different methods.

############################################################################################
# Author: Dr. Frank M. You
# Affiliation: Agricultural and Agri-Food Canada (AAFC), Winnipeg
# Contact: frank.you@agr.gc.ca
############################################################################################

#id=$year . "/" . $location . "/" . $trait
#avg id
#avg id and row
#avg id and col
#avg test plot (cp != 1)




use strict;
use Statistics::Distributions;       # Statistical test
use Statistics::Regression;          # Calculate regression coefficients
use Getopt::Std;
use vars qw ($opt_i);
getopts ('i:');


my $data_file = $opt_i;


if (!$data_file) {
    print "Usage: \n";
    print "perl $0 \n";
    print "     -i MAD design phynotypic data \n";
    exit;
}


# Files of MAD data analysis Step 1 from SAS data analysis 
my $subplot_anova_file = "unadjusted_subplot_anova_stats.txt";
my $plot_anova_file = "unadjusted_plot_anova_stats.txt";

################################################################################
# check if all input files are available
if (not -e $subplot_anova_file) {
	print "No $subplot_anova_file file exists! Please check the file name is correct.\n";
	exit;
}
if (not -e $plot_anova_file) {
	print "No $plot_anova_file file exists! Please check the file name is correct.\n";
	exit;
}
################################################################################

# Read in adjusted data by Method 1 (M1)and calculate regression coefficient and collect necessary data
my ($p_row_means, $p_col_means, $p_row_col_values, $p_row_col_testplot_means, $p_total_means, $row_ids, $col_ids, $control_genotypes)
	= &read_adjusted_data($data_file);

# Calculate regression coefficient of mean of four test genotypes or test genotypes+ subplot controls on plot control
my $p_reg_coes = calculate_regression ($p_row_col_values, $p_row_col_testplot_means, $row_ids, $col_ids);


# read subplot ANOVA results
my ($subplot_ERROR, $subplot_BLOCK) = &read_subplot_anova_stats($subplot_anova_file);

# read plot ANOVA results
my ($plot_ERROR, $plot_ROW, $plot_COLUMN) = &read_plot_anova_stats($plot_anova_file);

# Generate ANOVA tables
my $anova_output_file = "MAD_ANOVA_result_summary.txt";
my $adj_methods = &print_anova_table($anova_output_file);

#Adjust test lines using both Method 1 and Method 3 and choose the best one based on ANOVA resutls
# Output files
# (1) Output one file with all design information and raw values, Method 1 adjusted
#     values, Method 2 adjusted values and the best adjusted values.
# (2) Output one file with all design information and the best adjusted values
#     (same as raw file except values are adjusted based on ANOVA information). this file will be used in SAS
#     for further interaction analysis in the next step

# (3) For each trait, output a file with the best adjusted values in all years and locations

# output files
my $adjusted_data_all_suggested         = "MAD_adjusted_all_data_suggested.txt";
my $adjusted_data_all_m1                = "MAD_adjusted_all_data_Method1.txt";
my $adjusted_data_all_m3                = "MAD_adjusted_all_data_Method3.txt";
my $adjusted_data_three_methods         = "MAD_adjusted_data_by_three_methods.txt";
my $adjusted_data_for_multi_locations   = "MAD_adjusted_data_for_genotypes_at_multi_locations.txt";
my $adjusted_data_subplot_controls_file = "MAD_adjusted_subplot_controls_data.txt";

my ($control_values, $control_values_m1, $control_values_m3,
	$plot_control_values, $subplot_control_values,
	$plot_control_values_m1, $subplot_control_values_m1,
	$plot_control_values_m3, $subplot_control_values_m3) 
	= &adjusting_testlines_data_M1_M3($data_file, $adjusted_data_all_m1, $adjusted_data_all_m3, 
		$p_row_means, $p_col_means, $p_total_means, $p_row_col_values, $p_reg_coes);

# Read in adjusted data by Method 1 (M1)and calculate regression coefficient and collect necessary data
my ($m1_row_means, $m1_col_means, $m1_row_col_values, $m1_row_col_testplot_means, $m1_total_means, $m1_row_ids, $m1_col_ids, $m1_control_genotypes)
	= &read_adjusted_data($adjusted_data_all_m1);

# Calculate regression coefficient of mean of four test genotypes or test genotypes+ subplot controls on plot control
my $m1_reg_coes = calculate_regression ($m1_row_col_values, $m1_row_col_testplot_means, $m1_row_ids, $m1_col_ids);

# Adjust for soil variation by M1+M3
my $adjusted_data_all_m13        = "MAD_adjusted_all_data_Method13.txt";
my ($control_values_m13, $plot_control_values_m13, $subplot_control_values_m13) = &adjusting_testlines_data_M13($adjusted_data_all_m1, $adjusted_data_all_m13,
	$m1_row_means, $m1_col_means, $m1_total_means, $m1_row_col_values, $m1_reg_coes);

# Calculate mean squares of controls
my $ms   = calculate_control_error_MS($control_values, $control_genotypes);
my $ms1  = calculate_control_error_MS($control_values_m1, $m1_control_genotypes);
my $ms3  = calculate_control_error_MS($control_values_m3, $m1_control_genotypes);
my $ms13 = calculate_control_error_MS($control_values_m13, $m1_control_genotypes);

# Calculate adjustment relative efficiency RE of unadjusted and adjusted data
my $adjusting_efficiency_file = "MAD_adjusted_RE.txt";
my $best_methods = &calculate_RE($ms, $ms1, $ms3, $ms13, $adjusting_efficiency_file);

# Calculate CV of plots and subplots of unadjusted and adjusted data
my $adjusting_CV_file = "MAD_adjusted_CV.txt";
&calculate_all_CV($adjusting_CV_file,
	$plot_control_values, $subplot_control_values,
	$plot_control_values_m1, $subplot_control_values_m1,
	$plot_control_values_m3, $subplot_control_values_m3,
	$plot_control_values_m13, $subplot_control_values_m13);

# Export the best adjustment results
&export_best_adjusted_values(
	$best_methods,                     # best methods based on RE (M1, M3, M1+M3)
	$adj_methods,                      # best adjustment method based on ANOVA results (M1, M3 or unadj.)
	$data_file,
	$adjusted_data_all_m1,
	$adjusted_data_all_m3,
	$adjusted_data_all_m13,
	$adjusted_data_all_suggested,
	$adjusted_data_three_methods,
	$adjusted_data_for_multi_locations
	);

#############################  End of the main program #########################

############################### Subroutines ####################################

################################################################################
# Export the best adjusted values based on RE of adjustment methods
################################################################################
sub export_best_adjusted_values {
	my ($best_methods, $adj_methods, $data_file, $adjusted_data_all_m1, $adjusted_data_all_m3,
	$adjusted_data_all_m13, $adjusted_data_all_suggested, $adjusted_data_three_methods,
	$adjusted_data_for_multi_locations) = @_;
	
	print "Export the best adjusted values to files...\n";

	open (OUT1, ">$adjusted_data_all_suggested") or die ("Can not open the file $adjusted_data_all_suggested: $!\n");
	open (OUT2, ">$adjusted_data_three_methods") or die ("Can not open the file $adjusted_data_three_methods: $!\n");
	
	my ($data_hash,     $id_hash, $header)     = &read_data($data_file);
	my ($data_hash_m1,  $id_hash_m1)           = &read_data($adjusted_data_all_m1);
	my ($data_hash_m3,  $id_hash_m3)           = &read_data($adjusted_data_all_m3);
	my ($data_hash_m13, $id_hash_m13)          = &read_data($adjusted_data_all_m13);
	
	print OUT1 "$header\n";
	print OUT2 "$header\tRaw value\tMethod 1\tMethod 3\tMethod 1+3\tSuggested value\tSuggested method\n";
	
	foreach my $uniq_id (sort(keys(%$data_hash))) {
		print OUT1 "$uniq_id";
		print OUT2 "$uniq_id\t$data_hash->{$uniq_id}\t$data_hash_m1->{$uniq_id}\t$data_hash_m3->{$uniq_id}\t$data_hash_m13->{$uniq_id}\t";
		
		my $id = $id_hash->{$uniq_id};
		my $best_value = $data_hash->{$uniq_id};
		my $best_method = 'unadj.';
		if ($data_hash->{$uniq_id} eq '.' || $data_hash->{$uniq_id} == 0) {
			$best_method = 'undef';
		} else {
			if ($adj_methods->{$id}) {
				if ($adj_methods->{$id} eq '1' || $adj_methods->{$id} eq '3') {
					if ($best_methods->{$id} eq '1') {
						$best_value = $data_hash_m1->{$uniq_id};
						$best_method = 'M1';
					} elsif ($best_methods->{$id} eq '3') {
						$best_value = $data_hash_m3->{$uniq_id};
						$best_method = 'M3';
					} elsif ($best_methods->{$id} eq '13') {
						$best_value = $data_hash_m13->{$uniq_id};
					    $best_method = 'M1+M3';
					}
				} else {  # undef
					$best_value = $data_hash->{$uniq_id};
					$best_method = 'unadj.';
				}
			}
		}
		print OUT1 "$best_value\n";
		print OUT2 "$best_value\t$best_method\n";
		
	}
	close OUT1;
	close OUT2;

	open (OUT3, ">$adjusted_data_for_multi_locations") or die ("Can not open the file $adjusted_data_for_multi_locations: $!\n");
	open (IN, "<$adjusted_data_all_suggested") or die ("Can not open the file $adjusted_data_all_suggested: $!\n");

	my $count = 0;
	my %genotype_adjusted_data = ();
	my %years_locations = ();
	my %traits = ();

	while (my $line = <IN>) {
		$count++;
		chomp($line);
		$line = &trim($line);
		
		if ($count == 1) {
			next;
		}
		next if ($line =~ /^\s*$/);

		my @cols          = split(/\t/, $line);
		my $year          = $cols[7];
		my $location      = $cols[8];
		my $genotype      = $cols[9];
		my $trait         = $cols[10];
		my $adj_suggested = $cols[11];

		$genotype_adjusted_data{$genotype}{$trait}{$year . '/' . $location} = $adj_suggested;
		$years_locations{$year . '/' . $location} = '0';
		$traits{$trait} = '0';
	}

	
	# print adjusted values for genotypes
	my @year_location_arr = ();
	foreach my $id (sort(keys(%years_locations))) {
		push(@year_location_arr, $id);
	}
	my @trait_arr = ();
	foreach my $id (sort(keys(%traits))) {
		push(@trait_arr, $id);
	}

	# print header	
	print OUT3 "Genotype\tTrait\t";
	for (my $i = 0; $i < @year_location_arr; $i++) {
		print OUT3 $year_location_arr[$i] . "\t";
	}
	print OUT3 "Mean\tCV (%)\n";
	
	for (my $j = 0; $j < @trait_arr; $j++) {
		foreach my $genotype (sort(keys(%genotype_adjusted_data))) {
			print OUT3 "$genotype\t";
			print OUT3 "$trait_arr[$j]\t";

			my @row_arr = ();
			for (my $i = 0; $i < @year_location_arr; $i++) {
				if ($genotype_adjusted_data{$genotype}{$trait_arr[$j]}{$year_location_arr[$i]}) {
					
					if ($genotype_adjusted_data{$genotype}{$trait_arr[$j]}{$year_location_arr[$i]} ne '.') {
						print OUT3 int($genotype_adjusted_data{$genotype}{$trait_arr[$j]}{$year_location_arr[$i]}*100+0.5)/100 . "\t";
						push (@row_arr, $genotype_adjusted_data{$genotype}{$trait_arr[$j]}{$year_location_arr[$i]});
					}else {
						print OUT3 ".\t";
					}
					
				} else {
					print OUT3 ".\t";
				}
			}
			my ($mean, $cv) = ('.', '.');
			if (@row_arr > 0) {
				($mean, $cv) = &calculate_cv(\@row_arr);
			}
			if ($mean ne '.') {
				print OUT3 int($mean*100 + 0.5)/100 . "\t";  # mean
			} else {
				print OUT3 ".\t";
			}
			if ($cv ne '.') {
				print OUT3 int($cv*100 + 0.5)/100 . "\n";  # cv
			} else {
				print OUT3 ".\n";
			}
		}
	}

	close IN;
	close OUT3;
}

################################################################################
# Read data file of unadjusted and adjusted values
################################################################################
sub read_data {
	my $data_file = shift;
	open (IN, "<$data_file") or die ("Can not open the file $data_file: $!\n");
	my $count = 0;
	
	my %data_hash = ();
	my %id_hash = ();
	my $header;
	
	while (my $line = <IN>) {
		$count++;
		chomp($line);
		$line = &trim($line);
		
		if ($count == 1) { # skip
			$header = $line;
			next;
		}
		next if ($line =~ /^\s*$/);

		my @cols        = split(/\t/, $line);
		my $year        = $cols[7];
		my $location    = $cols[8];
		my $trait       = $cols[10];
		my $trait_value = $cols[11];
		my $id = $year . "/" . $location . "/" . $trait;

		my $unique_id = '';
		for (my $i = 0; $i <= 10;  $i++) {
			$unique_id .= "$cols[$i]\t";
		}
		$data_hash{$unique_id} = $trait_value;		
		$id_hash{$unique_id} = $id;
	}
	close IN;

	return (\%data_hash, \%id_hash, $header);	
}

################################################################################
# Calculate mean square of control error effect
################################################################################
sub calculate_all_CV {
	my ($adjusting_CV_file, $plot_control_values, $subplot_control_values,
	$plot_control_values_m1, $subplot_control_values_m1,
	$plot_control_values_m3, $subplot_control_values_m3,
	$plot_control_values_m13, $subplot_control_values_m13) = @_;
	print "Calculate CV...\n";
	
	open (OUT, ">$adjusting_CV_file") or die ("Can not open the file $adjusting_CV_file: $!\n");
	print OUT "Year\tLocation\tTrait\tUnadj.\t\tCV(Method1)\t\tCV(Method3)\t\tCV(Method1+3)\t\n";
	print OUT "\t\t\tPlot\tSubplot\tPlot\tSubplot\tPlot\tSubplot\tPlot\tSubplot\n";
	#my @tmp = keys($plot_control_values);
	foreach my $id (sort(keys(%$plot_control_values))) {
		my @tmp = split(/\//, $id);
		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t";
		my ($mean, $cv) = &calculate_cv($plot_control_values->{$id});
		print OUT "$cv\t";
		($mean, $cv) = &calculate_cv($subplot_control_values->{$id});
		print OUT "$cv\t";
		($mean, $cv) = &calculate_cv($plot_control_values_m1->{$id});
		print OUT "$cv\t";
		($mean, $cv) = &calculate_cv($subplot_control_values_m1->{$id});
		print OUT "$cv\t";
		($mean, $cv) = &calculate_cv($plot_control_values_m3->{$id});
		print OUT "$cv\t";
		($mean, $cv) = &calculate_cv($subplot_control_values_m3->{$id});
		print OUT "$cv\t";
		($mean, $cv) = &calculate_cv($plot_control_values_m13->{$id});
		print OUT "$cv\t";
		($mean, $cv) = &calculate_cv($subplot_control_values_m13->{$id});
		print OUT "$cv\n";
	}
	close OUT;	
}

################################################################################
# Calculate mean square of control error effect
################################################################################
sub calculate_control_error_MS {
	my ($control_values, $control_genotypes) = @_;
	my %control_mean_squares = ();
		
	foreach my $id (keys(%$control_values)) {
		my $genotype_count = 0;
		my $total_count = 0;
		my $SS1 = 0;
		my $SS2 = 0;	
	
		foreach my $genotype (keys(%$control_genotypes)) {
			$genotype_count++;
			my $array = $control_values->{$id}{$genotype};

			my $sum = 0;
			if ($array) {
				for (my $i = 0; $i < @$array; $i++) {
					$sum += @$array[$i];
					$SS1 += @$array[$i] * @$array[$i];
					$total_count++;
				}
				my $mean = $sum/@$array;
				$SS2 += scalar(@$array) * $mean * $mean;
			}
			
		}
		if ($total_count - $genotype_count > 0) {
			$control_mean_squares{$id} = ($SS1 - $SS2) / ($total_count - $genotype_count);
		} else {
			
		}
	}
	
	return \%control_mean_squares;
}


################################################################################
# Calculate regression coefficients
################################################################################
sub calculate_regression {
	my ($row_col_values, $row_col_testplot_means, $row_ids, $col_ids) = @_;
	
	my %reg_coes = ();
	
	foreach my $id (keys(%$row_col_values)) {
		my $reg = Statistics::Regression->new( "title", [ "intercept", "slope"] );
		foreach my $row (keys(%$row_ids)) {
			foreach my $col (keys(%$col_ids)) {
				# Add data points
				if ($row_col_testplot_means->{$id}{$row}{$col} && $row_col_values->{$id}{$row}{$col}) {
					$reg->include( $row_col_testplot_means->{$id}{$row}{$col}, [ 1.0, $row_col_values->{$id}{$row}{$col}]);
				}
			}
		}
		my @theta = $reg->theta();
		$reg_coes{$id} = $theta[1];
	}		
	
	return \%reg_coes;
}

################################################################################
# Read and calculate data required for adjustment from unadjusted or adjusted data file
################################################################################
sub read_adjusted_data {
	my $adjusted_data_all = shift;
	my %row_col_values = ();
	my %row_col_testplot_means = ();
	my %row_col_testplot_counts = ();
	my %total_means = ();
	my %total_means_count = ();
	my %row_means = ();
	my %col_means = (); 
	my %row_means_count = ();
	my %col_means_count = (); 
	my %control_genotypes = ();
	
	# tempory data
	my %row_ids = ();
	my %col_ids = ();
	
	open (IN, "<$adjusted_data_all") or die ("Can not open the file $adjusted_data_all: $!\n");
	my $count = 0;
	
	while (my $line = <IN>) {
		$count++;
		chomp($line);
		$line = &trim($line);
		
		next if ($count == 1);
		next if ($line =~ /^\s*$/);

		my @cols     = split(/\t/, $line);
		my $year     = $cols[7];
		my $location = $cols[8];
		my $genotype = $cols[9];
		my $trait    = $cols[10];
		my $row      = $cols[2];
		my $col      = $cols[3];
		my $cp       = $cols[4];
		my $csp      = $cols[5];
		my $id = $year . "/" . $location . "/" . $trait;
		my $trait_value = $cols[11];

#		next if (not exists($used_traits->{$trait}));
		next if ($trait_value eq '.');
		$row_ids{$row} = 1;
		$col_ids{$col} = 1;

		if ($cp == 1 || $csp > 0) {
			$control_genotypes{$genotype} = 1;
		}
		
		if ($cp == 1) {
			$row_col_values{$id}{$row}{$col} = $trait_value;
#robm
print "$id\t$row\t$col\t$row_col_values{$id}{$row}{$col}\n";
			$total_means{$id} += $trait_value;
			$total_means_count{$id}++;
			$row_means{$id}{$row} += $trait_value;
			$row_means_count{$id}{$row}++;
			$col_means{$id}{$col} += $trait_value;
			$col_means_count{$id}{$col}++;
		} else {
			$row_col_testplot_counts{$id}{$row}{$col}++;
			$row_col_testplot_means{$id}{$row}{$col} += $trait_value;
		}
		
	}
	close IN;

	foreach my $id (keys(%total_means)) {
		$total_means{$id} /= $total_means_count{$id};
	}

	foreach my $id (keys(%row_means)) {
		foreach my $row (keys(%row_ids)) {
			if ($row_means_count{$id}{$row}) {
				$row_means{$id}{$row} /= $row_means_count{$id}{$row};
			} else {
				$row_means{$id}{$row} = 0;
			}
		}
		foreach my $col (keys(%col_ids)) {
			if ($col_means_count{$id}{$col}) {
				$col_means{$id}{$col} /= $col_means_count{$id}{$col};
			} else{
				$col_means{$id}{$col} = 0;
			}
		}
	}
	
	foreach my $id (keys(%row_col_testplot_means)) {
		foreach my $row (keys(%row_ids)) {
			foreach my $col (keys(%col_ids)) {
				if ($row_col_testplot_counts{$id}{$row}{$col} && $row_col_testplot_means{$id}{$row}{$col}) {
					$row_col_testplot_means{$id}{$row}{$col} /= $row_col_testplot_counts{$id}{$row}{$col};
				}
			}
		}
	}
	
	return (\%row_means, \%col_means, \%row_col_values, \%row_col_testplot_means, \%total_means, \%row_ids, \%col_ids, \%control_genotypes);
}

sub calculate_RE {
	my ($ms, $ms1, $ms3, $ms13, $adjusting_efficiency_file) = @_;
	print "Calculate REs ...\n";

	my %best_methods = (); 
	open (OUT, ">$adjusting_efficiency_file") or die ("Can not open the file $adjusting_efficiency_file: $!\n");
	print OUT "Year\tLocation\tTrait\tRE(Method1)\tRE(Method3)\tRE(Method1+3)\n";
	
	foreach my $id (sort(keys(%$ms))) {
		my ($year, $location, $trait) = split(/\//, $id);
		print OUT "$year\t$location\t$trait\t";

		my $RE_method1 =  'undef';
		my $RE_method3 =  'undef';
		my $RE_method13 = 'undef';

		if ($ms1->{$id}) {
			$RE_method1 =  int($ms->{$id}/$ms1->{$id} * 100 * 100 + 0.5)/100;
		}
		if ($ms3->{$id}) {
			$RE_method3 =  int($ms->{$id}/$ms3->{$id} * 100 * 100 + 0.5)/100;
		}
		if ($ms13->{$id}) {
			$RE_method13 = int($ms->{$id}/$ms13->{$id} * 100 * 100 + 0.5)/100;
		}
		
		$best_methods{$id} = '1';
		
		if ($RE_method1 ne 'undef' && $RE_method3 ne 'undef' && $RE_method13  ne 'undef') {
			my @arr = ($RE_method1, $RE_method3, $RE_method13);
			my $max = $arr[0];
			my $ind = 1;
			for (my $i = 1; $i < @arr; $i++) {
				if ($max < $arr[$i]) {
					$max = $arr[$i];
					$ind = $i + 1;
				}
			}
			if($ind == 3) {
				$ind = '13';
			} elsif ($ind == 2) {
				$ind = '3';
			}
				
			$best_methods{$id} = $ind;
		} else {
			$best_methods{$id} = 'undef';
		}
		
		print OUT "$RE_method1\t$RE_method3\t$RE_method13\n";
	}
	close OUT;
	
	return \%best_methods;
}


sub print_check_total_means {
	my $check_total_mean_file = shift;
	print "Print check total means ...\n";

	open (OUT, ">$check_total_mean_file") or die ("Can not open the file $check_total_mean_file: $!\n");
	print OUT "Year\tLocation\tTrait\tMean\n";
	
	foreach my $id (sort(keys(%$p_total_means))) {
		my ($year, $location, $trait) = split(/\//, $id);
		print OUT "$year\t$location\t$trait\t";
		print OUT int($p_total_means->{$id}*100 + 0.5)/100 . "\n";
	}
	close OUT;
}


################################################################################
# Adjustfor soil heterogeneity using Method 1 (M1) and Method 3 (M3)
# (1) output adjusted data by M1 and M3
# (2) generate some data sets for calculating RE for unadjusted data, M1 and M3
################################################################################
sub adjusting_testlines_data_M1_M3 {
	my ($data_file, $adjusted_data_all_m1, $adjusted_data_all_m3, 
	$row_means, $column_means, $total_means, $row_col_values, $reg_coes) = @_;
	
	print "Adjust test lines and controls using both Method 1 and 3 ...\n";

	open (OUT1, ">$adjusted_data_all_m1") or die ("Can not open the file $adjusted_data_all_m1: $!\n");
	open (OUT2, ">$adjusted_data_all_m3") or die ("Can not open the file $adjusted_data_all_m3: $!\n");
	
	open (IN, "<$data_file") or die ("Can not open the file $data_file: $!\n");
	my $count = 0;
	my @header_arr = ();
	my %genotype_adjusted_data = ();
	my %years_locations = ();
	my %traits = ();
	my %control_values_m1 = ();
	my %control_values_m3 = ();
	my %control_values    = ();
	my %plot_control_values = ();
	my %subplot_control_values = ();
	my %plot_control_values_m1 = ();
	my %subplot_control_values_m1 =();
	my %plot_control_values_m3 = ();
	my %subplot_control_values_m3 = ();

	while (my $line = <IN>) {
		$count++;
		chomp($line);
		$line = &trim($line);
		
		if ($count == 1) {
			@header_arr = split(/\t/, $line);
			print OUT1 "$line\n";
			print OUT2 "$line\n";

			next;
		}
		next if ($line =~ /^\s*$/);

		my @cols     = split(/\t/, $line);
		my $year     = $cols[7];
		my $location = $cols[8];
		my $genotype = $cols[9];
		my $trait    = $cols[10];
		my $row      = $cols[2];
		my $col      = $cols[3];
		my $cp       = $cols[4];
		my $csp      = $cols[5];
		my $id = $year . "/" . $location . "/" . $trait;


		for (my $i = 0; $i <= 10;  $i++) {
			print OUT1 "$cols[$i]\t";
			print OUT2 "$cols[$i]\t";
		}
		
		my ($adj_method1, $adj_method3);
#		next if (not exists($used_traits->{$trait}));
		
		my $trait_raw_value = $cols[11];


		##################################################################
		#   Get adjusted values
		##################################################################
		
		if ($trait_raw_value eq '.' || $trait_raw_value == 0) {  # ignore missing value
			$adj_method1 = $trait_raw_value;
			$adj_method3 = $trait_raw_value;
		} else {
			($adj_method1, $adj_method3) =
			&find_adjusted_values($year, $location, $trait,  $row, $col, $trait_raw_value,
				$row_means, $column_means, $total_means, $row_col_values, $reg_coes);
		}

		if ($trait_raw_value ne '.') {
			$trait_raw_value = int($trait_raw_value*100+0.5)/100;
		}
			
		if ($adj_method1 ne '.') {
			$adj_method1 = int($adj_method1*100+0.5)/100;
		}
		if ($adj_method3 ne '.') {
			$adj_method3 = int($adj_method3*100+0.5)/100;
		}
	
		##################################################################
		#   Get control values
		##################################################################
		if ($trait_raw_value && $trait_raw_value ne '.' && ($cp == 1 || $csp > 0)) {
			if (exists($control_values{$id}{$genotype})) {
				my $array = $control_values{$id}{$genotype};
				push(@$array, $trait_raw_value);
				$control_values{$id}{$genotype} = $array;
			} else{
				my @array = ();
				push(@array, $trait_raw_value);
				$control_values{$id}{$genotype} = \@array;
			}
			if ($cp == 1) {
				if (exists($plot_control_values{$id})) {
					my $array = $plot_control_values{$id};
					push(@$array, $trait_raw_value);
					$plot_control_values{$id} = $array;
				} else {
					my @array = ();
					push(@array, $trait_raw_value);
					$plot_control_values{$id} = \@array;
				}
			}
			if ($csp > 0) {
				if (exists($subplot_control_values{$id})) {
					my $array = $subplot_control_values{$id};
					push(@$array, $trait_raw_value);
					$subplot_control_values{$id} = $array;
				} else {
					my @array = ();
					push(@array, $trait_raw_value);
					$subplot_control_values{$id} = \@array;
				}
			}			
		}
			
		if ($adj_method1 && $adj_method1 ne '.' && ($cp == 1 || $csp > 0)) {
			if (exists($control_values_m1{$id}{$genotype})) {
				my $array = $control_values_m1{$id}{$genotype};
				push(@$array, $adj_method1);
				$control_values_m1{$id}{$genotype} = $array;
			} else{
				my @array = ();
				push(@array, $adj_method1);
				$control_values_m1{$id}{$genotype} = \@array;
			}
			if ($cp == 1) {
				if (exists($plot_control_values_m1{$id})) {
					my $array = $plot_control_values_m1{$id};
					push(@$array, $adj_method1);
					$plot_control_values_m1{$id} = $array;
				} else {
					my @array = ();
					push(@array, $adj_method1);
					$plot_control_values_m1{$id} = \@array;
				}
			}
			if ($csp > 0) {
				if (exists($subplot_control_values_m1{$id})) {
					my $array = $subplot_control_values_m1{$id};
					push(@$array, $adj_method1);
					$subplot_control_values_m1{$id} = $array;
				} else {
					my @array = ();
					push(@array, $adj_method1);
					$subplot_control_values_m1{$id} = \@array;
				}
			}
		}
		
		if ($adj_method3 && $adj_method3 ne '.' && ($cp == 1 || $csp > 0)) {
			if (exists($control_values_m3{$id}{$genotype})) {
				my $array = $control_values_m3{$id}{$genotype};
				push(@$array, $adj_method3);
				$control_values_m3{$id}{$genotype} = $array;
			} else{
				my @array = ();
				push(@array, $adj_method3);
				$control_values_m3{$id}{$genotype} = \@array;
			}
			if ($cp == 1) {
				if (exists($plot_control_values_m3{$id})) {
					my $array = $plot_control_values_m3{$id};
					push(@$array, $adj_method3);
					$plot_control_values_m3{$id} = $array;
				} else {
					my @array = ();
					push(@array, $adj_method3);
					$plot_control_values_m3{$id} = \@array;
				}
			}
			if ($csp > 0) {
				if (exists($subplot_control_values_m3{$id})) {
					my $array = $subplot_control_values_m3{$id};
					push(@$array, $adj_method3);
					$subplot_control_values_m3{$id} = $array;
				} else {
					my @array = ();
					push(@array, $adj_method3);
					$subplot_control_values_m3{$id} = \@array;
				}
			}
		}
			
		#output
		print OUT1 "$adj_method1\t";
		print OUT2 "$adj_method3\t";
		
		print OUT1 "\n";
		print OUT2 "\n";
	}
	close IN;
	close OUT1;
	close OUT2;
	
	return (\%control_values, \%control_values_m1, \%control_values_m3,
			\%plot_control_values,    \%subplot_control_values,
			\%plot_control_values_m1, \%subplot_control_values_m1,
			\%plot_control_values_m3, \%subplot_control_values_m3);
	
}


sub adjusting_testlines_data_M13 {
	my ($adjusted_data_all_m1, $adjusted_data_all_m13, $m1_row_means, $m1_col_means, 
	$m1_total_means, $m1_row_col_values, $m1_reg_coes) = @_;
	
	print "Adjust test lines and controls using both Method 1 + Method 3 ...\n";
	open (OUT, ">$adjusted_data_all_m13") or die ("Can not open the file $adjusted_data_all_m13: $!\n");
	
	open (IN, "<$adjusted_data_all_m1") or die ("Can not open the file $adjusted_data_all_m1: $!\n");
	my $count = 0;
	my @header_arr = ();
	my %control_values = ();
	my %plot_control_values_m13 = ();
	my %subplot_control_values_m13 = ();
	
	while (my $line = <IN>) {
		$count++;
		chomp($line);
		$line = &trim($line);
		
		if ($count == 1) {
			@header_arr = split(/\t/, $line);
			print OUT "$line\n";

			next;
		}
		next if ($line =~ /^\s*$/);

		my @cols     = split(/\t/, $line);
		my $year     = $cols[7];
		my $location = $cols[8];
		my $genotype = $cols[9];
		my $trait    = $cols[10];
		my $row      = $cols[2];
		my $col      = $cols[3];
		my $cp       = $cols[4];
		my $csp      = $cols[5];
		my $id = $year . "/" . $location . "/" . $trait;
		my $trait_raw_value = $cols[11];
#		next if (not exists($used_traits->{$trait}));

		for (my $i = 0; $i <= 10;  $i++) {
			print OUT "$cols[$i]\t";
		}
		
		my $adj_method13;

	
		##################################################################
		#   Get adjusted values
		##################################################################
		
		if ($trait_raw_value eq '.' || $trait_raw_value == 0) {  # ignore missing value
			$adj_method13 = $trait_raw_value;
		} else {
			$adj_method13 =
			&find_adjusted_values_m3($year, $location, $trait,  $row, $col,
				$trait_raw_value, $m1_row_means, $m1_col_means, $m1_total_means,
				$m1_row_col_values, $m1_reg_coes);
		}
			
		if ($adj_method13 ne '.') {
			$adj_method13 = int($adj_method13*100+0.5)/100;
		}
				
		##################################################################
		#   Get control values including plot and subplot controls
		##################################################################
		if ($adj_method13 && $adj_method13 ne '.' && ($cp == 1 || $csp > 0)) {
			if (exists($control_values{$id}{$genotype})) {
				my $array = $control_values{$id}{$genotype};
				push(@$array, $adj_method13);
				$control_values{$id}{$genotype} = $array;
			} else{
				my @array = ();
				push(@array, $adj_method13);
				$control_values{$id}{$genotype} = \@array;
			}
			if ($cp == 1) {
				if (exists($plot_control_values_m13{$id})) {
					my $array = $plot_control_values_m13{$id};
					push(@$array, $adj_method13);
					$plot_control_values_m13{$id} = $array;
				} else {
					my @array = ();
					push(@array, $adj_method13);
					$plot_control_values_m13{$id} = \@array;
				}
			}
			if ($csp > 0) {
				if (exists($subplot_control_values_m13{$id})) {
					my $array = $subplot_control_values_m13{$id};
					push(@$array, $adj_method13);
					$subplot_control_values_m13{$id} = $array;
				} else {
					my @array = ();
					push(@array, $adj_method13);
					$subplot_control_values_m13{$id} = \@array;
				}
			}
		}
			
		#output
		print OUT "$adj_method13\n";
	}
	close OUT;
	close IN;
	
	return (\%control_values, \%plot_control_values_m13, \%subplot_control_values_m13);
}


sub find_adjusted_values_m3 {
	my ($year, $location, $trait, $row, $col, $trait_raw_value,
		$row_means, $col_means, 
		$total_means, $row_col_values, $reg_coes) = @_;
	my $adj_method3;
	my $id = $year . "/" . $location . "/" . $trait;

	#Method 3
	my $xij = $row_col_values->{$id}{$row}{$col};
	if (!$xij) {
		$xij = ($row_means->{$id}{$row} + $col_means->{$id}{$col})/2;
	}

	if ($total_means->{$id} && $reg_coes->{$id}) {
		$adj_method3 = $trait_raw_value - $reg_coes->{$id} * ($xij - $total_means->{$id});
	} else {
		$adj_method3 = $trait_raw_value;
	}
	$adj_method3 = 0 if ($adj_method3 < 0);
	
	return $adj_method3;
}

sub find_adjusted_values {
	my ($year, $location, $trait, $row, $col, $trait_raw_value, $row_means, $column_means, $total_means, $row_col_values, $reg_coes) = @_;
	my ($adj_method1, $adj_method3, $adj_best) = ($trait_raw_value, $trait_raw_value, $trait_raw_value);
	my $id = $year . "/" . $location . "/" . $trait;

	# Method 1
	if ($total_means->{$id}) {
		$adj_method1 = $trait_raw_value - $row_means->{$id}{$row} - $column_means->{$id}{$col} + 2* $total_means->{$id};
	}
	
	# Method 3
	my $xij = $row_col_values->{$id}{$row}{$col};
	if (!$xij) {
		$xij = ($row_means->{$id}{$row} + $column_means->{$id}{$col})/2;
	}

	if ($reg_coes->{$id} && $total_means->{$id}) {
		$adj_method3 = $trait_raw_value - $reg_coes->{$id} * ($xij - $total_means->{$id});
	}
	
	$adj_method1 = 0 if ($adj_method1 < 0);
	$adj_method3 = 0 if ($adj_method3 < 0);
	return ($adj_method1, $adj_method3);
}

sub read_regression_coes {
	my $regression_coe_file = shift;
	
	print "Read gression coefficients for data adjustion using Methond 3 ...\n";
	
	my %reg_coes = ();
	open (IN, "<$regression_coe_file") or die ("Can not open the file $regression_coe_file: $!\n");
	while (my $line = <IN>) {
		chomp($line);
		$line = &trim($line);
		next if ($line =~ /^\s*$/);

		my @cols = split(/\t/, $line);
		my $id = $cols[0] . "/" . $cols[1] . "/" . $cols[2];
		my $coe = $cols[3];
		$reg_coes{$id} = $coe;
	}

	close IN;
	
	return \%reg_coes;
}


sub read_check_row_column_means {
	my $plot_check_row_column_means_file = shift;
	
	print "Read plot control row and column LS means of plot control for data adjusting using Method 1 ...\n";
	
	my %row_means = ();
	my %column_means = ();
	my %sum = ();
	my %counts = ();
	my %total_means = ();
	my %traits = ();
	
	open (IN, "<$plot_check_row_column_means_file") or die ("Can not open the file $plot_check_row_column_means_file: $!\n");
	while (my $line = <IN>) {
		chomp($line);
		$line = &trim($line);
		next if ($line =~ /^\s*$/);

		my @cols = split(/\t/, $line);
		my $id = $cols[0] . "/" . $cols[1] . "/" . $cols[2];
		my $row = $cols[3];
		my $col = $cols[4];
		my $mean = $cols[5];
		$traits{$cols[2]} = 'y';   # trait name
		
		$sum{$id} += $mean;
		$counts{$id}++;
	
		if ($row ne '.') {
			$row_means{$id}{$row} = $mean;
		} else {
			$column_means{$id}{$col} = $mean;
		}
		
	}

	close IN;
	
	foreach my $id (keys(%sum)) {
		$total_means{$id} = $sum{$id}/$counts{$id};
	}
	
	return (\%row_means, \%column_means, \%total_means, \%traits);
	
}

sub read_check_row_column_values {
	my ($plot_check_row_column_values_file, $row_means, $column_means) = @_;
	
	print "Read row and column observed values of plot control for data adjusting using Method 3 ...\n";
	
	my %row_col_values = ();
	
	open (IN, "<$plot_check_row_column_values_file") or die ("Can not open the file $plot_check_row_column_values_file: $!\n");
	while (my $line = <IN>) {
		chomp($line);
		$line = &trim($line);
		next if ($line =~ /^\s*$/);

		my @cols  = split(/\t/, $line);
		my $id    = $cols[0] . "/" . $cols[1] . "/" . $cols[2];
		my $row   = $cols[3];
		my $col   = $cols[4];
		my $value = $cols[5];
		
		if ($value eq '.' && $row_means->{$id}{$row} && $column_means->{$id}{$col}) {
			$value = ($row_means->{$id}{$row} + $column_means->{$id}{$col})/2;
		}
	
		$row_col_values{$id}{$row}{$col} = $value;
	}

	close IN;
	
	return \%row_col_values;
}

sub read_subplot_anova_stats {
	my $subplot_anova_file = shift;
	
	print "Read Subplot ANOVA results to determine data adjusting methods ...\n";
	
	my %subplot_ERROR = ();
	my %subplot_BLOCK = ();
	open (IN, "<$subplot_anova_file") or die ("Can not open the file $subplot_anova_file: $!\n");
	while (my $line = <IN>) {
		chomp($line);
		$line = &trim($line);
		next if ($line =~ /^\s*$/);

		my @cols = split(/\t/, $line);
		my $id = $cols[0] . "/" . $cols[1] . "/" . $cols[2];
		my $source = $cols[3];
		my $df     = $cols[4];
		my $ss     = $cols[5];
		my $f      = $cols[6];
		my $p      = $cols[7];
		
		if ($source eq 'ERROR') {
			$subplot_ERROR{$id}{$source}{'DF'} = $df;
			$subplot_ERROR{$id}{$source}{'SS'} = $ss;
			$subplot_ERROR{$id}{$source}{'F'} = $f;
			$subplot_ERROR{$id}{$source}{'P'} = $p;
		}
		elsif ($source eq 'Whole_plot') {
			$subplot_BLOCK{$id}{$source}{'DF'} = $df;
			$subplot_BLOCK{$id}{$source}{'SS'} = $ss;
			$subplot_BLOCK{$id}{$source}{'F'} = $f;
			$subplot_BLOCK{$id}{$source}{'P'} = $p;
		}
		elsif ($source eq 'Genotype') {    # plot and subplot controls (nomaly 3 genotypes)
			$subplot_BLOCK{$id}{$source}{'DF'} = $df;
			$subplot_BLOCK{$id}{$source}{'SS'} = $ss;
			$subplot_BLOCK{$id}{$source}{'F'} = $f;
			$subplot_BLOCK{$id}{$source}{'P'} = $p;
		} else {
			print "Error: $line\n";
		}
	}

	close IN;
	
	return (\%subplot_ERROR, \%subplot_BLOCK);
}

# Print aNOVA table for both plot and subplot controls
sub print_anova_table {
	my $anova_output_file = shift;

	print "Print MAD ANOVA summary ...\n";
	
	my %adj_methods = ();
	
	open (OUT, ">$anova_output_file") or die ("Can not open the file $anova_output_file: $!\n");
	print OUT "ANOVA summary\n";
#	print OUT "===================================================================\n";
	print OUT "\n";
	print OUT "Year\tLocation\tTrait\tSource\tDF\tSS\tMS\tF\tProb\>F\tSignificance\tSuggested adjustment method by ANOVA\n";
#	print OUT "===================================================================\n";
	print OUT "\n";

	my @ids = sort(keys(%$subplot_ERROR));
	
	for (my $i = 0; $i < @ids; $i++) {
		my ($year, $location, $trait) = split(/\//, $ids[$i]);

		my $method = "Method 1";
		my $adj_method;
		
		# ignore experiments with all missing data
		next if ($subplot_ERROR->{$ids[$i]}{'ERROR'}{'SS'} == 0 && $plot_ROW->{$ids[$i]}{'Row'}{'SS'} == 0 &&
				 $plot_COLUMN->{$ids[$i]}{'Column'}{'SS'} && $plot_ERROR->{$ids[$i]}{'ERROR'}{'SS'} == 0);
		
		my $error_ms = 0;
		if ($subplot_ERROR->{$ids[$i]}{'ERROR'}{'DF'} > 0) {
			$error_ms = $subplot_ERROR->{$ids[$i]}{'ERROR'}{'SS'}/$subplot_ERROR->{$ids[$i]}{'ERROR'}{'DF'};
		}
		next if ($error_ms == 0);
		
		my $rxc_ms = $plot_ERROR->{$ids[$i]}{'ERROR'}{'SS'}/$plot_ERROR->{$ids[$i]}{'ERROR'}{'DF'};
		my $rxc_f = $rxc_ms/$error_ms;
		
		my $row_significance = "ns";
		print OUT "$year\t$location\t$trait\tRow\t";
		print OUT $plot_ROW->{$ids[$i]}{'Row'}{'DF'} . "\t";
		print OUT $plot_ROW->{$ids[$i]}{'Row'}{'SS'} . "\t";
		print OUT $plot_ROW->{$ids[$i]}{'Row'}{'SS'}/$plot_ROW->{$ids[$i]}{'Row'}{'DF'} . "\t";
		print OUT $plot_ROW->{$ids[$i]}{'Row'}{'F'} . "\t";
		print OUT $plot_ROW->{$ids[$i]}{'Row'}{'P'} . "\t";
		if ($plot_ROW->{$ids[$i]}{'Row'}{'P'} <= 0.01) {
			$row_significance = "**";
		} elsif ($plot_ROW->{$ids[$i]}{'Row'}{'P'} <= 0.05) {
			$row_significance = "*";
		} else {
			$row_significance = "ns";
		}
		print OUT "$row_significance\n";

		my $col_significance = "ns";
		print OUT "$year\t$location\t$trait\tColumn\t";
		print OUT $plot_COLUMN->{$ids[$i]}{'Column'}{'DF'} . "\t";
		print OUT $plot_COLUMN->{$ids[$i]}{'Column'}{'SS'} . "\t";
		print OUT $plot_COLUMN->{$ids[$i]}{'Column'}{'SS'}/$plot_COLUMN->{$ids[$i]}{'Column'}{'DF'} . "\t";
		print OUT $plot_COLUMN->{$ids[$i]}{'Column'}{'F'} . "\t";
		print OUT $plot_COLUMN->{$ids[$i]}{'Column'}{'P'} . "\t";
		if ($plot_COLUMN->{$ids[$i]}{'Column'}{'P'} <= 0.01) {
			$col_significance = "**";
		} elsif ($plot_COLUMN->{$ids[$i]}{'Column'}{'P'} <= 0.05) {
			$col_significance = "*";
		} else {
			$col_significance = "ns";
		}
		print OUT "$col_significance\n";

		my $rxc_significance = "ns";
		print OUT "$year\t$location\t$trait\tRow x Column(Whole plot error)\t";
		print OUT $plot_ERROR->{$ids[$i]}{'ERROR'}{'DF'} . "\t";
		print OUT $plot_ERROR->{$ids[$i]}{'ERROR'}{'SS'} . "\t";
		print OUT $rxc_ms . "\t";
		print OUT $rxc_f . "\t";

		# tempry code
		my $df1 = $plot_ERROR->{$ids[$i]}{'ERROR'}{'DF'};
		my $df2 = $subplot_ERROR->{$ids[$i]}{'ERROR'}{'DF'};

		my $fprob=Statistics::Distributions::fprob ($df1,$df2,$rxc_f);
		print OUT "$fprob\t";
		
		if ($fprob <= 0.01) {
			$rxc_significance = "**";
		} elsif ($fprob <= 0.05) {
			$rxc_significance = "*";
		} else {
			$rxc_significance = "ns";
		}
		print OUT "$rxc_significance\t";
		
		# determine method
		if (($row_significance =~ /\*/ || $col_significance =~ /\*/) && $rxc_significance =~ /ns/) {
			$method = 'Method 1';
			$adj_method = 1;
		} elsif (($row_significance =~ /ns/ && $col_significance =~ /ns/) && $rxc_significance =~ /\*/) {
			$method = 'Method 3';
			$adj_method = 3;
		} elsif (($row_significance =~ /ns/ && $col_significance =~ /ns/) && $rxc_significance =~ /ns/) {
			$method = 'Unnecessary';
			$adj_method = 0;
			# According to ANOVA results, even row and column effects are not significant, still use Method 1
#			$method = 'Method 1';
#			$adj_method = 1;
		} elsif (($row_significance =~ /\*/ || $col_significance =~ /\*/) && $rxc_significance =~ /\*/) {
			$method = 'Method 1';
			$adj_method = 1;
		} elsif (($row_significance =~ /ns/ && $col_significance =~ /ns/) && $rxc_significance =~ /\*/) {
			$method = 'Method 1 or Method 3 but';
			if ($rxc_f > $plot_COLUMN->{$ids[$i]}{'Column'}{'F'} && $rxc_f > $plot_ROW->{$ids[$i]}{'Column'}{'F'}) {
				$method .= ' Method 3 is better';
				$adj_method = 3;
			} else {
				$method .= ' Method 1 is better';
				$adj_method = 1;
			}
		}
		print OUT "$method\n";
		$adj_methods{$ids[$i]} = $adj_method;
	
		# print ANOVA of subplotc controls 
		my $significance = "ns";
		print OUT "$year\t$location\t$trait\tWhole_plot\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Whole_plot'}{'DF'} . "\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Whole_plot'}{'SS'} . "\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Whole_plot'}{'SS'}/$subplot_BLOCK->{$ids[$i]}{'Whole_plot'}{'DF'} . "\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Whole_plot'}{'F'} . "\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Whole_plot'}{'P'} . "\t";
		if ($subplot_BLOCK->{$ids[$i]}{'Whole_plot'}{'P'} <= 0.01) {
			$significance = "**";
		} elsif ($subplot_BLOCK->{$ids[$i]}{'Whole_plot'}{'P'} <= 0.05) {
			$significance = "*";
		} else {
			$significance = "ns";
		}
		print OUT "$significance\n";
	
		$significance = "ns";
		print OUT "$year\t$location\t$trait\tControl\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Genotype'}{'DF'} . "\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Genotype'}{'SS'} . "\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Genotype'}{'SS'}/$subplot_BLOCK->{$ids[$i]}{'Genotype'}{'DF'} . "\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Genotype'}{'F'} . "\t";
		print OUT $subplot_BLOCK->{$ids[$i]}{'Genotype'}{'P'} . "\t";
		if ($subplot_BLOCK->{$ids[$i]}{'Genotype'}{'P'} <= 0.01) {
			$significance = "**";
		} elsif ($subplot_BLOCK->{$ids[$i]}{'Genotype'}{'P'} <= 0.05) {
			$significance = "*";
		} else {
			$significance = "ns";
		}
		print OUT "$significance\n";
	
		
		print OUT "$year\t$location\t$trait\tSubplot error\t";
		print OUT $subplot_ERROR->{$ids[$i]}{'ERROR'}{'DF'} . "\t";
		print OUT $subplot_ERROR->{$ids[$i]}{'ERROR'}{'SS'} . "\t";
		print OUT $error_ms . "\n";
		
		if ($ids[$i] ne $ids[$i-1]) {
			print OUT  "\n";
		}
		
	}
	print OUT "\n\n";
	return \%adj_methods;
}


sub read_plot_anova_stats {
	my $plot_anova_file = shift;
	
	print "Read plot ANOVA results to determine data adjusting methods ...\n";
	
	my %plot_ERROR = ();
	my %plot_ROW = ();
	my %plot_COLUMN = ();
	open (IN, "<$plot_anova_file") or die ("Can not open the file $plot_anova_file: $!\n");
	while (my $line = <IN>) {
		chomp($line);
		$line = &trim($line);
		next if ($line =~ /^\s*$/);

		my @cols = split(/\t/, $line);
		my $id = $cols[0] . "/" . $cols[1] . "/" . $cols[2];
		my $source = $cols[3];
		my $df     = $cols[4];
		my $ss     = $cols[5];
		my $f      = $cols[6];
		my $p      = $cols[7];
		
		if ($source eq 'ERROR') {
			$plot_ERROR{$id}{$source}{'DF'} = $df;
			$plot_ERROR{$id}{$source}{'SS'} = $ss;
			$plot_ERROR{$id}{$source}{'F'} = $f;
			$plot_ERROR{$id}{$source}{'P'} = $p;
		}
		elsif ($source eq 'Row') {
			$plot_ROW{$id}{$source}{'DF'} = $df;
			$plot_ROW{$id}{$source}{'SS'} = $ss;
			$plot_ROW{$id}{$source}{'F'} = $f;
			$plot_ROW{$id}{$source}{'P'} = $p;
		}
		elsif ($source eq 'Column') {
			$plot_COLUMN{$id}{$source}{'DF'} = $df;
			$plot_COLUMN{$id}{$source}{'SS'} = $ss;
			$plot_COLUMN{$id}{$source}{'F'} = $f;
			$plot_COLUMN{$id}{$source}{'P'} = $p;
		} else {
			print "Error: $line\n";
		}
	}

	close IN;
	
	return (\%plot_ERROR, \%plot_ROW, \%plot_COLUMN);
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

sub calculate_cv {
	my $row_arr = shift;
	my ($mean, $cv) = ();
	
	my $sum = 0;
	my $sum2 = 0;
	
	if ($row_arr) {
		my $n = @$row_arr;

		for (my $i = 0; $i < $n; $i++) {
			$sum += @$row_arr[$i];
			$sum2 += @$row_arr[$i] * @$row_arr[$i];
		}
		$mean = $sum/$n;
		if ($n > 1) {
			if ($mean > 0) {
				$cv = sqrt(($sum2 - $sum*$sum/$n)/($n-1))/$mean * 100;  # cv(%)
				$cv = int($cv*100 +0.5)/100;
			} else {
				$cv = '.';
			}
		} else {
			$cv = '.';
		}
	} else {
		$cv = '.';
	}
	
	return ($mean, $cv);
}
