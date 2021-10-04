#!/usr/bin/perl -w
use strict;
############################################################################################
# Program name: transpose_MAD_data.pl
# Version:      1.0
# Copyright:    Copyright 2013 by Dr. Frank You. All rights reserved
# Modified:     March, 2013
# Function:     transpose the MDA data to a format which can be used for ANOVA in SAS
# Input:        raw MAD phenotypic data file with the following format:

#  Record	Plot	Row	Colunm	Cp	Csp	Entry	Year	Location	Genotype	OIL	IOD	PAL	STE	OLE	LIO	LIN	YIELD
#  75	4575	2	6	0	0	1	2009	MD	CN18973	45.51	191.2	4.64	3.76	20.07	14.93	56.61	.
#  169	4669	4	7	0	0	2	2009	MD	CN18979	47.35	191.27	4.96	3.68	19.71	14.83	56.82	.
#  261	4761	6	8	0	0	3	2009	MD	CN18980	44.09	193.32	5.43	3.06	19.22	13.95	58.35	.
#  269	4769	6	7	0	0	4	2009	MD	CN18981	.	191.18	5.05	2.78	22	12.77	57.4	.
#  229	4729	5	6	0	0	5	2009	MD	CN18982	40.54	185.8	4.64	4.4	19.88	19.52	51.56	.
#  26	4526	1	6	0	0	6	2009	MD	CN18983	39.36	184.22	5.06	5.51	20.63	15.26	53.54	.

# Output:      a data file with the following format:

#  Record	Plot	Row	Colunm	Cp	Csp	Entry	Year	Location	Genotype	Trait	Value
#  75	4575	2	6	0	0	1	2009	MD	CN18973	OIL	45.51
#  75	4575	2	6	0	0	1	2009	MD	CN18973	IOD	191.2
#  75	4575	2	6	0	0	1	2009	MD	CN18973	PAL	4.64
#  75	4575	2	6	0	0	1	2009	MD	CN18973	STE	3.76
#  75	4575	2	6	0	0	1	2009	MD	CN18973	OLE	20.07
#  75	4575	2	6	0	0	1	2009	MD	CN18973	LIO	14.93
#  75	4575	2	6	0	0	1	2009	MD	CN18973	LIN	56.61
#  169	4669	4	7	0	0	2	2009	MD	CN18979	OIL	47.35
#  169	4669	4	7	0	0	2	2009	MD	CN18979	IOD	191.27
#  169	4669	4	7	0	0	2	2009	MD	CN18979	PAL	4.96
#  169	4669	4	7	0	0	2	2009	MD	CN18979	STE	3.68
#  169	4669	4	7	0	0	2	2009	MD	CN18979	OLE	19.71
#  169	4669	4	7	0	0	2	2009	MD	CN18979	LIO	14.83
#  169	4669	4	7	0	0	2	2009	MD	CN18979	LIN	56.82
#  ...
############################################################################################
# Author: Dr. Frank M. You
# Affiliation: Agricultural and Agri-Food Canada (AAFC), Winnipeg
# Contact: frank.you@agr.gc.ca
############################################################################################

use Getopt::Std;
use vars qw ($opt_i);
getopts ('i:');

my $data_file = $opt_i;


if (!$data_file) {
    print "Usage: \n";
    print "perl $0 \n";
    print "     -i MAD phenotypic data \n";
    exit;
}

my $tmp_file = $data_file . "_tmp.txt";
open (OUT, ">$tmp_file") or die ("Can not open the file $tmp_file: $!\n");

my $count = 0;
my @header_cols;
open (IN, "<$data_file") or die ("Can not open the file $data_file: $!\n");

my %missing_value_counts = ();
my %total_counts = ();

while (my $line = <IN>) {
	chomp($line);
	$line = &trim($line);
	next if ($line =~ /^\s*$/);

	$count++;
	if ($count == 1) {
		@header_cols = split(/\t/, $line);
		for (my $i = 0; $i <=9; $i++) {
			print OUT "$header_cols[$i]\t";
		}
		print OUT "Trait\tValue\n";
		next;
	}
    
    	
	my @cols = split(/\t/, $line);
	
	for (my $i = 10; $i < @cols; $i++) {
		for (my $j = 0; $j <=9; $j++) {
			print OUT "$cols[$j]\t";
		}

		print OUT "$header_cols[$i]\t$cols[$i]\n";
	
		my $id = $cols[7] . '_' . $cols[8] . '_' . $header_cols[$i];
		if (not exists($missing_value_counts{$id})) {
			$missing_value_counts{$id} = 0;
		}
		if ($cols[$i] eq '.') {
			$missing_value_counts{$id}++;
		}
		$total_counts{$id}++;
	
	}
}
close IN;
close OUT;

# remove experiments in which all are missing values
my $out_file = $data_file . "_converted.txt";
open (OUT, ">$out_file") or die ("Can not open the file $out_file: $!\n");

$count = 0;
@header_cols = ();
open (IN, "<$tmp_file") or die ("Can not open the file $tmp_file: $!\n");

while (my $line = <IN>) {
	chomp($line);
	$line = &trim($line);
	next if ($line =~ /^\s*$/);

	$count++;
	if ($count == 1) {
		print OUT "$line\n";
		next;
	}
    	
	my @cols = split(/\t/, $line);
	my $id = $cols[7] . '_' . $cols[8] . '_' . $cols[10];
	my $missing_count = $missing_value_counts{$id};
	if (not exists($missing_value_counts{$id})) {
#		$missing_count  = 0;
		print "???$id???\n";
	}
	if (not exists($total_counts{$id})) {
		
	}
	if ($missing_value_counts{$id} < $total_counts{$id}* 0.5) {
		print OUT "$line\n";
	} 
}
close IN;
close OUT;

# remove the temporary file
unlink $tmp_file;
print "\nConversion is done!\n";
print "\nThe $data_file has been converted to $out_file\n\n";
		
################################ End of the main program ###########################
####################################################################################    


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
