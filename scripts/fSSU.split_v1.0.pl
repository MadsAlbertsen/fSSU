#!/usr/bin/env perl
###############################################################################
#
#    fSSU.split_v1.0.pl
#
#	Splits sequences based on HMM classification (16S, 18S, 23S, 28S, etc).
#    
#    Copyright (C) 2016 Mads Albertsen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#locally-written modules
BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params
my $global_options = checkParams();

my $in_classification;
my $in_sequences;
my $outF;
my $trim;
my $minlength;

$in_classification = &overrideDefault("classification.txt",'in_classification');
$in_sequences = &overrideDefault("sequences.fa",'in_sequences');
$outF = &overrideDefault("out",'out_folder');
$trim = &overrideDefault(0,'trim');
$minlength = &overrideDefault(1200,'minlength');

my $header = "";
my $seq;
my $seqnr = 0;
my %h_seq;
my %h_classification;
my %h_evalue;
my %h_strand;
my %h_stats;
my %h_start;
my %h_end;
my %h_cat;
my $potential = 0;
my $total = 0;
 
######################################################################
# CODE HERE
######################################################################

### Read all sequences and store
	
open(IN_sequences, $in_sequences) or die("Cannot read file: $in_sequences\n");

while (my $line = <IN_sequences>)  {	
	chomp $line;
	if ($line =~ m/^>/) {
		$seqnr++;
		if($seqnr != 1){$h_seq{$header} = $seq;}
		$header = $line;
		$seq = "";
	}
	else{$seq = $seq.$line;}
	$h_seq{$header} = $seq;
}

close IN_sequences;

### Read classification and store best hit 
open(IN_classification, $in_classification) or die("Cannot read file: $in_classification\n");

while (my $line = <IN_classification>)  {	
	chomp $line;
	next if ($line =~ m/^#/);
	my @split = split(/\t/, $line);
	my @classification = split(/;/, $split[8]);
	$classification[0] =~ s/Name=//; 
	my $seqname = ">".$split[0];
	if (exists($h_classification{$seqname})){
		if($split[5] < $h_evalue{$seqname}){
			$h_classification{$seqname} = $classification[0];
			$h_evalue{$seqname} = $split[5];
			$h_strand{$seqname} = $split[6];
			$h_start{$seqname} = $split[3];
			$h_end{$seqname} = $split[4];		
		}	

	} 
	else{
		$h_classification{$seqname} = $classification[0];
		$h_evalue{$seqname} = $split[5];
		$h_strand{$seqname} = $split[6];
		$h_start{$seqname} = $split[3];
		$h_end{$seqname} = $split[4];
	}
}

close IN_classification;

foreach my $hit (keys %h_classification){
	$h_cat{$h_classification{$hit}} = 1; 
}

foreach my $hit (keys %h_cat){
	open(OUT, ">$outF/$hit.fa") or die("Cannot create file: $hit\n");
	close OUT;
}

open(OUT, ">$outF/unclassified.fa") or die("Cannot create file: unclassified.fa\n");	
close OUT;	

foreach my $id (keys %h_seq){
	$total++;
	my $len = length($h_seq{$id});
	if(($len < $minlength)){next;}
	$potential++;
	if (!exists($h_classification{$id})){		
		open(OUT, ">>$outF/unclassified.fa") or die("Cannot create file: unclassified.fa\n");	
		print OUT "$id\n$h_seq{$id}\n";
		close OUT;	
		$h_stats{"Unclassified"}++;
	} else{
		open(OUT, ">>$outF/$h_classification{$id}.fa") or die("Cannot create file: $h_classification{$id}\n");
		my $len = $h_end{$id}-$h_start{$id}+1;
		if($trim == 1){$h_seq{$id} = substr($h_seq{$id}, $h_start{$id}-1, $len)}
		if($h_strand{$id} eq "-"){
			$h_seq{$id} =~ tr/ACGTacgt/TGCAtgca/;  
			$h_seq{$id} = reverse($h_seq{$id}); 
		}
		print OUT "$id\n$h_seq{$id}\n";
		close OUT;	
		$h_stats{$h_classification{$id}}++;
	}
}

open(OUT_stats, ">$outF/stats_split.txt") or die("Cannot write file: stats_split.txt\n");

print OUT_stats "\n--- Stats ---:\n";

print OUT_stats "$potential sequences over $minlength bp (Total: $total)\n\n";

foreach my $class (sort {$h_stats{$b} <=> $h_stats{$a}} keys %h_stats) { 
	if (!exists($h_stats{$class})){$h_stats{$class} = 0;}
	print OUT_stats "$class:\t$h_stats{$class}\n";
}

exit;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "in_classification|c:s", "in_sequences|i:s", "out_folder|o:s", "trim|t:+", "minlength|m:s");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );
    
	#if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    #if(!exists $options{'infile'} ) { print "**ERROR: $0 : \n"; exec("pod2usage $0"); }

    return \%options;
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

__DATA__

=head1 NAME

    F16S.mock.stats.pl

=head1 COPYRIGHT

   copyright (C) 2016 Mads Albertsen

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

	Splits sequences based on HMM classification (16S, 18S, 23S, 28S, etc).

=head1 SYNOPSIS

fSSU.split.pl -c -i -o [-h -t]

 [-help -h]                    Displays this basic usage information
 [-in_sequences -i]            Input sequences in fasta format.
 [-in_classification -c]       Input classification.
 [-out_folder -o]              Output folder.
 [-trim -t]                    Trim sequence based on BARNAP positions (flag).
 [-minlength -m]               Remove sequences < X length based on BARNAP positions (default: 1200).

=cut
