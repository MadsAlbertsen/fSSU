#!/usr/bin/env perl
###############################################################################
#
#    F16S.cluster.split.pl
#
#	Splits a clustered file into sub-fasta files.
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
use Cwd;

#locally-written modules
BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params
my $global_options = checkParams();

my $in_cluster;
my $in_sequences;
my $outF;
my $minsize;
my $skip_fw;
my $skip_rw;

$in_cluster = &overrideDefault("cluster.uc",'in_cluster');
$in_sequences = &overrideDefault("sequences.fa",'in_sequences');
$outF = &overrideDefault("results.fa",'outF');
$minsize = &overrideDefault(3,'minsize');
$skip_fw = &overrideDefault(50,'skip_fw');
$skip_rw = &overrideDefault(40,'skip_rw');

my %h_cluster;
my %h_seq;
my $header = "";
my $good_seq = 0;
my $temp_seq = "";


######################################################################
# CODE HERE
######################################################################

### Read all sequences and store
	
open(IN_cluster, $in_cluster) or die("Cannot read file: $in_cluster\n");

while (my $line = <IN_cluster>)  {	
	chomp $line;
	my @sl = split(/\t/,$line);
	next if($sl[0] eq "C");
	$h_cluster{$sl[1]}++;
	$h_seq{$sl[8]} = $sl[1];
}

close IN_cluster;

#foreach my $key (sort {$h_cluster{$b} <=> $h_cluster{$a}} keys %h_cluster) { 
#	if ($h_cluster{$key} >= $minsize){
#		print "$key $h_cluster{$key}\n";
#	}
#}

open(IN_sequences, $in_sequences) or die("Cannot read file: $in_sequences\n");

my $dir = cwd();
my $sub_folder = "$dir/clusters";
if (-d $sub_folder){
	die("Aborting: The output folder already exists!\n");
	} else{
	mkdir $sub_folder;
	}

while (my $line = <IN_sequences>)  {	
	chomp $line;
	if ($line =~ m/^>/) {
		$header = $line;
		my @sl = split(/>/,$line);
		$temp_seq = $sl[1];
		if(exists($h_seq{$temp_seq})){
			if($h_cluster{$h_seq{$temp_seq}} >= $minsize){
				$good_seq = 1;
			}
		}
	} else {
		if ($good_seq == 1){
			$good_seq = 0;
			open(OUT, ">>$dir/clusters/$h_seq{$temp_seq}.fa") or die("Cannot create file: $h_seq{$temp_seq}\n");
			my $seqlen = length($line);
			my $seq = substr($line, $skip_fw, $seqlen-$skip_fw-$skip_rw);
			print OUT "$header\n$seq\n";
			close OUT;
		}
	}
	
}

close IN_sequences;


exit;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "in_cluster|c:s", "in_sequences|i:s", "minsize|m:s", "skip_fw|f:s", "skip_rw|r:s");
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

    F16S.cluster.stats.pl

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

F16S.cluster.split.pl -c -i [-h -m -r -f]

 [-help -h]                    Displays this basic usage information
 [-in_sequences -i]            Input sequences in fasta format.
 [-in_cluster -c]              Input cluster uc file.
 [-minsize -m]                 Remove clusters < X size (default: 3).
 [-skip_fw -f]                 Remove X bases from the start of the reads (default: 50)
 [-skip_rw -r]                 Remove X bases from the end of the reads (default: 40)

=cut

