#!/usr/bin/env perl

=head1 NAME

ambigram.pl - ambigram script by Jitendra Narayan

=head1 SYNOPSIS

ambigram.pl --reconstruct/-r

	ambigram.pl --reconstruct/-r --conf/-c <configuration file>

ambigram.pl --overlap/-m

	ambigram.pl --overlap/-m --conf/-c <configuration file>

ambigram.pl --plot/-p

	ambigram.pl --plot/-p --conf/-c <configuration file>

ambigram.pl --annot/-a

	ambigram.pl --annot/-a --conf/-c <configuration file>

ambigram.pl --full/-f

	ambigram.pl --full/-f --conf/-c <configuration file>

#NOTE: If you used extension approach in EBA, this script might not work
=cut

use strict;
use warnings;
use 5.010;

use Bio::SeqIO;
use Cwd;
use File::chdir;
use File::Copy;
use POSIX;
use File::Temp qw(tempfile);
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin;
use FindBin '$Bin';
use lib::abs qw{./utils};
use File::Remove;
use File::Path qw(make_path remove_tree);
use Capture::Tiny ':all';
use Getopt::Long;
use Tie::File;
use Try::Tiny;
use Data::Dumper;
#use Statistics::R;
use Math::Round;
use File::Find;
use Pod::Usage;
use Parallel::ForkManager;
use String::ProgressBar;
use lib "$FindBin::Bin/.";
require 'ambigram_module.pm';

#Basics mandatory ambigram variables
my (
$outfile, 	# Name for ambigram's main configuration file
$conf,
$reconstruct,
$overlap,
$plot,
$full,
$ancestral,
$annot,
$help,		# If you are asking for help
);

# Default settings here for ambigram
my $current_version = "0.1";	#ambigram version
my %opts; my $nam;

print <<'WELCOME';

 ---------- -----
    -/-A-\-
 ------ ---- ----
 >-ambigram v0.1-<

Citation - ambigram: ambigram a tool to check if chromosomal breakpoints retain meaning when viewed or interpreted from a different direction, perspective, or orientation.
License: Creative Commons Licence
Bug-reports and requests to: jnarayan81ATgmail.com

WELCOME

$|++;
#Get options for ambigram
GetOptions(
	\%opts,
	"conf|c=s" 	=> \$conf,
	"reconstruct|r" => \$reconstruct,
	"overlap|o" 	=> \$overlap,
	"plot|p" 	=> \$plot,
	"full|f" 	=> \$full,
	"ancestral|a" 	=> \$ancestral,
	"annot|x" 	=> \$annot,
	"help|h" 	=> \$help,
);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please check manual.') unless ($reconstruct or $overlap or $plot or $annot or $full or $ancestral);
reconstructHelp($current_version) if (($reconstruct) and (!$conf));
overlapHelp($current_version) if (($overlap) and (!$conf));
annotHelp($current_version) if (($annot) and (!$conf));
fullHelp($current_version) if (($full) and (!$conf));
ancestralHelp($current_version) if (($ancestral) and (!$conf));
plotHelp($current_version) if (($plot) and (!$conf));


pod2usage(-msg => 'Please supply a valid filename.') unless ($conf && -s $conf);

# used to measure total execution time
my $start_time = time();

#Store thr opted name
if ($reconstruct) {$nam = 'reconstruct';} elsif ($overlap) {$nam = 'overlap';} elsif ($plot) {$nam = 'plot';} elsif ($annot) {$nam = 'annot';} elsif ($full) {$nam = 'full';} elsif ($ancestral) {$nam = 'ancestral';} else { print "Missing parameters\n"; exit(1);}

my $project_config_file = $conf;

# Absolute path of the current working directory
my $ambigram_path = dirname(rel2abs($0)); #print " Path of the dir\n$ambigram_path --------\n";

# Parameters_ref - stores all user-defined parameters, such as file locations and program parameters
my $param_ref = read_config_files(\$project_config_file, \$ambigram_path);  # stores the configuration in a hash reference

# Check all the parameters for their correctness
parameters_validator($param_ref); #checkin if user set the parameters right

my ($SpsArray_ref, $SpsNumber, $SpsArrayTabed) = storeSPS("$param_ref->{ref1}/sps.txt","$param_ref->{ref2}/sps.txt", $param_ref);

#---------------------------------------
if ( ($reconstruct) or ($full) ){
# Delete the directory if already exist
if ((-e $param_ref->{out_dir}) and ($param_ref->{force} == 1)){ remove_tree( $param_ref->{out_dir});}

# Creating the needed directories if they don't exist
if (!-e $param_ref->{out_dir}) { mkdir ($param_ref->{out_dir}) || die ("Couldn't create the directory specified as '$param_ref->{out_dir}', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n"); }
else { die("Directory $param_ref->{out_dir} already exists.\n"); }

if (!-e "$param_ref->{out_dir}/results") { mkdir ("$param_ref->{out_dir}/results") || die ("Couldn't create the directory with the results of ambigram's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{out_dir}/intermediate_files") { mkdir ("$param_ref->{out_dir}/intermediate_files/") || die ("Couldn't create the directory with the steps of ambigram's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{out_dir}/intermediate_files/ambi") { mkdir ("$param_ref->{out_dir}/intermediate_files/ambi") || die ("Couldn't create the directory with the steps of ambigram's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{out_dir}/intermediate_files/stat") { mkdir ("$param_ref->{out_dir}/intermediate_files/stat") || die ("Couldn't create the directory with the steps of ambigram's analysis.\nDetails: $!\n"); }

#Copy file to the locations
copy($project_config_file, "$param_ref->{out_dir}/project_config");
#Create an intermediate folder
#Write the log files for all steps
open (LOG, ">", "$param_ref->{out_dir}/log.$nam") || die ('Could not create log file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
open (LOG_ERR, ">", "$param_ref->{out_dir}/log.err.$nam") || die ('Could not create log.err file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
open (SUMMARY, ">", "$param_ref->{out_dir}/results/$param_ref->{summary}.$nam") || die ('Could not create summary file. Please check writing permission in your current directory', "\n");

open (INTERMEDIATE, ">", "$param_ref->{out_dir}/intermediate_files/ambi/$param_ref->{intermediate}") || die ('Could not create all_ambi file. Please check writing permission in your current directory', "\n");

concatAll("$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "$param_ref->{ref1}/$param_ref->{resolution}");

concatAll("$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "$param_ref->{ref2}/$param_ref->{resolution}");

print "WORKS WELL $param_ref->{resolution} \n";
reconstructTar("$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "$param_ref->{threshold}", "$param_ref->{len}", "$param_ref->{ref1}/sps.txt", $SpsNumber, "$param_ref->{ref1Name}", $param_ref);

reconstructTar("$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "$param_ref->{threshold}", "$param_ref->{len}", "$param_ref->{ref2}/sps.txt", $SpsNumber, "$param_ref->{ref2Name}", $param_ref); # concated both species name in $SpsArray_ref

}

if ( ($overlap) or ($full) ){
#---------------------------------------
# overlap all TAR to TAR breakpoints

#Write the log files for all steps
openLog($nam);
openLogErr($nam);

my ($sequence_data_ref, $id2tmp_id_ref, $tmp_id2id_ref) = parse_genome($param_ref);

foreach my $spsName (@$SpsArray_ref) {
	print "Working on species: $spsName\n";
	#Lets create the reference based breakpoints stats
	if ($param_ref->{snameRef1} eq $spsName) {
		tar2ref("$param_ref->{out_dir}/output_$param_ref->{ref2Name}/$param_ref->{snameRef1}"."_brk_$param_ref->{ref2Name}.tar3", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{out_dir}/intermediate_files/OUT_$param_ref->{snameRef1}"."_brk_$param_ref->{ref1Name}", "$param_ref->{snameRef1}", $SpsNumber, "$param_ref->{out_dir}/intermediate_files/STAT_$param_ref->{snameRef1}"."_brk_$param_ref->{ref1Name}");
	next;
	}
	elsif ($param_ref->{snameRef2} eq $spsName) {
		tar2ref("$param_ref->{out_dir}/output_$param_ref->{ref1Name}/$param_ref->{snameRef2}"."_brk_$param_ref->{ref1Name}.tar3", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{out_dir}/intermediate_files/OUT_$param_ref->{snameRef2}"."_brk_$param_ref->{ref2Name}", "$param_ref->{snameRef2}", $SpsNumber, "$param_ref->{out_dir}/intermediate_files/STAT_$param_ref->{snameRef2}"."_brk_$param_ref->{ref2Name}");
		#tar2ref("gallus_gallus_brk_finch.tar3", "final_classify_chicken.eba7", "see", "gallus_gallus", 11, "aaa");
	next;
	}

 	checkOverlapsTAR ("$param_ref->{out_dir}/output_$param_ref->{ref1Name}/$spsName"."_brk_$param_ref->{ref1Name}.tar3", "$param_ref->{out_dir}/output_$param_ref->{ref2Name}/$spsName"."_brk_$param_ref->{ref2Name}.tar3", "$param_ref->{out_dir}/intermediate_files/OUT_$spsName"."_brk_$param_ref->{ref1Name}", $spsName, $SpsNumber, "$param_ref->{extend}", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/all_all.eba00", "$param_ref->{ref1}/sps.txt", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "aaaa", "$param_ref->{out_dir}/intermediate_files/STAT_$spsName"."_brk_$param_ref->{ref1Name}", $spsName, "$param_ref->{out_dir}/intermediate_files/stat/finalOut_$param_ref->{ref1Name}","$param_ref->{out_dir}/intermediate_files/stat/missedOut_$param_ref->{ref1Name}", "$param_ref->{ref2}/classification.eba", "$param_ref->{chekerExtend}", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{out_dir}/intermediate_files/stat/countSTAT_$param_ref->{ref1Name}");

	checkOverlapsTAR ("$param_ref->{out_dir}/output_$param_ref->{ref2Name}/$spsName"."_brk_$param_ref->{ref2Name}.tar3", "$param_ref->{out_dir}/output_$param_ref->{ref1Name}/$spsName"."_brk_$param_ref->{ref1Name}.tar3", "$param_ref->{out_dir}/intermediate_files/OUT_$spsName"."_brk_$param_ref->{ref2Name}", $spsName, $SpsNumber, "$param_ref->{extend}", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/all_all.eba00", "$param_ref->{ref2}/sps.txt", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "aaaa", "$param_ref->{out_dir}/intermediate_files/STAT_$spsName"."_brk_$param_ref->{ref2Name}", $spsName, "$param_ref->{out_dir}/intermediate_files/stat/finalOut_$param_ref->{ref2Name}","$param_ref->{out_dir}/intermediate_files/stat/missedOut_$param_ref->{ref2Name}", "$param_ref->{ref1}/classification.eba", "$param_ref->{chekerExtend}", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{out_dir}/intermediate_files/stat/countSTAT_$param_ref->{ref2Name}");


print LOG "Printing the $spsName STAT!\n" if $param_ref->{verbose};
}
}

if ($annot) { # add it later or ($full)
#---------------------------------------
# Annotate the results
#Write the log files for all steps
openLog($nam);
openLogErr($nam);
print "Thanks for your interest, still working on it :) :) \n\n"; exit;
my ($sequence_data_ref, $id2tmp_id_ref, $tmp_id2id_ref) = parse_genome($param_ref);

if (!$param_ref->{check}) { print LOG_ERR "You set NOT to check option in config line\n Change in config file if needed\n"; exit(0);}

if ($param_ref->{check}==1) { # Set zero in conf, if not interested
	print LOG "Checking in user provided GFF file for G4 associated features\n";
	if (!"$param_ref->{out_dir}/intermediate_files/ambi/final.q4") { print LOG_ERR "Missing final.q4 file; You might need to run --reconstruct first\n";}
 	checkInGFF("$param_ref->{out_dir}/intermediate_files/ambi/final.q4", "$param_ref->{reference_genome_gff}", "$param_ref->{out_dir}/intermediate_files/ambi/final_annotated.q4", $param_ref->{extend} );

}

analyticStat("$param_ref->{out_dir}/intermediate_files/ambi/final_annotated.q4", "$param_ref->{out_dir}/intermediate_files/stat/final_annotated.q4.stat", $param_ref->{score}, $sequence_data_ref);
}

#---------------------------------------
if ( ($plot) or ($full) ){
# Plot the results
#Write the log files for all steps
openLog($nam);
openLogErr($nam);
print "Creating final table and plotting\n";
#"perl $scriptBase/reformatTable.pl -i Result_$ref1 > barData_$ref1"
createResultTable("$param_ref->{out_dir}/intermediate_files/stat/countSTAT_$param_ref->{ref1Name}", "$param_ref->{out_dir}/intermediate_files/stat/finalOut_$param_ref->{ref1Name}","$param_ref->{out_dir}/intermediate_files/stat/missedOut_$param_ref->{ref1Name}", "$param_ref->{out_dir}/results/results_$param_ref->{ref1Name}");
createResultTable("$param_ref->{out_dir}/intermediate_files/stat/countSTAT_$param_ref->{ref2Name}", "$param_ref->{out_dir}/intermediate_files/stat/finalOut_$param_ref->{ref2Name}","$param_ref->{out_dir}/intermediate_files/stat/missedOut_$param_ref->{ref2Name}", "$param_ref->{out_dir}/results/results_$param_ref->{ref2Name}");

reformatTable("$param_ref->{out_dir}/results/results_$param_ref->{ref1Name}", "$param_ref->{out_dir}/results/results2plot_$param_ref->{ref1Name}");
reformatTable("$param_ref->{out_dir}/results/results_$param_ref->{ref2Name}", "$param_ref->{out_dir}/results/results2plot_$param_ref->{ref2Name}");

system ("Rscript $Bin/utils/plotResult.R $param_ref->{out_dir}/results/results_$param_ref->{ref1Name} $param_ref->{out_dir}/results/results_$param_ref->{ref2Name} $param_ref->{out_dir}/results/results2plot_$param_ref->{ref1Name} $param_ref->{out_dir}/results/results2plot_$param_ref->{ref2Name}");
system ("mv Rplots.pdf $param_ref->{out_dir}/results");
}


if ( ($ancestral) or ($full) ){
#---------------------------------------
# overlap all TAR to TAR breakpoints

#Write the log files for all steps
openLog($nam);
openLogErr($nam);

my ($sequence_data_ref, $id2tmp_id_ref, $tmp_id2id_ref) = parse_genome($param_ref);
#Store all the group names
my @ancestralRef1 = sort {$a cmp $b} (split /,/ , (storeClasses("$param_ref->{ref1}/classification.eba", "$param_ref->{ref1}/sps.txt", $param_ref, $param_ref->{snameRef1}, $param_ref->{ref1})));
my @ancestralRef2 = sort {$a cmp $b} (split /,/ , (storeClasses("$param_ref->{ref2}/classification.eba", "$param_ref->{ref2}/sps.txt", $param_ref, $param_ref->{snameRef2}, $param_ref->{ref2})));

#I expect they have the same classification groups in both references
if (@ancestralRef1 ~~ @ancestralRef2) {
  say "Great ancestral EBRs classifcation groups matches in both references";
}
else {
  say "It seems both the reference have different classifcation groups !!";
  exit(1);
};

reconstructAncestral("$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "$param_ref->{threshold}", "$param_ref->{len}", \@ancestralRef1, $SpsNumber, "$param_ref->{ref1Name}", $param_ref, "$param_ref->{snameRef2}", "$param_ref->{ref1}/classification.eba");

reconstructAncestral("$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "$param_ref->{threshold}", "$param_ref->{len}", \@ancestralRef1, $SpsNumber, "$param_ref->{ref2Name}", $param_ref, "$param_ref->{snameRef1}", "$param_ref->{ref2}/classification.eba");

foreach my $grpName (@ancestralRef1) {
	print "Working on ancestral group: $grpName\n";
	#Lets create the reference based breakpoints stats
	#if ($param_ref->{snameRef1} eq $spsName) {
		#tar2ref("$param_ref->{out_dir}/output_$param_ref->{ref2Name}/$param_ref->{snameRef1}"."_brk_$param_ref->{ref2Name}.tar3", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{out_dir}/intermediate_files/OUT_$param_ref->{snameRef1}"."_brk_$param_ref->{ref1Name}", "$param_ref->{snameRef1}", $SpsNumber, "$param_ref->{out_dir}/intermediate_files/STAT_$param_ref->{snameRef1}"."_brk_$param_ref->{ref1Name}");
	#next;
	#}
	#elsif ($param_ref->{snameRef2} eq $spsName) {
		#tar2ref("$param_ref->{out_dir}/output_$param_ref->{ref1Name}/$param_ref->{snameRef2}"."_brk_$param_ref->{ref1Name}.tar3", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{out_dir}/intermediate_files/OUT_$param_ref->{snameRef2}"."_brk_$param_ref->{ref2Name}", "$param_ref->{snameRef2}", $SpsNumber, "$param_ref->{out_dir}/intermediate_files/STAT_$param_ref->{snameRef2}"."_brk_$param_ref->{ref2Name}");
		#tar2ref("gallus_gallus_brk_finch.tar3", "final_classify_chicken.eba7", "see", "gallus_gallus", 11, "aaa");
	#next;
	#}

 	#checkOverlapsTAR ("$param_ref->{out_dir}/output_$param_ref->{ref1Name}/$spsName"."_brk_$param_ref->{ref1Name}.tar3", "$param_ref->{out_dir}/output_$param_ref->{ref2Name}/$spsName"."_brk_$param_ref->{ref2Name}.tar3", "$param_ref->{out_dir}/intermediate_files/OUT_$spsName"."_brk_$param_ref->{ref1Name}", $spsName, $SpsNumber, "$param_ref->{extend}", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/all_all.eba00", "$param_ref->{ref1}/sps.txt", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "aaaa", "$param_ref->{out_dir}/intermediate_files/STAT_$spsName"."_brk_$param_ref->{ref1Name}", $spsName, "$param_ref->{out_dir}/intermediate_files/stat/finalOut_$param_ref->{ref1Name}","$param_ref->{out_dir}/intermediate_files/stat/missedOut_$param_ref->{ref1Name}", "$param_ref->{ref2}/classification.eba", "$param_ref->{chekerExtend}", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{out_dir}/intermediate_files/stat/countSTAT_$param_ref->{ref1Name}");

	#checkOverlapsTAR ("$param_ref->{out_dir}/output_$param_ref->{ref2Name}/$spsName"."_brk_$param_ref->{ref2Name}.tar3", "$param_ref->{out_dir}/output_$param_ref->{ref1Name}/$spsName"."_brk_$param_ref->{ref1Name}.tar3", "$param_ref->{out_dir}/intermediate_files/OUT_$spsName"."_brk_$param_ref->{ref2Name}", $spsName, $SpsNumber, "$param_ref->{extend}", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/all_all.eba00", "$param_ref->{ref2}/sps.txt", "$param_ref->{ref2}/$param_ref->{resolution}/EBA_OutFiles/all.hsb", "aaaa", "$param_ref->{out_dir}/intermediate_files/STAT_$spsName"."_brk_$param_ref->{ref2Name}", $spsName, "$param_ref->{out_dir}/intermediate_files/stat/finalOut_$param_ref->{ref2Name}","$param_ref->{out_dir}/intermediate_files/stat/missedOut_$param_ref->{ref2Name}", "$param_ref->{ref1}/classification.eba", "$param_ref->{chekerExtend}", "$param_ref->{ref1}/$param_ref->{resolution}/EBA_OutFiles/final_classify.eba7", "$param_ref->{out_dir}/intermediate_files/stat/countSTAT_$param_ref->{ref2Name}");


print LOG "Printing the ancestral group $grpName STAT!\n" if $param_ref->{verbose};
}
}

#-------------------------------------------------------------------------------
print LOG "Analysis finished, closing ambigram.  \nGood bye :) \n" if $param_ref->{verbose};

close(SUMMARY);
close(INTERMEDIATE);
close(LOG_ERR);
close(LOG);

############## Subs ###################
sub openLog {
my $nam =shift;
open (LOG, ">", "$param_ref->{out_dir}/log.$nam") || die ('Could not create log file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
}

sub openLogErr {
my $nam =shift;
open (LOG_ERR, ">", "$param_ref->{out_dir}/log.err.$nam") || die ('Could not create log.err file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
}

sub openSummary {
my $nam =shift;
open (SUMMARY, ">", "$param_ref->{out_dir}/results/$param_ref->{summary}.$nam") || die ('Could not create summary file. Please check writing permission in your current directory', "\n");
}
__DATA__
