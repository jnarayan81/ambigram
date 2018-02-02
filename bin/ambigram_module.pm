use strict;
#use warnings;
use Bio::Root::Root;

# Turning off BioPerl warnings
$Bio::Root::Root::DEBUG = -1;
use re qw(eval);
use vars qw($matchStart);

sub overlapsHelp {
  my $ver = $_[0];
  print "\n  ambigram --overlaps $ver \n\n";
  print "    Usage: ambigram.pl --overlaps/-o --conf/-c <configuration file>\n\n";
  print "    To overlaps reads sequences into sanjivani reads\n\n";
  print "    The path to a valid ambigram configuration file. This file contains all parameters needed to execute  ambigram.\n\n";

exit(1);
}

sub reconstructHelp {
  my $ver = $_[0];
  print "\n  ambigram --reconstruct $ver \n\n";
  print "    Usage: ambigram.pl --reconstruct/-r --conf/-c <configuration file>\n\n";
  print "    To overlap the direct overllaping strings/reads sequences\n\n";
  print "    The path to a valid ambigram configuration file. This file contains all parameters needed to execute  ambigram.\n\n";

exit(1);
}

sub annotHelp {
  my $ver = $_[0];
  print "\n  ambigram --annot $ver \n\n";
  print "    Usage: ambigram.pl --annot/-x --conf/-c <configuration file>\n\n";
  print "    To annotate the ambigram results \n\n";
  print "    The path to a valid ambigram configuration file. This file contains all parameters needed to execute  ambigram.\n\n";

exit(1);
}

sub fullHelp {
  my $ver = $_[0];
  print "\n  ambigram --full $ver \n\n";
  print "    Usage: ambigram.pl --full/-f --conf/-c <configuration file>\n\n";
  print "    To run complete in one GO \n\n";
  print "    The path to a valid ambigram configuration file. This file contains all parameters needed to execute  ambigram.\n\n";

exit(1);
}

sub ancestralHelp {
  my $ver = $_[0];
  print "\n  ambigram --ancestral $ver \n\n";
  print "    Usage: ambigram.pl --ancestral/-a --conf/-c <configuration file>\n\n";
  print "    To run ancestral breakpoint checks \n\n";
  print "    The path to a valid ambigram configuration file. This file contains all parameters needed to execute  ambigram.\n\n";

exit(1);
}

sub plotHelp {
  my $ver = $_[0];
  print "\n  ambigram --plot $ver \n\n";
  print "    Usage: ambigram.pl --plot/-p --conf/-c <configuration file>\n\n";
  print "    To plot the ambigram results \n\n";
  print "    The path to a valid ambigram configuration file. This file contains all parameters needed to execute  ambigram.\n\n";

exit(1);
}

############################################################################################
# Read config files in the form element = value #comment --->
sub read_config_files {
  my $project_config_file = shift;
  my $ambigram_path = shift;
  my %param;
  #There is two config files, one for general settings and other for third party software
  open(my $user_config, "<", "$$project_config_file") || die ("Couldn't open the project configuration file: $!\n");
  open(my $ambigram_config, "<", "$$ambigram_path/../config_files/ambigram_config") || die ("The configuration file 'ambigram_config' couldn't be read, please check if the file exists and if its permissions are properly set.\n");

# BASIC PARAMETER FOR LOCATION AND FILE NAME --->
  $param{reads_dir} = read_config('reads_dir', '', $ambigram_config);

# INPUT FILES --->
  $param{data_dir} = read_config('data_dir', $param{reads_dir}, $user_config);

# PROJECT NAME --->
  $param{out_dir} = read_config('out_dir', $param{reads_dir}, $user_config);

# PROJECT CONFIGURATION --->
  $param{reference_genome_gff} = read_config('reference_genome_gff', '', $user_config);
  $param{verbose} = read_config('verbose', '', $user_config);
  $param{force} = read_config('force', '', $user_config);
  $param{mismatch} = read_config('mismatch', '', $user_config);
  $param{reverse} = read_config('reverse', '', $user_config);
  $param{overlaps} = read_config('overlaps', '', $user_config);
  $param{extend} = read_config('extend', '', $user_config);
  $param{chekerExtend} = read_config('chekerExtend', '', $user_config);
  $param{palsize} = read_config('palsize', '', $user_config);
  $param{check} = read_config('check', '', $user_config);
  $param{kmer} = read_config('kmer', '', $user_config);
  $param{strict} = read_config('strict', '', $user_config);

  $param{resolution} = read_config('resolution', '', $user_config);
  $param{ref1} = read_config('ref1', '', $user_config);
  $param{ref1Name} = read_config('ref1Name', '', $user_config);
  $param{snameRef1} = read_config('snameRef1', '', $user_config);
  $param{ref2} = read_config('ref2', '', $user_config);
  $param{ref2Name} = read_config('ref2Name', '', $user_config);
  $param{snameRef2} = read_config('snameRef2', '', $user_config);

  $param{threshold} = read_config('threshold', '', $user_config);
  $param{len} = read_config('len', '', $user_config);

# PROJECT PENALTY --->
  $param{mutation} = read_config('mutation', '', $user_config);
  $param{score} = read_config('score', '', $user_config);
  $param{consider} = read_config('consider', '', $user_config);


#GENERAL SETTINGS
  $param{mode} = read_config('mode', '', $user_config);

# QUALITY AND PERFORMANCE --->
  $param{max_processors} = read_config('max_processors', '', $user_config);

  close($user_config);

# PATH TO EXTERNAL PROGRAMS --->

  #$param{bedtools_path} = read_config('bedtools', $param{reads_dir}, $ambigram_config);

# OUTPUT NAMES --->
  $param{result_table} = read_config('result_table', '', $ambigram_config);
  $param{result_uncertain} = read_config('result_uncertain', '', $ambigram_config);
  $param{result_positive} = read_config('result_positive', '', $ambigram_config);
  $param{summary} = read_config('summary', '', $ambigram_config);
  $param{intermediate} = read_config('intermediate', '', $ambigram_config);
  $param{result_recombinants} = "recombinants"; #hard coded, must go to config_file

# EXTERNAL ERROR HANDLING --->
  $param{tries} = read_config('tries', '', $ambigram_config);
  close($ambigram_config);
  return \%param;
}

############################################################################################""
sub read_config { # file format element = value
  my ($parameter, $reads_dir, $config_file) = @_;

  seek($config_file, 0, 0);              # added to allow any order of parameters in the config files, preventing unfriendly error messages if the user changes the order
  while (my $line = <$config_file>){
    if ($line =~ /^\s*$parameter\s*=/) {    # the string to be searched in the file
      chomp ($line);
      $line =~ s/^\s*$parameter\s*=\s*//;   # removing what comes before the user input
      $line =~ s/#.*$//;                    # removing what comes after the user input (commentaries)
      $line =~ s/\s*$//;                    # removing what comes after the user input (space caracteres)
      $line =~ s/\$reads_dir/$reads_dir/;     # allows the use of "$reads_dir" in the config file as a reference to the said parameter
      if ($line eq 'undef' || $line eq '') { return; }
      else { return $line; }
    }
  }
  return;
}

############################################################################################""
# function to identify errors in the configuration files and direct the user to the needed adjustments
sub parameters_validator { #check for all parameters,
  my $param = shift;

  my $config_path = getcwd();
  $config_path =~ s/\/\w+$/\/config_files/;

# BASIC PARAMETER FOR LOCATION AND FILE NAME --->
  if (!defined $param->{reads_dir}) { die ("No path to ambigram was specified in ambigram_config at $config_path, please open this file and fill the parameter 'reads_dir'.\n"); }
  if (!-d $param->{reads_dir}) { die ("The path to ambigram isn't a valid directory, please check if the path in 'reads_dir' is correct: $param->{reads_dir}\n"); }
  if (!-w $param->{reads_dir}) { die ("You don't have permission to write in the ambigram directory, please redefine your permissions for this directory.\n"); }

# INPUT FILES --->
  if (!defined $param->{data_dir}) { die ("No path to the nucleotide files was specified in your project's configuration file, please fill the parameter 'data_dir'.\n"); }
  if (!-d $param->{data_dir}) { die ("The path to your project's nucleotide files isn't a valid directory, please check if the path in 'data_dir' is correct: $param->{data_dir}\n"); }
  if (!-r $param->{data_dir}) { die ("You don't have permission to read in your project's nucleotide directory, please redefine your permissions.\n"); }


# PROJECT CONFIGURATION --->

  if (!defined $param->{verbose}) { $param->{verbose} = 0; } # default value
  if (!defined $param->{force}) { $param->{force} = 0; } # default value
  if (!defined $param->{mismatch}) { $param->{mismatch} = 0; } # default value
  if (!defined $param->{reverse}) { $param->{reverse} = 0; } # default value
  if (!defined $param->{overlaps}) { $param->{overlaps} = 0; } # default value 0 to keep all
  if (!defined $param->{extend}) { $param->{extend} = 0; } # default value 0

# EXTERNAL ERROR HANDLING --->
  if (!defined $param->{tries} || $param->{tries} !~ /^\d+$/) { $param->{tries} = 3; } # must be a number, and not a negative one; also must be an integer

  if (!defined $param->{out_dir}) {die "Project directory not configured. Please set out_dir element in configuration file\n";}

}

# now the script loads all nucleotide sequence files to a hash structure,
# checks their validity and translates them to protein sequence
sub parse_genome {
  my ($param) = shift;
  opendir (my $nt_files_dir, $param->{data_dir}) || die ("Path to asseembly fasta files not found: $!\n");
  my (%sequence_data);
  print LOG ('Parsing overlapsd genome/contigs/scaffolds files', "\n") if $param->{verbose};
  my $id_numeric_component = 1;  # Generate unique IDs with each sequence later on
  while (my $file = readdir ($nt_files_dir)) {
    if (($file eq '.') || ($file eq '..') || ($file =~ /^\./) || ($file =~ /~$/)) { next; }  # Prevents from reading hidden or backup files
    my $file_content = new Bio::SeqIO(-format => 'fasta',-file => "$param->{data_dir}/$file");
    print LOG ('Reading file ', $file, "\n") if $param->{verbose};
    while (my $gene_info = $file_content->next_seq()) {
      my $sequence = $gene_info->seq();
      my $len = length ($sequence);
      my $accession_number = $gene_info->display_id;
      $sequence_data{$accession_number}{status} = "OK"; #everybody starts fine
      $sequence_data{$accession_number}{problem_desc} = "-"; #everybody starts fine
      if ($sequence_data{$accession_number}{status} eq "OK") { # Add check points here <<<<<<
        $sequence_data{$accession_number}{nuc_seq} = $sequence;
	$sequence_data{$accession_number}{len} = $len;
      }
    }
  }
  print LOG ('Done', "\n") if $param->{verbose};
  closedir ($nt_files_dir);
  return (\%sequence_data);
}

sub parse_gene_id {
  my @aux = split (/\(/, $_[0]);
  my $specie = $aux[1];
  $specie =~ s/\)//g;
  return ($aux[0], $specie);  # aux[0] has the if o the gene
}

sub mean {
  my @tmp = @{$_[0]};
  my $soma = 0;
  foreach my $value(@tmp) {
    $soma = $soma + $value;
  }
  my $mean = ($soma/($#tmp+1));
  return $mean;
}

############################################################################################""
# Prints the time taken by the tasks of the group before the codeml runs in a file. Used for summary
sub print_task_time {
  my ($ortholog_dir, $task_time, $name) = @_;
  my $f_tree_time = time() - $$task_time;
  open(my $fh_time_write, ">", "$$ortholog_dir/$name");
  print $fh_time_write ('Time taken by task: ', $f_tree_time, "\n");
#  print "$$fh_time_write\n";
  close ($fh_time_write);
  return;
}

############################################################################################""
sub write_summary {
  my ($param, $clusters_ref, $start_time) = @_;
  my ($f_tree, $model1, $model2, $model7, $model8) = (0,0,0,0,0);


  # printing time spent by the program to run
  my $total_time = time() - $$start_time;
  my ($hours, $minutes, $seconds) = divide_time(\$total_time);

  print SUMMARY ("Time spent: $hours:$minutes:$seconds ($total_time seconds)\n");
  foreach my $ortholog_group (keys %{$clusters_ref}) {
    if (-s "$param->{out_dir}/intermediate_files/$ortholog_group/time_id_rec") {
      open (my $fh_time_read, "<", "$param->{out_dir}/intermediate_files/$ortholog_group/time_id_rec");
      my $line = <$fh_time_read>;
      $f_tree += $1 if ($line =~ /:\s(\d+)/);
      close($fh_time_read);
    }
  }

  my $sequential_time = $f_tree+$model1+$model2+$model7+$model8;

  ($hours, $minutes, $seconds) = divide_time(\$sequential_time);
  print SUMMARY ("Total time (sequential run): $hours:$minutes:$seconds ($sequential_time seconds)\n");

  ($hours, $minutes, $seconds) = divide_time(\$f_tree);
  print SUMMARY (" - total time on building phylogenetic trees: $hours:$minutes:$seconds ($f_tree seconds)\n");

  return;
}

############################################################################################""
#Time check
sub divide_time {
  my $total_time = shift;

  my $hours = POSIX::floor( $$total_time / 3600 );
  my $minutes = POSIX::floor(($$total_time % 3600) / 60);
  if ($minutes < 10) { $minutes = '0' . $minutes; }
  my $seconds = $$total_time % 60;
  if ($seconds < 10) { $seconds = '0' . $seconds; }

  return ($hours, $minutes, $seconds);
}

############################################################################################""!!
#Function to move the files with wild
sub moveFiles {
    my ( $source_ref, $arc_dir ) = @_;
    my @old_files = @$source_ref;
    foreach my $old_file (@old_files)
         {
    #my ($short_file_name) = $old_file =~ m~/(.*?\.dat)$~;
    #my $new_file = $arc_dir . $short_file_name;
    move($old_file, $arc_dir) or die "Could not move $old_file to $arc_dir: $!\n";
   }
}

=pod

Comment here

=cut



########################################################################################"
#Store fasta file to hash
sub fastafile2hash {
    my $fastafile = shift @_;
    my %sequences;
    my $fh = &read_fh($fastafile);
    my $seqid;
    while (my $line = <$fh>) {
        if ($line =~ /^>(\S+)(.*)/) {
            $seqid = $1;
            $sequences{$seqid}{desc} = $2;
        }
        else {
            chomp $line;
            $sequences{$seqid}{seq}     .= $line;
            $sequences{$seqid}{len}     += length $line;
            $sequences{$seqid}{gc}      += ($line =~ tr/gcGC/gcGC/);
            $line =~ s/[^atgc]/N/ig;
            $sequences{$seqid}{nonatgc} += ($line =~ tr/N/N/);
        }
    }
    close $fh;
    return \%sequences;
}


########################################################################################"
#Open and Read a file
sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
}


sub extractSeq_new {
my ($name, $st, $end, $db)=@_;
my $seq = $db->seq($name, $st => $end);
return $seq;
}


########################################################################################"
#Process the GCAT sequence
sub processGCAT {
    my $sequence = shift;
    my @letters = split(//, $sequence);
    my $gccount = 0; my $totalcount = 0; my $gccontent = 0;
    my $acount = 0; my $tcount = 0; my $gcount = 0; my $ccount = 0; my $atcontent =0;
    foreach my $i (@letters) {
	if (lc($i) =~ /[a-z]/) { $totalcount++;}
	if (lc($i) eq "g" || lc($i) eq "c") { $gccount++; }
	if (lc($i) eq "a") { $acount++;}
	if (lc($i) eq "t") { $tcount++;}
	if (lc($i) eq "g") { $gcount++;}
	if (lc($i) eq "c") { $ccount++;}
    }
    if ($totalcount > 0) {
	$gccontent = (100 * $gccount) / $totalcount;
    }
    else {
	$gccontent = 0;
    }
    my $others=($totalcount-($acount+$tcount+$gcount+$ccount));
    return ($gccontent,$others,$totalcount,$gcount,$ccount,$acount,$tcount);

}


########################################################################################"
#Print the hash values
sub print_hash_final {
    my ($href,$fhandler)  = @_;
    while( my( $key, $val ) = each %{$href} ) {
        print $fhandler "$key\n";
	#print $fhandler "$key\t=>$val\n";
    }
}


########################################################################################"
# Returns 1 if present else 0
sub isInList {
   my $needle = shift;
   my @haystack = @_;
   foreach my $hay (@haystack) {
	if ( $needle eq $hay ) {
	return 1;
	}
   }
   return 0;
}


########################################################################################"
#Mean of GCAT
use List::Util qw(sum);
sub meanGCAT { return @_ ? sum(@_) / @_ : 0 }
#sub meanGCAT { return sum(@_)/@_; }


########################################################################################"
#  Sigmoid function
sub sigmoidFun {
my $val = shift;
   my ($h) = @_;
   return 1.0 / ( 1.0 + exp(-$val) );
}


########################################################################################"

sub sumArray {
    return( defined $_[0] ? $_[0] + sumArray(@_[1..$#_]) : 0 );
}

########################################################################################"

sub spacer { return (" " x 20000);}


########################################################################################"


sub geneBased {
my ($accession_number, $sequence)=@_;

# next if length($sequence) <= 10000; # Limit the minimum length

	print LOG "$accession_number\n";

	my $tmpf = new File::Temp( UNLINK => 1 );
	print $tmpf ">$accession_number\n$sequence\n";

	my $SStat = &processGCAT($sequence);
	my $seqLen = length($sequence);
	my $gcSStat = (split /\t/, "$SStat")[0];


}

########################################################################################"

=pod
sub round {

    my ($nr,$decimals) = @_;
    return (-1)*(int(abs($nr)*(10**$decimals) +.5 ) / (10**$decimals)) if $nr<0;
    return int( $nr*(10**$decimals) +.5 ) / (10**$decimals);

}
=cut


########################################################################################

sub findreads {

my ($name, $sequence, $increase, $reference, $mismatch, $param) = @_;
my $defLine = $name; my $rev;
my $SeqLength=length($sequence); # Print the sequence length and a newline
#my $pattern = "G{3,5}.{1,7}G{3,5}.{1,7}G{3,5}.{1,7}G{3,5}";
my $pattern = "G{2,5}.{1,7}G{2,5}.{1,7}G{2,5}.{1,7}G{2,5}";

if ($param->{reverse}) {$rev=1;} else {$rev=0}
for (my $aa=0; $aa<=$rev; $aa++) {  ## This loop to check the G4 twice, both for negative and positive.
if ($aa == 1) { #Instead of reverse completement i reverse the pattern
	#$pattern = "T{2,5}.{1,7}T{2,5}.{1,7}T{2,5}.{1,7}T{2,5}";
        $sequence=rcomplement($sequence);
	}

# Find and print out all the exact matches
my $epattern = fpattern($pattern, 0);
my @ematches = mpositions($epattern, $sequence);
if (@ematches) {
   #print "Exact matches:\n";
   pmatches(\@ematches, 0, $sequence, $defLine, $SeqLength, $aa, $reference, $param);
}

# Now do the same, but allow one mismatch
#
if ($mismatch) {
my $omismatch = fpattern($pattern, $mismatch);
my @amatches = mpositions($omismatch, $sequence);
if (@amatches) {
   #print "Matches with one possible mismatch\n";
   pmatches(\@amatches, $mismatch, $sequence, $defLine, $SeqLength, $aa, $reference, $param);
   }
 }
}

########################################################################################"

sub pmatches {
   my ($mref, $diff, $sequence, $defLine, $SeqLength, $aa, $reference, $param) = @_;
   my @matches=@$mref;
   foreach my $match (@matches) {
     my @matchCor = split(/\.\./, $match);

     $matchCor[0] =~ s/[^0-9]//g; $matchCor[1] =~ s/[^0-9]//g;
     my $len=$matchCor[1]-$matchCor[0];

     my $mstring=substr ($sequence, $matchCor[0],$len );

	if ($aa == 0) {
	my $seqFreqG = checkFreq($mstring,"G",$diff, $param);
	my $palRes = checkPal($mstring, $param); my $palDecision;
	if ($palRes) { $palDecision=1;} else {$palDecision=0;}

     #R(everse) and P(lush)
     print INTERMEDIATE "$mstring\t$defLine\t$reference\t$matchCor[0]\t$matchCor[1]\t$diff\t$len\tP\t$seqFreqG\t$palDecision\n";
	}
	elsif ($aa ==1) {
	my $seqFreqT = checkFreq($mstring,"G", $diff, $param);
	my $palRes = checkPal($mstring, $param); my $palDecision;
	if ($palRes) { $palDecision=1;} else {$palDecision=0;}

	print INTERMEDIATE "$mstring\t$defLine\t$reference\t$matchCor[0]\t$matchCor[1]\t$diff\t$len\tR\t$seqFreqT\t$palDecision\n";
	}
   }
}


sub checkFreq {
my ($str, $letter, $mutation, $param) = @_;
my @all; my $finalScore=0; my $longString=0;
while ($str =~ /($letter+)/g) {
push @all, length($1);
	@all = grep { $_ != 1 } @all; #Keep all except 1 (i.e delete 1)

	if (@all) {
	my $max = (sort { $b <=> $a } @all)[0];
	if ($max > $param->{consider}) { $longString=$max/$param->{consider}}
	}

	my $scoreGen = scalar @all/4; #4 is reads
	if ($scoreGen == 1) { $mutation = 0;}
	my $scorePen = $scoreGen + ($mutation * $param->{mutation}) + $longString;
	my $scoreFin = $scorePen;
	$finalScore=sigmoidFun ($scoreFin);
}
$str=join ',', @all;
return "$str\t$finalScore";
}

########################################################################################"
sub get_genome_sequence {
   my ($fpath) = @_;
   open GENOME, "<$fpath" or die "Can't open $fpath: $!\n";
   $_ = <GENOME>; # discard first line
   my $sequence = "";
   while (<GENOME>) {
      chomp($_);                  # Remove line break
      s/\r//;                     # And carriage return
      $_ = uc $_;                 # Make the line upper-case;
      $sequence = $sequence . $_; # Append the line to the sequence
   }
   return $sequence;
}



########################################################################################"
sub mpositions {
   my $pattern;
   local $_;
   ($pattern, $_) = @_;
   my @results;
   local $matchStart;
   my $iPattern = qr/(?{ $matchStart = pos() })$pattern/;

   while (/$iPattern/g) {
      my $nextStart = pos();
      push @results, "[$matchStart..$nextStart)";
      pos() = $matchStart+1;
   }
   return @results;
}


########################################################################################"
sub fpattern {
   my ($opattern, $misallowed) = @_;
   $misallowed >= 0
      or die "Number of mismatches must be greater than or equal to zero\n";
   my $npattern = mapproximate($opattern, $misallowed);
   return qr/$npattern/;
}


########################################################################################"
sub mapproximate {
   my ($pattern, $misallowed) = @_;
   if ($misallowed == 0) { return $pattern }
   elsif (length($pattern) <= $misallowed)
      { $pattern =~ tr/ACTG/./; return $pattern }
   else {
      my ($first, $rest) = $pattern =~ /^(.)(.*)/;
      my $amatch = mapproximate($rest, $misallowed);
      if ($first =~ /[ACGT]/) {
         my $amiss = mapproximate($rest, $misallowed-1);
         return "(?:$first$amatch|.$amiss)";
      }
      else { return "$first$amatch" }
   }
}
}


########################################################################################"
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^[\t\s]+//;
	$string =~ s/[\t\s]+$//;
	$string =~ s/[\r\n]+$//; ## remove odd or bad newline ...
	return $string;
}


########################################################################################"
#Get reverse complement
sub rcomplement {
        my $dna = shift;
	# reverse the DNA sequence
	my $revcomp = reverse($dna);
	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}
#print the lines

########################################################################################"
use Carp qw/croak/;
sub file_write {
        my $file = shift;
        open IO, ">$file" or croak "Cannot open $file for output: $!\n";
        print IO @_;
        close IO;
}

########################################################################################"
#Open and write a file
sub write_fh {
    my $filename = shift @_;
    my $filehandle;
    open $filehandle, ">$filename" or die $!;
    return $filehandle;
}


########################################################################################"
sub pLine {
my $msg = shift;
print LOG "$msg" x 80 . "\n";
#print ($msg x 20);
}


########################################################################################"

sub sorter {
	$a->[1] cmp $b->[1] ||
     $a->[3] <=> $b->[3]
  || $b->[4] <=> $a->[4]
}

########################################################################################"

sub detectStat {
my ($file, $outfile, $score, $sequence_data_ref)=@_;
my @terms;
my %genome=%{$sequence_data_ref};

my $fh = &read_fh($file);
my $out =&write_fh($outfile);
while (<$fh>) { chomp; push @terms, [split /\t/]; }
my $Ids_ref=extractIds(\@terms);
foreach my $id (@$Ids_ref) {
my ($all, %pal, %str);
my $len = $genome{$id}{len};
my $palN=0; my $palY=0; my $strP=0; my $strR=0;

	for my $item (@terms) {
    		if ($item->[1] eq $id) {
			if ($item->[9] >= $score) { # 0 score means all the count
			$pal{$item->[10]}++;
			$str{$item->[7]}++;
			}
		$all++;
    		}
	}
if ($pal{0}) { $palN=$pal{0};}
if ($pal{1}) { $palY=$pal{1};}
if ($str{P}) { $strP=$str{P};}
if ($str{R}) { $strR=$str{R};}
print $out "$id\t$all\t$palN\t$palY\t$strP\t$strR\n";
}

}

########################################################################################"

sub analyticStat {
my ($file, $outfile, $score, $sequence_data_ref)=@_;
my @terms;
my %genome=%{$sequence_data_ref};

my $fh = &read_fh($file);
my $out =&write_fh($outfile);
while (<$fh>) { chomp; push @terms, [split /\t/]; }
my $Ids_ref=extractIds(\@terms);
my $feature_ref=extractFeature(\@terms);

foreach my $id (@$Ids_ref) {
my %feature;
	for my $item (@terms) {
	if ($item->[11] eq $id) { if ($item->[9] >= $score) { $feature{$item->[13]}++; } } #$item->[9] >= 0 for score filter
	}
print $out "$id\t";
my $len = $genome{$id}{len};

foreach my $f(@$feature_ref) {
	if ($feature{$f}) { print $out "$feature{$f}/$len\t"; }
	else { print $out "0\t";} }
print $out "\n";
}

}

########################################################################################"
sub extractFeature {
my ($terms_ref)=@_;
my @allIds;
for my $item (@$terms_ref) {
     if ($item->[13]) { push @allIds , $item->[13];}
}
my @allIds_uniq=uniq(@allIds);
my @allIds_uniq_sorted = sort { lc($a) cmp lc($b) } @allIds_uniq;
return \@allIds_uniq_sorted;
}

########################################################################################"
sub extractIds {
my ($terms_ref)=@_;
my @allIds;
for my $item (@$terms_ref) {
     if ($item->[1]) { push @allIds , $item->[1];}
}
my @allIds_uniq=uniq(@allIds);
return \@allIds_uniq;
}

########################################################################################"
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

########################################################################################"
sub blockoverlapr {
my ($file, $outfile)=@_;
my @terms;

my $fh = &read_fh($file);
my $out =&write_fh($outfile);
while (<$fh>) { chomp; push @terms, [split /\t/]; }

my $biggest = 0;
my $id = '';
for my $term (sort sorter @terms) {
	$biggest = 0 if $id ne $term->[1];
    if ($term->[4] > $biggest) {
        print $out join "\t", @$term; print $out "\n";
        $biggest = $term->[4];
    }
    $id = $term->[1];
}

}

########################################################################################"

sub checkInGFF {
my ($readsFile, $gffFile, $outFile, $exSize)=@_;
my $fh = &read_fh($readsFile);
my $out =&write_fh($outFile);
while (my $row = <$fh>) {
  chomp $row; #print "$row\n";
  my @data = split /\t/, $row;
  my $breakST = $data[3]-$exSize; my $breakED=$data[4]+$exSize;
  my $blkSize = $data[4]-$data[3];
  my $gffIn = checkGFF ($gffFile, $data[1], $breakST, $breakED);

  #my $exSeq = extractSeq($ARGV[4],$data[1], $data[2],$data[3]);
  #my $palRes = checkPal($exSeq);

  my @gffIn=@$gffIn;
	if (@gffIn) { foreach my $line (@gffIn) { print $out "$row\t$line\n"; } }
	else { print $out "$row\tNA\n"; }

	#check the repeats in blocks then
	#my $trfInblk=checkTRF ($ARGV[1], $data[1], $data[2], $data[3]);
	#my @trfInblk=@$trfInblk; if (@trfInblk) {$trfInBlkStr = join ',', @trfInblk;} else {$trfInBlkStr='NA';}
	#my $exSeq = extractSeq($ARGV[4],$data[1], $data[2],$data[3]);
	#my $palRes = checkPal($exSeq);
	#if ($palRes) { $palDecision="Palindromic:$palRes";} else {$palDecision='NotPalindromic';}
	#$color='#181009';

	#my $coverage=checkCov ($ARGV[2], $data[1], $breakST, $breakED);
}
}

# Checks if a provided two coordinates overlaps or not it return 1 if overlaps
sub checkCorOverlaps {
my ($x1, $x2, $y1, $y2)=@_;
return $x1 <= $y2 && $y1 <= $x2;
}

sub checkGFF {
my ($file, $name, $cor1, $cor2)=@_;
my @gffData;
my $fh = &read_fh($file);
while (my $row = <$fh>) {
  	chomp $row;
	my @data = split /\t/, $row;
	next if $name ne $data[0];
	my $res=checkCorOverlaps ($data[3], $data[4], $cor1, $cor2);
	if ($res) {push @gffData, $row}

}
return \@gffData;
}

sub checkCov {
my ($file, $name, $cor1, $cor2)=@_;
my @covData; my $sum=0; my $avCov=0;
open(my $fh, '<:encoding(UTF-8)', $file) or die "Could not open file '$file' $!";
while (my $row = <$fh>) {
  	chomp $row;
	my @data = split /\t/, $row;
	next if $name ne $data[0];
	my $res=checkCorOverlaps ($data[1], $data[2], $cor1, $cor2);
	if ($res) {push @covData, $data[3]}

}
$sum += $_ for @covData;
$avCov=$sum/scalar(@covData);
return $avCov;
}

sub extractSeq {
my ($file, $chr, $st, $ed)=@_;
use Bio::DB::Fasta;
my $db = Bio::DB::Fasta->new($file);
my $seq = $db->seq($chr, $st => $ed);
return $seq;
}

sub checkPal {
my ($seq, $param)=@_;
my $pp = qr/(?: (\w) (?1) \g{-1} | \w? )/ix;
    while ($seq =~ /(?=($pp))/g) {
        return "$-[0] - $1" if length($1) > $param->{palsize};
    }
}


#Plot the DBJ graph
sub dbjGraph {
    my ( $k, $name, $text, $outFile) = @_;
    my %graph;
    my $out =&write_fh($outFile);
    foreach my $i ( 0 .. ( length $text ) - $k ) {
        my $kmer = substr $text, $i, $k;
        my $prefix = substr $kmer, 0, -1;
        my $suffix = substr $kmer, 1;
        push @{ $graph{$prefix} }, $suffix;
    }
    #return \%graph;
    foreach my $prefix ( sort keys %graph ) {
	my $mers=join ',' , sort @{$graph{$prefix}};
	print $out "$prefix\t$mers\n";
    	#printf "%s -> %s\n", $prefix, join q{,}, sort $graph{$prefix};
    }
}

#Store SPS detail
sub storeSPS {
my ($spsFile, $spsFile2, $param) = @_;
my @SpsArray; my $SpsNumber; my @SpsArray2; my $SpsNumber2;
open SPSFILE, "$spsFile" or die "cant open .sps file $!";
while (<SPSFILE>) {
	my $SpsLine=$_; chomp $SpsLine; @SpsArray=split /,/, lc($SpsLine);  $SpsNumber = scalar (@SpsArray); }
	my $SpsArrayTabed=join("\t", @SpsArray);
close SPSFILE or die "could not close file: $!\n";

open SPSFILE2, "$spsFile2" or die "cant open .sps file $!";
while (<SPSFILE2>) {
	my $SpsLine2=$_; chomp $SpsLine2; @SpsArray2=split /,/, lc($SpsLine2);  $SpsNumber2 = scalar (@SpsArray2); }
	my $SpsArrayTabed2=join("\t", @SpsArray2);
close SPSFILE2 or die "could not close file: $!\n";

my @outArray = keys %{{map {($_ => 1)} (@SpsArray, @SpsArray2)}}; #Combine both array by removing duplicates

return (\@outArray, $SpsNumber, $SpsArrayTabed);
}

##Recontruct the breakpoints
sub reconstructTar {
my ($finalEBA, $allHSB, $threshold, $length, $spsFile, $SpsNumber, $refName, $param)=@_;

my $version=0.1;

if (-d "$param->{out_dir}/output_$refName") {
deldir("$param->{out_dir}/output_$refName"); # or deldir($ARGV[0]) to make it commandline
} else { mkdir "$param->{out_dir}/output_$refName"; }

if (-f "Reconstruction_$refName.stats") { unlink "Reconstruction_$refName.stats";}
my $InFile=$finalEBA; #final_classify.eba file
print "$InFile\n";
my $threshold=$threshold; # threhold value to filter

#Store species
my @SpsArray; my $SpsNumber;
open SPSFILE, "$spsFile" or die $!;
if (-f "Reconstruction_$refName.stats") { unlink "Reconstruction_$refName.stats";}
while (<SPSFILE>) {
	my $SpsLine=$_; chomp $SpsLine; @SpsArray=split /,/, lc($SpsLine);  $SpsNumber = scalar (@SpsArray); } ## It read the species names from sps.txt file ... need to improve !!!
	my $SpsArrayTabed=join("\t", @SpsArray);

close SPSFILE or die "could not close file: $!\n";


foreach my $spsName(@SpsArray) {
my $outFile1="$param->{out_dir}/output_$refName/$spsName"."_brk_$refName.tar1";
my $outfile2="$param->{out_dir}/output_$refName/$spsName"."_brk_$refName.tar2";
my $outfile3="$param->{out_dir}/output_$refName/$spsName"."_brk_$refName.tar3";

open (OUTFILE1, ">$outFile1") or die "$0: open $outFile1: $!";
open (OUTSTAT, ">>$param->{out_dir}/output_$refName/Reconstruction_$refName.stats") or die "$0: open $param->{out_dir}/output_$refName/Reconstruction_$refName.stats: $!";

	open INFILE,  $InFile or die "$0: open $InFile: $!";
	my @array; my @index; my @nameArray; my $in; my $countReal; my @done; my $countBrk; my $countGap; my $total;
	while (<INFILE>) {
		chomp;
		my $line=trim($_);
		my @tmp = split /\t/, lc ($line);
		if ($. == 1) {
			@nameArray = split /\t/, $line;
			(@index)= grep { $nameArray[$_] eq "$spsName" } 0..$#nameArray;
			$in=$index[0];
			next;
		}
		next if $_ =~ /^\s*#/;
		next if !$tmp[$in];

		my @val = split /\,/, $tmp[$in];

		# In case more than one breakpoints ( separated with comma)
		foreach my $l(@val) {
		my $line=trim($l);

		my $seenIn=isInList ($l, @done);
		if (!$seenIn) { $countReal++; }
		if (index($l, "breakpoints") != -1) { $countBrk++; } elsif (index($l, "gap") != -1) { $countGap++; }
		$total++;

		print OUTFILE1 "$l\t$tmp[$SpsNumber+1]\t$tmp[$SpsNumber+7]\t$tmp[$SpsNumber+9]\t$tmp[$SpsNumber+10]\t$tmp[$SpsNumber+11]\t$tmp[$SpsNumber+12]\t$tmp[$SpsNumber+13]\n";

		push @done, $l;
		#push (@array,$_);
		#72412203--72417432=Breakpoints+0.925	Breakpoints	chicken:6.22711733452034e-07	237.8835341	0.05	20	1	18
		}
	}
	if (-z "Reconstruction_$refName.stats") { print OUTSTAT "spsName\tcountReal\tcountBrk\tcountGap\ttotal\n";}
	print OUTSTAT "$spsName\t$countReal\t$countBrk\t$countGap\t$total\n";
	undef @done;
	close INFILE or die "could not close $InFile file: $!\n";
close OUTSTAT or die "could not close OUTSTAT file: $!\n";
close OUTFILE1 or die "could not close $outFile1 file: $!\n";

reconstructBrk($outFile1, $outfile2, $spsName, $outfile3, $threshold, $allHSB, $param);
undef @index;
}

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub isInList {
	my $needle = shift;
	my @haystack = @_;
	foreach my $hay (@haystack) {
		if ( $needle eq $hay ) {
			return 1;
		}
	}
	return 0;
}


sub deldir {
  my $dirtodel = pop;
  my $sep = '\\'; #change this line to "/" on linux.
  opendir(DIR, $dirtodel);
  my @files = readdir(DIR);
  closedir(DIR);

  @files = grep { !/^\.{1,2}/ } @files;
  @files = map { $_ = "$dirtodel$sep$_"} @files;

  @files = map { (-d $_)?deldir($_):unlink($_) } @files;

  rmdir($dirtodel);
}

}

# Reconstruct the breaks
sub reconstructBrk {
my ($InFile, $OutFile, $spsName, $out, $threshold, $HSBInFile, $param) = @_;

open INFILE,  $InFile or die "$0: open $InFile: $!";
open (OUTFILE, ">$OutFile") or die "$0: open $OutFile: $!";
	my %inArray;
	while (<INFILE>) {
		chomp;
		if ($_ =~ /^\s*#/) { next; }
		$_=trim($_);
		my @array = split /\t/, lc($_);
		for (@array) { s/^\s+//; s/\s+$//;}
		my @array2 = split /\=/, lc($array[0]);
		my @array3 = split /\-\-/, lc($array2[0]);
		my $pahala=trim($array3[0]);
		my $dusara=trim($array3[1]);
		#next if $array[8] ne $species;
		#print "$spsName:$array[1]:$pahala\n";
		$inArray{"$spsName:$array[1]:$pahala"}=$_; #Here still the reference chr an coords
		$inArray{"$spsName:$array[1]:$dusara"}=$_;

	}
	close INFILE or die "could not tar2 file: $!\n";

	open HSBFILE,  $HSBInFile or die "$0: open HSB file $HSBInFile: $!";
	my @HSBArray;
	while (<HSBFILE>) {
		chomp;
		if ($_ =~ /^\s*#/) { next; }
		#chicken:100K	12	375235	955085	14	14785	688127	-	Meleagris_gallopavo	Chromosome
		my @array = split /\t/, lc($_);
		for (@array) { s/^\s+//; s/\s+$//;}
		next if lc($array[8]) ne $spsName; ## To print only required species
		#print "$array[8]:$array[1]:$array[2]\tdusarawala\n";
		#print "WARNING: If you have used extended breakpoint approach in EBA, this might not work\n";
		my $a1=$inArray{"$array[8]:$array[1]:$array[2]"}; #first HSB coord
		my $a2=$inArray{"$array[8]:$array[1]:$array[3]"}; #second HSB coord
		#print $a1;
		if(!$a1) { $a1="Tel\t\t\t\t\t\t\t";} if(!$a2) { $a2="Tel\t\t\t\t\t\t\t";}
my @aa = split /\=/, lc($a1); #get the coord and brk info
my @bb = split /\=/, lc($a2);
my @aa2 = split /\+/, $aa[1]; #get the brk info
my @bb2 = split /\+/, $bb[1];

		#Relying on the HSB orientation
		if($array[7] eq "-") { print OUTFILE "$_\t$a2\t$bb2[0]\t$a1\t$aa2[0]\n";} else { print OUTFILE "$_\t$a1\t$aa2[0]\t$a2\t$bb2[0]\n";}
	}
	close HSBFILE or die "could not close $HSBInFile file: $!\n";
close OUTFILE or die "could not close $OutFile file: $!\n";

reconstructTarget($OutFile, $out, $spsName, $threshold, $param);
}



# Reconstruct the breaks
sub reconstructBrk2 {
my ($InFile, $OutFile, $spsName, $out, $threshold, $HSBInFile, $param, $spsNum) = @_;
open (OUTFILE, ">>$OutFile") or die "$0: open $OutFile: $!";
open INFILE,  $InFile or die "$0: open $InFile: $!";
	while (<INFILE>) {
		chomp;
		if ($_ =~ /^\s*#/) { next; }
		my $line=trim($_);
		my @array = split /\t/, lc($line);
		for (@array) { s/^\s+//; s/\s+$//;}
		my @array2 = split /\=/, lc($array[0]);
		my @array3 = split /\-\-/, lc($array2[0]);
		my $pahala=trim($array3[0]);
		my $dusara=trim($array3[1]);

		my $firstC = findAncestralCor($HSBInFile, $spsName, $array[1], $pahala, 1);
		my $secondC = findAncestralCor($HSBInFile, $spsName, $array[1], $dusara, 2);
		
		#In case not overlaps in HSB
		if (!$firstC) { $firstC = "NA:NA:0\tNA:NA:0"; }
		elsif (!$secondC) { $secondC = "NA:NA:0\tNA:NA:0"; }

		print OUTFILE "$firstC\t$line\n";
		print OUTFILE "$secondC\t$line\n";
		
	}
	close INFILE or die "could not tar2 file $InFile : $!\n";
	close OUTFILE or die "could not close hhh $OutFile file: $!\n";
checkOverlapsInRef("$param->{ref2}/$param->{resolution}/EBA_OutFiles/final_classify.eba7", $OutFile, $out, $param, $spsNum);
}

sub checkOverlapsInRef {
my ($finaRefFile, $inF, $OutFile, $param, $spsNum) = @_;
open INFILE1,  $inF or die "$0: open $inF: $!";
open (OUTFILE, ">>$OutFile") or die "$0: open $OutFile: $!";
my $extend=1000;
while (<INFILE1>) {
chomp;
if ($_ =~ /^\s*#/) { next; }
my $line1=trim($_);
my @array1 = split /\t/, lc($line1);
my @fDetail = split /\:/, $array1[0];
my @sDetail = split /\:/, $array1[1];

	open INFILE2,  $finaRefFile or die "$0: open $finaRefFile: $!";
	while (<INFILE2>) {
	chomp;
	if ($_ =~ /^\s*#/) { next; }
	my $line=trim($_);
	my @array = split /\t/, lc($line);
	next if $array[12] ne $fDetail[1];
	my @finalDetailCor = split /\--/, $array[$spsNum+2];
	my $extST= $fDetail[2]-$extend; my $extEd= $sDetail[2]+$extend;
	my $ORes = checkOverlaps($finalDetailCor[0], $finalDetailCor[1], $extST, $extEd);
	if ($ORes) {
	print OUTFILE "$fDetail[1]\t$fDetail[2]\t$sDetail[2]\t$extST\t$extEd\t$line1\t$_\n";
	} }
close INFILE2 or die "could not close hmm $finaRefFile file: $!\n";
}
close OUTFILE or die "could not close hhh $OutFile file: $!\n";
close INFILE1 or die "could not close hmm $finaRefFile file: $!\n";
}

#Find ancestral coordiantes
sub findAncestralCor {
my ($HSBInFile, $spsName, $spsChr, $spsSt, $nn) = @_;
open HSBFILE,  $HSBInFile or die "$0: open HSB file $HSBInFile: $!";
my $corVal="NA\tNA"; my $flag=0;
	while (<HSBFILE>) {
		chomp;
		if ($_ =~ /^\s*#/) { next; }
		#chicken:100K	12	375235	955085	14	14785	688127	-	Meleagris_gallopavo	Chromosome
		my @array = split /\t/, lc($_);
		for (@array) { s/^\s+//; s/\s+$//;}
		next if lc($array[8]) ne $spsName; ## To print only required species
		#Loop over to check the nearby ... as these ancetral breakpoint coorindates 
		if (( $spsSt >= $array[2]  and   $spsSt <= $array[3] ) and  ($array[1] eq $spsChr)) { 
		$flag=1;
		
		#print "$spsChr\t$spsSt ---- $_\n";
		#Sometime the size go beyond HSB blocks
		my $localDist=$spsSt - $array[2];
		my $localCor = $array[5] + $localDist;

		my $extCor = 20000;
		my $secC= $localCor + $extCor;
		my $firC= $localCor - $extCor;

		my $cName='NA'; my $sName='NA';
		$sName=$array[8];
		$cName=$array[4];
		if ($nn == 1) {
		$corVal = "$sName:$cName:$localCor\t$sName:$cName:$secC";
		}
		elsif ($nn == 2) {
		$corVal = "$sName:$cName:$firC\t$sName:$cName:$localCor";
		}
		return $corVal;
		}
	}
	close HSBFILE or die "could not close hmm $HSBInFile file: $!\n";
return 0 if !$flag;
}

sub reconstructTarget {
my ($InFile, $OutFile, $spsName, $threshold, $param)=@_;
#threashold is for first and second best classification score differences

open (OUTFILE, ">$OutFile") or die "$0: open $OutFile: $!";
	open INFILE,  $InFile or die "$0: open $InFile: $!";
	my @array;
	while (<INFILE>) {
		chomp;
		#next if $. == 1; ## to next the header -- no header in tar2
		if ($_ =~ /^\s*#/) { next; }
		$_=trim($_);
		push (@array,$_); #Keep all  the tar2 in an array

	}
	close INFILE or die "could not close $InFile file: $!\n";
	#Sort this tar2 with target cordinate start and end
	my @sorted_array = sort { (split "\t", $a)[4] cmp (split "\t", $b)[4] || (split "\t", $a)[5] <=> (split "\t", $b)[5] && (split "\t", $a)[6] <=> (split "\t", $b)[6] } @array;

print OUTFILE "chr\tstart\tend\tfirstBestratio\tsecondBestRatio\tPercentageUsed\tBrk\tclass\tsameordiffClass\tsameordiffBrk\tborderUsed\tREFchr\tREFstart\tREFchr\tREFend\tsecondClass\n";

if (!$param->{strict}) {

	for (my $val=0; $val <= $#sorted_array; $val++ ) {
		my @array_first = split /\t/, $sorted_array[$val];
		my @array_next = split /\t/, $sorted_array[$val+1];

		if ($array_first[4] eq $array_next[4]) { #same chr
			my $brkDecision; my $newScore; my $classNum; my $borNum; my $secClass;

	#NOTE: first and last HSB coordinate will be ignored.
	#taegut:100k	1a	61261345	62559310	1	68978907	70838383	+	anas_platyrhynchos	chromosomes	61258206--61261345=breakpoints+0.909090909090909	1a	taeniopygia_guttata:0.903687300801826	7221.5850129769	0.1818181818	11	2	9	breakpoints	62559310--62561208=breakpoints+0.863636363636364	1a	taeniopygia_guttata:0.904088631832987	3649.1507838237	0.2727272727	11	3	8	breakpoints


	#taegut:100k	19	3015453	3143368	20	5478	207928	+	anas_platyrhynchos	chromosomes	3008057--3015453=breakpoints+0.863636363636364	19	taeniopygia_guttata:0.91030772058792	398.1409945019	0.2727272727	11	3	8	breakpoints	3143368--3164683=breakpoints+0.681818181818182	19	both_finches_crow:0.92773464582499	1	0.2727272727	11	3	6	breakpoints
	#taegut:100k	19	87441	3008057	20	220685	2931079	+	anas_platyrhynchos	chromosomes	Tel									3008057--3015453=breakpoints+0.863636363636364	19	taeniopygia_guttata:0.91030772058792	398.1409945019	0.2727272727	11	3	8	breakpoints


			my @className1=split /\:/, $array_first[21];
			my @className2=split /\:/, $array_next[12];
			my $secBor; my $secBor2;
			$classNum="singleClass"; $brkDecision="singleBrk"; $borNum=1; $secBor='NA';

			if (($array_first[27] eq "breakpoints") and ($array_next[18] eq "breakpoints")) { #If both border are breakpoints
			$brkDecision="sameBrk"; $borNum=2;
			if (($array_first[22] >= $array_next[13]) and ($array_first[22] >= $threshold) ) {
					$secBor=$array_next[13]; $secClass=$array_next[12];
					if($className1[0] eq $className2[0]) { $classNum="sameClass";}  else {  $classNum="diffClass";}
					$newScore="$array_first[22]\t$secBor\t$array_first[26]\t$array_first[27]\t$array_first[21]";
				}
			elsif (($array_first[22] <= $array_next[13]) and ($array_next[13] >= $threshold) ) {
					$secBor=$array_first[22]; $secClass=$array_first[21];
					if($className1[0] eq $className2[0]) { $classNum="sameClass";}  else {  $classNum="diffClass";}
					$newScore="$array_next[13]\t$secBor\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}
			elsif (($array_first[22] == $array_next[13]) and ($array_next[13] <= $threshold) ) { #since both are equal ... do not nned to check other brk cor
					$secBor=$array_next[13]; $secClass=$array_next[12];
					if($className1[0] eq $className2[0]) { $classNum="sameClass";}  else {  $classNum="diffClass";}
					$newScore="$array_first[22]\t$secBor\t$array_first[26]\t$array_first[27]\t$array_first[21]";
				}
			}


			elsif (($array_first[27] eq "breakpoints") and ($array_next[18] eq "gap")) { #If sec is gap
			if ($array_first[22] >= $threshold) {
					$secBor=$array_next[13]; $secClass=$array_next[12];
					$classNum="diffClass";
					$newScore="$array_first[22]\t$secBor\t$array_first[26]\t$array_first[27]\t$array_first[21]";
				}
			elsif ($array_first[22] <= $threshold) {
					$secBor=$array_next[13]; $secClass=$array_next[12];
					$classNum="diffClass";
					$newScore="$array_first[22]\t$secBor\t$array_first[26]\t$array_first[27]\t$array_first[21]";
				}
			}


			elsif (($array_first[27] eq "breakpoints") and (!$array_next[18])) { #If sec is TEL empty
			if ($array_first[22] >= $threshold) {
					$secBor=$array_next[13]; $secClass="NA";
					$classNum="diffClass";
					$newScore="$array_first[22]\t$secBor\t$array_first[26]\t$array_first[27]\t$array_first[21]";
				}
			elsif ($array_first[22] <= $threshold) {
					$secBor=$array_next[13]; $secClass="NA";
					$classNum="diffClass";
					$newScore="$array_first[22]\t$secBor\t$array_first[26]\t$array_first[27]\t$array_first[21]";
				}
			}

			elsif (($array_first[27] eq "gap") and (!$array_next[18])) { #If sec is TEL empty
			if ($array_first[22] >= $threshold) {
					$secBor=$array_next[13]; $secClass="NA";
					$classNum="diffClass";
					$newScore="$array_first[22]\t$secBor\t$array_first[26]\t$array_first[27]\t$array_first[21]";
				}
			elsif ($array_first[22] <= $threshold) {
					$secBor=$array_next[13]; $secClass="NA";
					$classNum="diffClass";
					$newScore="$array_first[22]\t$secBor\t$array_first[26]\t$array_first[27]\t$array_first[21]";
				}
			}



			elsif (($array_first[27] eq "gap") and ($array_next[18] eq "breakpoints")) { #If first is gap
			if ($array_next[13] >= $threshold) {
					$secBor=$array_first[22]; $secClass=$array_first[21];
					$classNum="diffClass";
					$newScore="$array_next[13]\t$secBor\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}
			elsif ($array_next[13] <= $threshold) {
					$secBor=$array_first[22]; $secClass=$array_first[21];
					$classNum="diffClass";
					$newScore="$array_next[13]\t$secBor\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}
			}

			elsif ((!$array_first[27]) and ($array_next[18] eq "breakpoints")) { #If first is TEL empty
			if ($array_next[13] >= $threshold) {
					$secBor=$array_first[22]; $secClass="NA";
					$classNum="diffClass";
					$newScore="$array_next[13]\t$secBor\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}
			elsif ($array_next[13] <= $threshold) {
					$secBor=$array_first[22]; $secClass="NA";
					$classNum="diffClass";
					$newScore="$array_next[13]\t$secBor\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}
			}

			elsif ((!$array_first[27]) and ($array_next[18] eq "gap")) { #If first is TEL empty
			if ($array_next[13] >= $threshold) {
					$secBor=$array_first[22]; $secClass="NA";
					$classNum="diffClass";
					$newScore="$array_next[13]\t$secBor\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}
			elsif ($array_next[13] <= $threshold) {
					$secBor=$array_first[22]; $secClass="NA";
					$classNum="diffClass";
					$newScore="$array_next[13]\t$secBor\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}
			}

			elsif (($array_first[27] eq "gap") and ($array_next[18] eq "gap")) { #both brk cords are GAP

			#Dont need to check it at all

			if ($array_next[13] >= $threshold) {
					$secBor=$array_first[22]; $secClass=$array_first[21];
					$classNum="diffClass";
					$newScore="$array_next[13]\t$secBor\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}
			elsif ($array_next[13] <= $threshold) {
					$secBor=$array_first[22]; $secClass=$array_first[21];
					$classNum="diffClass";
					$newScore="$array_next[13]\t$secBor\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}

			}

		else {
		#print "$array_first[27]  $array_next[18] -->>\n";
		$newScore="\t\t\t\t"; $classNum="NA"; $brkDecision="NA"; $borNum=0; # If this means both are empty
		undef $array_first[1]; undef $array_first[3]; undef $array_next[1]; undef $array_next[2]; undef $secClass;
		next;
		}


		print OUTFILE "$array_first[4]\t$array_first[6]\t$array_next[5]\t$newScore\t$classNum\t$brkDecision\t$borNum\t$array_first[1]\t$array_first[3]\t$array_next[1]\t$array_next[2]\t$secClass\n";

		}

		else	{ #the end of the chromosomes are here -- mostly one brk coordinate is missing
	#print "$array_first[27]  $array_next[18] <<--\n";
				if ($array_first[27]) {
				#print OUTFILE "$array_first[4]\t$array_first[6]\tEND\t$array_first[22]\t$array_first[13]\t$array_first[26]\t$array_first[27]\t$array_first[21]\tNA\tNA\tNA\tNA\t$array_first[1]\t$array_first[3]\tNA\tNA\n";
				}

				if ($array_next[18]) {
				#print OUTFILE "$array_next[4]\tEND\t$array_next[5]\t$array_next[13]\t$array_next[22]\t$array_next[17]\t$array_next[18]\t$array_next[12]\tNA\tNA\tNA\tNA\tNA\tNA\t$array_next[1]\t$array_next[2]\n";
				}
			}
		}

	#tar3 print this
	#chr	start	end	firstBestratio	secondBestRatio	PercentageUsed	Brk	class	sameordiffClass	sameordiffBrk	borderUsed	REFchr	REFstart	REFchr	REFend
	#1	24752130	24752143	3.8496540738	2.3664500767	2	breakpoints	anolis_carolinensis:0.151936462935248	diffClass	sameBrk	2	5	36755557	1a	23507243

	close OUTFILE or die "could not close $OutFile file: $!\n";
	}


#------

elsif ($param->{strict}) { # this is old strict function

	for (my $val=0; $val <= $#sorted_array; $val++ ) {
		my @array_first = split /\t/, $sorted_array[$val];
		my @array_next = split /\t/, $sorted_array[$val+1];

	if($array_first[4] eq $array_next[4]) {
			my $brkDecision; my $newScore; my $classNum; my $borNum;

			my @className1=split /\:/, $array_first[21];
			my @className2=split /\:/, $array_next[12];
			my $secBor; my $secBor2;

			if(($array_first[22] >= $array_next[13]) and ($array_first[27] eq "breakpoints") and ($array_first[22] >= $threshold) and ($array_next[18] eq "breakpoints")) {

					if (($array_first[22] >= $threshold) and ($array_next[13] >= $threshold)){
						if($className1[0] eq $className2[0]) { $classNum="sameClass";}  else {  $classNum="diffClass";}
						$brkDecision="sameBrk"; $borNum=2;
						$secBor=$array_next[13];
						}
					else { $classNum="singleClass"; $brkDecision="singleBrk"; $borNum=1; $secBor='NA';}

					$newScore="$array_first[22]\t$secBor\t$array_first[26]\t$array_first[27]\t$array_first[21]";

				}
			elsif(($array_first[22] <= $array_next[13]) and ($array_next[18] eq "breakpoints") and ($array_next[13] >= $threshold) and ($array_first[27] eq "breakpoints")) {

					if (($array_first[22] >= $threshold) and ($array_next[13] >= $threshold)){
						if($array_next[13] >= $threshold) { $classNum="sameClass";}  else {  $classNum="diffClass";}
						$brkDecision="sameBrk"; $borNum=2;
						$secBor2=$array_first[22];
						}
					else { $classNum="singleClass"; $brkDecision="singleBrk"; $borNum=1; $secBor2='NA';}

					$newScore="$array_next[13]\t$secBor2\t$array_next[17]\t$array_next[18]\t$array_next[12]";
				}
			elsif(($array_first[22] >= $array_next[13]) and ($array_first[27] eq "breakpoints") and ($array_first[22] >= $threshold) and ($array_next[18] eq "gap")) {
					$newScore="$array_first[22]\tNA\t$array_first[26]\t$array_first[27]\t$array_first[21]";
					$classNum="singleClass";
					$brkDecision="singleBrk";
					$borNum=1;

				}
			elsif(($array_first[22] <= $array_next[13]) and ($array_next[18] eq "breakpoints") and ($array_next[13] >= $threshold) and ($array_first[27] eq "gap")) {
					$newScore="$array_next[13]\tNA\t$array_next[17]\t$array_next[18]\t$array_next[12]";
					$classNum="singleClass";
					$brkDecision="singleBrk";
					$borNum=1;
				}
			elsif(($array_first[22] >= $array_next[13]) and ($array_first[27] eq "breakpoints") and ($array_first[22] >= $threshold) and (!$array_next[18])) {
					$newScore="$array_first[22]\tNA\t$array_first[26]\t$array_first[27]\t$array_first[21]";
					$classNum="singleClass";
					$brkDecision="singleBrk";
					$borNum=1;
				}

			elsif(($array_first[22] <= $array_next[13]) and ($array_next[18] eq "breakpoints") and ($array_next[13] >= $threshold) and (!$array_first[27])) {
					$newScore="$array_next[13]\tNA\t$array_next[17]\t$array_next[18]\t$array_next[12]";
					$classNum="singleClass";
					$brkDecision="singleBrk";
					$borNum=1;
				}
			elsif(($array_first[22] >= $array_next[13]) and ($array_next[18] eq "breakpoints") and ($array_next[13] >= $threshold) and ($array_first[27] eq "gap")) {
					$newScore="$array_next[13]\tNA\t$array_next[17]\t$array_next[18]\t$array_next[12]";
					$classNum="singleClass";
					$brkDecision="singleBrk";
					$borNum=1;
				}
			elsif(($array_first[22] <= $array_next[13]) and ($array_first[27] eq "breakpoints") and ($array_first[22] >= $threshold) and ($array_next[18] eq "gap")) {
					$newScore="$array_first[22]\tNA\t$array_first[26]\t$array_first[27]\t$array_first[21]";
					$classNum="singleClass";
					$brkDecision="singleBrk";
					$borNum=1;
				}

			else {
				$newScore="\t\t\t\t";
				$classNum="NA";
				$brkDecision="NA";
				$borNum=0;
				}
			#}

			print OUTFILE "$array_first[4]\t$array_first[6]\t$array_next[5]\t$newScore\t$classNum\t$brkDecision\t$borNum\t$array_first[1]\t$array_first[3]\t$array_next[1]\t$array_next[2]\n";
			}
		else	{
				#if(($array_first[27] eq "breakpoints") and ($array_first[22] >= $threshold)){
				#print "$array_first[4]\t$array_first[6]\tTel\t$array_first[22]\t$array_first[13]\t$array_first[26]\t$array_first[27]\t$array_first[21]\tNA\tNA\tNA\n";
				#}
				#else {
				#print "$array_first[4]\t$array_first[6]\tTel\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
				#}

			}
		}

	close OUTFILE or die "could not close $OutFile file: $!\n";
	}
	else { print "strict parameter is wrong\n Check config file\n"; }

} #---- sub end


#Concat files
sub concatAll {
my ($outfile, $location) = @_;
open my $OUT, '>', "$outfile" or die "$outfile: $!";

for my $file (glob "$location/*.txt") {
    open my $FH, '<', $file or die "$file: $!";
    while (<$FH>) {
        print {$OUT} $_ unless 0 == $.;
    }
}
close $OUT or die $!;

}


#check the TAR2TAR overlaps
sub checkOverlapsTAR {

my ($file_ref, $file_tar, $outfile, $refName, $spsNum, $length, $allall, $spsFile, $allHSB, $statfile, $outfile2, $identify, $dataFinal, $missedFinal, $groupClass, $chekerExtend, $finalEBA, $cntSTATFile) = @_;

open cntSTAT, ">>$cntSTATFile" or die $!;

open(my $infh, '<:encoding(UTF-8)', $groupClass) or die "Could not open file '$groupClass' $!";
my %classHash;
while (<$infh>) {
	chomp $_; my $line=$_;
	$line=trim($line);
	if (index($_,"#") == 0) { next; } # Lines starting with a hash mark are comments
	next if $line =~ /^\s*$/;
	my @tmpLine=split /\=/, $line;
	#next if $tmpLine[1] eq "";
	$classHash{$tmpLine[0]}=$tmpLine[1];
    }
close $infh;

# Store species
open SPSFILE, "$spsFile" or die $!;
my @SpsArray; my $SpsNumber;
while (<SPSFILE>) {
	my $SpsLine=$_; chomp $SpsLine; @SpsArray=split /,/, lc($SpsLine);  $SpsNumber = scalar (@SpsArray); } ## It read the species names from sps.txt file ... need to improve !!!
	my $SpsArrayTabed=join("\t", @SpsArray);

close SPSFILE or die "could not close file: $!\n";

# Use need to provide the location of all.hsb (of reference used in commandline).
my $allHSBFile = $allHSB;
#Use need to provide the location of all_all.eba00 (of target used in commandline).
my $allALLeba00 = $allall;

# User can find the above files in EBA tools output dir

my $filename = "$file_ref";
my $filename3 = "$outfile";
open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";
open(my $fh3, '>:encoding(UTF-8)', $filename3) or die "Could not open file '$filename3' $!";

my %count; my %countONE; my %countOTHER;my %countTarClass; my $allCount=0;my %allClass; my $flag=0; my @allCor; my @tmpAll;
my $sameClass_sameBrk=0; my $diffClass_sameBrk=0; my $singleClass_singleBrk=0; my $NA=0;
my @missinBrk;
my %done;
my @remain;
my $head1;
my $head2;

while (my $row = <$fh>) {
chomp $row;
if ($. == 1) { $head1= $row; next;}
my $lineNumR1=$.;
my @tmp1 = split('\t', $row);
# no warnings; # To avoid warnings for blank spaces
my @cla;
if ($tmp1[7]) { @cla =split('\:', $tmp1[7]);}
next if !$cla[0];
next if $cla[0] ne $refName; # group or species name provided by user
$countTarClass{$cla[0]}++;
push @allCor, "$tmp1[0]:$tmp1[1]:$tmp1[2]";

if ($tmp1[8] eq "sameClass" and $tmp1[9] eq "sameBrk") { $sameClass_sameBrk++; }
elsif ($tmp1[8] eq "diffClass" and $tmp1[9] eq "sameBrk") { $diffClass_sameBrk++; }
elsif ($tmp1[8] eq "singleClass" and $tmp1[9] eq "singleBrk") { $singleClass_singleBrk++; }
else { $NA++; }

my $filename2 = "$file_tar";
my $doneLine=0;
open(my $fh2, '<:encoding(UTF-8)', $filename2) or die "Could not open file '$filename2' $!";

 	while (my $row2 = <$fh2>) {
	chomp $row2;
	if ($. == 1) { $head2= $row2; next;}
	my $lineNumR2=$.;
	my @values2 = split('\t', $row2);
	my $spsNum=$spsNum;
	my $corStr="$values2[1]--$values2[2]";
	my @cor =split('\--', $corStr);
	next if !$values2[7];
	my @class =split('\:', $values2[7]);
	if ($flag==0) {$allClass{$class[0]}++;}

	my $changeSize=$length;

	my $newTarSt=$cor[0]-$changeSize;
	my $newTarEd=$cor[1]+$changeSize;
	my $newRefSt=$tmp1[1]-$changeSize;
	my $newRefEd=$tmp1[2]+$changeSize;

#print "$newTarSt,$newTarEd,$newRefSt,$newRefEd\n";

	if ($values2[0] eq $tmp1[0]) {
		my $OverRes = checkOverlaps($newTarSt,$newTarEd,$newRefSt,$newRefEd);
		if (($OverRes) and (!isInList("$tmp1[0]:$tmp1[1]:$tmp1[2]", @tmpAll))) { #To remove the duplicate hits
  		 	$count{$class[0]}++;
			$allCount++;
			#push @tmpAll, "$values2[0]:$values2[1]:$values2[2]";
			push @tmpAll, "$tmp1[0]:$tmp1[1]:$tmp1[2]"; ####
			#If anyone of the ratio is one ... count
			no warnings;
			if ( -z "$filename3" ) { print $fh3 "$head1\t\t$head2\n"; }
			if (($values2[3] == 1) or ($values2[4] == 1)) {$countONE{$class[0]}++;} else {$countOTHER{$class[0]}++;}
			print $fh3 "$row\t\t$row2\n";
			$doneLine=1;
			$done{$row2}=$lineNumR2;
     			}
		}
	if ($class[0] eq $refName) { push @remain, $row2;}
	}
$flag=1; # To check the final_classify.eba7 once

close $fh2;
if (!$doneLine) { push @missinBrk, $row;}
}

close $fh;
close $fh3;

my $multiHits="$countTarClass{$refName}:$allCount";
my $uniqueRecoEBRs= uniq (@allCor);

#Only values
print cntSTAT "$identify\t$countTarClass{$refName}\t$allCount\n";

#Create a file for detail stats
open(my $outfh4, '>:encoding(UTF-8)', $outfile2) or die "Could not open file '$outfile2' $!";

my $message = <<"END_MSG";

..........
Hello USER, ...........$refName analysis ..................... SEE NEXT TABLE BELOW?
..........

END_MSG


## Printing the TAR reconstructed Info
print $outfh4 "Overview of reconstructed table detail\n";
print $outfh4 "Species/Group Name \t Total Number of reconstructed EBRs \t Overlaps Count \t Unique reconstructed EBRs \t Duplicates reconstructed EBRs\t Tar\t sameClass_sameBrk \t diffClass_sameBrk \t singleClass_singleBrk \t Others\n";

foreach my $class (sort keys %countTarClass) {  # This will be only one as user check for indivisual
my $difference=($countTarClass{$class}-$uniqueRecoEBRs);
print $outfh4 "$class\t$countTarClass{$class}\t$multiHits\t$uniqueRecoEBRs\t$difference\tTarRecon\t$sameClass_sameBrk\t$diffClass_sameBrk\t $singleClass_singleBrk \t $NA\n"; }

print $outfh4 $message;

print $outfh4 "Reconstructed EBRs overlapping details ( for $file_ref : Used as reference)\n";
print $outfh4 "Species/Group Name \t Total Number of Overlapping EBRs \t Ratio_ONE \t Ratio_OTHER\n";
foreach my $str (sort keys %count) {
	my ($vv1, $vv2);
	if (%countONE) {$vv1=$countONE{$str}} else {$vv1='NA';}
	if (%countOTHER) {$vv2=$countOTHER{$str}} else {$vv2='NA';}
	print $outfh4 "$str\t$count{$str}\t$vv1\t$vv2\n";
	createResTable ($str,$count{$str},$vv1,$vv2, \%classHash, $dataFinal,$refName);
}


if (@missinBrk){
print $outfh4 "\n\nMISSED (Reference used EBRs from $file_ref ) OVERLAPS EBRs------------>>>\n";
print $outfh4 "chr\tstart\tend\tfirstBestratio\tsecondBestRatio\tPercentageUsed\tBrk\tclass\tsameordiffClass\tsameordiffBrk\tborderUsed\tChrREF\tStREF\tChrREF\tStREF\tfinalREFchrTAR\tfinalCor1TAR\tfinalREFchrTAR\tfinalCor2TAR\tdeciChrTAR1\tdeciCordiTAR1\tdeciBrkTAR1\tdeciChrTAR2\tdeciCordiTAR2\tdeciBrkTAR2\tbothBrkGap\tbrkDeci1\tpresence1\tbrkDeci2\tpresence2\tfinalDecision\n";
foreach my $v (@missinBrk) {
	my @vals=split('\t', $v);
	my ($refChr1, $corSt) = extractRefCors($vals[0],$vals[1],$allHSB, $chekerExtend, $refName);
	my ($refChr2, $corEd) = extractRefCors($vals[0],$vals[2],$allHSB, $chekerExtend, $refName);

	if (!$corSt) {$refChr1='0';$corSt='0';} if (!$corEd) {$refChr2='0'; $corEd='0';}

	my ($refChrVal1, $corVal1, $DeciVal1) = assignRefCorsVal($refChr1,$corSt, $chekerExtend, $refName, $allall);
        my ($refChrVal2, $corVal2, $DeciVal2) = assignRefCorsVal($refChr2,$corEd, $chekerExtend, $refName, $allall);

	my ($fclass1,$fc1) = assignClass($refChrVal1,$corVal1, $chekerExtend, $finalEBA, \%classHash, $refName);
        my ($fclass2,$fc2) = assignClass($refChrVal2,$corVal2, $chekerExtend, $finalEBA, \%classHash, $refName);

	if (!$corVal1) {$refChrVal1='NA';$corVal1='NA'; $DeciVal1='NA';} if (!$corVal2) {$refChrVal2='NA';$corVal2='NA'; $DeciVal2='NA';}
	my $finClass="$DeciVal1,$DeciVal2";
	#if (index($finClass, 'NA') != -1) { $finClass =~ s/NA//g;}
	my $finalD='NFD';
	if ((index($fc1, $refName) != -1) or (index($fc2, $refName) != -1)) { $finalD = "$refName"; }
	my $finalLine="$v\t$refChr1\t$corSt\t$refChr2\t$corEd\t$refChrVal1\t$corVal1\t$DeciVal1\t$refChrVal2\t$corVal2\t$DeciVal2\t$finClass\t$fclass1\t$fc1\t$fclass2\t$fc2\t$finalD\n";
	print $outfh4 "$finalLine";
	createMissedTable ($finalLine, $missedFinal, $refName);
	}
}

=pod
#Print the same classification 'missed' in TARGET genome
print $outfh4 $message;
print $outfh4 "MISSED/NO OVERLAP Breakpoint with same classification in Target ( for $file_tar : Used as target) \n";
print $outfh4 "chr\tstart\tend\tfirstBestratio\tsecondBestRatio\tPercentageUsed\tBrk\tclass\tsameordiffClass\tsameordiffBrk\tborderUsed\tChrREF\tStREF\tChrREF\tStREF\n";
@remain = do { my %seen; grep { !$seen{$_}++ } @remain };
@remain = grep { not exists $done{$_} } @remain;
foreach my $remain (@remain) { print $outfh4 "$remain\n";}
=cut


} #checkOverklaps sub end here


#sum all hits values
#use List::Util 'sum';
#my $value_count = sum values %count;

sub uniq {
my %seen;
return grep { !$seen{$_}++ } @_;
}

# Checks if a provided two coordinates overlaps or not
#my $OverRes = EBALib::CommonSubs::checkCorOverlaps ($val1[0],$val1[1],$val2[0],$val2[1]);
#if ($OverRes) { $j1++;} else {$j2++; }
sub checkOverlaps {
my ($x1, $x2, $y1, $y2)=@_;
return $x1 <= $y2 && $y1 <= $x2;
}


sub createMissedTable {
my ($line, $missedFinal, $substr) = @_;
open(my $fhOut, '>>',$missedFinal) or die "Could not open file '$missedFinal' $!";
print $fhOut "$line";
close $fhOut;
}


sub createResTable {
my ($str, $count, $value1, $value2, $classHash_ref, $dataFinal, $substr) = @_;
open(my $fhOut, '>>',$dataFinal) or die "Could not open file '$dataFinal' $!";
my $sps='NA';
my %classHash=%$classHash_ref;

if ($classHash{$str}) {
	if (index($classHash{$str}, $substr) != -1) { # -1 if found
		$sps=$substr;
		print $fhOut "$substr\t$str\t$count\t$value1\t$value2\t$sps\n";
		#print $fhOut "$substr\t$str\t$count\t$value1\t$value2\tcontains:$sps\n";
	}
	else {
		print $fhOut "$substr\t$str\t$count\t$value1\t$value2\t$sps\n";
	}
}
else {
	if ($substr eq $str) {
		$sps=$substr;
		print $fhOut "$substr\t$str\t$count\t$value1\t$value2\t$sps\n";
	}
	else {
		print $fhOut "$substr\t$str\t$count\t$value1\t$value2\t$sps\n";
	}
}
close $fhOut;
}


sub extractRefCors {
my ($chr, $cor, $file, $chekerExtend, $refName) = @_;
my $finalRefCor='NA';
	open(my $fh, '<:encoding(UTF-8)', $file) or die "Could not open file '$file' $!";
	while (my $row = <$fh>) {
		chomp $row;
		my @vals = split('\t', $row);
		next if $refName ne $vals[8]; # Interested in a species only
		my $st=$vals[5]-$chekerExtend; # If extension needed
		my $ed=$vals[6]+$chekerExtend;
		if (($cor >= $st && $cor <= $ed ) and ($vals[4] eq $chr)) {
			my $diff=$cor-$vals[5];
			my $finalRefCor=$vals[2]+$diff;
			return ($vals[1], $finalRefCor);
   			# $x is in the range
		}
	}
	close $fh;
}


sub assignRefCorsVal {
my ($chr, $cor, $chekerExtend, $refName, $file) = @_;
#my $file = "$allALLeba00"; # Provide target all_all.eba00 file
my $finalRefCorVal='NA';
	open(my $fh, '<:encoding(UTF-8)', $file) or die "Could not open file '$file' $!";
	while (my $row = <$fh>) {
		chomp $row;
		my @vals = split('\t', $row);
		next if $refName ne $vals[0]; # Interested in a species only
		my @cordi=split('\--', $vals[2]);
		my $st=$cordi[0]-$chekerExtend; # If extension needed
		my $ed=$cordi[1]+$chekerExtend;
		if (($cor >= $st && $cor <= $ed ) and ($vals[1] eq $chr)) {
			return ($vals[1], $vals[2], $vals[3]);
   			# $x is in the range
		}
	}
	close $fh;
}


sub assignClass {
my ($chr, $cor, $chekerExtend, $finalEBA, $classHash_ref, $refName) = @_;
my  @allClass;
my $classRes;
my %classHash=%$classHash_ref;
open(my $fh, '<:encoding(UTF-8)', $finalEBA) or die "Could not open file '$finalEBA' $!";
	while (my $row = <$fh>) {
		chomp $row;
		my @vals = split('\t', $row);
		next if !$cor;
		next if $vals[13] eq "Brk_Point";
		my @cordi=split('\--', $vals[13]);
		my @qcor= split('\--', $cor);
		my $st=$cordi[0] - $chekerExtend; # If extension needed
		my $ed=$cordi[1] + $chekerExtend;
		next if $vals[12] ne $chr;
		my $OverRes = checkOverlaps ($qcor[0],$qcor[1],$st,$ed);
		if ($OverRes) { push @allClass, $vals[18]; }
	}
	my $aClass= join ',', @allClass;
	if ($aClass) { $classRes = checkPresence (\@allClass, \%classHash, $refName); } else {$aClass='NC';}
	if (!$classRes) {$classRes= "NCA"}
	return ($aClass , $classRes);
	close $fh;
}

sub checkPresence {
my ($allClass_ref, $classHash_ref, $refName) = @_;
my %classHash=%$classHash_ref;
my @allClass=@$allClass_ref;
foreach my $val (@allClass) {
	my @cLine=split ('\:', $val);
	if ($cLine[0] eq $refName) {return $refName;}
	elsif ($classHash{$cLine[0]}) {
		if (index($classHash{$cLine[0]}, $refName) != -1) { return $refName; }
		}
	}
}


#Store the classification file
sub trim($)
{
my $string = shift;
$string =~ s/^[\t\s]+//;
$string =~ s/[\t\s]+$//;
$string =~ s/[\r\n]+$//; ## remove odd or bad newline ...
return $string;
}

#Create table for the TAR2TAR data
sub createResultTable {
my ($stat, $infile1, $infile2, $outfile) = @_;

open(my $rout, ">$outfile") or die "Could not open file '$outfile' $!"; #countSTAT_<SPS> file

open(my $infh, $stat) or die "Could not open file '$stat' $!"; #countSTAT_<SPS> file
while (<$infh>) {
	chomp $_; my $line=$_;
	$line=trim($line);
	next if $line =~ /^\s*$/;
	my @tmpLine = split /\t/, $line;
	my $spsName = $tmpLine[0];
	my $missedNum=$tmpLine[1]-$tmpLine[2];
	my ($realCnt, $etc) = countOccurance( $spsName, $infile1);
	my ($spsCnt, $others, $gap) = findInMissed($spsName, $infile2);
	my $totalDetected=$tmpLine[2]+$spsCnt+$gap;
	my $totalReal=$realCnt+$spsCnt+$gap;
	print $rout "$line($realCnt:$etc),$missedNum($spsCnt:$others),$others($gap)\t$totalDetected\t$totalReal\n";
	undef $realCnt; undef $etc; undef $spsCnt; undef $others; undef $gap;
    }
close $infh;

} #create table sub ends here

sub countOccurance {
my ($name, $infile)= @_;
my $cnt=0; my $others=0;
open(my $fh1, '<:encoding(UTF-8)', $infile) or die "Could not open file '$infile' $!"; #finalOut_<sps>
while (<$fh1>) {
	chomp $_; my $line=$_;
	$line=trim($line);
	next if $line =~ /^\s*$/;
	my @tmpLine = split /\t/, $line;
	my $spsName = $tmpLine[0];
	if (index($tmpLine[0], $name) != -1) {
		if (index($tmpLine[5], $name) != -1) {
			$cnt=$cnt+$tmpLine[4];
		}
		else { $others=$others+$tmpLine[4]; }
	}
    }
close $fh1;
return ($cnt, $others);
}



sub findInMissed {
my ($spsName, $infile)= @_;
my $cnt=0; my $others=0; my $gap=0;
open(my $fh2, '<:encoding(UTF-8)', $infile) or die "Could not open file '$infile' $!"; #missedOut_<sps>
while (<$fh2>) {
	chomp $_; my $line= lc ($_);
	$line=trim($line);
	next if $line =~ /^\s*$/;
	my @tmpLine = split /\t/, $line;
	if (index($tmpLine[7], $spsName) != -1) { #print "$tmpLine[7], $spsName --\n";
		if (index($tmpLine[30], $spsName) != -1) {
			$cnt++;
		}
		else {
			if (index($tmpLine[25], 'gap') != -1) { $gap++;}
			$others++;
		}
	}
    }
close $fh2;
return ($cnt, $others, $gap);
}

sub help {

print "OH";

}


sub reformatTable {
my ( $infile1, $finalOut) = @_;

open(my $frout, ">$finalOut") or die "Could not open file '$finalOut' $!"; #countSTAT_<SPS> file

open(my $infh, '<:encoding(UTF-8)', $infile1) or die "Could not open file '$infile1' $!"; #countSTAT_<SPS> file
while (<$infh>) {
	chomp $_; my $line=$_;
	$line=trim($line);
	next if $line =~ /^\s*$/;
	my @tmpLine = split /\t/, $line;
	for (my $val1=0; $val1<=scalar(@tmpLine); $val1++) {
		if ($val1 == 1) {
		print $frout "$tmpLine[0]\tBreak\t$tmpLine[1]\n";
		}
		elsif ($val1 == 2) {
		print $frout "$tmpLine[0]\tDetected\t$tmpLine[3]\n";
		}
		elsif ($val1 == 3) {
		print $frout "$tmpLine[0]\tReal\t$tmpLine[4]\n";
		}
		else {

		}
	}
    }
close $infh;

}


sub tar2ref {
my ($filename, $filename2, $filename3, $spsName, $spsNum, $statOut)= @_;
#User need to provide the TAR reconstructed.e and final_classify.eba7 file ... see the commandline usage below for more detail
#USAGE perl checkOverlaps.pl <filenameReconstructed.e> <final_classify.eba7> <OutFileName> <NameORgroupLookinFor> <spsNum>

#my $filename = "$ARGV[0]";
#my $filename3 = "$ARGV[2]";
open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";
open(my $fh3, '>:encoding(UTF-8)', $filename3) or die "Could not open file '$filename3' $!";

my %count; my %countONE; my %countOTHER;my %countTarClass; my $allCount=0;my %allClass; my $flag=0; my @allCor;
my $sameClass_sameBrk=0; my $diffClass_sameBrk=0; my $singleClass_singleBrk=0; my $NA=0;
my @missinBrk;

while (my $row = <$fh>) {
chomp $row;
next if $. == 1;
my @tmp1 = split('\t', $row);
no warnings; # To avoid warnings for blank spaces
my @cla =split('\:', $tmp1[7]);
next if $cla[0] ne $spsName; # group or species name provided by user
$countTarClass{$cla[0]}++;
push @allCor, "$tmp1[0]:$tmp1[1]:$tmp1[2]";

if ($tmp1[8] eq "sameClass" and $tmp1[9] eq "sameBrk") { $sameClass_sameBrk++; }
elsif ($tmp1[8] eq "diffClass" and $tmp1[9] eq "sameBrk") { $diffClass_sameBrk++; }
elsif ($tmp1[8] eq "singleClass" and $tmp1[9] eq "singleBrk") { $singleClass_singleBrk++; }
else { $NA++; }

#my $filename2 = "$ARGV[1]";
my $done=0;
open(my $fh2, '<:encoding(UTF-8)', $filename2) or die "Could not open file '$filename2' $!";

 	while (my $row2 = <$fh2>) {
	chomp $row2;
	next if $.==1;
	my @values2 = split('\t', $row2);
	#my $spsNum=$ARGV[4];
	my @cor =split('\--', $values2[$spsNum+2]);
	my @class =split('\:', $values2[$spsNum+7]);
	if ($flag==0) {$allClass{$class[0]}++;}


	if ($values2[$spsNum+1] eq $tmp1[0]) {
		my $OverRes = checkOverlaps($cor[0],$cor[1],$tmp1[1],$tmp1[2]);
		if ($OverRes) {
  		 	$count{$class[0]}++;
			$allCount++;
			if ($values2[$spsNum+9] == 1) {$countONE{$class[0]}++;} else {$countOTHER{$class[0]}++;}
			print $fh3 "$row\t\t$row2\n";
			$done=1;
     			}
		}
	}
$flag=1; # To check the final_classify.eba7 once

close $fh2;
if ($done != 1) { push @missinBrk, $row;}
}
close $fh;
close $fh3;

my $multiHits="$countTarClass{$spsName}:$allCount";
my $uniqueRecoEBRs= uniq (@allCor);


#Create a file for detail stats
open(my $outSTAT, '>:encoding(UTF-8)', $statOut) or die "Could not open file '$statOut' $!";

my $message = <<"END_MSG";

..........
Hello USER, ...........$spsName analysis ..................... SEE NEXT TABLE BELOW?
..........

END_MSG


## Printing the TAR reconstructed Info
print $outSTAT "Overview of reconstructed table detail\n";
print $outSTAT "Species/Group Name \t Total Number of reconstructed EBRs \t Multiple Hits in reference table \t Unique reconstructed EBRs \t Duplicates reconstructed EBRs\t Tar\t sameClass_sameBrk \t diffClass_sameBrk \t singleClass_singleBrk \t Others\n";

foreach my $class (sort keys %countTarClass) {  # This will be only one as user check for indivisual
my $difference=($countTarClass{$class}-$uniqueRecoEBRs);
print $outSTAT "$class\t$countTarClass{$class}\t$multiHits\t$uniqueRecoEBRs\t$difference\tTarRecon\t$sameClass_sameBrk\t$diffClass_sameBrk\t $singleClass_singleBrk \t $NA\n"; }

print $outSTAT $message;

print $outSTAT "Reconstructed EBRs overlapping details\n";
print $outSTAT "Species/Group Name \t Total Number of Overlapping EBRs \t Ratio_ONE \t Ratio_OTHER\n";
foreach my $str (sort keys %count) { print $outSTAT "$str\t$count{$str}\t$countONE{$str}\t$countOTHER{$str}\n"; }


if (@missinBrk){
print $outSTAT "\n\nMISSED OVERLAPS EBRs------------>>>\n";
print $outSTAT "chr\tstart\tend\tfirstBestratio\tsecondBestRatio\tPercentageUsed\tBrk\tclass\tsameordiffClass\tsameordiffBrk\tborderUsed\n";
foreach my $v (@missinBrk) { print $outSTAT "$v\n";}
}


#sum all hits values
use List::Util 'sum';
my $value_count = sum values %count;

}

#Flip the HSB coordinates
sub flipHSB {
my ($inHSB, $outHSB)=@_;
open (OUTHSB, ">$outHSB") or die "$0: open $outHSB: $!";
open INHSB,  $inHSB or die "$0: open $inHSB: $!";
while (<INHSB>) {
	chomp;
	my $line=trim($_);
	my @tmp = split /\t/, lc ($line);
	print OUTHSB "$tmp[0]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[7]\t$tmp[8]\t$tmp[9]\n";
	}
close OUTHSB or die "could not close OUTHSB file: $!\n";
close INHSB or die "could not close INHSB file: $!\n";
}

##Recontruct the breakpoints
sub reconstructAncestral {
my ($finalEBA, $allHSB, $threshold, $length, $AncestralNames_ref, $SpsNumber, $refName, $param, $refSciName2)=@_;

#if (-d "$param->{out_dir}/output_$refName") {
#deldir("$param->{out_dir}/output_$refName"); # or deldir($ARGV[0]) to make it commandline
#} else { mkdir "$param->{out_dir}/output_$refName"; }
#Flip the HSB
#flipHSB ($allHSB, "$param->{out_dir}/intermediate_files/all_flipped_$refName.hsb");
#my $newHSB="$param->{out_dir}/intermediate_files/all_flipped_$refName.hsb";

if (-f "ReconstructionAncestral_$refName.stats") { unlink "ReconstructionAncestral_$refName.stats";}
#my @SpsArray = @$AncestralNames_ref;
foreach my $spsName(@$AncestralNames_ref) {
  $spsName =lc $spsName;
  print "Checking $spsName --- $refName -- $refSciName2\n";
  my $outAncestralfile1="$param->{out_dir}/output_$refName/$spsName"."_brk_$refName.anc1";
  my $outAncestralfile2="$param->{out_dir}/output_$refName/$spsName"."_brk_$refName.anc2";
  my $outAncestralfile3="$param->{out_dir}/output_$refName/$spsName"."_brk_$refName.anc3";

  open (OUTFILE1, ">$outAncestralfile1") or die "$0: open $outAncestralfile1: $!";
  open (OUTSTAT, ">>$param->{out_dir}/output_$refName/ReconstructionAncestral_$refName.stats") or die "$0: open $param->{out_dir}/output_$refName/ReconstructionAncestral_$refName.stats: $!";

  open INFILE,  $finalEBA or die "$0: open $finalEBA: $!";
	my @array; my @index; my @nameArray; my $in; my $countReal=0; my @done; my $countBrk=0; my $countGap=0; my $total=0;
	while (<INFILE>) {
		chomp;
		my $line=trim($_);
		my @tmp = split /\t/, lc ($line);
    my @classification = split /\:/, $tmp[$SpsNumber+7];
    my @coordinate = split /\<-->/, $tmp[$SpsNumber+4]; #broad is +2
    next if $classification[0] ne $spsName;

    #ignore if below threashold
    next if $tmp[$SpsNumber+9] <= $threshold;
    #next if $tmp[$SpsNumber+12] >= ($SpsNumber/3);
    next if $tmp[$SpsNumber+12] >= $SpsNumber/5;
		# In case more than one breakpoints ( separated with comma)
    #1	8577687	8577826	12046.0447048244	947.055705950168	10	breakpoints	gallus_gallus:0.793690571736545	sameClass	sameBrk	2	1	8882608	1	9223943

		next if $_ =~ /^\s*#/;

    $countReal++; $countBrk++; $countGap=0;
		$total++;

    #Create a similar file format to reconstructTar creator
    #45377370--45382081=breakpoints+0.954545454545455	1	taeniopygia_guttata:0.884990322164547	2634.19668710963	0.0909090909090909	11	1	10
    print OUTFILE1 "$tmp[$SpsNumber+2]=$tmp[$SpsNumber]+$classification[1]\t$tmp[$SpsNumber+1]\t$tmp[$SpsNumber+7]\t$tmp[$SpsNumber+9]\t$tmp[$SpsNumber+10]\t$tmp[$SpsNumber+11]\t$tmp[$SpsNumber+12]\t$tmp[$SpsNumber+13]\n";
    #print OUTFILE1 "$line\n";

		#72412203--72417432=Breakpoints+0.925	Breakpoints	chicken:6.22711733452034e-07	237.8835341	0.05	20	1	18
	}
	if (-z "ReconstructionAncestral_$refName.stats") { print OUTSTAT "spsName\tcountReal\tcountBrk\tcountGap\ttotal\n";}
	print OUTSTAT "$spsName\t$countReal\t$countBrk\t$countGap\t$total\n";
	close INFILE or die "could not close $finalEBA file: $!\n";
close OUTSTAT or die "could not close OUTSTAT file: $!\n";
close OUTFILE1 or die "could not close $outAncestralfile1 file: $!\n";

reconstructBrk2($outAncestralfile1, $outAncestralfile2, $refSciName2, $outAncestralfile3, $threshold, $allHSB, $param, $SpsNumber);
}}



#store ancestral breakpoints detail
sub storeClasses {
my ($inClass, $inSps, $param, $name, $location) = @_;

open(my $infh, '<:encoding(UTF-8)', $inClass) or die "Could not open file '$inClass' $!";
my @allClass;
while (<$infh>) {
	chomp $_; my $line=$_;
	$line=trim($line);
	next if $line =~ /^\s*$/;
	my @tmpLine = split /\=/, $line;
	if (($tmpLine[0] ne "lineage") and ( $tmpLine[0] ne $name)) { push @allClass, $tmpLine[0]; }
}
my $classString=join ',', @allClass;
copy($inSps,"$param->{out_dir}/sps_$name.txt") or die "Copy failed: $!";

open(my $fh, '>>', "$param->{out_dir}/sps_$name.txt") or die "Could not open file '$param->{out_dir}/sps_$name.txt' $!";
say $fh "$classString";
close $fh;

return "$classString";

}

1;

__END__
