#!/usr/bin/perl
use strict;
#use warnings;
use Term::ANSIColor;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

## Contact jnarayan81@gmail.com
## perl RecontructTargetEBRs_Final.pl <final_classify.eba7> <Threshold>
## Need to files all.hsb and sps.txt
## Note: Generate the all.hsb by concatinating all HSB files. ( Should be of same resolutions)
## sps.txt from EBA current place.

# Updated
# Added commandline threshold options
# The final *.tar3 file contain the correspondig reference coordinates

my ($finalEBA, $threshold, $spsFile, $allHSB, $length, $help, $refName);
my $version=0.1;
GetOptions(
    'finalEBA|f=s' => \$finalEBA,
    'spsFile|s=s' => \$spsFile,
    'allHSB|a=s' => \$allHSB, # all HSB from same file
    'threshold|t=i' => \$threshold,
    'refName|r=s' => \$refName,
    'length|l=i' => \$length,
    'help|h' => \$help
) or  die("Error in command line arguments\n");

#&help($version) if $help;
#pod2usage("$0: \nI am afraid, no files given.")  if ((@ARGV == 0) && (-t STDIN));
if (!$finalEBA or !$spsFile or !$allHSB or !$threshold or !$refName) { help($version) }

if (-d "output_$refName") {
deldir("test"); # or deldir($ARGV[0]) to make it commandline
} else { mkdir "output_$refName"; }

#Store species
open SPSFILE, "$spsFile" or die colored ["bold red on_black"], $!;
my @SpsArray; my $SpsNumber;
if (-f "Reconstruction_$refName.stats") { unlink "Reconstruction_$refName.stats";}
while (<SPSFILE>) { 
	my $SpsLine=$_; chomp $SpsLine; @SpsArray=split /,/, lc($SpsLine);  $SpsNumber = scalar (@SpsArray); } ## It read the species names from sps.txt file ... need to improve !!!
	my $SpsArrayTabed=join("\t", @SpsArray); 

close SPSFILE or die colored ["bold red on_black"], "could not close file: $!\n";

my $InFile=$finalEBA; #final_classify.eba file
my $threshold=$threshold; # threhold value to filter
foreach my $spsName(@SpsArray) {
my $outFile1="output_$refName/$spsName"."_brk_$refName.tar1";
my $outfile2="output_$refName/$spsName"."_brk_$refName.tar2";
my $outfile3="output_$refName/$spsName"."_brk_$refName.tar3";

open (OUTFILE1, ">$outFile1") or die "$0: open $outFile1: $!";
open (OUTSTAT, ">>Reconstruction_$refName.stats") or die "$0: open Reconstruction_$refName.stats: $!";

	open INFILE,  $InFile or die "$0: open $InFile: $!";
	my @array; my @index; my @nameArray; my $in; my $countReal; my @done; my $countBrk; my $countGap; my $total;
	while (<INFILE>) { 
		chomp; 
		my $line=trim($_);
		my @tmp = split /\t/, lc ($line);
		if ($. == 1) { @nameArray = split /\t/, $line; (@index)= grep { $nameArray[$_] eq "$spsName" } 0..$#nameArray; $in=$index[0]; next;}
		if ($_ =~ /^\s*#/) { next; } 
		if (!$tmp[$in]) { next;} 

		my @val = split /\,/, $tmp[$in];
		
		# In case more than one breakpoints ( separated with comma)
		foreach my $l(@val) {
		my $line=trim($l);

		my $seenIn=isInList($l,@done);
		if (!$seenIn) {$countReal++;}
		if (index($l, "breakpoints") != -1) { $countBrk++; } elsif (index($l, "gap") != -1) {$countGap++;}
		$total++;

		print OUTFILE1 "$l\t$tmp[$SpsNumber+1]\t$tmp[$SpsNumber+7]\t$tmp[$SpsNumber+9]\t$tmp[$SpsNumber+10]\t$tmp[$SpsNumber+11]\t$tmp[$SpsNumber+12]\t$tmp[$SpsNumber+13]\n";

		push @done, $l;
		#push (@array,$_);
		#72412203--72417432=Breakpoints+0.925	Breakpoints	chicken:6.22711733452034e-07	237.8835341	0.05	20	1	18
		}
		
	}
	if (-z "Reconstruction_$refName.stats") {print OUTSTAT "spsName\tcountReal\tcountBrk\tcountGap\ttotal\n";}
	print OUTSTAT "$spsName\t$countReal\t$countBrk\t$countGap\t$total\n";
	undef @done;
	close INFILE or die "could not close $InFile file: $!\n";
close OUTSTAT or die "could not close OUTSTAT file: $!\n";
close OUTFILE1 or die "could not close $outFile1 file: $!\n";

reconstructBrk($outFile1, $outfile2, $spsName, $outfile3, $threshold, $allHSB);
	
undef @index;
}

sub reconstructBrk {
my ($InFile, $OutFile, $spsName, $out, $threshold, $allHSB) = @_;

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
		$inArray{"$spsName:$array[1]:$pahala"}=$_; 
		$inArray{"$spsName:$array[1]:$dusara"}=$_;

	}
	close INFILE or die "could not tar2 file: $!\n";

		
	my $HSBInFile=$allHSB;
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
		my $a1=$inArray{"$array[8]:$array[1]:$array[2]"};
		my $a2=$inArray{"$array[8]:$array[1]:$array[3]"};
		#print $a1;
		if(!$a1) { $a1="Tel\t\t\t\t\t\t\t";} if(!$a2) { $a2="Tel\t\t\t\t\t\t\t";}
my @aa = split /\=/, lc($a1);
my @bb = split /\=/, lc($a2);
my @aa2 = split /\+/, $aa[1];
my @bb2 = split /\+/, $bb[1];
		if($array[7] eq "-") { print OUTFILE "$_\t$a2\t$bb2[0]\t$a1\t$aa2[0]\n";} else { print OUTFILE "$_\t$a1\t$aa2[0]\t$a2\t$bb2[0]\n";}
	}
	close HSBFILE or die "could not close $HSBInFile file: $!\n";
close OUTFILE or die "could not close $OutFile file: $!\n";

reconstructTarget($OutFile, $out, $spsName, $threshold);
}



sub reconstructTarget {

use strict;
#use warnings;
use Term::ANSIColor;

my ($InFile, $OutFile, $spsName, $threshold)=@_;
open (OUTFILE, ">$OutFile") or die "$0: open $OutFile: $!";
	open INFILE,  $InFile or die "$0: open $InFile: $!";
	my @array;
	while (<INFILE>) { 
		chomp; 
		next if $. == 1; ## to nex the header
		if ($_ =~ /^\s*#/) { next; }  
		$_=trim($_);
		push (@array,$_);
		
	}
	close INFILE or die "could not close $InFile file: $!\n";
	my @sorted_array = sort { (split "\t", $a)[4] cmp (split "\t", $b)[4] || (split "\t", $a)[5] <=> (split "\t", $b)[5] && (split "\t", $a)[6] <=> (split "\t", $b)[6] } @array;

#my $threshold=45; ## Need to chage if required		
print OUTFILE "chr\tstart\tend\tfirstBestratio\tsecondBestRatio\tPercentageUsed\tBrk\tclass\tsameordiffClass\tsameordiffBrk\tborderUsed\tREFchr\tREFstart\tREFchr\tREFend\n";

for (my $val=0; $val <= $#sorted_array; $val++ ) {
	next if $val==0;
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

=pod
# Returns:
# 1 - empty
# 0 - not empty
# -1 - doesn't exist
# Definition of "empty" -- no files/folders/links except . and ..
sub isEmpty{
  my ($dir) = @_;
  my $file; my $dfh;
  if (opendir $dfh, $dir){
     while (defined($file = readdir $dfh)){
        next if $file eq '.' or $file eq '..';
        closedir $dhf;
        return 0;
     }
     closedir $dfh;
     return 1;
  }else{
     return -1;
  }
}
=cut

#Help section
sub help {
  my $ver = $_[0];
  print "\n ReconstructTargetEBRs_Final $ver\n\n";

  print "Usage: $0 -f final_classify.eba7 -a all.hsb -t 2 -l 0 -s sps.txt \n\n";
  print	"Options:\n";
  print "  'finalEBA|f=s' => finalEBA,\n";
  print "  'spsFile|s=s' => spsFile,\n";
  print "  'allHSB|a=s' => allHSB, # all HSB from same file\n";
  print "  'threshold|t=i' => threshold,\n";
  print "  'refName|r=s' => refName,\n";
  print "  'length|l=i' => length,\n";
  print "  'help|h' => help\n";

exit;
}
