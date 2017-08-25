use strict;
use warnings;
use Term::ANSIColor;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

# Version v1.0
# Contact jnarayan81@gmail.com
# User need to provide the REF reconstructed.tar3 and TAR secondReconstructed.tar3 file ... see the commandline usage below for more detail
# USAGE perl checkOverlaps.pl <firstReconstructed.tar3> <secondReconstructed.tar3> <OutFileName> <NameORgroupLookinFor> <spsNum> <numberToincreaseBreksSize>

#perl checkOverlaps_TAR2TAR.Pl anas_platyrhynchos_brk_chicken_ref.tar3 anas_platyrhynchos_brk_finch_ref.tar3 OUT anas_platyrhynchos 11 0 > STAT

my (
	$file_ref, 
	$file_tar, 
	$spsFile, 
	$allall, 
	$allHSB, 
	$help, 
	$refName, 
	$length, 
	$outfile, 
	$spsNum, 
	$statfile, 
	$outfile2, 
	$identify, 
	$dataFinal,
	$missedFinal, 
	$groupClass,
	$finalEBA, 
	$chekerExtend
);
my $version = 0.1;

GetOptions(
    'file_ref|r=s' => \$file_ref,
    'file_tar|t=s' => \$file_tar,
    'allall|a=s' => \$allall,
    'spsFile|s=s' => \$spsFile,
    'groupClass|g=s' => \$groupClass,
    'spsNum|c=i' => \$spsNum,
    'allHSB|b=s' => \$allHSB, # all HSB from same file
    'refName|n=s' => \$refName,
    'outfile|o=s' => \$outfile,
    'dataFinal|d=s' => \$dataFinal,
    'missedFinal|m=s' => \$missedFinal,
    'outfile2|e=s' => \$outfile2,
    'identify|i=s' => \$identify,
    'finalEBA|f=s' => \$finalEBA, #The final_out.eba7 file for the refernce
    'chekerExtend|x=i' => \$chekerExtend, #To extend the size (in bp) of missed breaks before check in final.eba file
    'statfile|z=s' => \$statfile,
    'length|l=i' => \$length, #Length to extend the size of the each breakpoints before check
    'help|h' => \$help
) or  die("Error in command line arguments\n");


#&help($version) if $help;
#pod2usage("$0: \nI am afraid, no files given.")  if ((@ARGV == 0) && (-t STDIN));
if (!$file_ref or !$file_tar or !$spsFile or !$allHSB or !$allall or !$spsNum or !$refName or !$outfile) { help($version) }

#if (-d "output_$refName") {
#deldir("test"); # or deldir($file_ref) to make it commandline
#} else { mkdir "output_$refName"; }

#Store the classification file
sub trim($)
{
my $string = shift;
$string =~ s/^[\t\s]+//;
$string =~ s/[\t\s]+$//;
$string =~ s/[\r\n]+$//; ## remove odd or bad newline ...
return $string;
}

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
open SPSFILE, "$spsFile" or die colored ["bold red on_black"], $!;
my @SpsArray; my $SpsNumber;
while (<SPSFILE>) { 
	my $SpsLine=$_; chomp $SpsLine; @SpsArray=split /,/, lc($SpsLine);  $SpsNumber = scalar (@SpsArray); } ## It read the species names from sps.txt file ... need to improve !!!
	my $SpsArrayTabed=join("\t", @SpsArray); 

close SPSFILE or die colored ["bold red on_black"], "could not close file: $!\n";


# Use need to provide the location of all.hsb (of reference used in commandline).
my $allHSBFile = $allHSB;
#Use need to provide the location of all_all.eba00 (of target used in commandline).
my $allALLeba00 = $allall;

# User can find the above files in EBA tools output dir

my $filename = "$file_ref";
my $filename3 = "$outfile";
open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";
open(my $fh3, '>:encoding(UTF-8)', $filename3) or die "Could not open file '$filename3' $!";

my %count; my %countONE; my %countOTHER;my %countTarClass; my $allCount=0;my %allClass; my $flag=0; my @allCor;
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
my $done=0;
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
		if ($OverRes) {
  		 	$count{$class[0]}++;
			$allCount++;
			#If anyone of the ratio is one ... count
			no warnings;
			if ( -z "$filename3" ) { print $fh3 "$head1\t\t$head2\n"; }
			if (($values2[3] == 1) or ($values2[4] == 1)) {$countONE{$class[0]}++;} else {$countOTHER{$class[0]}++;}
			print $fh3 "$row\t\t$row2\n"; 
			$done=1;
			$done{$row2}=$lineNumR2;
     			}
		}
	if ($class[0] eq $refName) { push @remain, $row2;} 
	}
$flag=1; # To check the final_classify.eba7 once

close $fh2;

if ($done != 1) { push @missinBrk, $row;}
}

close $fh;
close $fh3;

my $multiHits="$countTarClass{$refName}:$allCount";
my $uniqueRecoEBRs= uniq (@allCor);

#Only values
print "$identify\t$countTarClass{$refName}\t$allCount\n";

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
	my ($refChr1, $corSt) = extractRefCors($vals[0],$vals[1],$allHSB, $chekerExtend); 
	my ($refChr2, $corEd) = extractRefCors($vals[0],$vals[2],$allHSB, $chekerExtend);

	if (!$corSt) {$refChr1='0';$corSt='0';} if (!$corEd) {$refChr2='0'; $corEd='0';}

	my ($refChrVal1, $corVal1, $DeciVal1) = assignRefCorsVal($refChr1,$corSt, $chekerExtend);
        my ($refChrVal2, $corVal2, $DeciVal2) = assignRefCorsVal($refChr2,$corEd, $chekerExtend);

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

#sum all hits values
use List::Util 'sum';
my $value_count = sum values %count;

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
my ($chr, $cor, $file, $chekerExtend) = @_;
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
my ($chr, $cor, $chekerExtend) = @_;
my $file = "$allALLeba00"; # Provide target all_all.eba00 file
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
