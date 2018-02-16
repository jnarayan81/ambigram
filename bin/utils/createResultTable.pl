use strict;
use warnings;
use Term::ANSIColor;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

# Version v1.0
# Contact jnarayan81@gmail.com


my (
	$stat, 
	$infile1,
	$infile2,  
	$outfile,
	$help
);
my $version = 0.1;

GetOptions(
    'stat|s=s' => \$stat,
    'infile|i=s' => \$infile1,
    'missedfile|m=s' => \$infile2,
    'out|o=s' => \$outfile,
    'help|h' => \$help
) or  die("Error in command line arguments\n");


#&help($version) if $help;
#pod2usage("$0: \nI am afraid, no files given.")  if ((@ARGV == 0) && (-t STDIN));
if (!$stat or !$infile1 or !$infile2 ) { help($version) }

#Store the classification file
sub trim($)
{
my $string = shift;
$string =~ s/^[\t\s]+//;
$string =~ s/[\t\s]+$//;
$string =~ s/[\r\n]+$//; ## remove odd or bad newline ...
return $string;
}

open(my $infh, '<:encoding(UTF-8)', $stat) or die "Could not open file '$stat' $!"; #countSTAT_<SPS> file
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
	print "$line($realCnt:$etc),$missedNum($spsCnt:$others),$others($gap)\t$totalDetected\t$totalReal\n";
	#undef $realCnt; undef $etc; undef $spsCnt; undef $others; undef $gap;
    }
close $infh;


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
