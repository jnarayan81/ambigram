use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

# Version v1.0
# Contact jnarayan81@gmail.com


my (
	$infile1,
	$help
);
my $version = 0.1;

GetOptions(
    'infile|i=s' => \$infile1,
    'help|h' => \$help
) or  die("Error in command line arguments\n");


#&help($version) if $help;
#pod2usage("$0: \nI am afraid, no files given.")  if ((@ARGV == 0) && (-t STDIN));
if (!$infile1 ) { help($version) }

#Store the classification file
sub trim($)
{
my $string = shift;
$string =~ s/^[\t\s]+//;
$string =~ s/[\t\s]+$//;
$string =~ s/[\r\n]+$//; ## remove odd or bad newline ...
return $string;
}

open(my $infh, '<:encoding(UTF-8)', $infile1) or die "Could not open file '$infile1' $!"; #countSTAT_<SPS> file
while (<$infh>) {
	chomp $_; my $line=$_;
	$line=trim($line);
	next if $line =~ /^\s*$/;
	my @tmpLine = split /\t/, $line;
	for (my $val1=0; $val1<=scalar(@tmpLine); $val1++) {
		if ($val1 == 1) {
		print "$tmpLine[0]\tBreak\t$tmpLine[1]\n";
		}
		elsif ($val1 == 2) {
		print "$tmpLine[0]\tDetected\t$tmpLine[3]\n";
		}
		elsif ($val1 == 3) {
		print "$tmpLine[0]\tReal\t$tmpLine[4]\n";
		}
		else {

		}
	}
    }
close $infh;

sub help {

print "OH";

}
