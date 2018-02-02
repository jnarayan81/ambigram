
---. .-. .---
  --\'A'/--
     \ /
     " "
---ambigram v0.1---

# ambigram

ambigram is a bioinformatics tool to check if chromosomal breakpoints retain meaning when viewed or interpreted from a different direction, perspective, or orientation.

# ambigram
ambigram: evolutionary chromosomes breakpoints validator

Copyright 2018 Jitendra Narayan <jitendra.narayan@unamur.be> or jnarayan81@gmail.com;

## LICENSE

ambigram is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

ambigram is distributed with the hope that it will be useful for reseachers worldwide, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with ambigram, in a file called COPYING. If not, see <http://www.gnu.org/licenses/>.

## INSTALLATION

See INSTALLATION steps below:

1) Place the package into your environment, add
call perl ambigram.pl

USING PERL 5.x

ambigram does support Perl 5.x, and no plans exist to provide any support to Perl6. For a strong biological analysis package for perl 5, see perl https://www.perl.org/ and bioperl: http://bioperl.org/

Perl modules

You also need to have Perl, BioPerl, and some other modules installed in your
machine:
```
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
```
You can check if above mentioned modules installed in your machine with
'perl -M<module> -e 1'. It will return an error message if it isn't installed.

E.g. "perl -MBio::SeqIO -e 1"

To install these modules, you can do either through the CPAN or manually downloading
(http://search.cpan.org/) and compiling them. To use CPAN, you can do by
writing:

> perl -MCPAN -e 'install <module>'

To install manually, search for the most recent version of these modules and,
for each, download and type the following (should work most of the time, except
for modules like BioPerl, of course):

> tar -zxvf <module.tar.gz>
> perl Makefile.PL
> make
> make test
> make install

For more detail, you could visit: http://bioinformaticsonline.com/blog/view/710/how-to-install-perl-modules-manually-using-cpan-command-and-other-quick-ways

## DOCUMENTATION

For documentation, see ambigram_manual.pdf in the base project folder.

This documentation may be out of date depending on whether or not the developers did their job and re-generated the documentation before the release. If you suspect that the documentation is out of date, or if you are using code from the repository (and not from a release), you can re-generate the documentation or contact the authors jnarayan81ATgmail.com.

## RELEASE HISTORY

0.1.O - 30 Jan 2018

## OUTPUT FORMAT

ambigram outfile columns: [ NEED UPDATE ]

* CHROM       Reference entry where ambigram occurs
* OUTERSTART  The 5' most boundary estimate of where ambigram begins
* START       Best guess as the the exact start of the ambigram
* INNERSTART  The 3' most boundary estimate of where ambigram begins
* INNEREND    The 5' most boundary estimate of where ambigram ends
* End         Best guess as the the exact end of the ambigram
* OUTEREND    The 3' most boundary estimate of where the ambigram ends
* TYPE        Variant type. One of MIS, INS, INSZ, or DEL
* SIZE        Estimate size of the ambigram
* INFO        More information associated with the calls

ambigram result columns:

* ID        Unique identifier of the call
* CHROM     Reference entry where q4 occurs
* START     Start point
* END       End point
* EXSEQ     Average amount of sequence left between
* GFF       GFF entry

ambigram tabfile column

* SEQ       ambigram sequence
* NAME      Name of the ambigram
* REFNAME   Reference name
* START     Start coordinate of ambigram
* END       End coordinate of ambigram
* MUTATION  Number of random mutation
* SIZE      Size of the ambigram
* STRAND    Strand of of the ambigram
* ISLAND    Number of G island
* SCORE     Score of the ambigram


## ANNOTATION DESCRIPTIONS
Coming Soon

## EXTRA
Coming Soon

## FAQ

Can I report bugs  or ask questions?
Please report your issues to ticketing system.

## CONTRIBUTION

Feel free to clone this repository and use it under the licensing terms.

Additionally, as the project is on github, you may submit patches, ticket requests, edit the wiki, send pull requests - anything you like and have the permissions to do. I will encourage any suggestions from followers :)

As always, you can contact the authors at <jitendra.narayan@unamur.be> or jnarayan81@gmail.com ;
