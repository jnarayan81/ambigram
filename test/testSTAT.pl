
analyticStat("/home/jitendra/ETC/quadraG/out/intermediate_files/quadra/final_annotated.q4", "final_annotated.q4.stat", 0);


########################################################################################"

sub analyticStat {
my ($file, $outfile, $score)=@_;
my @terms;

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
#CDS\tgene\tmRNA\tUTR
foreach my $f(@$feature_ref) { if ($feature{$f}) { print $out "$feature{$f}\t"; } else { print $out "0\t";} }
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
#Open and write a file
sub write_fh {
    my $filename = shift @_;
    my $filehandle;
    open $filehandle, ">$filename" or die $!;
    return $filehandle;
}

