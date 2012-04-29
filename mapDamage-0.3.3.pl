#!/usr/bin/perl -w

# mapDamage.pl version 0.3.3, 120119

# Add option -c to chained the three commands in a row with standard plot
# added the plot options to be passed from map ->merge -> plot
# read bam files
# solved read start and end + $around fall apart the reference extremi
# add writing tables for freq C>T and G>A reverse

# This program is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation

# BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
# FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
# OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
# PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
# OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
# TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
# PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
# REPAIR OR CORRECTION.

# Authors: Aurelien Ginolhac & Ludovic Orlando
# Center of Geogenetics
# University of Copenhagen, Denmark
# contact: aginolhac@snm.ku.dk

# If you use this script, please cite the associated publication:
# Ginolhac A, Rasmussen M, Gilbert MT, Willerslev E, Orlando L. 2011
# mapDamage: testing for damage patterns in ancient DNA sequences. Bioinformatics. 2011 27(15):2153-5

use strict;
use warnings;
use Getopt::Std;

my $version = "0.3.3";

&printUsage if (@ARGV < 1);

my $command = shift(@ARGV);
my %func = (map=>\&map, merge=>\&merge, plot=>\&plot);

die("Unknown command \"$command\", expected: map | merge | plot\n") if (!defined($func{$command}));

&{$func{$command}};
exit(0);

# Main map subroutine ###

sub map {

my %opts = (i=>undef, d=>undef, l=>70, a=>10, f=>undef, r=>undef, j=>undef, u=>undef, k=>undef, c=>undef, m=>undef, b=>undef, t=>undef,  );

getopts('i:d:r:fl:a:j:ukcm:b:y:t:',\%opts);

if(defined($opts{i}) and defined($opts{r})) { # mandatory options
    
    
my $inputAncient = $opts{i};

my $folder;
my $dbFasta = $opts{r};
if(defined($opts{a})){
	if($opts{a} !~/^[0-9]+$/ ) {die("around is expected to be an integrer\n");}
	if($opts{a} > 100) {print"Warning: the around option was set to $opts{a}, are you sure to retrieve such large regions? continuing anyway\n";}
	if($opts{a} < 1) {die("The around position must be > 1 nucleotide\n");}
}
if(defined($opts{l})){
	if($opts{l} !~/^[0-9]+$/) {die("length is expected to be an integrer\n");}
	if($opts{l} < 1) {die("The length must be > 1 nucleotide\n");}
	if($opts{l} > 99998) {die("The read length choosen, $opts{l}, is too large. Must be < 99998 nucleotides\n");}
}


my $cpus = $opts{t};
my $around = $opts{a};
my $seqLength = $opts{l};
my $read = {};
my $read2 = {};
my $lengthDistrib = {};
my $pattern = 'nuclMis';
my $chromo	= 'hitPerChrom';
my $summary = 'dnaFrag';
my $classN = 'readComp';
my $misinc  = 'nuclSumMis';
my $types	= 'nuclComp';
my $typesRev	= 'nuclCompReverse';
my $max=0;
my $cpt = 0;
my $i = 0;
my $j=0;


&checkExecutable;

# exit if files are empty
&checkFileEmpty($inputAncient);
&checkFileEmpty($dbFasta);	# inspect reference file, size
my $reference = &checkFastaRef($dbFasta); # inspect reference file, format


# First, exclude hits which matched both on ancient & second genome
# write no second hits in a new file
# first will be indexed by the read name 

my $fileOUT = $inputAncient.'.filtered';
$fileOUT =~ s/\.bam/.sam/;
my $unique = 0;
my $first; 
if( defined($opts{u})) { # if user want to remove non-unique hits ie., presence of suboptimal hits in X1 tag
	$unique = 1;
	$first = readSAMfile($inputAncient,$unique);
	&writeSamFile ($first,$fileOUT);
}
else { # if nothing to remove, just use the file
	$fileOUT = $inputAncient;
	if($fileOUT =~ /\.bam$/) {$fileOUT = "samtools view $inputAncient | "; }
}


if(defined($opts{j})) {
	if(!defined($opts{u})) {
		$first = readSAMfile($inputAncient,$unique);
		$fileOUT = $inputAncient.".filtered";
		$fileOUT =~ s/\.bam/.sam/;
	}
	
	my $inputFile  = $opts{j};
	&checkFileEmpty($inputFile);
	# second will be indexed by the read name 
	my $second  = readSAMfile($inputFile,$unique);
	# now compare 
	&compare2SamFiles ($first,$second,$fileOUT);
	undef($second); # free memory
}

# free memory
undef($first);

# filtered file written, trimmed for multiple hits or / and second file
open(F_READ,"$fileOUT") or die ("Cannot open main file: $fileOUT\n");

$cpt=0; 

# initializing lengthDistrib & opening file to write distrib
$lengthDistrib = initializeHash(0,800,'+','-',%$lengthDistrib);

while ( my $line = <F_READ>  ){ # store in RAM read information
		
	chomp($line);
	my @temp = split(/\t/,$line);
	next if($temp[1] & 0x0004) ; # if no filtering was performed, skip at least unmapped read
	### reading SAM information per read, 0-offset
	# 0 read name, 
	# 1 strand
	# 2 genomic reference
	# 3 position, strand forward 1-based leftmost
	# 4 score alignement MAQ
	# 5 CIGAR, alignment code
	# 6 mate reference, forget about it 
	# 7 1-based position of mate, forget about it
	# 8 inferred insert size, forget about it
	# 9 read sequence
	# 10 read quality phred score
	$cpt++;	# increment counter, read will be indexed by 1-offset
	
	if($temp[1] & 0x0010) {	$read->{$cpt}{strand} = 16; } # at this step, unmapped read were removed
	else { $read->{$cpt}{strand} = 0; }
	$temp[2] =~ s/[:\|\s]//g;  # character '|' or others must be removed because of their interpretation by the shell 
	# the reference name must be same
	die("The reference $temp[2] doesn't correspond to any chromosomes in the reference, be sure to map and run mapDamage with the same reference file\n") if (!defined( $reference->{$temp[2]}));
	$read->{$cpt}{ref} = $temp[2];	
	$read->{$cpt}{name} = $temp[0]; # keep track of read name

	if ( $temp[3] =~ m/([0-9]+)/ ) {
		$read->{$cpt}{startPos} = $1; 
	}
	else { print "Problem with SAM format(position field)! $temp[3]\n";exit 1 }
	if(defined($temp[5] )){
		if ( $temp[5] =~ m/([0-9]+[MIDNSHP]+)/ ) { $read->{$cpt}{cigar} = $temp[5]; }
	}
	else { print "Problem with SAM format(CIGAR field)! $line\n";exit 1 }
	if ( $temp[9] =~ m/([acgtnACGTN]+)/ ) { $read->{$cpt}{seq} = uc($1); } # put sequence in upper case
	else { print "Problem with SAM format(sequence field)! $temp[9]\n";exit 1 }
	my $posEnd = $temp[3] + length($temp[9]); 
	$read->{$cpt}{endPos} = $posEnd;

	# compute size length 
	my $lg = length($read->{$cpt}{seq});
	if($max < $lg) {
		$max = $lg;
	}
	if ( $read->{$cpt}{strand} & 0x0010) {
		
		$lengthDistrib->{'-'}{$lg}++;	
	}	else {
		$lengthDistrib->{'+'}{$lg}++;
	}	
}
close F_READ;
print "$cpt reads put into RAM\n";

# before writing results, create folder
# with user name if provided via -d option
if (defined($opts{d})) {$folder = $opts{d};}
elsif(defined($opts{u})){ $folder = "results_".$opts{i}.".filtered"; $folder =~ s/\.bam/.sam/;}
else {$folder  = "results_".$opts{i}; $folder =~ s/\.bam/.sam/;}
mkdir($folder);

# write Distriblength one file per strand
open(PLUS,">$folder/lengthDistribStrd+.txt") or die "Cannot open file: $folder/lengthDistribStrd+.txt\n";
print PLUS "Length\tOccurences\n";
foreach my $lg (0..$max) {
	print PLUS " $lg\t".$lengthDistrib->{'+'}{$lg}."\n";
}
close PLUS;
open(MINUS,">$folder/lengthDistribStrd-.txt") or die "Cannot open file: $folder/lengthDistribStrd-.txt\n";
print MINUS "Length\tOccurences\n";
foreach my $lg (0..$max) {
	print MINUS " $lg\t".$lengthDistrib->{'-'}{$lg}."\n";
}
close MINUS;

# aware user if computation doesn't encompass the whole data, ie reads length > 70 bp
if ( $max > $seqLength ) {	
	print "#######\nWARNING! some reads are larger than the choosen length for computation ($max > $seqLength)\n#######\n";
}
print "recovery of genomic regions (-$around to +$around nt around match alignment)\n";
open(F_FASTA,"$dbFasta") or die "Cannot open file: $dbFasta\n";

my $bedFile = "$inputAncient.bed";
open(BED, ">$bedFile") or die ("Can't write the bed file $bedFile\n");

foreach my $key (sort {$a <=> $b} keys(%$read))	{
	
	# first round of looking at the CIGAR field to get the exact genomic coordinates to be retrieved
	my (@CIG) =  $read->{$key}{cigar} =~ m/\d+\w/g;
   
	my $nbDel = 0;
	my $nbIns =0;
  
	foreach my $var (@CIG)   {  # foreach group of integrer[MDI] ie, Match, Deletion or Insertion
	
		if (my ($I) = $var =~ m/(\d+)[IS]/)  {$nbIns  += $I;  } # if an insertion or a soft clipping, count how many
   
		elsif (my ($D) = $var =~ m/(\d+)D/)  { $nbDel  += $D;  } # if an deletion, count how many
	}
	# then define les region limits to retrieve
	my $start = ($read->{$key}{startPos}-$around-1); # bedtools is 0-offset based
	if( $start < 0 ) {
		#print "Problem with starting position, set to zero!($start for $key ".$read->{$key}{startPos}.")\n";
		$start = 0;
	#	print "for $read->{$key}{name} left! $read->{$key}{startPos} flag $read->{$key}{strand}\n ";
		$read->{$key}{shortLeft} = $read->{$key}{startPos} - 1; # keep track of how many bases we actually get
	}
	# for end position, adjust for gaps, ie., minus gap number in reference, add gap number in read
	my $end = ($read->{$key}{endPos}+$around+$nbDel-$nbIns-1);
	my $ref = $read->{$key}{ref};
	# if the reach the reference end, stop when needed
	if($end > $reference->{$ref}) {
		$end = $reference->{$ref}-1 ;
	#	print "for $read->{$key}{name} right! $end and endpos $read->{$key}{endPos} pr ref $ref flag $read->{$key}{strand}\n";
		$read->{$key}{shortRight} = abs($end - $read->{$key}{endPos} + 1); # keep track of how many bases we actually get
	}
	# soft clipping at right position could be missing
	if($read->{$key}{endPos} >  $reference->{$ref} ) {
		if ( $read->{$key}{strand} & 0x0010 ) {
			$read->{$key}{softClip5} = $read->{$key}{endPos} -  $reference->{$ref};
		}
		else {
			$read->{$key}{softClip3} = $read->{$key}{endPos} -  $reference->{$ref};
		}
		$read->{$key}{shortRight} = 0 ;
	#	print "for $read->{$key}{name} $read->{$key}{endPos} sup $read->{$key}{ref},  $read->{$key}{softClip3} ###################\n"; 
	}
	print BED "$ref\t$start\t$end\n";
}
close BED;
print "Genomic regions were written in $bedFile\n";
print "Retrieving genomic regions using fastaFromBed, could take a while...\n";
&checkFileEmpty($bedFile);
open(BED, "$bedFile") or die("Cannot read file $bedFile\n");
my $bedFasta = "${inputAncient}_ref.tab";
system("fastaFromBed -fi $dbFasta -bed $bedFile -fo $bedFasta -tab");
&checkFileEmpty($bedFasta);
print "All genomic regions retrieved and written in $bedFasta, reading this file...\n";
open(BFASTA, "$bedFasta") or die ("Cannot read file $bedFasta\n");
$cpt = 0;
my $reference2 = {};
while ( my $line = <BFASTA>) {
	
	$cpt++;
	chomp($line);
	my ($ref_coord, $seq) = split(/\t/, $line);
	my @ref = split(/:/, $ref_coord);
	if ( $seq eq "") { 
		die("Problem reading genomic sequence! ($dbFasta $ref[0])\n"); # exit if sequence empty
	}	
	if ( $read->{$cpt}{ref} ne $ref[0]) { # QC is we read in a same order
		die("Problem reading genomic sequence! ($ref[0] doesn't match ".$read->{$cpt}{ref}." for ".$read->{$cpt}{name}.")\n"); 
	}
	$reference2->{$ref[0]} = $reference->{$ref[0]}; # keep only the references we need if many
	$seq = uc($seq); # MUST BE IN UPPERCASE FOR COMPARISON READ<->REFERENCE
	$seq =~ s/[\n\r\e\t]//g;	 # erase carrier return or any space characters if any
	my $seq2 = $read->{$cpt}{seq};
	my (@CIG) =  $read->{$cpt}{cigar} =~ m/\d+\w/g;
   
	my $nbDel = 0;
	my $nbIns =0;	
	my $curr_pos = 0; 
	foreach my $var (@CIG) { 
			if (my ($M) = $var =~ m/(\d+)[MX=]/) { $curr_pos += $M; } # match, only update the current position
		elsif (my ($I) = $var =~ m/(\d+)I/)  {
			substr($seq,($curr_pos+$around),0,"-" x $I);  #insert x gap in the REFERENCE since insertion. Must correct with $around
			$curr_pos  += $I;
			$nbIns += $I;
		}
		elsif (my ($D) = $var =~ m/(\d+)D/) {
			substr($seq2,$curr_pos,0,"-" x $D); #insert x gap in the READ since deletion
			$curr_pos  += $D;
			$nbDel += $D;
		}  
		elsif( my ($S) = $var =~ m/(\d+)S/) {
			
			if( $curr_pos == 0 ) {
				if($read->{$cpt}{strand} & 0x0010 ) { $read->{$cpt}{softClip3} = $S; }
				else {  $read->{$cpt}{softClip5} = $S;}
			}
			else {	if($read->{$cpt}{strand} & 0x0010 ) { $read->{$cpt}{softClip5} = $S; }
				else {  $read->{$cpt}{softClip3} = $S;}
			}
			$curr_pos  += $S;
		}
	}
	my $before = $around;
	my $after = $around;
	# if read matches complementary reverse seq
	# same for read sequence which was reverse complemented
	if ( $read->{$cpt}{strand} & 0x0010 ) {
		# reverse REFERENCE sequence
		$seq =~ tr/ACGT/TGCA/;
		$read->{$cpt}{seqRef} = reverse( $seq );
		# reverse READ sequence
		$seq2 =~ tr/ACGT/TGCA/;
		$read->{$cpt}{seq} = reverse( $seq2 );
	}
	else {
		$read->{$cpt}{seqRef} = $seq;
		$read->{$cpt}{seq} = $seq2;
	}
	
	if(defined($read->{$cpt}{shortRight})) { 
		if ( $read->{$cpt}{strand} & 0x0010 ) {
			$before = $read->{$cpt}{shortRight};
		}
		else {
			$after = $read->{$cpt}{shortRight};
		}
	#	print "for $read->{$cpt}{name} revcomp $before $after\n ";
	}
	
	if(defined($read->{$cpt}{shortLeft})) { 
		if ( $read->{$cpt}{strand} & 0x0010 ) {
			$after = $read->{$cpt}{shortLeft};
		}
		else {
			$before = $read->{$cpt}{shortLeft};
		}
	#	print "for $read->{$cpt}{name} $before $after\n ";
	}	
	
	# create after and before regions for handy computation
	# {before} will contain the $around nucleotides just before the read
	# and {after} will contain the $around nucleotides just after the read
	$read->{$cpt}{before} = substr($read->{$cpt}{seqRef},0,$before); 
	my $summ = (length($read->{$cpt}{seqRef})-$after);
	$read->{$cpt}{after} = substr($read->{$cpt}{seqRef},$summ,$after);
	if(($cpt%50000) == 0) {printf("%8d genomic regions computed\n",$cpt); }
}
close BFASTA;

if(!defined($opts{k})) { 
	system("rm $bedFasta $bedFile");
	print"$bedFasta $bedFile deleted (use -k to keep these files)\n";
}
print "All genomic regions put in memory, computing and writing results in $folder\n";
close (F_FASTA); # can close the reference genome file
 
if( defined($opts{f})) {
	# if -f option, open fasta file handle
	open(FILEFASTA,">$inputAncient.fasta") or die "Cannot open file: $dbFasta\n";
	foreach my $key (sort {$a <=> $b} keys(%$read))	{
		die("The ref sequence for ".$read->{$key}{name}." wasn't retrieved\n") if (!defined($read->{$key}{before}));
		my $clip5p = 0;
		my $clip3p = 0;
		if(defined($read->{$key}{softClip5})) { # adjust the alignment with clipped bases
			$clip5p = $read->{$key}{softClip5};
		}
		if(defined($read->{$key}{softClip3})) { 
			$clip3p = $read->{$key}{softClip3};
		}
	#	print "".$read->{$key}{name}." CIG $read->{$key}{cigar} std $read->{$key}{strand} lg ".length($read->{$key}{seq})." clip5 $clip5p clip3 $clip3p\n";
		print FILEFASTA ">Ref_before\n".$read->{$key}{before}."\n";
		print FILEFASTA ">Ref_".$read->{$key}{ref}."\n".$read->{$key}{seqRef}."\n";
		# strand is appended to the read name
		print FILEFASTA ">".$read->{$key}{name}."_".$read->{$key}{strand}."\n";
		my $gap = 0;
		while ($gap < length($read->{$key}{before}) - $clip5p){
			print FILEFASTA "-";
			$gap++
		}
		print FILEFASTA "".$read->{$key}{seq}."\n>Ref_after\n";
		$gap = 0;
		while ($gap <(length($read->{$key}{before})+length($read->{$key}{seq})-$clip3p-$clip5p)){
			print FILEFASTA "-";
			$gap++
		}
		print FILEFASTA "".$read->{$key}{after}."\n";
	}
	close (FILEFASTA);
	print"Alignments were written in the fasta file: $inputAncient.fasta\n";
	
}
$cpt = 0;
undef($reference); # useless now


### Compute damage tracking ###

# Compute per chromosome
foreach my $chr (sort {$a cmp $b } keys(%$reference2))	
{
	$cpt++;
	### Composition for each read
	my $file_OUT	= $classN.'_'.$chr.'.txt';
	open(FILE_OUT,">$folder/$file_OUT");
	print FILE_OUT "Seq\tLength\tAl\tCl\tGl\tTl\tdell\tinsl\tAr\tCr\tGr\tTr\tdelr\tinsr\n";

	### Hits per chromosome
	my $filehit = $folder.'/'.$chromo.'_'.$chr.'_'.$around.'_'.$seqLength.'.txt';
	open(FILEHIT,">$filehit");
	print FILEHIT "Chr\tStrd\tNb\tLength\tMax\n";
	my $read2 = {}; 
	# new hash readRev which contains reads and ref in the reverse WITHOUT complementation
	# in order to align left and avoid stair-like effects with reads of different lengths
	my $readRev = {};
	my $cptMinus=0;
	my $cptPlus=0;
	my $totLengthMinus=0;
	my $totLengthPlus=0;
	my $chrLength= 0;
	# fetch chr length
	$chrLength = $reference2->{$chr};
	
	foreach my $key (keys(%$read)) {
			
		# look for read which have mapped to the current chromosome
		if(defined($chr) and defined($read->{$key}{ref})){
			if ( $read->{$key}{ref} eq $chr ){
				
				if(defined($read->{$key}{softClip5})) {$read2->{$key}{softClip5}= $read->{$key}{softClip5} ; }
				if(defined($read->{$key}{softClip3})) {$read2->{$key}{softClip3}= $read->{$key}{softClip3} ; }
				$read2->{$key}{before} = $read->{$key}{before};
				$read2->{$key}{after} = $read->{$key}{after};
				$readRev->{$key}{seq} = reverse($read->{$key}{seq});
				$readRev->{$key}{seqRef} = reverse($read->{$key}{seqRef});
				# count hit per chr on both strands
				if($read->{$key}{strand} == 0){
					$cptPlus++;
					$totLengthPlus = $totLengthPlus + length($read->{$key}{seq});
				}
				else { 
					$cptMinus++;
					$totLengthMinus = $totLengthMinus + length($read->{$key}{seq});
				}
			}
		}
		
	}
	print FILEHIT "$chr\t+\t$cptPlus\t$totLengthPlus\t$chrLength\n"; # write hitPerChrom data
	print FILEHIT "$chr\t-\t$cptMinus\t$totLengthMinus\t$chrLength\n";
	close (FILEHIT);
	
	# hash table which will contain the composition of read and around 
	# nuclCompRev will only contain reversed read to align from 3'
	my %nuclComp;
	my %nuclCompRev;
	foreach my  $strd (0,16){ ### Initializing both hash tables
		foreach my $i (1..$around){
			foreach my $base ('A','C','G','T','-','+','sum'){ # '-' means deletion in the read, '+' means insertion in the read (= gap in the reference sequence)
				$j = 99999+$i;
				${${$nuclComp{-$i}}{$strd}}{$base} = 0;
				${${$nuclComp{$j}}{$strd}}{$base} = 0;	
			}
		}
	}	
	foreach my $base ('A','C','G','T','-','+','sum'){
	
		foreach my $i (1..$seqLength){
		
			foreach my $strd (0,16){
			
				${${$nuclComp{$i}}{$strd}}{$base} = 0;
				# same for reverted
				${${$nuclCompRev{$i}}{$strd}}{$base} = 0;
			
			}
		}
	}	
	

	foreach my $key (keys(%$read2)) # do the damage tracking for the read of current chromosome
	{
		### nuclComp, nuclCompReverse files and readComp files
		
		my %localACGT;	# hash indexed by position, per strand and per base for the read for readComp file
		my %regionACGT;	# idem for region, ie before + read + after for readComp file
		
		foreach my $letter ('A','C','G','T','-','+'){ # initialize Hash
			$localACGT{$letter} = 0;
			$regionACGT{$letter} = 0;
		}
		my @temp = split(//,$read->{$key}{before}); # compute base composition of region BEFORE
		my $strand = $read->{$key}{strand};
		foreach my $i (1..scalar(@temp))	{
			${${$nuclComp{-scalar(@temp)+$i-1}}{$strand}}{$temp[$i-1]}++; 
			$regionACGT{$temp[$i-1]}++;
		}
		@temp = split(//,$read->{$key}{after}); # compute base composition of region AFTER
		foreach my $i (1..scalar(@temp))	{
			${${$nuclComp{99999+$i}}{$strand}}{$temp[$i-1]}++; 
			$regionACGT{$temp[$i-1]}++;
		}
		@temp = split(//,$read->{$key}{seq}); # compute base composition for read sequence
		my @tempRev = reverse(@temp);	# reverse sequence for nuclCompRev
		
		my $newSeqLength = $seqLength; # avoid access to not allocated table, adjust to the read length
		if ( length($read->{$key}{seq}) < $seqLength ) {
			$newSeqLength = length($read->{$key}{seq});
		}
		my $nbDel=0;
		foreach my $i (1..$newSeqLength){
			${${$nuclComp{$i}}{$strand}}{$temp[$i-1]}++;
			${${$nuclCompRev{$i}}{$strand}}{$tempRev[$i-1]}++;
			$localACGT{$temp[$i-1]}++;
			if( $temp[$i-1] eq '-'){
				$nbDel++;	# count deletion to get number of letter in sequence, without gaps
			}
		}
		# Counting of A, C, T, G and deletion ('-') done, now count insertion number ('+') into read sequences
		# which means count the gap number into the REFERENCE sequence adjusted to the read sequence limits
		die("The ref sequence for ".$read->{$key}{name}." wasn't retrieved\n") if (!defined($read->{$key}{seqRef}));
	#	my $seqREF = substr($read->{$key}{seqRef},$around,length($read->{$key}{seq}));
		my $seqREF = substr($read->{$key}{seqRef},length($read->{$key}{before}),length($read->{$key}{seq}));
	#	print "name ".$read->{$key}{name}." lg before ".length($read->{$key}{before})." lg seq ".length($read->{$key}{seq})." lg after ".length($read->{$key}{after})." newseqlg ".$newSeqLength."\n";
		my @tempREF = split(//,$seqREF);
		my $nbIns=0; # reset counter
		foreach my $i (1..$newSeqLength){
			
			$regionACGT{$tempREF[$i-1]}++; # compute reference matches read
			if( $tempREF[$i-1] eq '-'){
				${${$nuclComp{$i}}{$strand}}{'+'}++;
				$localACGT{'+'}++;
				
				$nbIns++;	# count insertion to get number of letter in reference, without gaps
			}
		}	
		# do not count deletion in the length
		my $long = (length($read->{$key}{seq})-$nbDel);
		
		foreach my $letter ('A','C','G','T','-','+'){
			
			if($long != 0){ # compute base computation
	
				$regionACGT{$letter} = $regionACGT{$letter}  / ( $long + (2 * $around) - $nbIns) ;
				$localACGT{$letter}  = $localACGT{$letter}  / $long;
			}
		}
		print FILE_OUT "".$read->{$key}{name}."\t$long";
		foreach my $letter ('A','C','G','T','-','+')	{ 
			if(defined($localACGT{$letter})){printf FILE_OUT "\t%.13f",$localACGT{$letter}; }
		}
		foreach my $letter ('A','C','G','T','-','+')	{ 
			if(defined($regionACGT{$letter})){ printf FILE_OUT "\t%.13f",$regionACGT{$letter}; }
		}
		print FILE_OUT "\n";
	}
	close(FILE_OUT);	

	#   nucleotide composition
	my $fileOUT = $summary.'_'.$chr.'_'.$around.'_'.$seqLength.'.txt';

	open(FILEOUT,">$folder/$fileOUT");
	print FILEOUT "Strd\tPos\tA\tC\tG\tT\tdel\tins\tTot";
	
	foreach my $i (1..$around){ # do the sum for base composition per strand for before and after
		foreach my $base ('A','C','G','T','-','+'){
			foreach my $strd (0,16){
																					  
				if(defined(${${$nuclComp{-$i}}{$strd}}{sum}) and defined(${${$nuclComp{-$i}}{$strd}}{$base})){ # sum for before region
					${${$nuclComp{-$i}}{$strd}}{sum} = ${${$nuclComp{-$i}}{$strd}}{sum} + ${${$nuclComp{-$i}}{$strd}}{$base};
				}
				if(defined(${${$nuclComp{99999+$i}}{$strd}}{sum}) and defined(${${$nuclComp{99999+$i}}{$strd}}{$base})){ # sum for after region
					${${$nuclComp{99999+$i}}{$strd}}{sum} = ${${$nuclComp{99999+$i}}{$strd}}{sum} + ${${$nuclComp{99999+$i}}{$strd}}{$base};
				}
			}
		}
	}
	# now write result in the file
	foreach my $strd (0,16){
		foreach my $i (1..$around)	{
			print FILEOUT "\n$strd\t-$i"; 
			foreach my $base ('A','C','G','T','-','+','sum')	{ # for before region
				if(defined(${${$nuclComp{-$i}}{$strd}}{$base})){
					print FILEOUT "\t${${$nuclComp{-$i}}{$strd}}{$base}"; 
				}
			}
			$j = 99999+$i;
			print FILEOUT "\n$strd\t$j"; 
			foreach my $base ('A','C','G','T','-','+','sum')	{ # for after region
				if(defined(${${$nuclComp{$j}}{$strd}}{$base})){
					print FILEOUT "\t${${$nuclComp{$j}}{$strd}}{$base}"; 
				}
			}
		}
	}
	# compute sum for sequence
	foreach my $i (1..$seqLength)
	{
		foreach my $base ('A','C','G','T','-','+')
		{
			foreach my $strd (0,16)
			{
				if(defined(${${$nuclComp{$i}}{$strd}}{sum}) and defined( ${${$nuclComp{$i}}{$strd}}{$base})) {
					${${$nuclComp{$i}}{$strd}}{sum} = ${${$nuclComp{$i}}{$strd}}{sum} + ${${$nuclComp{$i}}{$strd}}{$base};
				}
				if(defined(${${$nuclCompRev{$i}}{$strd}}{sum}) and defined( ${${$nuclCompRev{$i}}{$strd}}{$base})) {
					${${$nuclCompRev{$i}}{$strd}}{sum} = ${${$nuclCompRev{$i}}{$strd}}{sum} + ${${$nuclCompRev{$i}}{$strd}}{$base};
				}
			}
		}
	}
	# write results for nucleotide composition by strand
	foreach my $strd (0,16)
	{
		foreach my $i (1..$seqLength)	
		{
			print FILEOUT "\n$strd\t$i"; 
			foreach my $base ('A','C','G','T','-','+','sum')	{ 
				if(defined(${${$nuclComp{$i}}{$strd}}{$base})){
					print FILEOUT "\t${${$nuclComp{$i}}{$strd}}{$base}"; 
				}
			}
		}
	}
	# write results for nucleotide composition by strand for reverse sequence
	# arbitratry choose that this reversed sequence is from position 200000 to 200000 + $seqLength
	foreach my $strd (0,16)
	{
		foreach my $i (1..$seqLength)	
		{
			print FILEOUT "\n$strd\t".($i+199999).""; 
			foreach my $base ('A','C','G','T','-','+','sum')	{ 
				if(defined(${${$nuclCompRev{$i}}{$strd}}{$base})){
					print FILEOUT "\t${${$nuclCompRev{$i}}{$strd}}{$base}"; 
				}
			}
		}
	}
	close(FILEOUT);

	# Misincorporations
	my %nuclmis;
	my %basesG;
	my %nuclmisRev;
	my %basesGRev;
	my $fileMismOUT	= $types.'_'.$chr.'.txt';
	my $fileMismRevOUT	= $typesRev.'_'.$chr.'.txt';
	open(FILEMISMOUT,">$folder/$fileMismOUT");
	open(FILEMISMREVOUT,">$folder/$fileMismRevOUT");
	# initialising
	foreach my $key  ('A','C','G','T')	{ 
		foreach my $position (0..$seqLength-1)	{ 
			${$basesG{$position}}{$key} = 0;
			${$basesGRev{$position}}{$key} = 0;
		} 
	}
	foreach my  $key2 ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G','T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S') {
		foreach my $position (0..$seqLength-1)	{ 
			${$nuclmis{$key2}}{$position} = 0 ; 
			${$nuclmisRev{$key2}}{$position} = 0 ; 
		}
	} 
	
	# file open for mismatch pattern
	$fileOUT = $folder.'/'.$pattern.'_'.$chr.'_'.$around.'_'.$seqLength.'.txt';
	open(FILEOUT,">$fileOUT");
	my %mismPos;
	my %total;
	my $fileTOT = $folder.'/'.$misinc.'_'.$chr.'_'.$around.'_'.$seqLength.'.txt';
	open(FILETOT,">$fileTOT");
	
	foreach my $strand(0,16){
		foreach my  $key2 ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G','T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S') {
			foreach my $position (0..$seqLength-1)	{ 
				${${$mismPos{$key2}}{$position}}{$strand} = 0;
	
			}
		}
	}
	# initializing hash tables to avoid missing values
	foreach my $strand (0,16) {
		foreach my $key2 ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G','T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S','total')
		{
			${$total{$key2}}{$strand}=0;
		}
	}
		
	foreach my $key (keys(%$read2)){
		
		my $strd = $read->{$key}{strand};	
		my $seqREF = substr($read->{$key}{seqRef},length($read->{$key}{before}),length($read->{$key}{seq}));
		my $seqREFRev = substr($readRev->{$key}{seqRef},length($read->{$key}{before}),length($readRev->{$key}{seq}));
		my @ref = split(//,$seqREF);
		my @seq = split(//,$read->{$key}{seq});
		my @refRev = split(//,$seqREFRev);
		my @seqRev = split(//,$readRev->{$key}{seq});
		my $clip5p =0; 
		my $clip3p=0; 
		if(defined($read->{$key}{softClip5})) { 
			$clip5p = $read->{$key}{softClip5} ;
			for (my $i=0; $i < $clip5p; $i++) {
				${${$mismPos{'S'}}{$i}}{0}++;	
				${$nuclmis{'S'}}{$i}++;
				${$total{'S'}}{0}++;	
				${$total{total}}{0}++;					
			}
		}
		if(defined($read->{$key}{softClip3})) { 
			$clip3p = $read->{$key}{softClip3} ;
			for (my $i=0; $i < $clip3p; $i++) {
				${${$mismPos{'S'}}{$i}}{16}++;	
				${$total{'S'}}{16}++;		
				${$nuclmisRev{'S'}}{$i}++;
				${$total{total}}{16}++;	
			}			
			
		}
		my $newSeqLength = $seqLength; # avoid access to not allocated table, adjust to the read length
		if ( length($read->{$key}{seq}) < $seqLength ) {
			$newSeqLength = length($read->{$key}{seq});
		}
	#	print "ref @ref $clip5p\nseq @seq $clip3p\n";
		foreach my $position (0..$newSeqLength-1- $clip5p - $clip3p){
			if(defined($ref[$position]) and defined($refRev[$position]) and defined($seq[$position+$clip5p] ) and defined($seqRev[$position+$clip3p] ) ){
				${$basesG{$position}}{$ref[$position]}++;	# compute base composition of reference
				${$basesGRev{$position}}{$refRev[$position]}++;
				# if base are different AND comparable = a letter in read sequence
				if (( $ref[$position]  ne $seq[$position+$clip5p] ) and ($seq[$position+$clip5p] ne "")) {
					my $mis = $ref[$position] .">".$seq[$position+$clip5p];	# hash is indexed by the mismatches ie letter>other letter. Letter are A, C, G, T and -. N are not taken into account
					${${$mismPos{$mis}}{$position}}{$strd}++;	
					${$nuclmis{$mis}}{$position}++;
					${$total{$mis}}{$strd}++;						
					${$total{total}}{$strd}++;	
				}
				# same for reverse nuclComp
				if (( $refRev[$position]  ne $seqRev[$position+$clip3p] ) and ($seqRev[$position+$clip3p] ne "")) {
					my $mis = $refRev[$position] .">".$seqRev[$position+$clip3p];				
					${$nuclmisRev{$mis}}{$position}++;
				}
			}
			else {
				die("ref $ref[$position] refrev $refRev[$position] seq $seq[$position+$clip5p] revseq $seqRev[$position+$clip3p] not defined 5s $clip5p 3s $clip3p $read->{$key}{name}");
			}
		}		
	}
	print FILEMISMOUT "Pos\tA\tC\tG\tT\tTot\tG>A\tC>T\tA>G\tT>C\tA>C\tA>T\tC>G\tC>A\tT>G\tT>A\tG>C\tG>T\tA>-\tT>-\tC>-\tG>-\t->A\t->T\t->C\t->G\tS\n";
	print FILEMISMREVOUT "Pos\tA\tC\tG\tT\tTot\tG>A\tC>T\tA>G\tT>C\tA>C\tA>T\tC>G\tC>A\tT>G\tT>A\tG>C\tG>T\tA>-\tT>-\tC>-\tG>-\t->A\t->T\t->C\t->G\tS\n";

	foreach my $position (sort { $a <=> $b } (keys(%basesG)))
	{
		print FILEMISMOUT $position+1; 
		print FILEMISMREVOUT $position+1; 
		my $sumACGT = 0;
		my $sumRevACGT = 0;
		
		foreach my $key ('A','C','G','T')	{ 
			$sumACGT = $sumACGT + ${$basesG{$position}}{$key}; 
			$sumRevACGT = $sumRevACGT + ${$basesGRev{$position}}{$key}; 
			print FILEMISMOUT "\t${$basesG{$position}}{$key}"; 
			print FILEMISMREVOUT "\t${$basesGRev{$position}}{$key}"; 
		}
		print FILEMISMOUT "\t$sumACGT";
		print FILEMISMREVOUT "\t$sumRevACGT";
		foreach my $key2 ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G','T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S') {
			print FILEMISMOUT "\t${$nuclmis{$key2}}{$position}";
			print FILEMISMREVOUT "\t${$nuclmisRev{$key2}}{$position}";
		}
		print FILEMISMOUT "\n";
		print FILEMISMREVOUT "\n";
	}
	
	close(FILEMISMOUT);	
	close(FILEMISMREVOUT);
	
	# write results
	print  FILEOUT "Strd\tPos";
	print FILETOT "Strd";
	foreach my $key  ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G','T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S') {
		print FILEOUT "\t$key";
		print FILETOT "\t$key";

	}
	print  FILEOUT "\n";
	print FILETOT "\tTotal\n";
	foreach my $strd (0,16){
			
		print FILETOT "$strd";
		foreach my $position (0..$seqLength-1){
			my $position2 = $position + 1; # NB. 0-offset shifted in 1-offset for final output 
			print FILEOUT "$strd\t$position2";
			foreach my $key ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G','T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S') {
			
				print FILEOUT "\t".${${$mismPos{$key}}{$position}}{$strd}."";
			}	
			print FILEOUT "\n";
		}
		foreach my $key2 ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G','T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S') {
					
			print FILETOT "\t".${$total{$key2}}{$strd}."";
		}
		print FILETOT "\t".${$total{total}}{$strd}."\n";
	}
	close (FILETOT);
	close (FILEOUT);
	# undef to be sure that read2 will be empty for the next analyzed chromosome
	undef(%$read2);
	printf ("Computation done for reference %5s [%2d/%2d]\n", $chr, $cpt,scalar(keys(%$reference2)) );
}

if(defined($opts{c}) ) {
		
	my $plotl; my $plota; my $ploty; my $plott;
	if( defined($opts{m})) { $plotl = $opts{m};} else {$plotl = 25;}
	if( defined($opts{b})) { $plota = $opts{b};} else {$plota = 10;}
	if( defined($opts{y})) { $ploty = $opts{y};} else {$ploty = 0.3;}
	if( defined($opts{t})) { $plott = $opts{t};} else {$plotl = "plot";}
	&merge($folder, $plotl, $plota, $ploty, $plott);
}
else {
	print "\nmap terminated, you may use $0 merge -d $folder\n";
}
exit 0;
}
else{
	&printUsageMap;
}
}

## subroutine merge

sub merge {

my %opts = (d=>undef, c=>undef);

if(defined($_[0]) ) {
	$opts{d}= $_[0];
	$opts{c}= 1;
}

getopts('d:',\%opts);

if(defined($opts{d}) ) {

my $folder 	= $opts{d};
print "Merging files from $folder\n";

my $fileNC = 'nuclComp.all.txt';
my $fileNCR = 'nuclCompReverse.all.txt';
my $fileDF = 'dnaFrag.all.txt';
my $fileHIT = 'hitPerChrom.all.txt';
my $fileRCOM = 'readComp.all.txt';
my $fileSUM = 'nuclSumMis.all.txt';
my $fileMIS = 'nuclMis.all.txt';

### Reading sequences in $folder
opendir(FILES, "$folder") or die "Can't open folder $folder\n";
my @files = readdir(FILES); 
closedir(FILES);

my $hit=(); my $dFrag=(); my $mis =(); my $sumMis=(); my $nComp =(); my $nCompR =(); my $rComp =();

foreach my $file (@files) { 

	&checkFileEmpty($file);
	# regex modified to accept more diverse reference names
	if ( $file =~ m/^nuclComp_([0-9a-zA-Z]+).*\.txt/ ){		
		$nComp = ReadFile($folder,$file,$nComp,$1); # add chromosome name 
	}
	if ( $file =~ m/^nuclCompReverse_([0-9a-zA-Z]+).*\.txt/ ){		
		$nCompR = ReadFile($folder,$file,$nCompR,$1);
	}
	if ( $file =~ m/^dnaFrag_(.+)_[0-9]+_[0-9]+\.txt/ ){		
		$dFrag = ReadFile($folder,$file,$dFrag,$1);
	}
	if ( $file =~ m/^hitPerChrom_/ ){		
		$hit = ReadFile($folder,$file,$hit,"");
	}
	if ( $file =~ m/^readComp_(.+)\.txt/ ){		
		$rComp = ReadFile($folder,$file,$rComp,$1);
	}
	if ( $file =~ m/^nuclMis_(.+)_[0-9]+_[0-9]+\.txt/ ){		
		$mis = ReadFile($folder,$file,$mis,$1);
	}
	if ( $file =~ m/^nuclSumMis_(.+)_[0-9]+_[0-9]+\.txt/){		
		$sumMis = ReadFile($folder,$file,$sumMis,$1);
	}			
}

# nuclComp
open(FILEOUT,">$folder/$fileNC");
print FILEOUT "Chr\tPos\tA\tC\tG\tT\tTot\tG>A\tC>T\tA>G\tT>C\tA>C\tA>T\tC>G\tC>A\tT>G\tT>A\tG>C\tG>T\tA>-\tT>-\tC>-\tG>-\t->A\t->T\t->C\t->G\tS\n";
foreach my $line (@$nComp)	{ print FILEOUT $line; }
close(FILEOUT);
print "Printing $fileNC done\n";

#nuclCompReverse
open(FILEOUT,">$folder/$fileNCR");
print FILEOUT "Chr\tPos\tA\tC\tG\tT\tTot\tG>A\tC>T\tA>G\tT>C\tA>C\tA>T\tC>G\tC>A\tT>G\tT>A\tG>C\tG>T\tA>-\tT>-\tC>-\tG>-\t->A\t->T\t->C\t->G\tS\n";
foreach my $line (@$nCompR)	{ print FILEOUT $line; }
close(FILEOUT);
print "Printing $fileNCR done\n";

# hitPerChromosome
open(FILEOUT,">$folder/$fileHIT");
print FILEOUT "Chr\tStrd\tNb\tLength\tMax\n";
foreach my $line (@$hit)	{ print FILEOUT $line; }
close(FILEOUT);
print "Printing $fileHIT done\n";

# dnaFragmentation
open(FILEOUT,">$folder/$fileDF");
print FILEOUT "Chr\tStrd\tPos\tA\tC\tG\tT\tdel\tins\tTot\n";
foreach my $line (@$dFrag)	{ print FILEOUT $line; }
close(FILEOUT);
print "Printing $fileDF done\n";

# read Composition
open(FILEOUT,">$folder/$fileRCOM");
print FILEOUT "Chr\tSeq\tLength\tAl\tCl\tGl\tTl\tdell\tinsl\tAr\tCr\tGr\tTr\tdelr\tinsr\n";
foreach my $line (@$rComp)	{ print FILEOUT $line; }
close(FILEOUT);
print "Printing $fileRCOM done\n";

# Misincorporation
open(FILEOUT,">$folder/$fileMIS");
print FILEOUT "Chr\tStrd\tPos\tG>A\tC>T\tA>G\tT>C\tA>C\tA>T\tC>G\tC>A\tT>G\tT>A\tG>C\tG>T\tA>-\tT>-\tC>-\tG>-\t->A\t->T\t->C\t->G\tS\n";
foreach my $line (@$mis)	{  print FILEOUT $line; }
close(FILEOUT);

# Sum Misincorporation
open(FILEOUT,">$folder/$fileSUM");
print "Printing $fileSUM ...\n";
print FILEOUT "Chr\tStrd\tG>A\tC>T\tA>G\tT>C\tA>C\tA>T\tC>G\tC>A\tT>G\tT>A\tG>C\tG>T\tA>-\tT>-\tC>-\tG>-\t->A\t->T\t->C\t->G\tS\tTotal\n";
foreach my $line (@$sumMis)	{ print FILEOUT $line; }
close(FILEOUT);
print "Printing $fileMIS done\n";

if(defined($opts{c}) ) {
	
	my $plotl; my $plota; my $ploty; my $plott;
	if( defined($_[1])) { $plotl = $_[1];} else {$plotl = 25;}
	if( defined($_[2])) { $plota = $_[2];} else {$plota = 10;}
	if( defined($_[3])) { $ploty = $_[3];} else {$ploty = 0.3;}
	if( defined($_[4])) { $plott = $_[4]; } else {$plott = "plot";}
	&plot($folder, $plotl, $plota, $ploty, $plott);
}
else {
	print "\nmerge terminated, you may use $0 plot -d $folder\n"
}
exit 0;
}
else{
	&printUsageMerge;
}
}

## subroutine merge

sub plot {

my $output = `R --version`;
die ("R is not present in your \$PATH, please install it on your system\n") if (!defined($output) or $output  !~ /version/) ;

my %opts = (d=>undef, l=>25, a=>10, m=>0.3, t=>undef);

getopts('d:l:a:m:t:',\%opts);

if(defined($_[0]) ) {
	$opts{d}= $_[0];
	if( defined($_[1])) { $opts{l} = $_[1];}
	if( defined($_[2])) { $opts{a} = $_[2];}
	if( defined($_[3])) { $opts{m} = $_[3];}
	if( defined($_[4])) { $opts{t} = $_[4];}
}
#print "length $opts{l} title set with $opts{t} and y-limit $opts{m}\n";
if(defined($opts{d}) and -d $opts{d}) {
    
    
my $folder = $opts{d};
my $title ;
if( $opts{d} =~ m/results_(.+)\.[sb]am/) {
	$title = $1;
}
else {
	$title = "plot";
}

if(defined($opts{t})) {
	$opts{t} =~  s/[:\|\s]//g;
	$title = $opts{t};
}

if($opts{a} != 20){
	if($opts{a} !~/^[0-9]+$/ ) {die("around is expected to be an integrer\n");}
	if($opts{a} > 100) {print"Warning: the around option was set to $opts{a}, are you sure to plot such large regions? continuing anyway\n";}
	if($opts{a} < 1) {die("The around position must be > 1 nucleotide\n");}
}
if($opts{l} != 70){
	if($opts{l} !~/^[0-9]+$/) {die("length is expected to be an integrer\n");}
	if($opts{l} < 1) {die("The length must be > 1 nucleotide\n");}
	if($opts{l} > 99998) {die("The read length choosen, $opts{l}, is too large. Must be < 99998 nucleotides\n");}
}
if($opts{m} != 0.3){
	if($opts{m} !~/^[0-9\.]+$/) {die("y plot limit is expected to be a decimal number\n");}
	if($opts{m} <= 0 or $opts{m} > 1) {die("y plot limit must be >= 0 and < 1\n");}
}

my $file = "mapDamage.R";
open(FILEF,">$file") or die("Cannot open file $file\n");
print FILEF "com<-read.csv(file=\"$folder/dnaFrag.all.txt\",sep=\"\\t\")\n";
print FILEF "idx <- (com\$Pos > 99999) & (com\$Pos < 200000)\n";
print FILEF "com\$Pos[idx] <- com\$Pos[idx] - 299999\n";
print FILEF "idx <- (com\$Pos > 199999)\n";
print FILEF "com\$Pos[idx] = -com\$Pos[idx] - 1\n";
print FILEF "pdf(file=\"FragMisincorporation_$title.pdf\", title=\"mapDamage-$version $title\")\npar(oma=c(4,2,2,2),mar=c(1,2,1,1))\nlayout(matrix(c(1,2,3,4,5,6,7,8,9,9,10,10), 3, 4, byrow=TRUE))\n";
# for letter A, before - read
print FILEF "n<-numeric\nplot(com\$Pos,com\$A/com\$Tot,pch='.',xlim=c(-$opts{a},$opts{a}),ylim=c(0,0.5),col=\"green\",main=\"A\",";
print FILEF "cex.axis=0.8,las=2,xlab=\"\",ylab=\"\",lab=c(".(2*$opts{a}).",6,0.2), axes=FALSE)\n";
print FILEF "axis(side=2,labels=TRUE,line=0,las=2,cex.axis=0.8)\naxis(side=1,labels=FALSE)\n";
print FILEF "mtext(\"Frequency\",side=2,line=2.5,cex=0.7)\n";
print FILEF "rect(0.5,0,$opts{a}.5,0.5,border=\"darkgrey\")\nsegments($opts{a}.5,0,$opts{a}.5,0.5,col=\"white\",lwd=2)\n";
print FILEF "for (i in c(-$opts{a}:-1,1:$opts{a})) { n<-c(n,mean(com\$A[(com\$Pos==i)]/com\$Tot[(com\$Pos==i)],na.rm=T)) }\n";
print FILEF "points(c(-$opts{a}:-1,1:$opts{a}),n[2:".(2*$opts{a}+1)."],pch=20,col=\"green\",type=\"b\")\n";
# for letter A, read - after
print FILEF "n<-numeric\nplot(com\$Pos,com\$A/com\$Tot,pch='.',xlim=c(".(-200000-$opts{a}).",".(-200000+$opts{a})."),ylim=c(0,0.5),col=\"green\",main=\"A\",";
print FILEF "cex.axis=0.8,las=2,xlab=\"\",ylab=\"\",lab=c(".(2*$opts{a}).",6,0.2), axes=FALSE)\n";
print FILEF "axis(side=4,labels=FALSE,line=0,las=2,cex.axis=0.8)\naxis(side=1,labels=FALSE)\n";
print FILEF "rect(".(-200000-$opts{a}).".5,0,".(-200000).".5,0.5,border=\"darkgrey\")\nsegments(".(-200000-$opts{a}).".5,0,".(-200000-$opts{a}).".5,0.5,col=\"white\",lwd=2)\n";
print FILEF "for (i in c(".(-200000-$opts{a}).":-200001,-199999:".(-200000+$opts{a}).")) { n<-c(n,mean(com\$A[(com\$Pos==i)]/com\$Tot[(com\$Pos==i)],na.rm=T)) }\n";
print FILEF "points(c(".(-200000-$opts{a}).":-200001,-199999:".(-200000+$opts{a})."),n[2:".(2*$opts{a}+1)."],pch=20,col=\"green\",type=\"b\")\n";
print FILEF "mtext(\"$title\", side=3, line=1.2, cex=0.8)\n";
#print FILEF "text(labels=c(\"$title\"), c(2,2),c(3,3))\n";
# for letter C, before - read
print FILEF "n<-numeric\nplot(com\$Pos,com\$C/com\$Tot,pch='.',xlim=c(-$opts{a},$opts{a}),ylim=c(0,0.5),col=\"blue\",main=\"C\",";
print FILEF "cex.axis=0.8,las=2,xlab=\"\",ylab=\"\",lab=c(".(2*$opts{a}).",6,0.2), axes=FALSE)\n";
print FILEF "axis(side=2,labels=FALSE,line=0,las=2,cex.axis=0.8)\naxis(side=1,labels=FALSE)\n";
print FILEF "rect(0.5,0,$opts{a}.5,0.5,border=\"darkgrey\")\nsegments($opts{a}.5,0,$opts{a}.5,0.5,col=\"white\",lwd=2)\n";
print FILEF "for (i in c(-$opts{a}:-1,1:$opts{a})) { n<-c(n,mean(com\$C[(com\$Pos==i)]/com\$Tot[(com\$Pos==i)],na.rm=T)) }\n";
print FILEF "points(c(-$opts{a}:-1,1:$opts{a}),n[2:".(2*$opts{a}+1)."],pch=20,col=\"blue\",type=\"b\")\n";
# for letter C, read - after
print FILEF "n<-numeric\nplot(com\$Pos,com\$C/com\$Tot,pch='.',xlim=c(".(-200000-$opts{a}).",".(-200000+$opts{a})."),ylim=c(0,0.5),col=\"blue\",main=\"C\",";
print FILEF "cex.axis=0.8,las=2,xlab=\"\",ylab=\"\",lab=c(".(2*$opts{a}).",6,0.2), axes=FALSE)\n";
print FILEF "axis(side=4,labels=TRUE,line=0,las=2,cex.axis=0.8)\naxis(side=1,labels=FALSE)\n";
print FILEF "rect(".(-200000-$opts{a}).".5,0,".(-200000).".5,0.5,border=\"darkgrey\")\nsegments(".(-200000-$opts{a}).".5,0,".(-200000-$opts{a}).".5,0.5,col=\"white\",lwd=2)\n";
print FILEF "for (i in c(".(-200000-$opts{a}).":-200001,-199999:".(-200000+$opts{a}).")) { n<-c(n,mean(com\$C[(com\$Pos==i)]/com\$Tot[(com\$Pos==i)],na.rm=T)) }\n";
print FILEF "points(c(".(-200000-$opts{a}).":-200001,-199999:".(-200000+$opts{a})."),n[2:".(2*$opts{a}+1)."],pch=20,col=\"blue\",type=\"b\")\n";
# for letter G, before - read
print FILEF "n<-numeric\nplot(com\$Pos,com\$G/com\$Tot,pch='.',xlim=c(-$opts{a},$opts{a}),ylim=c(0,0.5),col=\"black\",main=\"G\",";
print FILEF "cex.axis=0.8,las=2,xlab=\"\",ylab=\"\",lab=c(".(2*$opts{a}).",6,0.2),axes=FALSE)\n";
print FILEF "axis(side=2,labels=TRUE,line=0,las=2,cex.axis=0.8)\naxis(side=1,labels=TRUE,las=2,cex.axis=0.8)\n";
print FILEF "rect(0.5,0,$opts{a}.5,0.5,border=\"darkgrey\")\nsegments($opts{a}.5,0,$opts{a}.5,0.5,col=\"white\",lwd=2)\n";
print FILEF "mtext(\"Frequency\",side=2,line=2.5,cex=0.7)\n";
print FILEF "for (i in c(-$opts{a}:-1,1:$opts{a})) { n<-c(n,mean(com\$G[(com\$Pos==i)]/com\$Tot[(com\$Pos==i)],na.rm=T)) }\n";
print FILEF "points(c(-$opts{a}:-1,1:$opts{a}),n[2:".(2*$opts{a}+1)."],pch=20,col=\"black\",type=\"b\")\n";
# for letter G, read - after
print FILEF "n<-numeric\nplot(com\$Pos,com\$G/com\$Tot,pch='.',xlim=c(".(-200000-$opts{a}).",".(-200000+$opts{a})."),ylim=c(0,0.5),col=\"black\",main=\"G\",";
print FILEF "cex.axis=0.8,las=2,xlab=\"\",ylab=\"\",lab=c(".(2*$opts{a}).",6,0.2), axes=FALSE)\n";
print FILEF "axis(side=4,labels=FALSE)\n";
print FILEF "axis(side=1,labels=seq(from=-$opts{a},to=$opts{a},by=1),at=seq(from=".(-200000-$opts{a}).",to=".(-200000+$opts{a}).",by=1),las=2,cex.axis=0.8)\n";
print FILEF "rect(".(-200000-$opts{a}).".5,0,".(-200000).".5,0.5,border=\"darkgrey\")\nsegments(".(-200000-$opts{a}).".5,0,".(-200000-$opts{a}).".5,0.5,col=\"white\",lwd=2)\n";
print FILEF "for (i in c(".(-200000-$opts{a}).":-200001,-199999:".(-200000+$opts{a}).")) { n<-c(n,mean(com\$G[(com\$Pos==i)]/com\$Tot[(com\$Pos==i)],na.rm=T)) }\n";
print FILEF "points(c(".(-200000-$opts{a}).":-200001,-199999:".(-200000+$opts{a})."),n[2:".(2*$opts{a}+1)."],pch=20,col=\"black\",type=\"b\")\n";
# for letter T, before - read
print FILEF "n<-numeric\nplot(com\$Pos,com\$T/com\$Tot,pch='.',xlim=c(-$opts{a},$opts{a}),ylim=c(0,0.5),col=\"red\",main=\"T\",";
print FILEF "cex.axis=0.8,las=2,xlab=\"\",ylab=\"\",lab=c(".(2*$opts{a}).",6,0.2), axes=FALSE)\n";
print FILEF "axis(side=2,labels=FALSE)\naxis(side=1,labels=TRUE,las=2,cex.axis=0.8)\n";
print FILEF "rect(0.5,0,$opts{a}.5,0.5,border=\"darkgrey\")\nsegments($opts{a}.5,0,$opts{a}.5,0.5,col=\"white\",lwd=2)\n";
print FILEF "for (i in c(-$opts{a}:-1,1:$opts{a})) { n<-c(n,mean(com\$T[(com\$Pos==i)]/com\$Tot[(com\$Pos==i)],na.rm=T)) }\n";
print FILEF "points(c(-$opts{a}:-1,1:$opts{a}),n[2:".(2*$opts{a}+1)."],pch=20,col=\"red\",type=\"b\")\n";
# for letter T, read - after
print FILEF "n<-numeric\nplot(com\$Pos,com\$T/com\$Tot,pch='.',xlim=c(".(-200000-$opts{a}).",".(-200000+$opts{a})."),ylim=c(0,0.5),col=\"red\",main=\"T\",";
print FILEF "cex.axis=0.8,las=2,xlab=\"\",ylab=\"\",lab=c(".(2*$opts{a}).",6,0.2),axes=FALSE)\n";
print FILEF "axis(side=4,labels=TRUE,line=0,las=2,cex.axis=0.8)\n";
print FILEF "axis(side=1,labels=seq(from=-$opts{a},to=$opts{a},by=1),at=seq(from=".(-200000-$opts{a}).",to=".(-200000+$opts{a}).",by=1),las=2,cex.axis=0.8)\n";
print FILEF "rect(".(-200000-$opts{a}).".5,0,".(-200000).".5,0.5,border=\"darkgrey\")\nsegments(".(-200000-$opts{a}).".5,0,".(-200000-$opts{a}).".5,0.5,col=\"white\",lwd=2)\n";
print FILEF "for (i in c(".(-200000-$opts{a}).":-200001,-199999:".(-200000+$opts{a}).")) { n<-c(n,mean(com\$T[(com\$Pos==i)]/com\$Tot[(com\$Pos==i)],na.rm=T)) }\n";
print FILEF "points(c(".(-200000-$opts{a}).":-200001,-199999:".(-200000+$opts{a})."),n[2:".(2*$opts{a}+1)."],pch=20,col=\"red\",type=\"b\")\n\n\n";


### Misincorporation patterns
print FILEF "nucl<-read.csv(file=\"$folder/nuclComp.all.txt\",sep=\"\\t\",check.names = FALSE)\n";
print FILEF "vec<-(1:$opts{l})\nsumhit<-data.frame(cbind(vec,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))\n";
#print FILEF "idx <- nucl\$A==0\n";
#print FILEF "nucl\$A[idx] <- 10^10\n";
#print FILEF "idx <- nucl\$C==0\n";
#print FILEF "nucl\$C[idx] <- 10^10\n";
#print FILEF "idx <- nucl\$G==0\n";
#print FILEF "nucl\$G[idx] <- 10^10\n";
#print FILEF "idx <- nucl\$T==0\n";
#print FILEF "nucl\$T[idx] <- 10^10\n";
#print FILEF "idx <- nucl\$Tot==0\n";
#print FILEF "nucl\$Tot[idx] <- 10^10\n";
print FILEF "insertion<-0\n";
print FILEF "for (i in vec)	{\n";
print FILEF "sumhit[i,2]<-sum(nucl[(nucl\$Pos == i), \"G>A\"])/sum(nucl[(nucl\$Pos == i), \"G\"])\n";
print FILEF "sumhit[i,3]<-sum(nucl[(nucl\$Pos == i), \"C>T\"])/sum(nucl[(nucl\$Pos == i), \"C\"])\n";
print FILEF "sumhit[i,4]<-sum(nucl[(nucl\$Pos == i), \"A>G\"])/sum(nucl[(nucl\$Pos == i), \"A\"])\n";
print FILEF "sumhit[i,5]<-sum(nucl[(nucl\$Pos == i), \"T>C\"])/sum(nucl[(nucl\$Pos == i), \"T\"])\n";
print FILEF "sumhit[i,6]<-sum(nucl[(nucl\$Pos == i), \"A>C\"])/sum(nucl[(nucl\$Pos == i), \"A\"])\n";
print FILEF "sumhit[i,7]<-sum(nucl[(nucl\$Pos == i), \"A>T\"])/sum(nucl[(nucl\$Pos == i), \"A\"])\n";
print FILEF "sumhit[i,8]<-sum(nucl[(nucl\$Pos == i), \"C>G\"])/sum(nucl[(nucl\$Pos == i), \"C\"])\n";
print FILEF "sumhit[i,9]<-sum(nucl[(nucl\$Pos == i), \"C>A\"])/sum(nucl[(nucl\$Pos == i), \"C\"])\n";
print FILEF "sumhit[i,10]<-sum(nucl[(nucl\$Pos == i), \"T>G\"])/sum(nucl[(nucl\$Pos == i), \"T\"])\n";
print FILEF "sumhit[i,11]<-sum(nucl[(nucl\$Pos == i), \"T>A\"])/sum(nucl[(nucl\$Pos == i), \"T\"])\n";
print FILEF "sumhit[i,12]<-sum(nucl[(nucl\$Pos == i), \"G>C\"])/sum(nucl[(nucl\$Pos == i), \"G\"])\n";
print FILEF "sumhit[i,13]<-sum(nucl[(nucl\$Pos == i), \"G>T\"])/sum(nucl[(nucl\$Pos == i), \"G\"])\n";
print FILEF "sumhit[i,14]<-sum(nucl[(nucl\$Pos == i), \"A>-\"])/sum(nucl[(nucl\$Pos == i), \"A\"])\n";
print FILEF "sumhit[i,15]<-sum(nucl[(nucl\$Pos == i), \"T>-\"])/sum(nucl[(nucl\$Pos == i), \"T\"])\n";
print FILEF "sumhit[i,16]<-sum(nucl[(nucl\$Pos == i), \"C>-\"])/sum(nucl[(nucl\$Pos == i), \"C\"])\n";
print FILEF "sumhit[i,17]<-sum(nucl[(nucl\$Pos == i), \"G>-\"])/sum(nucl[(nucl\$Pos == i), \"G\"])\n";
print FILEF "sumhit[i,22]<-sum(nucl[(nucl\$Pos == i), \"S\"])/sum(nucl[(nucl\$Pos == i), \"Tot\"])\n";
print FILEF "NTs <- c(\"A\", \"C\", \"G\", \"T\")\n";
print FILEF "for (j in seq(NTs)) {\n";
print FILEF "insertion[j]<-sprintf(\"->%s\", NTs[j])\n";
print FILEF "sumhit[i,17 + j] <-sum(nucl[(nucl\$Pos == i), insertion[j]]) /sum(nucl[(nucl\$Pos == i), NTs])\n";
print FILEF "}\n}\n";
# write a table for C>T frequencies
print FILEF "write.table(sumhit[,3], file=\"$folder/5pCtoT_freq.txt\", sep=\"\\t\", quote=FALSE, col.names=c(\"pos\\t5pC>T\"))\n";
print FILEF "plot(sumhit\$vec,sumhit\$V4,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1,type=\"l\",xlab=\"\",ylab=\"\",cex.lab=0.8,cex.axis=0.6,las=2,lab=c($opts{l},11,0.1),axes=FALSE)\n";
print FILEF "axis(side=1,labels=TRUE,las=2,cex.axis=0.8)\n";
if( $opts{m} > 0.1) { print FILEF "axis(side=2,labels=TRUE,las=2,cex.axis=0.8)\n"; }
else { print FILEF "axis(side=2,labels=seq(from=0,to=$opts{m},by=0.01),at=seq(from=0,to=$opts{m},by=0.01),las=2,cex.axis=0.8)\n"; }
print FILEF "rect(0.5,-0.01,$opts{l}.5,$opts{m},border=\"darkgrey\")\nsegments($opts{l}.5,-0.01,$opts{l}.5,$opts{m},col=\"white\",lwd=2)\n";
print FILEF "mtext(\"Frequency\",side=2,line=2.8,cex=0.7)\n";
print FILEF"lines(sumhit\$vec,sumhit\$V5,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF"lines(sumhit\$vec,sumhit\$V6,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF"lines(sumhit\$vec,sumhit\$V7,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF"lines(sumhit\$vec,sumhit\$V8,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF"lines(sumhit\$vec,sumhit\$V9,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF"lines(sumhit\$vec,sumhit\$V10,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF"lines(sumhit\$vec,sumhit\$V11,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF"lines(sumhit\$vec,sumhit\$V12,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF"lines(sumhit\$vec,sumhit\$V13,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
# soft clipping in orange
print FILEF"lines(sumhit\$vec,(sumhit\$V22),xlim=c(1,5),ylim=c(0,$opts{m}),col=\"orange\",lwd=1)\n";
# deletions in green
print FILEF"lines(sumhit\$vec,(sumhit\$V14+sumhit\$V15+sumhit\$V16+sumhit\$V17),xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"green\",lwd=1)\n";
# insertions in purple
print FILEF"lines(sumhit\$vec,(sumhit\$V18+sumhit\$V19+sumhit\$V20+sumhit\$V21),xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"purple\",lwd=1)\n";
# G>A in blue
print FILEF"lines(sumhit\$vec,sumhit\$V2,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"blue\",lwd=2)\n";
# C>T in red
print FILEF"lines(sumhit\$vec,sumhit\$V3,xlim=c(1,$opts{l}),ylim=c(0,$opts{m}),col=\"red\",lwd=2)\n";

# do the same at 3'-end
print FILEF "nucl<-read.csv(file=\"$folder/nuclCompReverse.all.txt\",sep=\"\\t\",check.names = FALSE)\n";
print FILEF "vec<-(1:$opts{l})\nsumhit<-data.frame(cbind(vec,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))\n";
#print FILEF "idx <- nucl\$A==0\n";
#print FILEF "nucl\$A[idx] <- 10^10\n";
#print FILEF "idx <- nucl\$C==0\n";
#print FILEF "nucl\$C[idx] <- 10^10\n";
#print FILEF "idx <- nucl\$G==0\n";
#print FILEF "nucl\$G[idx] <- 10^10\n";
#print FILEF "idx <- nucl\$T==0\n";
#print FILEF "nucl\$T[idx] <- 10^10\n";
#print FILEF "idx <- nucl\$Tot==0\n";
#print FILEF "nucl\$Tot[idx] <- 10^10\n";
print FILEF "for (i in vec)	{\n";
print FILEF "sumhit[i,2]<-sum(nucl[(nucl\$Pos == i), \"G>A\"])/sum(nucl[(nucl\$Pos == i), \"G\"])\n";
print FILEF "sumhit[i,3]<-sum(nucl[(nucl\$Pos == i), \"C>T\"])/sum(nucl[(nucl\$Pos == i), \"C\"])\n";
print FILEF "sumhit[i,4]<-sum(nucl[(nucl\$Pos == i), \"A>G\"])/sum(nucl[(nucl\$Pos == i), \"A\"])\n";
print FILEF "sumhit[i,5]<-sum(nucl[(nucl\$Pos == i), \"T>C\"])/sum(nucl[(nucl\$Pos == i), \"T\"])\n";
print FILEF "sumhit[i,6]<-sum(nucl[(nucl\$Pos == i), \"A>C\"])/sum(nucl[(nucl\$Pos == i), \"A\"])\n";
print FILEF "sumhit[i,7]<-sum(nucl[(nucl\$Pos == i), \"A>T\"])/sum(nucl[(nucl\$Pos == i), \"A\"])\n";
print FILEF "sumhit[i,8]<-sum(nucl[(nucl\$Pos == i), \"C>G\"])/sum(nucl[(nucl\$Pos == i), \"C\"])\n";
print FILEF "sumhit[i,9]<-sum(nucl[(nucl\$Pos == i), \"C>A\"])/sum(nucl[(nucl\$Pos == i), \"C\"])\n";
print FILEF "sumhit[i,10]<-sum(nucl[(nucl\$Pos == i), \"T>G\"])/sum(nucl[(nucl\$Pos == i), \"T\"])\n";
print FILEF "sumhit[i,11]<-sum(nucl[(nucl\$Pos == i), \"T>A\"])/sum(nucl[(nucl\$Pos == i), \"T\"])\n";
print FILEF "sumhit[i,12]<-sum(nucl[(nucl\$Pos == i), \"G>C\"])/sum(nucl[(nucl\$Pos == i), \"G\"])\n";
print FILEF "sumhit[i,13]<-sum(nucl[(nucl\$Pos == i), \"G>T\"])/sum(nucl[(nucl\$Pos == i), \"G\"])\n";
print FILEF "sumhit[i,14]<-sum(nucl[(nucl\$Pos == i), \"A>-\"])/sum(nucl[(nucl\$Pos == i), \"A\"])\n";
print FILEF "sumhit[i,15]<-sum(nucl[(nucl\$Pos == i), \"T>-\"])/sum(nucl[(nucl\$Pos == i), \"T\"])\n";
print FILEF "sumhit[i,16]<-sum(nucl[(nucl\$Pos == i), \"C>-\"])/sum(nucl[(nucl\$Pos == i), \"C\"])\n";
print FILEF "sumhit[i,17]<-sum(nucl[(nucl\$Pos == i), \"G>-\"])/sum(nucl[(nucl\$Pos == i), \"G\"])\n";
print FILEF "sumhit[i,22]<-sum(nucl[(nucl\$Pos == i), \"S\"])/sum(nucl[(nucl\$Pos == i), \"Tot\"])\n";
print FILEF "for (j in seq(NTs)) {\n";
print FILEF "insertion[j]<-sprintf(\"->%s\", NTs[j])\n";
print FILEF "sumhit[i,17 + j] <-sum(nucl[(nucl\$Pos == i), insertion[j]]) /sum(nucl[(nucl\$Pos == i), NTs])\n";
print FILEF "}\n}\n";
# write a table for G>A frequencies
print FILEF "write.table(sumhit[,2], file=\"$folder/3pGtoA_freq.txt\", sep=\"\\t\", quote=FALSE, col.names=c(\"pos\\t3pG>A\"))\n";
print FILEF "plot(-sumhit\$vec,sumhit\$V4,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1,type=\"l\",xlab=\"\",ylab=\"\",las=2,lab=c($opts{l},11,0.1),axes=FALSE)\n";
print FILEF "axis(side=1,labels=TRUE,las=2,cex.axis=0.8)\n";
if( $opts{m} > 0.1) { print FILEF "axis(side=4,labels=TRUE,las=2,cex.axis=0.8)\n"; }
else { print FILEF "axis(side=4,labels=seq(from=0,to=$opts{m},by=0.01),at=seq(from=0,to=$opts{m},by=0.01),las=2,cex.axis=0.8)\n"; }
print FILEF "rect(-$opts{l}.5,-0.01,-0.5,$opts{m},border=\"darkgrey\")\nsegments(-$opts{l}.5,-0.01,-$opts{l}.5,$opts{m},col=\"white\",lwd=2)\n";
print FILEF "lines(-sumhit\$vec,sumhit\$V5,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF "lines(-sumhit\$vec,sumhit\$V6,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF "lines(-sumhit\$vec,sumhit\$V7,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF "lines(-sumhit\$vec,sumhit\$V8,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF "lines(-sumhit\$vec,sumhit\$V9,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF "lines(-sumhit\$vec,sumhit\$V10,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF "lines(-sumhit\$vec,sumhit\$V11,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF "lines(-sumhit\$vec,sumhit\$V12,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
print FILEF "lines(-sumhit\$vec,sumhit\$V13,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"grey\",lwd=1)\n";
# soft clipping in orange
print FILEF"lines(-sumhit\$vec,(sumhit\$V22),xlim=c(-5,-1),ylim=c(0,$opts{m}),col=\"orange\",lwd=1)\n";
# insertions
print FILEF "lines(-sumhit\$vec,(sumhit\$V14+sumhit\$V15+sumhit\$V16+sumhit\$V17),xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"green\",lwd=1)\n";
# insertions
print FILEF "lines(-sumhit\$vec,(sumhit\$V18+sumhit\$V19+sumhit\$V20+sumhit\$V21),xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"purple\",lwd=1)\n\n";
# C>T in red
print FILEF "lines(-sumhit\$vec,sumhit\$V3,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"red\",lwd=2)\n";
# G>A in blue
print FILEF "lines(-sumhit\$vec,sumhit\$V2,xlim=c(-$opts{l},-1),ylim=c(0,$opts{m}),col=\"blue\",lwd=2)\n";
print FILEF "dev.off()\n";
close (FILEF);
print "plotting using R, it may take a while...\n";
system("R CMD BATCH mapDamage.R");
print "pdf FragMisincorporation_$title.pdf generated\n";
exit 0;
}
else{
	&printUsagePlot;
}

}
######################################
#               SECONDARY SUBROUTINES DEFINITION
######################################

sub checkFileEmpty {
	
	my $file = $_[0];
	die("file $file empty! -> Exit\n") if ( -z $file ) ; # exit if file is empty
}

sub readSAMfile {  # for map

	my $file = $_[0];
	my $uniq = $_[1];
	if ( $uniq == 1 ) { print "Reading $file, remove unmapped reads and reads with alternative hits according to BWA (XT and X1 tags)\n"; }
	my $read_ref = {};
	my $line = ();
	my $input = $file;
	my $cpt=0;
	my $i=0;
	my$j=0;
	my $flag = 0;
	my $pair=0;
		 
	if($file =~ /\.bam$/) {$file = "samtools view $file | "; } 
	
	open(FILESAM,"$file") or die("Cannot open file $input\n");
	while (my $line = <FILESAM>){
 
		chomp($line);
		next if ( $line  =~ m/^@/); 	# if SAM file provided with hearder, skip it
		my @field = split("\t",$line);
		
		if ( defined($field[5]) and $field[5] =~ m/([0-9]+[MIDNSHPX=])+|\*/ ) { # check if it looks like a SAM file using the CIGAR field
			
			my $nbBestHit = 0; # intialize like unique hits if BWA wasn't used for mapping.
			my $xt = 'U';
			if($line =~ m/XT:A:(\w+)\s/) { $xt = $1;} 
			if($line =~ m/X1:i:(\d+)\s/) { $nbBestHit = $1;} 
			
			
			if ( $field[1] & 0x0004 ) {	# remove unmapped reads 
					$j++;
			}
			elsif ( $uniq == 1 and ($nbBestHit > 0 or $xt ne 'U') ) { # removed non unique hits ie., reads with suboptimal hits or declared non unique by BWA
					$j++
			}
			else{
				if(defined($read_ref->{$field[0]}{sam})) {
					# paired-end run, change name, hash is indexed by read name, must be an unique key
					$field[0].='/2';
					$pair++;
				}
				$read_ref->{$field[0]}{sam} = $line; 	
				$i++;
			}
		$cpt++;
		}
		else { die("Problem with SAM format (CIGAR field $field[5])!\n"); exit; }
		if ($cpt % 500000 == 0 && $cpt !=0) { printf"%8d reads put in RAM ($j trashed, $pair pairs found)/%8d analysed\n",$i,$cpt ;}	
		
	}
	close FILESAM;
	printf("OK, %8d reads put in RAM ($j trashed, $pair pairs found)\n",$i) ;	
	return($read_ref);
}

sub compare2SamFiles { # for map
	
	my $first = $_[0];
	my $second = $_[1];
	my $fileOUT = $_[2];
	my $cpt=0;
	my $i=0;
	my$j=0;
	
	open(FILEOUT,">$fileOUT") or die("Cannot write file $fileOUT\n");
	foreach  my $key (keys(%$first)) {
	
		if( exists( $second->{$key} ) ) {
			# read present in both files, do nothing
			$j++;
		}
		else{ 
			# if read doesn't exist in second, write whole SAM line in a separate file
			$i++;
			if(defined($first->{$key}{sam})){print FILEOUT"".$first->{$key}{sam} ."\n";}
		}
		if ($cpt % 10000 == 0 && $cpt !=0) { printf"%8d reads specific to file $first /%8d analysed\n",$i,$cpt ;}	
		$cpt++;
	}
	print "$i reads were written to $fileOUT ($j reads in both files removed)\n";
	# close file pointers
	close (FILEOUT);
}

sub writeSamFile { # for map
	
	my $read = $_[0];
	my $file =  $_[1];
	my $i = 0;
	open(FILEOUT,">$file") or die("Cannot write output file $file\n");
	foreach  my $key (keys(%$read)) {
		
		if(defined($read->{$key}{sam})){
			print FILEOUT"".$read->{$key}{sam} ."\n";
			$i++;
		}
	}
	close (FILEOUT);
	print"$i reads written to file $file\n";
}

sub initializeHash { # for map
	my $a = $_[0];
	my $b = $_[1];
	my $c = $_[2];
	my $d = $_[3];
	my $hash = $_[4];
	
	foreach my $pos ($a..$b) {
		foreach my $str ($c,$d){
			$hash->{$str}{$pos}=0;
		}
	}
	return($hash);
}	
sub checkExecutable {
     
	my $output = `samtools 2>&1 1>/dev/null`;
	die ("$output samtools is not present in your \$PATH\n") if (!defined($output) or $output  !~ /Version/) ;
}

sub checkFastaRef { # for map
	
	print "@_\n";

	my $ref =  $_[0];
	my $h = {};
	my $cpt = 0;
	my $sum = 0;
	# built reference index if not present
	if (!-e $ref.".fai") { system("samtools faidx $ref");}
	open(F,$ref.".fai") or die ("Cannot open reference index $ref.fai");
	while (my $line = <F>) {
		chomp($line);
		my @field = split("\t",$line,);
		# first field is the name of  one or several chromosomes in reference
		# second field is the length in nucleotide 
		if(defined($field[0]) and defined($field[1]) and ($field[1] =~ m/^[0-9]+$/)){  # does it look like a good fai file
			  
			$field[0] =~ s/[:\|\s]//g;  # character '|' or others must be removed because of their interpretation by the shell and samtools command
			$h->{$field[0]} = $field[1];
			$cpt++;
			$sum = $sum +  $field[1];
		} else{ die ("The format of reference index file $ref.fai is not correct, you should erase it");}
	}
	print"The reference is assumed to contain $cpt chromosome(s) and $sum nucleotides in total\n";
	return($h);
}

sub ReadFile	{ # for merge
	
	my $d = $_[0];
	my $file = $_[1];
	my $result = $_[2];
	my $chr = $_[3];
	open(FILE,"<$d/$file") or die("Cannot open file $d/$file\n");
	my $lineNB = 1;
	while ( my $line = <FILE> )	{
		
		if ( $lineNB != 1 ) {# skip the header
			if(defined($line) and defined($chr)) {
				chomp($line);
				if($chr eq "") {
					$line = $line."\n"; # just add the newline character (for hitPerChrom files)
				}
				else {
					$line = $chr."\t".$line."\n"; # add the chromosome information
				}
				push(@$result, $line);
			}
		}
		$lineNB++;
	}
	close(FILE);
	return($result);
}

sub printUsageMap{
    print "
USAGE: $0 map -i input_SAM/BAM -d directory -r reference_fasta -l read_length -a reference_around
	-i : input SAM/BAM-file
	-r : input reference fasta file, with headers identical to the ones in the input_SAM/BAM file
		
	OPTIONAL
	-l : maximal read length, in nucleotides to consider [default: 70]
	-a : size, in nucleotides, of the genomic region to be retrieved before and after reads [default: 10]
	-d : folder for writing output-files, if not already present, this folder will be automatically created [default name: results_input_SAM.filtered]
	-j : second input SAM-file without header; reads present in both SAM files will not be processed
	-u : removes all non-unique reads based on the X1 and XT tags from the BWA mapper
	-k : keep bed file and fasta file containing all the genomic regions [default: delete these files]
	-f : outputs all read aligned against the reference genome, in a fasta file format [default: not active]
	-c : complete analysis, map will be automatically followed by both merge and plot steps [default: not active]
	when -c is active, the following are available to control the plot step:
	-b : the number of reference nucleotides to consider for ploting base composition in the region located upstream and downstream of every read [default: 10]
	-y : graphical y-axis limit for nucleotide misincorporation frequencies [default: 0.30]
	-m : read length, in nucleotides, considered for plotting nucleotide misincorporations [default: 25]
	\n";

}

sub printUsageMerge{
    print "
USAGE: $0 merge -d directory

	-d :folder with output-files generated by the map command
	
	OPTIONAL
	-c : complete analysis, merge will be automatically followed by the plot step [default: not active]
	\n";
}

sub printUsagePlot{
    print "
USAGE: $0 plot -d directory

	-d : folder with output-files generated by the merge command
	
	OPTIONAL
	-l : read length, in nucleotides, considered for plotting nucleotide misincorporations [default: 25]
	-m : graphical y-axis limit for nucleotide misincorporation frequencies [default: 0.30]
	-a : the number of reference nucleotides to consider for ploting base composition in the region located upstream and downstream of every read [default: 10]
	-t : title used for both graph and filename [default: folder name without extension]\n";
}

 sub printUsage{
    print "
USAGE: $0 <command> <arguments>

Commands:
	map	main function, perform damage patterns [samtools and bedtools required] 
	merge	merge all chromosomes
	plot	plot the graphics [R environment required]\n";
	exit 0;
}
