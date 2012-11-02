#!/usr/bin/perl
#Mappability.pl
use warnings;
use strict;
use File::Basename;
use File::Copy;
use Bio::SeqIO;

#===========
my $mappability_yml_file = shift @ARGV;
my $config_file = YAML::Tiny->new;
$config_file= YAML::Tiny->read($mappability_yml_file) or die "Cannot read config file: $!\n";


my $fasta_file_path = $config_file->[0]->{'genome_fasta'};
my $avg_read_length = $config_file->[0]->{'avg_read_length'};
my $frag_len = $config_file->[0]->{'frag_len:'};
my $bin_size = $config_file->[0]->{'bin_sizes'};

#===========

#check fragment length, bin size and read length are numbers
unless (($frag_len && $bin_size && $avg_read_length) ~~ /^[0-9]+$/) {
	die "average_read_length, fragment length and bin size must be integers";
}

#===========================

#get Mappability code and extract it in Mappability directory
!system 'wget','http://archive.gersteinlab.org/proj/PeakSeq/Mappability_Map/Code/Mappability_Map.tar.gz' or die "Fail downloading Mappability_Map.tar.gz : $!" ;
mkdir 'Mappability', 0755 or warn "Cannot make Mappability directory: $!";
!system 'tar','xvzf','Mappability_Map.tar.gz', '-C', 'Mappability' or die "tar unexpected fail : $!";
unlink 'Mappability_Map.tar.gz' or warn "deleting Mappability_Map.tar.gz fail : $!";
#===========================

#parse filename, dir path and suffix from the input genome
my ($filename,$directories,$suffix) = fileparse($fasta_file_path,qr/\.[^.]*/);
print "$filename, $directories, $suffix\n";

#===========================
my @Seq_IDs; #array where the different identifiers of the mutifasta are going to be saved
my %ID_length; #where lengths of each seqID are going to be saved
my $seqio = Bio::SeqIO->new(-file => "$fasta_file_path", '-format' => 'Fasta'); #read the fasta file
 
#foreach sequence in the file 
while ( my $seq = $seqio -> next_seq) {
	my $id = lc($seq->id); #get id (ex. Chr) in lowercase
	push (@Seq_IDs,$id); #save id
	my $length = $seq->length; #get sequence length
	$ID_length {$id} = $length;#save length indexed by seqid
	#create a fasta foreach sequence
	my $out = Bio::SeqIO->new(-file => ">$id.fa", '-format' => 'Fasta');
	$out->write_seq($seq);
}

#create fa_chr directory and move into it the .fa files
mkdir 'fa_chr', 0755 or warn "Cannot make fa_chr directory: $!";
my @fa_files = glob '*.fa';
foreach (@fa_files) {
	if ($_ ne $fasta_file_path) {
		rename ($_,"fa_chr/".$_) or warn "Cannot move file: $!";
	}
}

#go to /Mappability directory / make and include the folder in the system PATH
chdir $ENV{'PWD'}.'/Mappability' or die "Cannot chdir to Mappability :$!";
!system 'make' or warn "Apperently it's sth wrong with the make output: $!";
$ENV{'PATH'} = $ENV{'PWD'}."/Mappability".":$ENV{'PATH'}";

#go to /fa_chr and execute compyle.py to generate mappability files
chdir $ENV{'PWD'}.'/fa_chr' or die "Cannot chdir to chr_fa :$!";
!system 'python', "$ENV{PWD}\/Mappability/compile.py", "$avg_read_length" or die "Error with compile.py: $!";

#create the directory mappability_files and store the mappabity relevat files
mkdir 'mappability_files', 0755 or warn "Cannot make mappability_files directory: $!";
my @mappability_files = glob '*b.out';
foreach (@mappability_files)  {
	rename ($_,"mappability_files/".$_) or warn "Cannot move mappability files:$!";
}

#get mosaics_preprocessing_scripts.tar.gz and extracted into /fa_chr
!system 'wget','http://www.stat.wisc.edu/~keles/Software/mosaics/mosaics_preprocessing_scripts.tar.gz' or die "Fail downloading mosaics_preprocessing_scripts.tar.gz : $!" ;
#mkdir 'mosaics_preprocessing_scripts', 0755 or warn "Cannot make mosaics_preprocessing_scripts': $!";
!system 'tar','xvzf','mosaics_preprocessing_scripts.tar.gz', '-C', 'mappability_files' or die "tar unexpected fail : $!";
unlink 'mosaics_preprocessing_scripts.tar.gz' or warn "deleting mosaics_preprocessing_scripts.tar.gz fail : $!";

#go to mappability_files directory
chdir 'mappability_files' or die "Cannot chdir to mappability_files: $!";

#process mappability files using mosaics scripts for preprocessing and at the end obtain a bin file
foreach (@Seq_IDs) {
	my ($chr_id)= /[^0-9]+([0-9]+$)/;
	!system ("python cal_binary_map_score.py $chr_id 1 $ID_length{$_} > $_\_map_binary.txt") or die "cal_binary_map_score fails: $!"; 
	!system ("perl process_score_java.pl $_\_map_binary.txt $chr_id\_map $avg_read_length $frag_len $bin_size") or die "Preprocess mappability binary files to bin-level files fail : $!";
}

#return to /fa_chr where chr*.fa are located and run mosaics preprocessing script cal_binary_GC_N_score.pl

chdir '../' or die "Cannot chdir :$!";
foreach (@fa_files) {
	copy ($_,"mappability_files/".$_) or die "Cannot copy chr*.fa files : $!";
}
chdir 'mappability_files' or die "Cannot chdir :$!";

foreach (@Seq_IDs) {
	my ($chr_id)= /[^0-9]+([0-9]+$)/;
	!system "perl cal_binary_GC_N_score.pl $_.fa $chr_id 1" or die "cal_binary_GC_N_score.pl fails :$!";
	!system "perl process_score.pl $_\_GC_binary.txt $chr_id\_GC $frag_len $bin_size" or die "process_score.pl fails for GC: $!";
	!system "perl process_score.pl $_\_N_binary.txt $chr_id\_N $frag_len $bin_size" or die "process_score.pl for N fails : $!";
}

foreach (glob '*fragL*_bin*[0-9].txt') {
	my ($fn,$dir,$suf) = fileparse($_,qr/\.[^.]*/);
	open IN,'<',"$_" or die "Cannot read bin files to add Chr column: $!";
	open OUT,'>',"$fn\_colChr$suf" or die "Cannot write file adding Chr column: $!";
	my ($chr_id) = /(chr[0-9]+)_/i;
	$chr_id =~ s/chr/Chr/;
	print OUT $chr_id."\t".$_ while (<IN>);
	close IN; close OUT;
}

my @types = ('map', 'GC', 'N');
join_class_files($_) foreach @types;


mkdir $ENV{'PWD'}.'/MOSAICS_REQUIRE_FILES', 0755 or warn "Cannot create directory :$!";
foreach (glob 'whole_*') {
	rename ($_, $ENV{'PWD'}.'/MOSAICS_REQUIRE_FILES/'.$_) or die "Cannot move final files:$!";
}

print "MOSAICS Bin generator, Done!!\n";


sub join_class_files {

	my $type = shift;

	my @files = glob "*_$type\_*colChr*";
	open OUT, '>', "whole_$type\_fragL$frag_len\_bin$bin_size.txt" or die "Cannot write file : $!";
	foreach (@files) {
		open IN, '<',"$_" or die  "Cannot read colChr file : $!";
		while (<IN>) {
			print OUT $_;
		}
		close IN;
	}
	close OUT;
}


