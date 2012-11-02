#!/usr/bin/perl
use warnings;
use strict;
use YAML::Tiny;
use File::Basename;
use File::Copy;

my $config_file = YAML::Tiny->new;
# Open the config
$config_file= YAML::Tiny->read( 'fastqc_config_file.yml' );
print "config_file:$config_file\n";
#Get YAML parameters
#input files
my $dif= 0;
my $count = 0;
my @in_files;
while ($dif == 0) {
	my $in_files = $config_file->[0]->{in_files}->[$count];
	if ($in_files) {
		#die "$in_files file doesn't exists: $!" unless ( -e $in_files && -d $in_files == 0);
		$count += 1;
		push (@in_files,$in_files);
	}
	else {
		$dif = 1;
		last;
	}
}
#output directory
my $fastqc_dir = $config_file->[0]->{fastqc_dir};
if ($fastqc_dir) {	
	$fastqc_dir =~ s/\/$//;
	#die "$fastqc_dir is not a valid directory: $!" unless (-d $fastqc_dir)
}
else {
	warn "Using current directory as fastQC executable location"; 
}

#get YAML parameters
my $outdir = $config_file->[0]->{outdir};
my $threads = $config_file ->[0]->{threads};
my $format = $config_file->[0]->{format};
my $casava = $config_file->[0]->{casava};
my $noextract = $config_file->[0]->{nonextract};
my $no_group = $config_file->[0]->{no_group};
my $contaminants = $config_file->[0]->{contaminants};
my $kmers =$config_file->[0]->{kmers};
#===============================================================

#FastQC command construction (cmd concatenation)
my $fastqc_cmd = "fastqc";
$fastqc_cmd = $fastqc_dir.'/'.$fastqc_cmd if ($fastqc_dir);
$fastqc_cmd = $fastqc_cmd." --outdir ".$outdir if ($outdir);
$fastqc_cmd = $fastqc_cmd." --threads ".$threads if ($threads && $threads =~ /^[0-9]+$/);
$fastqc_cmd = $fastqc_cmd." --format ".$format if ($format);
$fastqc_cmd = $fastqc_cmd." --casava " if ($casava);
$fastqc_cmd = $fastqc_cmd." --noextract" if ($noextract);
$fastqc_cmd = $fastqc_cmd." --nogroup" if ($no_group);
$fastqc_cmd = $fastqc_cmd." --contaminants ".$contaminants if ($contaminants);
$fastqc_cmd = $fastqc_cmd." --kmers".$kmers if ($kmers);

my $input_files = join(' ',@in_files);
print $fastqc_cmd." ".$input_files."\n";
!system $fastqc_cmd." ".$input_files or die "fastqc fail : $!";


chdir $outdir;
mkdir 'Results_fastqc';
foreach (@in_files) {
	my ($filename,$directories,$suffix) = fileparse($_,qr/\.[^.]*/);
	print "$filename$suffix\_fastqc\t"."$filename$suffix\_fastqc.zip\n";
	move("$filename$suffix\_fastqc","Results_fastqc/$filename$suffix\_fastqc");
	move("$filename$suffix\_fastqc.zip","Results_fastqc/$filename$suffix\_fastqc.zip");
}



=d

 The options for the program as as follows:
    
    -h --help       Print this help file and exit
    
    -v --version    Print the version of the program and exit
    
    -o --outdir     Create all output files in the specified output directory.
                    Please note that this directory must exist as the program
                    will not create it.  If this option is not set then the 
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.
                    
    --casava        Files come from raw casava output. Files in the same sample
                    group (differing only by the group number) will be analysed
                    as a set rather than individually. Sequences with the filter
                    flag set in the header will be excluded from the analysis.
                    Files must have the same names given to them by casava
                    (including being gzipped and ending with .gz) otherwise they
                    won't be grouped together correctly.
                   
    --extract       If set then the zipped output file will be uncomressed in
                    the same directory after it has been created.  By default
                    this option will be set if fastqc is run in non-interactive
                    mode.
                   
    --noextract     Do not uncompress the output file after creating it.  You
                    should set this option if you do not wish to uncompress
                    the output when running in non-interactive mode.
                    
    --nogroup       Disable grouping of bases for reads >50bp. All reports will
                    show data for every base in the read.  WARNING: Using this
                    option will cause fastqc to crash and burn if you use it on
                    really long reads, and your plots may end up a ridiculous size.
                    You have been warned!
                    
    -f --format     Bypasses the normal sequence file format detection and
                    forces the program to use the specified format.  Valid
                    formats are bam,sam,bam_mapped,sam_mapped and fastq
                    
    -t --threads    Specifies the number of files which can be processed
                    simultaneously.  Each thread will be allocated 250MB of
                    memory so you shouldn't run more threads than your
                    available memory will cope with, and not more than
                    6 threads on a 32 bit machine
                  
    -c              Specifies a non-default file which contains the list of
    --contaminants  contaminants to screen overrepresented sequences against.
                    The file must contain sets of named contaminants in the
                    form name[tab]sequence.  Lines prefixed with a hash will
                    be ignored.
                    
   -q --quiet       Supress all progress messages on stdout and only report errors.
=cut
