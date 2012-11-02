#!/usr/bin/perl
#bed2MySQL.pl;
use warnings;
use strict;
use File::Basename;
use YAML::Tiny;
use File::Copy;
use DBI;
use DBD::mysql;

my @bed_files = @ARGV;

#Database login info
my $dsn = "DBI:mysql:GFF3_store:localhost";
my $user = "du_GFF3_store";
my $passwd = "santafek";

my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #connection to the DB
 #All peaks info is going to be stored in a MySQL table
 #Auto-increment ID must use because there is no simple way to define the uniqueness of a peak

#========RESET AUTOINCREMENT=================
my $check_if_empty_sql = "SELECT COUNT(*) FROM BED_Peaks_Info";
my $value = check_existence($check_if_empty_sql);

unless(check_existence($check_if_empty_sql)) { #when table is empty -> reset AUTO_INCREMENT
	print "VALUE $value\n";
	my $reset_autoincrement_sql = "ALTER TABLE BED_Peaks_Info AUTO_INCREMENT=1";
	sth($reset_autoincrement_sql);
}
#=============================================

#insert bed_info into seq_id table
 #Fill BED_Peaks_Info Table
my $insert_bed_sql = "INSERT INTO BED_Peaks_Info(file_path,filename,suffix,chr_seqID,start,end,program,score) VALUES" ;
my $question_marks = '?'x8;
$question_marks = join(',',split('',$question_marks));
$insert_bed_sql .= "($question_marks)";
print "$insert_bed_sql\n";
my $insert_bed_sth = $dbh->prepare($insert_bed_sql);
#============================


#Take multiple bed files
foreach my $bed_file (@bed_files) {
	my ($filename,$directories,$suffix) = fileparse($bed_file,qr/\.[^.]*/);
	#Need to quote file, directory and suffix
	my $directories_quote = $dbh->quote($directories);
	my $filename_quote = $dbh->quote($filename);
	my $suffix_quote = $dbh->quote($suffix);
	#=======================================
	#check file have bed extension.
	die "File extension must be bed\n" unless ($suffix =~ /bed/i);
	open IN,'<',$bed_file or die "Cannot open $bed_file";

	while (<IN>) {
		chomp;
		next if $_ !~ /\t/;
		my ($chr,$start,$end,$program,$score) = split('\t',$_);
		my $chr_quote = $dbh->quote($chr);
		my $program_quote = $dbh->quote($program);

		#check unique by $chr $start $end $program combination
		#Info that must be saved in the DB $filename, $directories, $suffix (check if is bed), $chr, $start, $end, $program and $score
		
		my $check_geneid_new_sql = "SELECT COUNT(*) FROM SEQID_LANDMARK WHERE Name LIKE $chr_quote";#Check if id is already present in the table
		if (check_existence($check_geneid_new_sql)) {
			
			my $check_peak_exist_sql = "SELECT COUNT(*) FROM BED_Peaks_Info WHERE file_path = $directories_quote AND filename LIKE $filename_quote AND suffix LIKE $suffix_quote AND chr_seqID LIKE $chr_quote AND start = '$start' AND end = '$end'";
			 unless (check_existence($check_peak_exist_sql)) {
			 	#INSERT INTO BED_Peaks_Info(file_path,filename,suffix,chr_seqID,start,end,program,score) VALUES
			 	#print "($directories,$filename,$suffix,$chr,$start,$end,$program,$score)\n";
			 	$insert_bed_sth->execute($directories,$filename,$suffix,$chr,$start,$end,$program,$score);
			 }
			 else {
			 	warn "This peak from $filename is already in the database\n";
			 }
		}
		else  {
			warn "Looks like the Chromosomes/Scaffolds names from $filename not match the Chromosomes/Scaffolds names provided in the annotation file\n";
		}
	}
	close IN;
}


sub sth { #prepare, execute a sql
	my $insert_sth = $dbh->prepare($_[0]);#prepare , look for basic sql errors 
	$insert_sth->execute(); # action takes part here!!
	return 1;
}

sub check_existence { #send COUNT(*) to check if register is new
	my $check_if_new = $dbh->selectrow_array($_[0]); #if it's already in the table value = 1 
	#print "check_if_new = $check_if_new\n";
	if ($check_if_new) {
		return 1;
	}
	else {
		return 0;
	}
}
