#!/usr/bin/perl
#bed_analysis.pl
use warnings;
use strict;

	
use List::Util qw(min max sum);
use POSIX;
use File::Copy;
use File::Basename;

use DBI;
use DBD::mysql;
use YAML::Tiny;

#+++++++++Description+++++++++++++++
#Takes bed files loaded to DB and estimate the summit location for-each of them.
#To accomplish summit estimation .bam of the aligned reads is needed (provided in YAML config file)
#Filter size (min, max) limits must be provided, the criteria applies to all the bed files.
#Score filter, a cut-off for each peak-caller result could be provided

#++++++++Script Usage+++++++++++++++
#bed_analysis.pl <config_file>
#In the pipeline this script is called by : 


#======Config-file===========
my $bed_manipulation_config_file = shift @ARGV; #config file 
die "bed_manipulation_config_file ($bed_manipulation_config_file) doesn't exists: $!" unless (-e $bed_manipulation_config_file); #check if config_file exists
die "bed_manipulation_config_file ($bed_manipulation_config_file) doesn't end in .yml" unless ($bed_manipulation_config_file =~ /\.yml$/); #check file suffix .yml mandatory
#================================

#=====Open config_file===========
my $config_file = YAML::Tiny->new;
$config_file= YAML::Tiny->read("$bed_manipulation_config_file") or die "Cannot read YAML config file $bed_manipulation_config_file: $!";
#=================================

#=======Files input=================
my $bam_file = $config_file->[0]->{'bam_file_path'}; #From YAML config file
my @bed_files; #FROM MySQL
my $shuffle_tries = $config_file->[0]->{'shuffle_tries'}; #shuffle iterations	
#===================================


#=========Database login info================
my $dsn = "DBI:mysql:GFF3_store:localhost";
my $user = "du_GFF3_store";
my $passwd = "santafek";

my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #conection to the DB
#============================================

#=========Extract bed files names from DB======
my $select_files_sql = "SELECT file_path,filename,suffix FROM BED_Peaks_Info GROUP BY filename";
my $sth_select_files = sth($select_files_sql);
while (my $ref = $sth_select_files->fetchrow_arrayref()) {
	my ($file_path,$filename,$suffix) = ($ref->[0],$ref->[1],$ref->[2]);
	my $file = $file_path.$filename.$suffix;
	push(@bed_files, $file);
}
#=============================================

#=======Extract chromosomes/scaffolds names ===============
my $select_chr_sql  = "SELECT DISTINCT(chr_seqID) FROM BED_Peaks_Info";
my $sth_select_chr = sth($select_chr_sql);
my @peaks_chrs;
while(my $ref = $sth_select_chr->fetchrow_arrayref()) {
	my $chr = $ref->[0];
	push(@peaks_chrs,$chr);
}
#=======================================================

#==========================UPDATE DB with the compute summit location=============================
warn "Starting Summit estimation\n";
foreach my $file (@bed_files) {
	warn "Performing Summit estimation for $file...";
	my ($fn,$dir,$suf) = fileparse($file,qr/\.[^.]*/); 
	my $fn_quo = $dbh->quote($fn);
	my $dir_quo = $dbh->quote($dir);
	my $suf_quo = $dbh->quote($suf);
	my $check_summits_not_present_sql = "SELECT COUNT(*) FROM BED_Peaks_Info WHERE (file_path LIKE $dir_quo AND filename LIKE $fn_quo AND suffix LIKE $suf_quo) AND Calculate_Summit IS NULL";
	if (check_existence($check_summits_not_present_sql) != 0) {
		open IN , '<', $file  or die "Cannot open $file : $!"; # Open bed file
		#open OUT, '>', "$fn\_summits$suf" or die "Cannot write in $fn\_summits$suf: $!";#Output files
		my (@chrs,@starts,@ends,@sources,@scores);
		while (<IN>) {
			next unless $_ ~~ /\t/; #Important if bed files had headers e.g MOSAiCS output
			my ($chr,$start,$end,$source,$score) = /(^\S+)\t([0-9]+)\t([0-9]+)\t(\S+)\t(\S+)/; #line fields better with split
			push (@chrs,$chr); push(@starts,$start); push(@ends,$end); push(@sources,$source); push(@scores,$score); # save all in independent arrays, order conserved
		}
		close IN;
		
		my $count = 0;
		#This procedure take place foreach peak
		

		foreach my $start (@starts) {
			#print "START:$_\n";
			#Run mpileup to get the depth at each base of the peak region.
			#Save profiles in pileup.tmp , warnings in pileup2.tmp
			!system "samtools mpileup -r $chrs[$count]:$starts[$count]-$ends[$count] $bam_file > $fn\_pileup.tmp 2> $fn\_pileup_2.tmp" or die "Cannot execute samtools mpileup : $!";
			open IN_1, '<', "$fn\_pileup.tmp"  or die "Cannot open $fn\_pileup.tmp"; #take a look to mpileup file (file shape: /^[position]\t[count]$)

			#Extract the position (base) that have the higher coverage (tag count)
			my (@pos,@depth);
			my %pos_depth;
			while (<IN_1>) {
				my ($pos,$depth) = /^\S+\t([0-9]+)\t\S+\t([0-9]+)\t/ or die "Not match found :$!";
				push(@pos,$pos); push(@depth,$depth);
				$pos_depth{$pos} = $depth;
			}
			close IN_1;
			#Get the maximum from @summit_depth
			my $summit_depth = max (@depth);
			#Again open pileup.tmp file
			#TODO:Waste of time the same could be accomplish just calling %pos_depth
			#I remember, I need to look the file again becouse multiple lines could have the same coverage
			open IN_1, '<', "$fn\_pileup.tmp"  or die "Cannot open $fn\_pileup.tmp";
			my @summit_pos;
			while (<IN_1>) { 
				my ($pos,$depth) = /^\S+\t([0-9]+)\t\S+\t([0-9]+)\t/ or die "Not match found :$!";
				if ($depth == $summit_depth) {
					push (@summit_pos,$pos);#push the max set
				}
			}

			#====================================
			#Average of the max set , floor just aproximates to the smaller number if the original number is exactly in the middle x.5
			my $summit_pos = floor(sum (@summit_pos)/ @summit_pos);
			#print OUT "$chrs[$count]\t$summit_pos\t".($summit_pos+1)."\t$sources[$count]\t$scores[$count]\n"; #Print line following BED format
			my $chr_quo = $dbh->quote($chrs[$count]);

			my $update_summit_field_sql = "UPDATE BED_Peaks_Info SET Calculate_Summit = '$summit_pos' WHERE file_path LIKE $dir_quo AND filename LIKE $fn_quo AND suffix LIKE $suf_quo AND chr_seqID LIKE $chr_quo AND start = '$start' AND end = '$ends[$count]' AND Calculate_Summit IS NULL";
			my $update_summit_field_sth = $dbh->prepare($update_summit_field_sql);
			$update_summit_field_sth->execute() or warn "Nothing to update\n";
			#print "$update_summit_field_sql\n";
			close IN_1;
			#Remove .tmp files
			unlink "$fn\_pileup.tmp" or die "Cannot remove $fn\_pileup.tmp : $!";
			unlink "$fn\_pileup_2.tmp" or die "Cannot remove $fn\_pileup_2.tmp : $!";
			$count += 1;
		}
		#close OUT;
	}
	warn "Done\n";
}
warn "Summit estimation finished\n\n";
	
#=================Filter part=========================
warn "Filtering Files...\n";
my $size_filter = $config_file->[0]->{'filter_size_limit'};
my ($min_size,$max_size) = split(',',$size_filter) if $size_filter;
my @bed_files_FilteredPeaks;
my @bed_files_FilteredSummits;

foreach my $file (@bed_files) {
	warn "Filtering $file...";
	my ($fn,$dir,$suf) = fileparse($file,qr/\.[^.]*/);
	my $fn_quo = $dbh->quote($fn);
	my $dir_quo = $dbh->quote($dir);
	my $suf_quo = $dbh->quote($suf);
	my $select_filter_peaks_sql;
	open OUT_peak,'>',$fn."_".$min_size."_".$max_size."_peaks".$suf or die "Cannot write: $!";
	open OUT_summit,'>',$fn."_".$min_size."_".$max_size."_summits".$suf or die "Cannot write: $!";
	
	#Save output files names for later processes
	push(@bed_files_FilteredPeaks,$fn."_".$min_size."_".$max_size."_peaks".$suf);
	push(@bed_files_FilteredSummits,$fn."_".$min_size."_".$max_size."_summits".$suf);

	
	foreach my $program (keys %{$config_file->[0]->{'filter_score_more_than'}}) {
		
		
		my $filter_score = $config_file->[0]->{'filter_score_more_than'}->{$program};
		
		if ($filter_score and $fn eq $program) {
			$select_filter_peaks_sql = "SELECT chr_seqID, start, end, program, score, Calculate_Summit FROM BED_Peaks_Info WHERE (file_path LIKE $dir_quo AND filename LIKE $fn_quo AND suffix LIKE $suf_quo) AND score >= $filter_score ";
		}
		elsif (!($filter_score) and $fn eq $program) {
			$select_filter_peaks_sql = "SELECT chr_seqID, start, end, program, score, Calculate_Summit FROM BED_Peaks_Info WHERE (file_path LIKE $dir_quo AND filename LIKE $fn_quo AND suffix LIKE $suf_quo) ";
		}
	}
	
	$select_filter_peaks_sql .= "AND (end-start) >= '$min_size'" if $min_size;
	$select_filter_peaks_sql .= "AND (end-start)<= '$max_size'" if $max_size;

	my $sth_select_filter_peaks = sth($select_filter_peaks_sql);
	while (my $ref = $sth_select_filter_peaks->fetchrow_arrayref()) {
		my ($chr,$start,$end,$program,$score,$summit) = ($ref->[0],$ref->[1],$ref->[2],$ref->[3],$ref->[4],$ref->[5]);
		print OUT_peak "$chr\t$start\t$end\t$program\t$score\n";
		print OUT_summit "$chr\t$summit\t".($summit + 1)."\t$program\t$score\n";
	}
	close OUT_peak;
	close OUT_summit;
	warn "Done\n";
}
warn "BED_files filtered generated\n";
#=========================================================================================


#=============Fetch Annotation Info from the DB (specific features files)==================
warn "Starting Fetching Annotation info from DB...";
my $GFF3_fetch = 'GFF3_all_features'; #basemane of output file 
!system "perl fetch_GFF3_store_anntotation.pl $GFF3_fetch" or die "Cannot fetch the features info from the DB\n";
!system "sortBed -i $GFF3_fetch.bed > $GFF3_fetch\_sorted.bed" or die "Cannot execute sortBed: $!"; #Sort generated bed file
unlink "$GFF3_fetch.bed"; #Remove unsorted file 
open GFF3_ALL,'<',$GFF3_fetch."_sorted.bed" or die "Cannot read $GFF3_fetch\_sorted.bed: $!"; #open GFF3_ALL_sorted.bed
#Filehandles to split GFF3_ALL by features
open OUT_TSS,'>',"$GFF3_fetch\_TSS.bed" or die "Cannot write : $!";
open OUT_SC,'>',"$GFF3_fetch\_SC.bed" or die "Cannot write : $!";
open OUT_nonOverlap,'>',"$GFF3_fetch\_nonOverlapping.bed" or die "Cannot write : $!";
open OUT_EI,'>',"$GFF3_fetch\_exon_intron.bed" or die "Cannot write : $!";


while (<GFF3_ALL>) {
	chomp;
	my ($chr,$start,$end,$feature) = $_ =~ /^(\w+)\t(\d+)\t(\d+)\t(.+);.+;/;
	if ($chr ~~ @peaks_chrs) { #Check if Chr match
		if ($feature eq 'TSS') { #split the TSS info
			print OUT_TSS $_."\n";
		}
		elsif ($feature eq 'start_codon') { #split the SC info
			print OUT_SC $_."\n";
			
		}
		elsif ($feature !~ /exon|intron/) { #split non_overlaping features (UTR's, CDS, non-CDS)
			print OUT_nonOverlap $_."\n";
		}
		elsif ($feature ~~ /exon|intron/) { #split exon\intron info
			print OUT_EI $_."\n";
		}
	}

}
close OUT_TSS; close OUT_SC; close OUT_nonOverlap; close OUT_EI;
 #Annotation files that are going to be used for distance analysis
my @distance_analysis_files = ("$GFF3_fetch\_TSS.bed","$GFF3_fetch\_SC.bed");
#Annotation files that are going to be used to overlapping analysis
my @overlapping_analysis_files = ("$GFF3_fetch\_nonOverlapping.bed","$GFF3_fetch\_exon_intron.bed"); 
warn "Done\n";
#===================================================================

#========================Closest gene Analysis=======================
warn "Starting Closest Analysis\n";
my @closest_analysis_files; #closest .2bed files names saved in this array
foreach my $bed_summit_file (@bed_files_FilteredSummits) {
	warn "\tPerforming Closest Analysis for $bed_summit_file...";
	foreach my $annotation_distance_files (@distance_analysis_files) {
		my $feature;
		if ($annotation_distance_files =~ /_(TSS)\.bed/) {
			$feature = $1;
		}
		elsif ($annotation_distance_files =~ /_(SC)\.bed/) {
			$feature = $1;
		}
		my $closest_file = Closest_job($feature,$annotation_distance_files,$bed_summit_file);
		push(@closest_analysis_files,$closest_file);
	}
	warn "Done\n";
}
warn "Closest Analysis performed\n\n";


#======================Save closest information into the DB and graphics compute info===========

my @Graph_original_files; #tab files for distance analysis
foreach my $closest_file (@closest_analysis_files) {
	#compute distance graphics files
	my $Graph_Info_file = Graph_Info($closest_file,'original','-1'); #-1 to differenciate from the shuffle files
	push(@Graph_original_files,$Graph_Info_file);

	#DB store part
	open IN_CLOSEST,'<',"$closest_file" or die "Cannot open $closest_file:$!";
	my @closest_file_parts = split ('_',$closest_file);
	my $part_no =0;
	my @original_name_parts;

	warn "Saving closest information of $closest_file...";
	foreach my $closest_file_part (@closest_file_parts) {
		#extract the basename for the corresponding closest .2bed file, needed to serch into BED_files DB table
		last if ($closest_file_part ~~ $min_size && $closest_file_parts[$part_no+1] ~~ $max_size && $closest_file_parts[$part_no+2] ~~ 'summits' && $max_size && $closest_file_parts[$part_no+3] ~~ 'closest');
		push(@original_name_parts,$closest_file_part);
		$part_no += 1;
	}
	my $original_name = join('_',@original_name_parts);
	my $original_name_quo = $dbh->quote($original_name);
	#prepare insert sql (general, with placeholders)
	my $insert_Closest_Analysis_sql = "INSERT INTO Closest_Analysis (Peak_ID,RNA_ID,Distance,Feature) VALUES(?,?,?,?)";
	my $insert_Closest_Analysis_sth = $dbh->prepare($insert_Closest_Analysis_sql);

	while (<IN_CLOSEST>) {
		my @fields = split('\t',$_);
		my ($chr,$summit,$annotation_info,$distance) = ($fields[0],$fields[1],$fields[8],$fields[-1]);#extract info from the closest file
		my $chr_quo = $dbh->quote($chr);
		my ($feature,$RNA_ID) = split (';',$annotation_info);
		my $select_Peak_ID_sql = "SELECT Peak_ID FROM BED_Peaks_Info WHERE filename LIKE $original_name_quo AND chr_seqID LIKE $chr_quo AND Calculate_Summit = '$summit'";
		#print "$select_Peak_ID_sql\n";
		my $Peak_ID = recover_info($select_Peak_ID_sql);
		#print "$Peak_ID\n";
		my $check_RNA_ID_sql = "SELECT COUNT(*) FROM Parent_RNA WHERE  RNA_ID = '$RNA_ID'";
		my $check_closest_entry_sql  = "SELECT COUNT(*) FROM Closest_Analysis WHERE Peak_ID ='$Peak_ID' AND RNA_ID  LIKE '$RNA_ID' AND Distance = '$distance' AND Feature LIKE '$feature'";
		if (check_existence($check_RNA_ID_sql)) {
			 unless (check_existence($check_closest_entry_sql)) {
			 	$insert_Closest_Analysis_sth->execute($Peak_ID,$RNA_ID,$distance,$feature)
			 }
		}
		else {
			warn "\tRNA_ID not present in the DB: $!";
		}
	}
	close IN_CLOSEST;

	warn "Done\n";
}
warn "Closest Information saved into DB table Closest_Analysis\n\n";

#======================Overlap Analysis========================
warn "Starting Overlap Analysis\n";
my @intersect_analysis_files;
#print "bed_files_FilteredSummits: @bed_files_FilteredSummits\n";
foreach my $bed_summit_file (@bed_files_FilteredSummits) {
	warn "\tPerforming Overlap Analysis for $bed_summit_file ...";
	foreach my $annotation_overlap_file (@overlapping_analysis_files) {
		my $feature;
		if ($annotation_overlap_file =~ /_(nonOverlapping)\.bed/) {
			$feature = $1;
		}
		elsif ($annotation_overlap_file =~ /_(exon_intron)\.bed/) {
			$feature = $1;
		}
		#print "$bed_summit_file\n";
		my $Intersect_file = Overlap_job($feature,$annotation_overlap_file,$bed_summit_file);
		push(@intersect_analysis_files,$Intersect_file);
	}
	warn "Done\n";
}
warn "Overlap Analysis performed\n\n";

#==================Save Overlapping info into DB ====================
foreach my $overlap_file (@intersect_analysis_files) {
	
	warn "Saving overlap information of $overlap_file...";
	open IN_OVERLAP,'<',"$overlap_file" or die "Cannot open $overlap_file:$!";
	my @overlap_file_parts = split ('_',$overlap_file);
	my $part_no =0;
	my @original_name_parts;
	
	foreach my $overlap_file_part (@overlap_file_parts) {
		last if ($overlap_file_part ~~ $min_size && $overlap_file_parts[$part_no+1] ~~ $max_size && $overlap_file_parts[$part_no+2] ~~ 'summits' && $max_size && $overlap_file_parts[$part_no+3] ~~ 'vs');
		push(@original_name_parts,$overlap_file_part);
		$part_no += 1;
	}
	my $original_name = join('_',@original_name_parts);
	my $original_name_quo = $dbh->quote($original_name);

	my $insert_Overlap_Analysis_sql = "INSERT INTO Overlap_Analysis (Peak_ID,RNA_ID,Feature,Feature_No) VALUES(?,?,?,?)";
	my $insert_Overlap_Analysis_sth = $dbh->prepare($insert_Overlap_Analysis_sql);
	while (<IN_OVERLAP>) {
		my @fields = split('\t',$_);
		my ($chr,$summit,$annotation_info) = ($fields[0],$fields[1],$fields[8]);
		my $chr_quo = $dbh->quote($chr);
		my ($feature,$RNA_ID) = split(';',$annotation_info);

		my ($feature_type,$feature_index);
		if ($feature =~ /(five|three|UTR)/i) {
			($feature_type,$feature_index) = ($feature,'NULL');
		}
		else {
			my @feature_parts = split('_',$feature);
			$feature_index = pop(@feature_parts);
			warn "feature index must be an integer: $!" if ($feature_index !~ /[0-9]+/);
			$feature_type = join ('_',@feature_parts);
		}

		my $feature_type_quo = $dbh->quote($feature_type);
		my $select_Peak_ID_sql = "SELECT Peak_ID FROM BED_Peaks_Info WHERE filename LIKE $original_name_quo AND chr_seqID LIKE $chr_quo AND Calculate_Summit = '$summit'";
		
		my $Peak_ID = recover_info($select_Peak_ID_sql);
		my $check_RNA_ID_sql = "SELECT COUNT(*) FROM Parent_RNA WHERE  RNA_ID = '$RNA_ID'";
		my $check_overlap_entry = "SELECT COUNT(*) FROM Overlap_Analysis WHERE Peak_ID ='$Peak_ID' AND RNA_ID  LIKE '$RNA_ID'  AND Feature LIKE $feature_type_quo AND Feature_No = '$feature_index'";
		if (check_existence($check_RNA_ID_sql)) {
			unless(check_existence($check_overlap_entry)) {
				$insert_Overlap_Analysis_sth->execute($Peak_ID,$RNA_ID,$feature_type,$feature_index)
			}
		}
		else {
			warn "RNA_ID not present in the DB: $!";
		}

	}
	close IN_OVERLAP;
	warn "Done\n";
}
warn "Closest Information saved into DB table Closest_Analysis\n\n";

#=======================================SHUFFLE DATA=================================
#========Generate genome.tab file===========
warn "Generation genome.tab file with chr's lengths...";
my %chrID_length;
my  $genome_info_file = 'genome_lengths.tab';
open TAB_OUT,'>',$genome_info_file or die "Cannot write on $genome_info_file: $!";
my $select_chr_info_sql = "SELECT Name,Length FROM SEQID_LANDMARK";
my $sth_chr_info = sth($select_chr_info_sql);
while (my $ref = $sth_chr_info->fetchrow_arrayref()) {
	my ($chr,$length) = ($ref->[0],$ref->[1]);
	
	print TAB_OUT "$chr\t$length\n";
}
close TAB_OUT;
warn "Done\n";
#===========================================


#=========Generate- Shuffle bed datasets=======
warn "Starting Genereate-Shuffle bed data-sets...\n";
my @output_shuffle;
foreach my $bed_file_summits (@bed_files_FilteredSummits) {
	warn "Shuffling $bed_file_summits...";
	my $index = 0;
	while ($index < $shuffle_tries) {
		my $shuffle_file = ShuffleBed($bed_file_summits,$genome_info_file,$index);
		$index++;
		push(@output_shuffle,$shuffle_file);
	}
	warn "Done\n";
}
my $starting_dir = "$ENV{'PWD'}";
unless (-e 'ShuffleFiles') { 
	mkdir 'ShuffleFiles' or die "Cannot create ShuffleFiles directory: $!";
}
chdir ("ShuffleFiles") or die "Cannot change to ShuffleFiles directory: $!";
foreach (glob '*') {
	unlink $_ or die "Cannot remove files in ShuffleFiles: $!";
}
chdir $starting_dir or die "Cannot change to dir $starting_dir: $!";

print $_."\n" foreach (@output_shuffle);
move("$_","ShuffleFiles") foreach (@output_shuffle);
chdir ("ShuffleFiles") or die "Cannot change to ShuffleFiles directory\n";
warn "Genereate-Shuffle bed data-sets finished\n";
#===============================================


#========================Closest SHUFFLE gene Analysis=======================
warn "Starting Closest SHUFFLE Analysis\n";
my @closest_analysis_shuffle_files;
foreach my $bed_summit_file (@output_shuffle) {
	warn "\tPerforming Closest SHUFFLE Analysis for $bed_summit_file...";
	foreach my $annotation_distance_files (@distance_analysis_files) {
		my $feature;
		if ($annotation_distance_files =~ /_(TSS)\.bed/) {
			$feature = $1;
		}
		elsif ($annotation_distance_files =~ /_(SC)\.bed/) {
			$feature = $1;
		}
		my $closest_file = Closest_job($feature,$starting_dir."/".$annotation_distance_files,$bed_summit_file);
		push(@closest_analysis_shuffle_files,$closest_file);
	}
	warn "Done\n";
}
warn "Closest SHUFFLE Analysis performed\n\n";
#=======================================================================



#======================Overlap Analysis========================
warn "Starting Overlap SHUFFLE Analysis\n";
my @intersect_analysis_shuffle_files;
foreach my $bed_summit_file (@output_shuffle) {
	warn "\tPerforming Overlap SHUFFLE Analysis for $bed_summit_file...";
	foreach my $annotation_overlap_file (@overlapping_analysis_files) {
		
		my $feature;
		if ($annotation_overlap_file =~ /_(nonOverlapping)\.bed/) {
			$feature = $1;
		}
		elsif ($annotation_overlap_file =~ /_(exon_intron)\.bed/) {
			$feature = $1;
		}
		my $Intersect_file = Overlap_job($feature,$starting_dir."/".$annotation_overlap_file,$bed_summit_file);
		push(@intersect_analysis_shuffle_files,$Intersect_file);
	}
	warn "Done\n";
}
warn "Overlap SHUFFLE Analysis performed\n\n";
#======================================================================

#======================Save Closest info==================
warn "Concatenating graphics_info_files...";
my @Graph_shuffle_files; #tab files containig distace graphics information
foreach my $closest_shuffle_file (@closest_analysis_shuffle_files) {
	my ($number) = $closest_shuffle_file =~ /shuffle_(\d+)_closest_/;
	my $Graph_Info_file = Graph_Info($closest_shuffle_file,'shuffle',$number);
	push(@Graph_shuffle_files,$Graph_Info_file);
}

open GRAPH_CLOSEST_INFO,'>','Graph_closest_all.tab' or die "Cannot write: $!";

foreach my $graph_file (@Graph_shuffle_files) {
	open my $fh, $graph_file or die "Cannot read $_: $!";
	while (<$fh>) {
		print GRAPH_CLOSEST_INFO $_;
	}
	close $fh;
	unlink $graph_file or die "Cannot remove $graph_file: $!";
}

foreach my $graph_file (@Graph_original_files) {
	open my $fh, "../$graph_file" or die "Cannot read ../$_: $!";
	while (<$fh>) {
		print GRAPH_CLOSEST_INFO $_;
	}
	close $fh;
	unlink "../$graph_file" or die "Cannot remove $../$graph_file: $!";
}
close GRAPH_CLOSEST_INFO;
warn "Done\n\n";
#==============================================================
chdir $starting_dir or die "Cannot change to $starting_dir directory: $!";

#=====================Intersect Analysis========================
warn "Starting multiIntersectBed Analysis...";
open my $fh, "-|", "multiIntersectBed -header -g $genome_info_file -empty -i @bed_files_FilteredPeaks;"; #generate mutiIntersect file for bedfiles provided

my $insert_peaksList_sql = "INSERT INTO Peak_Combinations (list_peaks) VALUES (?)";
my $insert_peaksList_sth = $dbh->prepare($insert_peaksList_sql);

my $insert_intersectInfo_sql  = "INSERT INTO Intersect_Programs_Info (chr_seqID,start,end,num_peaks,list_peaks,group_index,Peak_ID) VALUES (?,?,?,?,?,?,?)";
my $insert_intersectInfo_sth = $dbh->prepare($insert_intersectInfo_sql);

my $group_index = 0;
while (<$fh>) {
	chomp;
	next if ~~ /chrom\tstart\tend/;
	my @fields = split ('\t',$_);
	my $chrm_quo = $dbh->quote($fields[0]);

	if ($fields[3] == 0) {
		$group_index += 1;
	}

	my $check_existence_peakList = "SELECT list_peaks FROM Peak_Combinations WHERE list_peaks = '$fields[4]'";
	unless (check_existence($check_existence_peakList)) {
		$insert_peaksList_sth->execute($fields[4]);
	}

	my $recover_peaks_ID = "SELECT Peak_ID FROM BED_Peaks_Info WHERE start <= $fields[1] AND end >= $fields[2] AND chr_seqID LIKE $chrm_quo";
	$recover_peaks_ID .= "AND end-start >= '$min_size'" if $min_size;
	$recover_peaks_ID .= "AND end-start <= '$max_size'" if $max_size;

	my $sth_recover_peaks_ID = sth($recover_peaks_ID);
	while (my $ref = $sth_recover_peaks_ID->fetchrow_arrayref()) {
		my $peak_ID = $ref->[0];

		my $check_existence_intersect_sql  = "SELECT COUNT(*) FROM Intersect_Programs_Info WHERE chr_seqID LIKE $chrm_quo AND start = '$fields[1]' AND Peak_ID = '$peak_ID'";
		unless (check_existence($check_existence_intersect_sql)) {
			
			#my @imp_fields = splice(@fields,0,4);
			#print "@imp_fields\n";
			#push(@imp_fields,$group_index,$peak_ID);
			$insert_intersectInfo_sth->execute($fields[0],$fields[1],$fields[2],$fields[3],$fields[4],$group_index,$peak_ID);
		}
	}

}
close $fh;
warn "Done\n\n";
#===============================================================

$dbh->disconnect;
#==========Subrutines===============

sub sth { #prepare,execute a sql
	my $insert_sth = $dbh->prepare($_[0]); #prepare , look for basic sql errors 
	$insert_sth->execute(); # action takes part here!!
	return $insert_sth;
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

sub recover_info { #recover scalar info from the DB
	my ($recover_info) = $dbh->selectrow_array($_[0]);
	return $recover_info;
}

sub Closest_job {
	my $feature = shift;
	my $annotation_file = shift;
	my $bed_file = shift;

	my ($fn,$dir,$suf) = fileparse($bed_file,qr/\.[^.]*/);
	my $cmd_closest = "closestBed -D b -a $bed_file -b $annotation_file > $fn\_closest_$feature.2bed";
	!system $cmd_closest or die "Cannot execute closestBed: $!";
	return "$fn\_closest_$feature.2bed";
}

sub Overlap_job {
	my $feature = shift;
	my $annotation_file = shift;
	my ($fn_ann,$dir_ann,$suf_ann) = fileparse($annotation_file,qr/\.[^.]*/);
	my $bed_file = shift;
	my ($fn,$dir,$suf) = fileparse($bed_file,qr/\.[^.]*/);
	my $cmd_overlapping =  "intersectBed -wo -a $bed_file -b $annotation_file > $fn\_vs_$fn_ann.2bed" or die "Cannot execute intersectBed : $!"; #Generate intersect files .2bed
	!system $cmd_overlapping  or die "Cannot execute intersectBed: $!";
	return "$fn\_vs_$fn_ann.2bed";
}


sub ShuffleBed {
    #my @output_shuffle;
    my $bed_file = shift;
	my $genome_info_file = shift;
	my $index = shift;
    my ($fn,$dir,$suf) = fileparse($bed_file,qr/\.[^.]*/);
    open my $fh,"shuffleBed -i $bed_file -g $genome_info_file | sortBed -i |"  or die "Cannot execute shuffleBed|sortBed: $!";
    open OUT, '>', "$fn\_shuffle\_$index$suf" or die "Cannot write in $fn\_shuffle\_$index$suf : $!";
    print OUT while (<$fh>);
    close OUT;
    close $fh;
    return "$fn\_shuffle\_$index$suf";
}

sub Graph_Info {
	my $closest_file = shift;
	my $type_file = shift;
	my $number = shift;
	my ($fn,$dir,$suf) = fileparse($closest_file,qr/\.[^.]*/);
	my $closest_info_file = "$fn\_graphInfo.tab";
	my %peak_distance;
	my @distances;

	#======original Name============
	my @closest_file_parts = split ('_',$closest_file);
	my $part_no =0;
	my @original_name_parts;
	foreach my $closest_file_part (@closest_file_parts) {
		last if ($closest_file_part ~~ $min_size && $closest_file_parts[$part_no+1] ~~ $max_size && $closest_file_parts[$part_no+2] ~~ 'summits');
		push(@original_name_parts,$closest_file_part);
		$part_no += 1;
	}
	my $original_name = join('_',@original_name_parts);
	#===========================================
	
	#======Feature=================
	my ($feature) = $closest_file =~ /_([A-Z]+).2bed$/i; 
	#=======Distance Graph Info==================

	open my $fh,'<',$closest_file or die "Cannot open $closest_file: $!";
	while (<$fh>) {
		chomp;
		my @fields = split("\t",$_);
		my ($chr,$start,$end,$distance) = ($fields[0],$fields[1],$fields[2],$fields[-1]);
		$peak_distance{$chr.$start.$end} = $distance;
	}
	close $fh;
	
	foreach my $distance (values %peak_distance) {
		push(@distances,$distance);
		
	}

	#my @sort_distances = sort {$a <=> $b} @distances;
	
	open my $fh_out,'>',$closest_info_file or die "Cannot write into $closest_info_file: $!"; 
	
	for (my $count = -10000; $count <=10000; $count += 50) {
		my @distances_range = grep { $_ >= $count && $_ < $count + 50} @distances;
		my $average_position = ($count+($count+50))/2;
		print $fh_out "$original_name\t$number\t$feature\t$type_file\t$average_position\t".@distances_range."\n";
	}
	close $fh_out;
	return $closest_info_file;
	#===========================================

}

