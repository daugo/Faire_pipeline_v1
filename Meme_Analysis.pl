#!/usr/bin/perl
#Meme_Analysis.pl

use warnings;
use strict;

use File::Basename;
use DBI;
use DBD::mysql;
use YAML::Tiny;
use Bio::DB::Fasta;
use XML::Twig;

#==================

#======Config-file===========
my $Motif_config_file = shift @ARGV; #config file 
die "Motif_config_file ($Motif_config_file) doesn't exists: $!" unless (-e $Motif_config_file); #check if config_file exists
die "Motif_config_file ($Motif_config_file) doesn't end in .yml" unless ($Motif_config_file =~ /\.yml$/); #check file suffix .yml mandatory
#================================

#=====Open config_file===========
my $config_file = YAML::Tiny->new;
$config_file= YAML::Tiny->read("$Motif_config_file") or die "Cannot read YAML config file $Motif_config_file: $!";
#=================================

#=========Database login info================
my $dsn = "DBI:mysql:GFF3_store:localhost";
my $user = "du_GFF3_store";
my $passwd = "santafek";

my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #conection to the DB
#============================================

my $up_interest_distance = $config_file->[0]->{'Distances'}->{'up_interest_distance'};
my $down_interest_distance = $config_file->[0]->{'Distances'}->{'down_interest_distance'};
my $reference_feature = $config_file->[0]->{'reference_feature'};
my $fasta_file_path = $config_file->[0]->{'fasta_file_path'};

my $flanking_size = $config_file->[0]->{'flanking_size'};
my $promoters_distance = $config_file->[0]->{'promoters_distance'};
my $TSS_file = "GFF3_all_features_TSS_mod_unique.bed";
my $SC_file = "GFF3_all_features_SC_mod_unique.bed";
my $meme_path = $config_file->[0]->{'meme_path'};
#============================================
my $extract_datasets_names_sql  = "SELECT DISTINCT(filename) FROM BED_Peaks_Info";
my $sth_extract_datasets_names = sth($extract_datasets_names_sql);
my @filenames;

while (my $ref = $sth_extract_datasets_names->fetchrow_arrayref()) {
	my $filename = $ref->[0];
	push(@filenames,$filename);
}

#=================Flaking Sequences of interesting genes ==================

warn "Extract flanking interesting sequences...\n";
my @flanking_interest_sequences;

foreach my $filename (@filenames) {
	my $select_interesting_peaks ="
	SELECT BPI.chr_seqID, BPI.Calculate_Summit, BPI.score, CA.Distance
	FROM BED_Peaks_Info BPI JOIN Closest_Analysis CA ON BPI.Peak_ID = CA.Peak_ID 
	JOIN Parent_RNA PR ON PR.RNA_ID = CA.RNA_ID 

	WHERE CA.Distance <= $up_interest_distance AND CA.Distance >= $down_interest_distance AND CA.Feature = '$reference_feature' AND filename = '$filename' 
	GROUP BY PR.Parent_ID 
	ORDER BY BPI.chr_seqID,BPI.Calculate_Summit";
	my ($fn,$dir,$suf) = fileparse($filename,qr/\.[^.]*/);
	my $selected_peaks_file = "$fn\_selected_peaks.txt";
	open my $fh,'>',$selected_peaks_file or die "Cannot write in $selected_peaks_file: $!";
	my $sth_select_interesting_peaks = sth($select_interesting_peaks);
	while (my $ref = $sth_select_interesting_peaks->fetchrow_arrayref()) {
		my ($chr,$summit,$score,$distance) = ($ref->[0],$ref->[1],$ref->[2],$ref->[3]);
		
		print $fh "$chr\t$summit\t$filename\t$score\t$distance\n";
	}
	close $fh;
	my $flanking_sequences_out = get_flanking_seq($fasta_file_path,$selected_peaks_file,$flanking_size);
	push(@flanking_interest_sequences,$flanking_sequences_out);
}

warn "Done\n";

# ==========================================================================

#================Background Sequences =======================
warn "Extractiong Background sequnces...\n";
my $promoters_seq_out_file;
$promoters_seq_out_file = get_flanking_TSS_SC($fasta_file_path,$TSS_file,$promoters_distance) if ($reference_feature ~~ /TSS/i);
$promoters_seq_out_file = get_flanking_TSS_SC($fasta_file_path,$SC_file,$promoters_distance) if ($reference_feature ~~ /SC/i);
warn "Done\n";

warn "Performing fasta-get-markov for background sequences...\n";
my ($fn_promoters,$dir_promoters,$suf_promoters) = fileparse($promoters_seq_out_file,qr/\.[^.]*/);
!system "$meme_path/fasta-get-markov -m 5 < $promoters_seq_out_file > $fn_promoters.markov" or die "Cannot execute fasta-get-markov : $!";
warn "Done\n";

#===============================================================

#==============Meme-chip command and parse of results (meme and mast) to database============


my $background = "$fn_promoters.markov";

my $db_jaspar = $config_file->[0]->{'Meme_Chip'}->{'db_jaspar'};
my $nmeme = $config_file->[0]->{'Meme_Chip'}->{'nmeme'};
my $ccut = $config_file->[0]->{'Meme_Chip'}->{'ccut'};
my $run_mast = $config_file->[0]->{'Meme_Chip'}->{'run_mast'};
my $run_ama = $config_file->[0]->{'Meme_Chip'}->{'run_ama'};

my $meme_mod = $config_file->[0]->{'Meme_Chip'}->{'meme_mod'};
my $meme_minw = $config_file->[0]->{'Meme_Chip'}->{'meme_minw'};
my $meme_maxw = $config_file->[0]->{'Meme_Chip'}->{'meme_maxw'};
my $meme_nmotifs = $config_file->[0]->{'Meme_Chip'}->{'meme_nmotifs'};
my $meme_minsites = $config_file->[0]->{'Meme_Chip'}->{'meme_minsites'};
my $meme_maxsites = $config_file->[0]->{'Meme_Chip'}->{'meme_maxsites'};
my $meme_time = $config_file->[0]->{'Meme_Chip'}->{'meme_time'};
my $meme_p = $config_file->[0]->{'Meme_Chip'}->{'meme_p'};
my $meme_norevcomp = $config_file->[0]->{'Meme_Chip'}->{'meme_norevcomp'};
my $meme_pal = $config_file->[0]->{'Meme_Chip'}->{'meme_pal'};


my $insert_motifs_sql = "INSERT INTO Meme_motifs (motif_ID,width,sites,ic,re,llr,e_value,bayes_thr,elapsed_time,Reg_ex,filename) VALUES (?,?,?,?,?,?,?,?,?,?,?)";
my $insert_motifs_sth = $dbh->prepare($insert_motifs_sql);
my $insert_sites_sql = "INSERT INTO Meme_Sites (motifID,seq_meme_id,Peak_ID,position,strand,p_value,left_flank,site,rigth_flank,filename) VALUES (?,?,?,?,?,?,?,?,?,?)";
my $insert_sites_sth = $dbh->prepare($insert_sites_sql);
my $insert_mast_meme_sql = "INSERT INTO Meme_Mast (Peak_ID, p_value,e_value) VALUES (?,?,?)";
my $insert_mast_meme_sth = $dbh->prepare($insert_mast_meme_sql);


foreach my $flanking_sequence_file (@flanking_interest_sequences) {
	my ($filename) = $flanking_sequence_file =~ /^(.+)_selected_peaks_Flanking/;
	my $filename_qu = $dbh->quote($filename);
	warn "Running MEME-CHIP for $filename (May take several time depending of the number of sequences to analyse)...\n";
	my $meme_chip_folder = meme_job($meme_path,$background,$flanking_sequence_file,$db_jaspar,$nmeme,$ccut,$run_mast,$run_ama,$meme_mod,$meme_minw,$meme_maxw,$meme_nmotifs,$meme_minsites,$meme_maxsites,$meme_time,$meme_p,$meme_norevcomp,$meme_pal);
	
	warn "Saving $filename motifs info into the DB...\n";
	chdir "$meme_chip_folder/meme_out/" or die "Cannot change directory: $!";
	my @out_files = parse_Meme_XML('meme.xml');
	foreach my $out_file (@out_files) {
		
		open IN_motif,'<',$out_file or die "Cannot open $out_file: $!";
		if ($out_file =~ /motifs\.csv$/) {
			while (<IN_motif>) {
				chomp;
				#MotifID\tWidth\tSites\tic\tre\tllr\te_value\tbayes_threshold\telapsed_time\tRegularExp\n";
				next if ($_ =~ /^MotifID\tWidth/);
				my ($motif_id,$width,$sites,$ic,$re,$llr,$e_value,$bayes_thr,$elapsed_time,$Reg_ex) = split ('\t',$_);
				my $motif_id_qu = $dbh->quote($motif_id);
				
				my $check_motif_sql = "SELECT COUNT(*) FROM Meme_motifs WHERE motif_ID = $motif_id_qu AND filename = $filename_qu";
				
				unless (check_existence($check_motif_sql)) {
					
					$insert_motifs_sth->execute($motif_id,$width,$sites,$ic,$re,$llr,$e_value,$bayes_thr,$elapsed_time,$Reg_ex,$filename);
				}
			}
		}
		elsif ($out_file =~ /sites\.csv$/) {
			while (<IN_motif>) {
				chomp;
				#"MotifID\tSeqID\tSeqName\tPosition\tStrand\tpvalue\tleft_flank\tsite\tright_flank\n";
				next if ($_ =~ /^MotifID\tSeqID/);
				
				my ($motif_id,$seq_meme_id,$seq_name,$position,$strand,$pvalue,$left_flank,$site,$right_flank)= split('\t',$_);
				
				my ($chr_id,$summit) = $seq_name =~ /^(.+)_(\d+)_/;

				if ($strand eq 'plus') {
					$strand = '+';
				}
				elsif ($strand eq 'minus') {
					$strand = '-';
				}

				my $chr_id_qu = $dbh->quote($chr_id);
				my $recover_peakID_sql  = "SELECT Peak_ID FROM BED_Peaks_Info WHERE filename = $filename_qu AND chr_seqID = $chr_id_qu AND Calculate_Summit = '$summit'";
				
				my $Peak_ID = recover_info($recover_peakID_sql);
				my $motif_id_qu = $dbh->quote($motif_id);
				my $check_motif_id_sql = "SELECT COUNT(*) FROM Meme_motifs WHERE motif_ID = $motif_id_qu";
				my $check_site_sql  = "SELECT COUNT(*) FROM Meme_Sites WHERE motifID = $motif_id_qu AND Peak_ID = $Peak_ID";
				if (check_existence($check_motif_id_sql) == 1 and check_existence($check_site_sql)==0) {
					$insert_sites_sth->execute($motif_id,$seq_meme_id,$Peak_ID,$position,$strand,$pvalue,$left_flank,$site,$right_flank,$filename);
				}
				
			}
		}
	close IN_motif;
	}
	chdir "../../";


	chdir "$meme_chip_folder/meme_mast_out/" or die "Cannot change directory: $!";
	my $out_mast_file = parse_Mast_XML('mast.xml');

	open IN_mast,'<',$out_mast_file or die "Cannot open $out_mast_file; $!";
	while (<IN_mast>) {
		chomp;
		#"ID\tName\tCombined p value\tE-value\n";
		next if ($_=~ /^ID\tName\tCombined/);
		my ($mast_id,$name,$combined_p_value,$e_value) = split('\t',$_);
		my ($chr_id,$summit) = $name =~ /^(.+)_(\d+)_/;
		my $chr_id_qu = $dbh->quote($chr_id);
		my $recover_peakID_sql  = "SELECT Peak_ID FROM BED_Peaks_Info WHERE filename = $filename_qu AND chr_seqID = $chr_id_qu AND Calculate_Summit = '$summit'";
		my $Peak_ID = recover_info($recover_peakID_sql);
		my $check_Peak_ID_sql = "SELECT COUNT(*) FROM  Meme_Sites WHERE Peak_ID = '$Peak_ID'";
		my $check_mast_sql  = "SELECT COUNT(*) FROM Meme_Mast WHERE Peak_ID = '$Peak_ID'";
		if (check_existence($check_Peak_ID_sql) > 0 and check_existence($check_mast_sql) == 0) {
			$insert_mast_meme_sth->execute($Peak_ID,$combined_p_value,$e_value);
		}
	}
	close IN_mast;
	chdir "../../";
	warn "Done\n";
}




#===========================================================
$dbh->disconnect;

sub sth { #prepare,execute a sql
	my $insert_sth = $dbh->prepare($_[0]); #prepare , look for basic sql errors 
	$insert_sth->execute(); # action takes part here!!
	return $insert_sth;
}

sub recover_info { #recover scalar info from the DB
	my ($recover_info) = $dbh->selectrow_array($_[0]);
	return $recover_info;
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

sub get_flanking_seq {
	my $fasta_file_path = shift;
	my $peaks_file = shift;
	my $flanking_size = shift;

	my $fasta_db = Bio::DB::Fasta->new($fasta_file_path) or die "Cannot create db from fasta file: $!"; #create the index of the fasta file

	my ($filename,$directories,$suffix) = fileparse($peaks_file,qr/\.[^.]*/);

	open IN,'<',$peaks_file or die "Cannot open $peaks_file : $!";
	my $flanking_seq_out = "$filename\_Flanking$flanking_size.fasta";
	open OUT,'>',$flanking_seq_out or die "Cannot create fasta output file: $!";
	while (<IN>) {
	        chomp;
	        my ($chr,$summit,$program) = split('\t',$_);
	        my $seqobj = $fasta_db-> subseq($chr,$summit-$flanking_size,$summit+$flanking_size);
	        print OUT ">$chr\_$summit\_$program\n" if $seqobj;
	        print OUT "$seqobj\n" if $seqobj;
	}
	close IN;
	close OUT;
	return $flanking_seq_out;
}


sub get_flanking_TSS_SC {
	my $fasta_file_path = shift;
	my $TSSorSC_BED = shift;
	my $flanking_size = shift;

	my $fasta_db = Bio::DB::Fasta->new($fasta_file_path) or die "Cannot create db from fasta file: $!"; #create the index of the fasta file

	my ($filename,$directories,$suffix) = fileparse($TSSorSC_BED,qr/\.[^.]*/);

	open IN,'<',"$TSSorSC_BED" or die "Cannot open TSS or Start codon file : $!";
	my $promoters_seq_out_file = "$filename\_promoters\_$flanking_size.fasta";
	open OUT,'>',$promoters_seq_out_file or die "Cannot create TSS or Start codon flanking sequences file : $!";

	while (<IN>) {
	        chomp;
	        #Extract promoter regions
	        my ($chr,$location,undef,$name,undef,$strand) = split('\t',$_);
	        my $fasta_label = ">$chr\_$location\_$strand\_$name";
	        #print "$chr\_$location\_$feature\_$strand\_$name\n";

	        my $seqobj;
	        #+++++++++++.==========> Strand +
	        $seqobj = $fasta_db-> subseq($chr,$location-$flanking_size,$location) if $strand ~~ /\+/;
	        #<==========.----------- Strand -
	        $seqobj = $fasta_db-> subseq($chr,$location,$location+$flanking_size) if $strand ~~ /-/;
	        print OUT $fasta_label."\n" if $seqobj;
	        print OUT $seqobj."\n" if $seqobj;
	        }
	close OUT;
	close IN;
	return $promoters_seq_out_file;
}


sub meme_job {
	my $meme_path = shift;
	my $background = shift;	my $sequences_file= shift; my $db_jaspar = shift;
	my $nmeme = shift; my $ccut = shift;
	my $run_mast = shift; my $run_ama = shift;

	my $meme_mod = shift;
	my $meme_minw = shift; my $meme_maxw = shift; 
	my $meme_nmotifs = shift; 
	my $meme_minsites = shift; my $meme_maxsites = shift;
	my $meme_time = shift;
	my $meme_p = shift;
	my $meme_norevcomp = shift; my $meme_pal = shift;

	my ($fn,$dir,$suf) = fileparse($sequences_file,qr/\.[^.]*/);

	my $meme_chip_cmd = "$meme_path/meme-chip ";
	$meme_chip_cmd .= "-o $fn ";
	$meme_chip_cmd .= "-db $db_jaspar " if $db_jaspar;
	$meme_chip_cmd .= "-bfile $background " if $background;
	$meme_chip_cmd .= "-ccut $ccut " if $ccut;
	$meme_chip_cmd .= "-nmeme $nmeme " if $nmeme;
	$meme_chip_cmd .= "-run-mast " if $run_mast;
	$meme_chip_cmd .= "-run-ama " if $run_ama;
	$meme_chip_cmd .= "-meme-mod $meme_mod " if $meme_mod;
	$meme_chip_cmd .= "-meme-minw $meme_minw " if $meme_minw;
	$meme_chip_cmd .= "-meme-maxw $meme_maxw " if $meme_maxw;
	$meme_chip_cmd .= "-meme-nmotifs $meme_nmotifs " if $meme_nmotifs;
	$meme_chip_cmd .= "-meme-minsites $meme_minsites " if $meme_minsites;
	$meme_chip_cmd .= "-meme-maxsites $meme_maxsites " if $meme_minsites;
	$meme_chip_cmd .= "-meme-time $meme_time " if $meme_time;
	$meme_chip_cmd .= "-meme-p $meme_p " if $meme_p;
	$meme_chip_cmd .= "-meme-norevcomp $meme_norevcomp " if $meme_norevcomp;
	$meme_chip_cmd .= "-meme-pal $meme_pal " if $meme_pal;
	$meme_chip_cmd .= "$sequences_file";
	print "MEME-CHIP-CMD: $meme_chip_cmd\n";
	#!system $meme_chip_cmd or die "Cannot execute meme-chip for $sequences_file : $!";
	return $fn;
}

sub parse_Meme_XML{
	my $file = shift;
	my ($fn,$dir,$suf) = fileparse($file,qr/\.[^.]*/);
	my $twig=XML::Twig->new();
	$twig->parsefile($file) or die  'cannot parse $infile';
	my %seqmap;
	my $out_motifs=$fn.'.motifs.csv';
	my $out_sites =$fn.'.sites.csv';
	my $out_motiffasta=$fn.'.motifs.fasta';
	open MOTIFS, ">$out_motifs" or die 'cannot open $out_motifs for writting';
	print MOTIFS "MotifID\tWidth\tSites\tic\tre\tllr\te_value\tbayes_threshold\telapsed_time\tRegularExp\n";
	open SITES,  ">$out_sites"  or die 'cannot open $out_sites  for writting';
	print SITES "MotifID\tSeqID\tSeqName\tPosition\tStrand\tpvalue\tleft_flank\tsite\tright_flank\n";
	open MOTIFSEQ, ">$out_motiffasta";
	my $root= $twig->root;
	my $trainingset=$root->first_child( 'training_set');
	my @sequencemapping=$trainingset->children( 'sequence');
	foreach my $seqmap(@sequencemapping){
		$seqmap{$seqmap->att('id')}=$seqmap->att('name');
	}
	my @motifs= $root->children( 'motifs');
	foreach my $motifs(@motifs){
		my @motif = $motifs->children('motif');
		foreach my $motif(@motif){
			my $motif_id = $motif->att('id');
			my $motif_width = $motif->att('width');
			my $motif_sites = $motif->att('sites');
			my $motif_ic = $motif->att('ic');
			my $motif_re = $motif->att('re');
			my $motif_llr = $motif->att('llr');
			my $motif_e_value = $motif->att('e_value');
			my $motif_bayes_threshold = $motif->att('bayes_threshold');
			my $motif_elapsed_time = $motif->att('elapsed_time');
			my $regexp = $motif->first_child('regular_expression')->text;
			$regexp=~s/\n//g;
			print MOTIFS "$motif_id\t$motif_width\t$motif_sites\t$motif_ic\t$motif_re\t$motif_llr\t$motif_e_value\t$motif_bayes_threshold\t$motif_elapsed_time\t$regexp\n";
			my @contributing_sites=$motif->children('contributing_sites');
			foreach my $contributing_sites(@contributing_sites){
				my @contributing_site=$contributing_sites->children('contributing_site');
				foreach my $contributing_site(@contributing_site){
					my $contributing_site_seq_id=$contributing_site->att('sequence_id');
					my $contribuing_site_seq_name = $seqmap{$contributing_site_seq_id};
					my $contributing_site_position=$contributing_site->att('position');
					my $contributing_site_strand=$contributing_site->att('strand');
					my $contributing_site_pvalue=$contributing_site->att('pvalue');
					my $contributing_site_left_flank=$contributing_site->first_child('left_flank')->text;
					my $contributing_site_right_flank=$contributing_site->first_child('right_flank')->text;
					my @site=$contributing_site->children('site');
					my @letters=$site[0]->children('letter_ref');
					my @hit_seq;
					my $poscounter=0;
					foreach my $letter(@letters){
						my $residue= $letter->att('letter_id');
						$residue=~s/\n//g;
						$residue=~s/letter_//g;
						$hit_seq[$poscounter]=$residue;
						$poscounter++;
					}
					my $hit_seq=join('',@hit_seq);
					print MOTIFSEQ ">".$motif_id."_".$seqmap{$contributing_site_seq_id}."\n$hit_seq\n";
					print SITES "$motif_id\t$contributing_site_seq_id\t$contribuing_site_seq_name\t$contributing_site_position\t$contributing_site_strand\t$contributing_site_pvalue\t$contributing_site_left_flank\t$hit_seq\t$contributing_site_right_flank\n";

				}
			}
		}
	}
	close MOTIFS; close SITES; close MOTIFSEQ;
	return ($out_motifs,$out_sites,$out_motiffasta);
#	$twig->print;
}


sub parse_Mast_XML{
	my $file = shift;
	my ($fn,$dir,$suf) = fileparse($file,qr/\.[^.]*/);
	my $twig=XML::Twig->new();
	$twig->parsefile($file) or die  'cannot parse $infile';
	my $out_hits=$fn.'.hits.csv';
	open HITS, ">$out_hits" or die 'cannot open $out_motifs for writting';
	print HITS "ID\tName\tCombined p value\tE-value\n";
	my $root= $twig->root;
	my @sequences= $root->children( 'sequences');
	foreach my $sequences(@sequences){
		my @sequence = $sequences->children('sequence');
		foreach my $sequence(@sequence){
			my $sequence_id = $sequence->att('id');
			my $sequence_name = $sequence->att('name');
			my @score = $sequence->children('score');
			foreach my $score(@score){
				my $score_combined_pvalue=$score->att('combined_pvalue');
				my $score_evalue=$score->att('evalue');
				print HITS "$sequence_id\t$sequence_name\t$score_combined_pvalue\t$score_evalue\n";
			}
		}
	}
	close HITS;
	return ($out_hits);
#	$twig->print;
}
