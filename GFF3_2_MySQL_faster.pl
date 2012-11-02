#!/usr/bin/perl
#GFF3_2_MySQL.pl
use warnings;
use strict;
use File::Basename;
use DBI;
use DBD::mysql;
use YAML::Tiny;
use Bio::SeqIO;
#use Bio::DB::Fasta;
#use 5.010001;

#===========Database login info=============
my $dsn = "DBI:mysql:GFF3_store2:localhost";
my $user = "du_GFF3_store2";
my $passwd = "santafek";

my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #conection to the DB
#===================================

#======Config-file===========
my $bed_manipulation_config_file = shift @ARGV; #config file 
die "Config_file ($bed_manipulation_config_file) doesn't exists: $!" unless (-e $bed_manipulation_config_file); #check if config_file exists
die "Config_file ($bed_manipulation_config_file) doesn't end in .yml" unless ($bed_manipulation_config_file =~ /\.yml$/);
#================================

#=====Open config_file===========
my $config_file = YAML::Tiny->new;
# Open the config
$config_file= YAML::Tiny->read("$bed_manipulation_config_file") or die "Cannot read YAML config file $bed_manipulation_config_file: $!";
#=================================


my $fasta_file = $config_file->[0]->{'fasta_file'};
my ($fn_fasta,$dir_fasta,$suf_fasta) = fileparse($fasta_file,qr/\.[^.]*/);
die "Fasta suf must match '.fa' or .fasta" unless ($suf_fasta =~ /fa$|fasta$/i);
#print $dir_fasta.$fn_fasta.$suf_fasta."\n";
die "Fasta file $fasta_file doesn't exists: $!" unless (-e $dir_fasta.$fn_fasta.$suf_fasta);
my $gff3_file =  $config_file->[0]->{'gff3_file'};
die "Gff3 file $gff3_file doesn't exists: $!" unless (-e $gff3_file);

#==========Get Chr ID_lengths========================
my @Seq_IDs; #array where the different identifiers of the mutifasta are going to be saved
my %ID_length; #where lengths of each seqID are going to be saved
my $seqio = Bio::SeqIO->new(-file => "$fasta_file", '-format' => 'Fasta'); #read the fasta file


#foreach sequence in the fasta file  -> get sequence lengths and stored
while ( my $seq = $seqio -> next_seq) {
        my $id = $seq->id; #get id (ex. Chr)
        push (@Seq_IDs,$id); #save id
        my $length = $seq->length; #get sequence length
        $ID_length {$id} = $length;#save length indexed by seqid
}

#============================

#Save all the seqids in @seqid
open IN, $gff3_file or die "Cannot open gff3 file provided in the config file #bed_manipulation_config_file"; #open file gff for parsing and fill the Genes table


#save seq_ids names
my @seqid = ();
while (<IN>) {
	next if $_~~ /(^#)|(^$)/; #skip the headers
	chomp $_;
	my @fields = split(/\t/,$_);
	my $seqid = $fields[0];
	if ($seqid ~~ @Seq_IDs) {
		push(@seqid,$seqid);
	}
	else  {
		warn "Chromosomes ids don't match with the seqids listed into GFF3 file
		Check that gff and fasta Chromosomes Ids match";
	}
	

}
close IN; 

#my insert sql queries with placeholders to increase the insert speed

#insert seq_id into seq_id table
my $insert_seqid_sql = "INSERT INTO SEQID_LANDMARK (Name,Length) VALUES (?,?)"; #Fill SEQID_LANDMARK Table
my $insert_seqid_sth = $dbh->prepare($insert_seqid_sql);

#insert gene info en Gene table
my $insert_geneinfo_sql = "INSERT INTO Gene (seqid,Source,Type,Start,End,Strand,Gene_ID,Description_GFF) VALUES (?,?,?,?,?,?,?,?)";
my $insert_geneinfo_sth = $dbh->prepare($insert_geneinfo_sql);

#insert RNAs info in Parent_RNA
my $insert_rnainfo_sql = "INSERT INTO Parent_RNA (seqid, Source,Type,Start,End,Strand,RNA_ID,Description_GFF,Parent_ID) VALUES (?,?,?,?,?,?,?,?,?)";
my $insert_rnainfo_sth = $dbh->prepare($insert_rnainfo_sql);

#inser Features_info into Features table
my $insert_featureinfo_sql = "INSERT INTO Features (seqid, Source, Feature, Start, End, Strand, Phase, Description_GFF, Parent_ID) VALUES (?,?,?,?,?,?,?,?,?)";
my $insert_featureinfo_sth = $dbh->prepare($insert_featureinfo_sql);


#find the unique seqids
my @uniques = uniq(@seqid); #list without duplicates
foreach my $seq_id (@uniques) {
	#$seq_id = $dbh->quote($seq_id);
	my $check_seqid_new_sql = "SELECT COUNT(*) FROM SEQID_LANDMARK WHERE Name = '$seq_id'"; #Check if name is already present in the table
	unless (check_existence($check_seqid_new_sql)) {
		#my $insert_seqid_sql = "INSERT INTO SEQID_LANDMARK (Name) VALUES ($seq_id)"; #Fill SEQID_LANDMARK Table
		#sth($insert_seqid_sql);
		$insert_seqid_sth->execute($seq_id,$ID_length{$seq_id});
	}
	else {
		warn "seqid $seq_id had already been saved in the SEQID_LANDMARK table\n";
	}
}

open IN, $gff3_file; #open file gff for parsing and fill the Genes table

while (my $line = <IN>) {
	next if $line ~~ /(^#)|(^$)/; #skip the headers , metainfo
	chomp $line;
	#print $line."\n";
	
	#recover fields spliting by tab
	my ($seqid, $source, $type, $start, $end, undef, $strand, $phase, $attributes) = split(/\t/,$line);
	
	#quoting the necesary fields
	my $seqid_quote = $dbh->quote($seqid);
	my $source_quote = $dbh->quote($source);
	my $type_quote = $dbh->quote($type);
	my $attributes_quote = $dbh->quote($attributes);
	
	if  ( $attributes !~ (/Parent=/) and $strand !~ (/\./) and $type !~ /protein/i) { #Check if the line contain a gene-like  type (excludes features, RNAs and Protein), strand check allows to  exclude landmark lines that could be present 
		my ($gene_id) = $attributes =~ /ID=([^;]+)/; #recover gene_id
		my $gene_id_quote = $dbh->quote($gene_id);
		
		my $check_geneid_new_sql = "SELECT COUNT(*) FROM Gene WHERE Gene_ID = $gene_id_quote"; #Check if id is already present in the table
		unless (check_existence($check_geneid_new_sql)) {
			$insert_geneinfo_sth->execute($seqid, $source,$type,$start,$end,$strand,$gene_id,$attributes);
		}
		else {
			warn "gene_id $gene_id had already been saved in the Gene table\n";
		}
	}
	if ( $attributes ~~ /Parent=([^;]+)/) { #select lines that have a Parent= tag (RNAs, Features)
		my $Parent = $1;
		#print "Parent\t$Parent\n";
		my @Parents = split(',',$Parent);
		foreach my $parent_id (@Parents) {
			#print "parentid : $line\t$_\n";
			my $parent_id_quote = $dbh->quote($parent_id);
			
			#check if parent IDs had already been saved in Parents tables
			my $match_parent_Gene_sql = "SELECT Gene_ID FROM Gene WHERE Gene_ID = $parent_id_quote"; #for RNAs lines
			my $Gene_ID = recover_info($match_parent_Gene_sql);
			my $match_parent_rna_sql = "SELECT RNA_ID FROM Parent_RNA WHERE RNA_ID = $parent_id_quote";#for Features lines
			my $RNA_ID = recover_info($match_parent_rna_sql);
			
			if ($Gene_ID) { # In case the current line contains RNA-level info
				#print "mRNAs : $line \n";
				my ($rna_id) = $attributes =~ /ID=([^;]+)/; #recover rna_id
				my $rna_id_quote = $dbh->quote($rna_id);
				
				my $check_rna_id_new_sql = "SELECT COUNT(*) FROM Parent_RNA WHERE RNA_ID =$rna_id_quote";
				unless (check_existence($check_rna_id_new_sql)) {
					#my $insert_rnainfo_sql = "INSERT INTO Parent_RNA (seqid, Source,Type,Start,End,Strand,RNA_ID,Description_GFF,Parent_ID) VALUES ($seqid_quote, $source_quote, $type_quote, '$start', '$end', '$strand', $rna_id_quote, $attributes_quote, $parent_id_quote)";
					#print $insert_rnainfo_sql."\n";
					#print $line."\t".$rna_id."\t".$parent_id."\n";
					#sth($insert_rnainfo_sql);
					$insert_rnainfo_sth->execute($seqid, $source, $type, $start, $end, $strand, $rna_id, $attributes, $parent_id);
				}
				else {
					warn "rna_id $rna_id had already been saved in the Gene table \n";
				}
			}

			elsif ($RNA_ID) { #In case the current line contains Features-level info
				#print "features : $line\n";
				my $check_feature_new_sql = "SELECT COUNT(*) FROM Features WHERE seqid = $seqid_quote AND Feature = $type_quote AND Start = '$start' AND End = '$end' AND Parent_ID = $parent_id_quote";
				#print $check_feature_new_sql."\n";
				unless (check_existence($check_feature_new_sql)) {
					my $check_if_empty_sql = "SELECT COUNT(*) FROM Features";
					unless(recover_info($check_if_empty_sql)) { #when table is empty -> reset AUTO_INCREMENT
						my $reset_autoincrement_sql = "ALTER TABLE Features AUTO_INCREMENT=1";
						sth($reset_autoincrement_sql);
					}
					#my $insert_featureinfo_sql = "INSERT INTO Features (seqid, Source, Feature, Start, End, Strand, Phase, Description_GFF, Parent_ID) VALUES ($seqid_quote, $source_quote, $type_quote, '$start', '$end', '$strand', '$phase', $attributes_quote, $parent_id_quote)";
					#sth($insert_featureinfo_sql);
					$insert_featureinfo_sth->execute($seqid, $source, $type, $start, $end, $strand, $phase, $attributes, $parent_id);
				}
				else {
					warn "This feature had already been saved in the Gene table \n";
				}
			}
		}
	}
}

close IN;
$dbh->disconnect;



sub uniq { #check the diffent posibilities that are contained in an array
	my %seen = ();
	my @uniques = ();
	foreach my $seqid (@_) {
		unless ($seen{$seqid}) {
			push @uniques, $seqid;
			#print "$seqid\n";
			$seen{$seqid} = 1;
		}
	}
	return @uniques;
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

sub sth { #prepare,execute a sql
	my $insert_sth = $dbh->prepare($_[0]);#prepare , look for basic sql errors 
	$insert_sth->execute(); # action takes part here!!
	return 1;
}

sub recover_info { #recover scalar info from the DB
	my ($recover_info) = $dbh->selectrow_array($_[0]);
	return $recover_info;
}
