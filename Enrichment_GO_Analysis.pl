#!/usr/bin/perl
#Extract_GO_of_Interest_Genes.pl

use warnings;
use strict;


use File::Basename;
use DBI;
use DBD::mysql;
use YAML::Tiny;

use Data::Dumper;
#==================

#======Config-file===========
my $GOExtract_config_file = shift @ARGV; #config file 
die "GOExtract_config_file ($GOExtract_config_file) doesn't exists: $!" unless (-e $GOExtract_config_file); #check if config_file exists
die "GOExtract_config_file ($GOExtract_config_file) doesn't end in .yml" unless ($GOExtract_config_file =~ /\.yml$/); #check file suffix .yml mandatory
#================================

#=====Open config_file===========
my $config_file = YAML::Tiny->new;
$config_file= YAML::Tiny->read("$GOExtract_config_file") or die "Cannot read YAML config file $GOExtract_config_file: $!";
#=================================

#=========Database login info================
my $dsn = "DBI:mysql:GFF3_store:localhost";
my $user = "du_GFF3_store";
my $passwd = "santafek";

my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #conection to the DB
#============================================

#=========Take selection criteria parameters========================
my $max_interest_distance = $config_file->[0]->{'Distances'}->{'max_interest_distance'};
my $min_interest_distance = $config_file->[0]->{'Distances'}->{'min_interest_distance'};
my $reference_feature = $config_file->[0]->{'reference_feature'};

#============================================



#===============Extract Peaks different filenames===============
my $extract_datasets_names_sql  = "SELECT DISTINCT(filename) FROM BED_Peaks_Info";
my $sth_extract_datasets_names = sth($extract_datasets_names_sql);
my @filenames;
while (my $ref = $sth_extract_datasets_names->fetchrow_arrayref()) {
	my $filename = $ref->[0];
	push(@filenames,$filename);
	}
#============================================


#==============Save Filename <-> Index relation===========
my @interest_files;
my %filename_intersect_combination;

my %filename_index;
my $peakCaller_index_sql ="SELECT BED_Peaks_Info.filename,Intersect_Programs_Info.list_peaks 
FROM 
BED_Peaks_Info JOIN Intersect_Programs_Info ON BED_Peaks_Info.Peak_ID = Intersect_Programs_Info.Peak_ID 
	WHERE Intersect_Programs_Info.num_peaks =1 
	GROUP BY BED_Peaks_Info.filename";
	my $sth_peakCaller_index = sth($peakCaller_index_sql);
	while (my $ref = $sth_peakCaller_index->fetchrow_arrayref()) {
		my ($filename,$index) = ($ref->[0],$ref->[1]);
		$filename_index{$filename} = $index;
	}


#===================Interest related Genes for each peak-caller results=============

warn "Interest related Genes  for individual peak-caller results...";
foreach my $filename (@filenames) {
	my $filename_quo = $dbh->quote($filename);
	#==================================
	my $extract_GO_related_sql  = "
	SELECT G.Gene_ID 
	FROM 
	Closest_Analysis CA 
	JOIN BED_Peaks_Info BPI ON (CA.Peak_ID = BPI.Peak_ID)
	JOIN Parent_RNA PR ON (PR.RNA_ID = CA.RNA_ID) 
	JOIN Gene G ON (G.Gene_ID = PR.Parent_ID) 

	WHERE CA.Distance <= $max_interest_distance AND CA.Distance >= $min_interest_distance AND BPI.filename = $filename_quo AND CA.Feature = '$reference_feature'
	GROUP BY G.Gene_ID,CA.Peak_ID,CA.Distance 
	ORDER BY G.Gene_ID ASC";
	#==================================
	my $sth_extract_GO_related_sql = sth($extract_GO_related_sql);
	open my $fh,'>',"Interest_genes\_$filename.txt" or die "Cannot write in GO_mapping\_$filename.tab: $!";
	push (@interest_files,"Interest_genes\_$filename.txt");
	$filename_intersect_combination{"Interest_genes\_$filename.txt"} = $filename_index{$filename};
	while (my $ref = $sth_extract_GO_related_sql->fetchrow_arrayref()) {
		my $Gene_ID = $ref->[0];
		print $fh "$Gene_ID\n";
	}
	close $fh;
}
warn "Done\n";

#=================================================================================

#=================Interest related Genes for peaks supported by more than one peak-calller======
warn "Starting Interest related Genes  for peaks supported by more than one peak-calller...";
my $get_max_intersect_peaks_sql = "SELECT max(num_peaks) FROM Intersect_Programs_Info";
my $max_num_intersect_peaks = recover_info($get_max_intersect_peaks_sql);

my $get_intersect_combinations_sql = "SELECT DISTINCT(list_peaks) FROM Intersect_Programs_Info WHERE num_peaks >1";
my $sth_get_intersect_combinations = sth($get_intersect_combinations_sql);


while (my $ref = $sth_get_intersect_combinations->fetchrow_arrayref()) {
	my $intersect_combination  = $ref->[0];
	my $number_intersect_peaks = split(',',$intersect_combination);
	my $get_related_filenames_sql = "SELECT DISTINCT(BPI.filename) FROM Intersect_Programs_Info IPI JOIN BED_Peaks_Info BPI ON IPI.Peak_ID = BPI.Peak_ID  WHERE list_peaks = '$intersect_combination'";
	my $sth_get_related_filenames = sth($get_related_filenames_sql);
	my @related_filenames;
	while (my $ref = $sth_get_related_filenames->fetchrow_arrayref()) {
		my $related_filename = $ref->[0];
		push (@related_filenames,$related_filename);
		#push ( @{$intersect_combination_filenames{$intersect_combination}}, $related_filename);	
	}
	my $extract_GO_intersect_regions_sql = "
		SELECT  G.Gene_ID 
		FROM
		Intersect_Programs_Info IPI 
		JOIN BED_Peaks_Info BPI ON IPI.Peak_ID = BPI.Peak_ID 
		JOIN Closest_Analysis CA ON BPI.Peak_ID = CA.Peak_ID 
		JOIN Parent_RNA PR ON CA.RNA_ID = PR.RNA_ID 
		JOIN Gene G ON PR.Parent_ID = G.Gene_ID  

		WHERE IPI.list_peaks = '$intersect_combination' AND IPI.num_peaks = $number_intersect_peaks AND CA.Distance >= $min_interest_distance AND CA.Distance <= $max_interest_distance AND CA.Feature = '$reference_feature'
		GROUP BY G.Gene_ID,IPI.group_index 
		ORDER BY IPI.group_index, G.Gene_ID";
	my $sth_extract_GO_intersect_regions = sth($extract_GO_intersect_regions_sql);
	my $combination_label = join('_Intersect_',@related_filenames);
	open my $fh,'>',"Interest_genes\_$combination_label.txt" or die "Cannot write in Interest_genes\_$combination_label.tab: $!";
	push (@interest_files,"Interest_genes\_$combination_label.txt");
	$filename_intersect_combination{"Interest_genes\_$combination_label.txt"} = $intersect_combination;
	while (my $ref = $sth_extract_GO_intersect_regions->fetchrow_arrayref()) {
		my ($Gene_ID) = $ref->[0];
		print $fh "$Gene_ID\n";
	}
	close $fh;
}
warn "Done\n";


#==============================================

#===============Get universe ==============================
warn "Starting GO Mapping of Universe...";
my $extract_universe_sql = "
SELECT  Gene.Gene_ID, GO_Annotation.GO_ID
FROM Features 
JOIN Parent_RNA ON Features.Parent_ID = Parent_RNA.RNA_ID 
JOIN Gene ON Gene.Gene_ID = Parent_RNA.Parent_ID 
JOIN GO_Annotation ON Gene.Gene_ID = GO_Annotation.locus_name  

WHERE Features.Feature = 'exon' AND Gene.seqid IN 
	(SELECT DISTINCT(chr_seqID) FROM BED_Peaks_Info) 
GROUP BY Gene.Gene_ID,GO_Annotation.GO_ID
ORDER BY Gene.Gene_ID ASC
";
my $sth_extract_universe = sth($extract_universe_sql);
my $universe_file = "GO_mapping\_universe.tab";
open my $fh,'>',"$universe_file" or die "Cannot write in GO_mapping\_universe.tab: $!";

my %GeneID_GOIDs;
while (my $ref = $sth_extract_universe->fetchrow_arrayref()) {
	my ($Gene_ID,$GO_ID) = ($ref->[0],$ref->[1]);
	push ( @{$GeneID_GOIDs{$Gene_ID}}, $GO_ID);
}
foreach my $gene_id (sort {$a cmp $b} keys %GeneID_GOIDs) {
	my $go_ids_list  = join(', ',@{$GeneID_GOIDs{$gene_id}});
	print $fh "$gene_id\t$go_ids_list\n"
}
close $fh;

warn "Done\n";

#===================Top GO Jobs====================

warn "Starting TopGO jobs...";

my $graph_No_GO_cat =$config_file->[0]->{'graph_No_GO_cat'};
my $Num_top_GO_save =$config_file->[0]->{'Num_top_GO_save'};

my @term_enrichment_out_files;
my%termEnrichmentOut_filename;
foreach my $interest_file (@interest_files) {
	my @out_files = top_GO_job($universe_file,$interest_file,$graph_No_GO_cat,$Num_top_GO_save);
	push(@term_enrichment_out_files,@out_files);
	$termEnrichmentOut_filename{$_} = $interest_file foreach (@out_files);
}
warn "Done\n\n";
#==============================Insert GO enrichment results on the database======================

warn "Saving GO_Enrichment results on the DB...";
my $insert_GO_Enrichment_sql = "INSERT INTO GO_Enrichment (Rank,GO_ID,Annotated,Significant,Expected,Rank_Classic,classic,weight,list_peaks) VALUES (?,?,?,?,?,?,?,?,?)";
my $insert_GO_Enrichment_sth = $dbh->prepare($insert_GO_Enrichment_sql);

foreach my $term_enrichment_file (@term_enrichment_out_files) {
	open my $fh,'<',$term_enrichment_file or die "Cannot open $term_enrichment_file:$!";
	my $intersect_combination;
	foreach (keys %termEnrichmentOut_filename) {
		if ($term_enrichment_file eq $_) {
			my $interest_file = $termEnrichmentOut_filename{$_};
			$intersect_combination = $filename_intersect_combination{$interest_file};
		}
	}
	while (<$fh>) {
		chomp;
		my @fields = split('\t',$_);
		#print "@fields\t$intersect_combination\n";
		next if ($_ ~~ /GO.ID\sTerm/);
		my $check_GO_ID = "SELECT COUNT(*) FROM GO_IDs WHERE GO_ID = '$fields[1]'";

		if (check_existence($check_GO_ID)) {
			my $check_GO_Enrichment_sql  = "SELECT COUNT(*) FROM GO_Enrichment WHERE
			GO_ID = '$fields[1]'
			AND list_peaks = '$intersect_combination'";

			unless (check_existence($check_GO_Enrichment_sql)) {
				
				$insert_GO_Enrichment_sth->execute($fields[0],$fields[1],$fields[3],$fields[4],$fields[5],$fields[6],$fields[7],$fields[8],$intersect_combination);
			}
		}

	}
	close $fh;
}
warn "Done\n\n";


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

sub top_GO_job {

	my $universe_file = shift;
	my $interest_genes_file =shift;
	my $graph_No_GO_cat =shift;
	my $Num_top_GO_save =shift;

	my ($fn,$dir,$suf) = fileparse($interest_genes_file,qr/\.[^.]*/);
	$fn =~ s/Interest_genes_//;
	my @ontology_groups = qw (MF BP CC);
	my @out_files;
	foreach my $ontology_group (@ontology_groups) {	
		print "FILE:$fn\_$ontology_group.R\n";
		my $top_GO_script = "
rm(list=ls())
library(topGO)
sink('$fn\_$ontology_group\_results.txt')
ontology_group <- '$ontology_group' #MF,BP,CC
map_file <- '$universe_file'
Interesting_genes <- '$interest_genes_file'
graph_show_size <- $graph_No_GO_cat
topNo <- $Num_top_GO_save

geneID2GO <- readMappings(file = map_file)
geneNames <- names(geneID2GO)
myInterestingGenes <- read.table(Interesting_genes)

geneList <- factor(as.integer(geneNames %in% myInterestingGenes\$V1))
names(geneList) <- geneNames

GOdata <- new('topGOdata', ontology = ontology_group, allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)

test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'Fisher test')
resultFis <- getSigGroups(GOdata, test.stat)

test.stat <- new('weightCount', testStatistic = GOFisherTest, name = 'Fisher test', sigRatio = 'ratio')
resultWeight <- getSigGroups(GOdata, test.stat)

pvalFis <- score(resultFis)
pvalWeight <- score(resultWeight, whichGO = names(pvalFis))

cat ('Correlation fisher vs weigths')
cor(pvalFis, pvalWeight)

allRes <- GenTable(GOdata, classic = resultFis, weight = resultWeight,orderBy = 'weight', ranksOf = 'classic', topNodes = topNo)
write.table(allRes,file='$fn\_$ontology_group\_results.tab',sep='\t',quote=F)

pdf('$fn\_$ontology_group\_graphics.pdf')
showSigOfNodes(GOdata, score(resultFis), firstSigNodes = graph_show_size , useInfo = 'all')
showSigOfNodes(GOdata, score(resultWeight), firstSigNodes = graph_show_size, useInfo = 'all')
hist(pvalFis, 50, xlab = 'p-values')
hist(pvalWeight, 50, xlab = 'p-values')
dev.off()
";
		open my $fh,'>',"$fn\_$ontology_group.R" or die "Cannot create $fn\_$ontology_group.R: ";
		print $fh "$top_GO_script";
		close $fh;
		#!system "Rscript $fn\_$ontology_group.R 2> $fn\_$ontology_group\_report.txt" or die "Cannot execute Rscript $fn\_$ontology_group.R";
		push(@out_files,"$fn\_$ontology_group\_results.tab");
	}

	return @out_files;

}


