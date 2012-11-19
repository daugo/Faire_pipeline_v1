#!/usr/bin/perl
#fetch_GFF3_store.pl
use warnings;
use strict;
use DBI;
use DBD::mysql;
use Data::Dumper;
use 5.010001;

my $basename = $ARGV[0];
#=========Database login info===============
my $dsn = "DBI:mysql:GFF3_store:localhost";
my $user = "du_GFF3_store";
my $passwd = "santafek";
my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #conection to the DB
#=========================================

my %Gene_RNA_Features;

my $match_parent_Gene_sql = "SELECT Type, Gene_ID FROM Gene ";
#my $fields = recover_arrayref(sth($match_parent_Gene_sql));
my $sth_match_Gene= sth($match_parent_Gene_sql);
#Output file
if (-e $basename.".bed") {
	die "File with the same name exists at $ENV{'PWD'}";
}
open BED, '>>', $basename.'.bed' or die "Cannot write in $basename.bed"; 

while (my $ref = $sth_match_Gene->fetchrow_arrayref()) {
	#print "$ref->[0]\t$ref->[1]\n";
	my ($gene_id,$type_gene) = ($ref->[1],$ref->[0]);
	#my %gene_parent= ($ref->[1] => $ref->[0]);
	#print $gene_parent{$ref->[1]}."\n";
	my $match_parent_Parent_RNA_sql = "SELECT Type, Strand, RNA_ID, seqid FROM Parent_RNA WHERE Parent_ID = '$ref->[1]'";
	my $sth_match_Parent_RNA = sth($match_parent_Parent_RNA_sql);
	
	while (my $ref = $sth_match_Parent_RNA->fetchrow_arrayref()) {
		#print "@{$ref}\n";
		#my %rna_parent =($ref->[1] => $ref->[0]);
		my ($rna_id, $type_rna,$strand, $seqid) = ($ref->[2],$ref->[0],$ref->[1],$ref->[3]);
		#print "$strand\n";
		my $match_parent_Features_sql = "SELECT ID, Feature, Start, End, Strand FROM Features WHERE Parent_ID = '$ref->[2]' ORDER BY Start ASC";
		
		#"print "$match_parent_Features_sql\n";
		my $sth_match_Features = sth($match_parent_Features_sql);
		#print "$sth_match_Features\n";
		my %exon;
		my @exons_index;
		my %CDS;
		my @CDS_index;
		while (my $ref = $sth_match_Features->fetchrow_arrayref()) {
			my ($ID, $feature,$start,$end,$strand) = ($ref->[0],$ref->[1],$ref->[2],$ref->[3],$ref->[4]);
			given ($feature) {
				when (/exon/) {
					#print "$ID\n";
					push (@exons_index,$ID);
					push @{$exon{$ID}},$feature,$start,$end;
					#open EXON, '>>', $basename.'_exon.bed' or die "Cannot write in $basename\_exon.bed";
					#print EXON "$ID\t$feature\t$start\t$end\n";
					#print "$ID\t$feature\t$start\t$end\n";
					#close EXON;
				}
				when (/CDS/) {
					push (@CDS_index,$ID);
					push @{$CDS{$ID}},$feature,$start,$end;
				}
				when (/five/) {
					print BED "$seqid\t$start\t$end\t$feature;$rna_id;$gene_id\t500\t$strand\n";
				}
				when (/three/) {
					print BED "$seqid\t$start\t$end\t$feature;$rna_id;$gene_id\t500\t$strand\n";
				}
			}
			#if ($feature ~~ /./) {
			#	push @{$Gene_RNA_Features{$gene_id}{$rna_id}{$ID}},$feature,$start,$end,$strand;
			#}
		}
		my $value=0;
		if (@exons_index) {
			foreach my $exon (@exons_index) {
				if ($strand eq '+') {
					my ($feature,$start,$end) = @{$exon{$exon}};
					print BED "$seqid\t".($start-1)."\t".($start)."\tTSS;$rna_id;$gene_id\t500\t$strand\n" if $value == 0;
					print BED "$seqid\t".($start-1)."\t".($end-1)."\t$feature\_".($value+1).";$rna_id;$gene_id\t500\t$strand\n";
					 unless ($value == @exons_index-1) {
						 my ($feature_future,$start_future,$end_future)=@{$exon{$exons_index[$value+1]}};
						print BED "$seqid\t".($end)."\t".($start_future-2)."\tintron_".($value+1).";$rna_id;$gene_id\t500\t$strand\n";
					}
				}
				elsif ($strand eq '-') {
					my $value_neg = @exons_index - $value -1;
					my ($feature,$start,$end) = @{$exon{$exons_index[$value_neg]}};
					print BED "$seqid\t".($end-2)."\t".($end-1)."\tTSS;$rna_id;$gene_id\t500\t$strand\n" if $value_neg == @exons_index -1;
					print BED "$seqid\t".($start-1)."\t".($end-1)."\t$feature\_".($value+1).";$rna_id;$gene_id\t500\t$strand\n";
					unless ($value_neg == 0) {
						my ($feature_past,$start_past,$end_past)=@{$exon{$exons_index[$value_neg-1]}};
						print BED "$seqid\t".($end_past)."\t".($start-2)."\tintron\_".($value+1).";$rna_id;$gene_id\t500\t$strand\n";
					}
				}
				$value += 1;
			}
		}
		$value =0;
		if (@CDS_index) {
			foreach my $CDS (@CDS_index) {
				if ($strand eq '+') {
					my ($feature,$start,$end) = @{$CDS{$CDS}};
					print BED "$seqid\t".($start-1)."\t".($start+2)."\tstart_codon;$rna_id;$gene_id\t500\t$strand\n" if $value == 0;
					print BED "$seqid\t".($start-1)."\t".($end-1)."\t$feature\_".($value+1).";$rna_id;$gene_id\t500\t$strand\n";
					unless ($value == @CDS_index-1) {
						my ($feature_future,$start_future,$end_future)=@{$CDS{$CDS_index[$value+1]}};
						print BED "$seqid\t".($end)."\t".($start_future-2)."\tnon-CDS_".($value+1).";$rna_id;$gene_id\t500\t$strand\n";
					}
				}
				elsif ($strand eq '-') {
					my $value_neg = @CDS_index - $value -1;
					my ($feature,$start,$end) = @{$CDS{$CDS_index[$value_neg]}};
					print BED "$seqid\t".($end-4)."\t".($end-1)."\tstart_codon;$rna_id;$gene_id\t500\t$strand\n" if $value_neg == @CDS_index -1;
					print BED "$seqid\t".($start-1)."\t".($end-1)."\t$feature\_".($value+1).";$rna_id;$gene_id\t500\t$strand\n";
					unless ($value_neg == 0) {
						my ($feature_past,$start_past,$end_past)=@{$CDS{$CDS_index[$value_neg-1]}};
						print BED "$seqid\t".($end_past)."\t".($start-2)."\tnon-CDS_".($value+1).";$rna_id;$gene_id\t500\t$strand\n";
					}
				}
				$value += 1;
			}
		}
	}
}
#print Dumper \%Gene_RNA_Features;
#print "$fields\n";
$dbh->disconnect;
close BED;

sub sth { #prepare,execute a sql
	my $insert_sth = $dbh->prepare($_[0]); #prepare , look for basic sql errors 
	$insert_sth->execute(); # action takes part here!!
	return $insert_sth;
}

