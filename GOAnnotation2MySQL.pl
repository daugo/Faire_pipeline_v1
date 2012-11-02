#!/usr/bin/perl
#GOAnnotation2MySQL.pl

use warnings;
use strict;

use DBI;
use DBD::mysql;
use YAML::Tiny;

#======Config-file===========
my $GOAnnotation_config_file = shift @ARGV; #config file 
die "GOAnnotation_config_file ($GOAnnotation_config_file) doesn't exists: $!" unless (-e $GOAnnotation_config_file); #check if config_file exists
die "GOAnnotation_config_file ($GOAnnotation_config_file) doesn't end in .yml" unless ($GOAnnotation_config_file =~ /\.yml$/); #check file suffix .yml mandatory
#================================

#=====Open config_file===========
my $config_file = YAML::Tiny->new;
$config_file= YAML::Tiny->read("$GOAnnotation_config_file") or die "Cannot read YAML config file $GOAnnotation_config_file: $!";
#=================================

#=========Database login info================
my $dsn = "DBI:mysql:GFF3_store:localhost";
my $user = "du_GFF3_store";
my $passwd = "santafek";

my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #conection to the DB
#============================================

my @fields_DB;
my @fields_indexes;
foreach my $field_index_label (sort {$config_file->[0]->{'columns_info_order'}->{$a} <=> $config_file->[0]->{'columns_info_order'}->{$b}} keys %{$config_file->[0]->{'columns_info_order'}}) {
	if ($config_file->[0]->{'columns_info_order'}->{$field_index_label}) {
		push (@fields_DB,$field_index_label);
		die 'Only the index of the column must be provided, entry provided is not a number between 1 to n' unless ($config_file->[0]->{'columns_info_order'}->{$field_index_label} =~ /^\d+$/);
		push (@fields_indexes,$config_file->[0]->{'columns_info_order'}->{$field_index_label});
	}
}

my $fields_list = join(', ',@fields_DB);
#print $fields_list."\n\n";

#==================
my $placeholders_list = '?'x scalar(@fields_DB);
my @placeholders = split('',$placeholders_list);
$placeholders_list = join(', ',@placeholders);
#===================

my $insert_GO_IDs_sql = "INSERT INTO GO_IDs (GO_ID) VALUES (?)";
my $insert_GO_IDs_sth = $dbh->prepare($insert_GO_IDs_sql);

my $insert_GOAnnotation_sql = "INSERT INTO GO_Annotation ($fields_list) VALUES ($placeholders_list)";
print "$insert_GOAnnotation_sql\n";
my $insert_GOAnnotation_sth = $dbh->prepare($insert_GOAnnotation_sql);

open my $fh,'<',$config_file->[0]->{'GOAnnotation_file'} or die "Cannot open GOAnnotation file: $!";

while (<$fh>) {
	chomp;
	my @imp_fields_info;
	my @imp_fields_info_quo;
	my @fields = split ('\t',$_);
	foreach my $imp_column (@fields_indexes) {
		push(@imp_fields_info,$fields[$imp_column-1]);
		push(@imp_fields_info_quo,$dbh->quote($fields[$imp_column-1]));
	}
	my $insert_values = join(', ',@imp_fields_info);
	my $count = 0;

	my $check_GO_Annotation_existence_sql = "SELECT COUNT(*) FROM GO_Annotation WHERE ";

	my $Gene_ID = $fields[$config_file->[0]->{'columns_info_order'}->{'locus_name'}-1];
	$Gene_ID = $dbh->quote($Gene_ID);

	my $check_gene_id_sql = "SELECT COUNT(*) FROM Gene WHERE Gene_ID LIKE $Gene_ID";
	if (check_existence($check_gene_id_sql)) {
		my  @checks_statments; 
		foreach (@fields_DB) {
			my $check_statment= "$_ LIKE $imp_fields_info_quo[$count]";
			push (@checks_statments,$check_statment);
			$count += 1;
		} 
		$check_GO_Annotation_existence_sql.= join(' AND ',@checks_statments);

		#=============Fill GO_IDs Table=======================

		my $GO_ID_index = ($config_file->[0]->{'columns_info_order'}->{'GO_ID'});
		my $GO_ID  = $fields[$GO_ID_index-1];
		my $check_GO_IDs_sql = "SELECT COUNT(*) FROM GO_IDs WHERE GO_ID = '$GO_ID'";
		unless (check_existence($check_GO_IDs_sql)) {
			$insert_GO_IDs_sth->execute($GO_ID);
		}
		#================================================

		#===========Fill GO Annotation Table =====================
		unless (check_existence($check_GO_Annotation_existence_sql)) {
			$insert_GOAnnotation_sth->execute(@imp_fields_info);
		}
	}
	else {
		warn "Gene ID: $Gene_ID is not saved in the database";
	}


}
close $fh;
$dbh->disconnect;



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
