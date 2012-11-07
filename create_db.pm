package create_db;
#Load_database;
use warnings;
use strict;
use YAML::Tiny;
use DBI;


sub Load_DB_Schema {
	#======Config-file===========
	my $bed_config_file = shift; #config file 
	#print "$bed_config_file\n";
	die "Config_file ($bed_config_file) doesn't exists: $!" unless (-e $bed_config_file); #check if config_file exists
	die "Config_file ($bed_config_file) doesn't end in .yml" unless ($bed_config_file =~ /\.yml$/);
	#================================

	#=====Open config_file===========
	my $config_file = YAML::Tiny->new;
	# Open the config
	$config_file= YAML::Tiny->read("$bed_config_file") or die "Cannot read YAML config file $bed_config_file: $!";
	#=================================
	my $dump_schema = shift;

	my $host = $config_file->[0]->{'MySQL'}->{'host'};
	my $user = $config_file->[0]->{'MySQL'}->{'user'};
	my $pass = $config_file->[0]->{'MySQL'}->{'password'};
	my $db_name = "Faire_pipeline_results";
	my $db_user = "Faire";
	my $db_pass = "faire-seq";

	my $dsn = "dbi:mysql::$host";
	my $dbh = DBI->connect($dsn, $user, $pass) or die "Unable to connect: $DBI::errstr\n";

	$dbh->do("CREATE DATABASE $db_name") or die "Cannot create database: $DBI::errstr\n";
	$dbh->do("CREATE USER $db_user\@$host") or die "Cannot create user: $DBI::errstr\n";
	$dbh->do("GRANT ALL ON $db_name.* TO '$db_user'\@'$host' IDENTIFIED BY '$db_pass'");

	$dbh->disconnect();
	

	!system "mysql -u$db_user -p$db_pass $db_name < $dump_schema" or die "Cannot load Database Schema\n";

	return 1;

}

1;