package change_format_fork;
#change_format_fork.pm
use warnings;
use strict;
use change_format_v3;


sub sam2bed_fork {
	my $safe_disk_space = shift;
	#print "SAFE DISK: $safe_disk_space\n";
	my @files = @_;
	print $files[0]."\n";
	foreach (@files) {
		if ($_ =~ /.sam/) {
			#TODO (better check)
			print "file check succesfull\n";
		}
		else {
			die "file(s) extension(s) incorrect";
		}
	}
	my @childs = ();
	for (my $i = 0; $i < scalar(@files) ; $i++) {
		my $pid = fork();
		if ($pid) { #parent 
			push(@childs, $pid);
		}
		elsif ($pid == 0) {
			#thing that has to be done;
			#print "$files[$i]\n";
			change_format_v3::sam2bedCoordinateSort($safe_disk_space,$files[$i]);
			#sleep(1);
			exit(0);
		}
		else {
			die "couldn’t fork: $!\n";
			}
		#print "BEFORE FOR BRACKET\n";
	}
	#print "AFTER FOR BRACKET\n";
	
	foreach (@childs) {
		waitpid($_, 0);
	}
	return 1;
}

1;
