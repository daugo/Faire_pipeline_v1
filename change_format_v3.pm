package change_format_v3;
#change_format_v3.pm
use warnings; use strict;
use File::Basename;

#sam2bedCoordinateSort(@ARGV);


sub sam2bedCoordinateSort {
	my $safe_disk_space = shift;
	#print "SAFE DISK: $safe_disk_space\n";
	my $path_in = shift;
	#print $path_in."\n";
	my $path_out = shift; #path where the user/module want the output files to be located
	#chomp($path_in);#,$path_out);
	#die "File $path_in not found" unless (-f $path_in); # check if the input file is a plain text
	my ($filename,$directories,$suffix) = fileparse($path_in,qr/\.\w+(\.gz|[^.]*)/); #save the diferent parts of the path, the regex for the suffix accepts .any and .any.gz
	die "Output entry must be a directory" unless (!$path_out or -d $path_out);
	my $dir_out = dirname($path_out) if $path_out; #save the directory path where data is going to be write
	$dir_out = $directories unless $path_out;
	
	#die "File $path_out not a valid directory" unless ( -d $path_out);
	
	if ($suffix ~~ /sam/) {
		my $in_sam = $path_in;
		my $out_bam = $dir_out.$filename.".bam";
		my $command_sam = "samtools view -bS -q 1 -o $out_bam $in_sam";
		my ($run)  = command($command_sam);
		if ($safe_disk_space) {
			unlink $in_sam or warn "Cannot remove $in_sam:$!";
		}
		if ($run) {
			print "sam2bam Succes!!\n";
			my $out_bam_sort = $dir_out.$filename."_sorted";
			my $command_bamsort = "samtools sort $out_bam $out_bam_sort";
			($run) = command($command_bamsort);
			if ($safe_disk_space) {
				unlink $out_bam or warn "Cannot remove $out_bam:$!";
			}
			my $command_bamindex = "samtools index $out_bam_sort.bam";
			($run) = command($command_bamindex);
			if ($run) {
				print "bam2bamsort Succes!!\n";
				print "bam_sort index Succes!!\n";
				$out_bam_sort = $dir_out.$filename."_sorted.bam";
				my $out_bed = $dir_out.$filename."_sorted.bed";
				my $command_bed = "bamToBed -i $out_bam_sort > $out_bed";
				($run) = command($command_bed);
				print "bamsort2bed Succes!!\n" if $run;
			}
		}
		
	}
}
=c
	elsif ($suffix ~~ /bam/) {
		my $in_bam = $path_in;
		my $out_bam_sort = $dir_out.$filename."_sorted";
		my $command_bamsort = "samtools sort $in_bam $out_bam_sort";
		my ($run) = command($command_bamsort);
		if ($run) {
			print "bam2bamsort Succes!!\n";
			$out_bam_sort = $dir_out.$filename."_sorted.bam";
			my $out_bed = $dir_out.$filename."_sorted.bed";
			my $command_bed = "bamToBed -i $out_bam_sort > $out_bed";
			($run)= command($command_bed);
			print "bamsort2bed Succes!!\n" if $run;
		}
	}
	
	elsif ($suffix ~~ /bed/) {
		print "file is already in bed format \n";
	}
	else { warn "file extension don't supported"};
=cut

sub command {
	my $command = shift;
	print "command: $command\n";
	my ($system_ID) = system("$command");
	unless ($system_ID) { #Not just one becouse BED_TOOLS ouput msj is distinct
		return 1;
	}
	else {
		return 0;
	}
}



1;
