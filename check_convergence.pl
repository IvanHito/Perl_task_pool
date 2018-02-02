#!/usr/bin/perl

use strict;
use Math::Trig;

###############################################################################
##                            Glogal Parameters                              ##
###############################################################################

our $gBasePath = "/home/ivan/git/Perl_task_pool";            ## should be same as in manage_pool.pl
##our $gBasePath = "/Users/ivanl/git/Perl_task_pool";        ## MakBook
our $gFldData = "data";                                      ## should be same as in manage_pool.pl
our $gFnPool = "pool";                                  ## should be same as in manage_pool.pl
our $gFnRun = "p_run";                                       ## should be same as in manage_pool.pl
our $gFnDone = "p_done";                                     ## should be same as in manage_pool.pl
our $gFnTmpPool = ".tmp_pool2";                              ## should be different!
our $gFnTmpRun = ".tmp_run2";                                ## should be different!
our $gFnTmpDone = ".tmp_done2";                              ## should be different!
##our $gFldDummy = "task_dummy_gre";
our $gFldDummy = "task_dummy";
our $gFldAllData = "all_data_files";

##our $TST_WAIT = "waiting";

###############################################################################
##                               The Program                                 ##
###############################################################################

my %hTSt;

print "\n";
my %hFilesFound = check_files();
print "Found : ";
prt_hash(\%hFilesFound);
print "\n";

print "\n";

if (1>0) {
	my @arrTn = (1953..2056);
	##my @arrTn = (1853..1857);
	my @arrBad1MD = check_firstStep_convergence(\@arrTn);
	print "Convergence check finished. Inappropriate convergence: \n";
	foreach my $tnBad (@arrBad1MD) {
		print "  $tnBad ";
		$hTSt{$tnBad} = "fail";
	}
	print "\n";
}

if (-1>0){
	#%hTSt = (1853=>"wait", 1856=>"wait", 1860=>"wait", 1866=>"wait", 1878=>"wait", 
	#1881=>"wait", 1883=>"wait", 1886=>"wait", 1906=>"wait", 1911=>"wait", 1916=>"wait", 
	#1917=>"wait", 1918=>"wait", 1921=>"wait", 1922=>"wait", 1923=>"wait", 1929=>"wait" );
	print "Change status: \n";
	#%hTSt = (1933=>"done", 1934=>"done", 1935=>"done", 1939=>"done", 1941=>"done");
	prt_hash(\%hTSt);
	print "\n";
	change_tasks_status(\%hTSt);
}
print "Finished \n";

print "\n";

###############################################################################
##                            sub declarations                               ##
###############################################################################

#sub pattern{
#  ##   <<<   Input parameters   >>>   ##
#  ##   <<<   ----------------   >>>   ##
#}

sub prt_hash{
  ##   <<<   Input parameters   >>>   ##
	my %ch = %{$_[0]};
  ##   <<<   ----------------   >>>   ##
	keys %ch; # reset the internal iterator so a prior each() doesn't affect the loop
	print "[";
	while(my($k, $v) = each %ch) {
		print " $k => $v ";
	}
	keys %ch;
	print "]";
}

## %hTNewSt = ( taskN => newState )
sub change_tasks_status {
  ##   <<<   Input parameters   >>>   ##
	my %hTNewSt = %{$_[0]};
  ##   <<<   ----------------   >>>   ##
	my $cl; my $clCp;
	my @arrLine; my @arrOutL;
	my $tState;
	my $tN;
	my $tPath;
	my $nParams;
	my %hFSt;
	my $fn;
	my @args;
	my $ftmp;
	my $poolChanged;
	my %hSTATUSfns = ("DONE"=>0, "WAIT"=>0, "SKIP"=>0, "RUN"=>0, "FAIL"=>0);
	open(FPOOL, "<$::gBasePath/$::gFnPool") or die "Could not open /$::gFnPool : $! \n\n";
	open($ftmp, ">$::gBasePath/$::gFnTmpPool") or die "Could not open /$::gFnTmpPool : $! \n\n";
	while($cl = <FPOOL>){
		if ($cl =~ /#/) {print $ftmp $cl; next;}
		$clCp = $cl;
	  chomp $cl; $cl =~ s/^\s+//;    ## get rid of space in the beginning and new line in the end
	  @arrLine = split /\s+/, $cl;
	  $nParams = @arrLine;
	  if ($nParams != 7){
			print "WARNING!!! Wrong pool string length: $nParams \n";
			print $ftmp $clCp;
			next;
		}
		$tState = $arrLine[4];
		$tN = $arrLine[0];
		keys %hTNewSt; # reset the internal iterator so a prior each() doesn't affect the loop
		$poolChanged = 0;
		while(my($taskNum, $taskNewState) = each %hTNewSt) {
			if ($taskNum == $tN) {
				if ($tState eq $taskNewState){ print "WARNING!!! New state of task $tN is same as old : $tState \n"; }
				elsif ($tState eq "run"){ print "WARNING!!! Task $tN is on the run! We will not touch it. \n"; }
				else{
					$tPath = $arrLine[3];
					%hFSt = check_arbir_files($tPath, \%hSTATUSfns);
					## Delete found STATE files
					foreach $fn (keys %hFSt) {
						if ($hFSt{$fn} == 1){
							$hFSt{$fn} = 0;
							@args = ("rm", "$tPath/$fn");
							system(@args) == 0 or die "system @args failed: $?\n\n";
							print "From $tPath removed : $fn\n";
						}
					}
					if ($taskNewState eq "wait"){
						@args = ("touch", "$tPath/WAIT");
						system(@args) == 0 or die "system @args failed: $?\n\n";
						print "In $tPath Created : WAIT\n";
						$arrLine[4] = "waiting";
						write_pool_str($ftmp, @arrLine);
						$poolChanged = 1;
						print "\n";
					}
					elsif ($taskNewState eq "skip"){
						@args = ("touch", "$tPath/SKIP");
						system(@args) == 0 or die "system @args failed: $?\n\n";
						print "In $tPath Created : SKIP\n";
						$arrLine[4] = "skip";
						write_pool_str($ftmp, @arrLine);
						$poolChanged = 1;
						print "\n";
					}
					elsif ($taskNewState eq "done"){
						@args = ("touch", "$tPath/DONE");
						system(@args) == 0 or die "system @args failed: $?\n\n";
						print "In $tPath Created : DONE\n";
						$arrLine[4] = "done";
						write_pool_str($ftmp, @arrLine);
						$poolChanged = 1;
						print "\n";
					}
					elsif ($taskNewState eq "fail"){
						@args = ("rm", "-r", "$tPath/result");
						system(@args) == 0 or die "system @args failed: $?\n\n";
						print "From $tPath removed : result/\n";
						@args = ("rm", "$tPath/VASP/CHGCAR", "$tPath/VASP/CHG", "$tPath/VASP/WAVECAR");
						system(@args) == 0 or die "system @args failed: $?\n\n";
						print "From $tPath removed : VASP/CHGCAR, VASP/CHG, VASP/WAVECAR\n";
						@args = ("touch", "$tPath/FAIL");
						system(@args) == 0 or die "system @args failed: $?\n\n";
						print "In $tPath Created : FAIL\n";
						$arrLine[4] = "!!!fail!!!";
						write_pool_str($ftmp, @arrLine);
						$poolChanged = 1;
						print "\n";
					}
					else {
						print "WARNING!!! Unknown new task status : $taskNewState\n";
					}
				}
				print "\n";
				last;
			}
		}
		if ($poolChanged == 0) {print $ftmp $clCp;}
	}
	close(FPOOL);
	close($ftmp);

	print "Copy $::gFnTmpPool > $::gFnPool\n";
	## copy from pool tmp to pool
	open(FPOOL, ">$::gBasePath/$::gFnPool") or die "Could not open $::gFnPool: $!\n\n";
	open($ftmp, "<$::gBasePath/$::gFnTmpPool") or die "Could not open $::gFnTmpPool: $!\n\n";
	while($cl = <$ftmp>){print FPOOL $cl;}
	close $ftmp;
	close FPOOL;

} ## end change_tasks_status

sub check_firstStep_convergence {
  ##   <<<   Input parameters   >>>   ##
	my @arrNTs = @{$_[0]};
  ##   <<<   ----------------   >>>   ##
	my $cl; my $clCp;
	my @arrLine; my @arrOutL;
	my $tState;
	my $pInUse = 0;
	my $nParams;
  my $nChecked = 0;
	my $tPath;
	my $dir;
	my @files;
	my $foundOut;
	my $dE;
	my @arrBad = ();
	open(FPOOL, "<$::gBasePath/$::gFnPool") or die "Could not open /$::gFnPool : $! \n\n";
	while($cl = <FPOOL>){
  	if ($cl =~ /#/) { next;}
		chomp $cl; $cl =~ s/^\s+//;    ## get rid of space in the beginning and new line in the end
  	@arrLine = split /\s+/, $cl;
	  $nParams = @arrLine;
  	if ($nParams != 7){print "WARNING!!! Wrong pool string length: <$cl> \n";next;}
		foreach my $tn (@arrNTs){ if ($tn == $arrLine[0]) {
		  $tState = $arrLine[4];
			if (($tState eq "run") or ($tState eq "waiting") or ($tState eq "skip")) {
				print "Watch the state of task $tn! It is <$tState> \n";
			}
			$nChecked++;
			$tPath = $arrLine[3];
	    opendir $dir, "$::gBasePath/$tPath/VASP" or die "Cannot open directory $::gBasePath/$tPath/VASP : $! \n\n";
	    my @files = readdir $dir;
	    closedir $dir;
			$foundOut = 0;
			foreach my $cfile (@files){if ($cfile eq "out_1"){$foundOut = 1; last}}
			if ($foundOut == 0) {
				print "ERROR For task $tn : out_1 is not presented in $::gBasePath/$tPath/VASP \n\n";
				last;
			}
			open(FOUT, "<$::gBasePath/$tPath/VASP/out_1");
			print "In out_1 of task $tn: \n";
			$dE = 999.666;
			while ($cl = <FOUT>){
				if ($cl =~ /Error/){ next; }
				if ($cl =~ /DAV:/){
					chomp $cl; $cl =~ s/^\s+//;
					@arrOutL = split /\s+/, $cl;
					##print "  $arrOutL[0] $arrOutL[1]  $arrOutL[3]\n";
					$dE = $arrOutL[3];
				}
			}
			close(FOUT);
			print "Final dE = $dE";
			if (abs($dE) > 1.0E-6){
				print "  BAD convergence\n";
				push @arrBad, $tn;
			}
			else {print "  good convergence\n"}
			print "\n";
			last;
		}}
		if ($nChecked == @arrNTs){last;}
	}
	close(FPOOL);
	return @arrBad;
}

sub check_files {
	my $foundPool = 0;
	my $foundData = 0;
	my $foundRun = 0;
	my $foundDone = 0;
	my $dir;
	my $cfile;
	opendir $dir, $::gBasePath or die "Cannot open directory $::gBasePath: $! \n\n";
	my @files = readdir $dir;
	closedir $dir;
	foreach $cfile (@files){
  	##print "$cfile \n";
	  if ($cfile eq $::gFldData){$foundData = 1}
  	if ($cfile eq $::gFnPool){$foundPool = 1}
	  if ($cfile eq $::gFnRun){$foundRun = 1}
  	if ($cfile eq $::gFnDone){$foundDone = 1}
	}
	if ($foundPool == 0){die "No POOL FILE! $::gFnPool is not in $::gBasePath \n\n";}
	if ($foundData == 0){die "No DATA FOLDER! $::gFldData is not in $::gBasePath \n\n";}
	my %hFilesFound = (
			"pool" => $foundPool,
			"data" => $foundData,
			"run"  => $foundRun,
			"done" => $foundDone,
	);
	return %hFilesFound;
}

sub check_arbir_files {
	##   <<<   Input parameters   >>>   ##
	my $dirPath      = $_[0];       ## dir path
	my %hFilesStatus = %{$_[1]};    ##
  ##   <<<   ----------------   >>>   ##
	my $cfile;
	my $dir;
	opendir $dir, $dirPath or die "Cannot open directory $dirPath: $! \n\n";
	my @files = readdir $dir;
	closedir $dir;
	my($fn, $fs);
	foreach $fn (keys %hFilesStatus) { $hFilesStatus{$fn} = 0;}
	foreach $cfile (@files){
  	##print "$cfile \n";
		##keys %hFilesStatus; # reset the internal iterator so a prior each() doesn't affect the loop
		foreach $fn (keys %hFilesStatus) {
			if ($fn eq $cfile){ $hFilesStatus{$fn} = 1; last; }
		}
	}
	return %hFilesStatus
}

sub write_pool_str{
  my $fh = $_[0];
  my $pStr = sprintf("  %05d     %25s  %2d  %20s  %8s  %8d  %8d\n",
    $_[1],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7]);
  print $fh $pStr;
}



















##
