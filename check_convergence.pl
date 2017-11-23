#!/usr/bin/perl

use strict;
use Math::Trig;

###############################################################################
##                            Glogal Parameters                              ##
###############################################################################

our $gBasePath = "/home/ivan/git/Perl_task_pool";            ## should be same as in manage_pool.pl
##our $gBasePath = "/Users/ivanl/git/Perl_task_pool";        ## MakBook
our $gFldData = "data";                                      ## should be same as in manage_pool.pl
our $gFnPool = "pool_copy";                                  ## should be same as in manage_pool.pl
our $gFnRun = "p_run";                                       ## should be same as in manage_pool.pl
our $gFnDone = "p_done";                                     ## should be same as in manage_pool.pl
our $gFnTmpPool = ".tmp_pool2";                              ## should be different!
our $gFnTmpRun = ".tmp_run2";                                ## should be different!
our $gFnTmpDone = ".tmp_done2";                              ## should be different!
##our $gFldDummy = "task_dummy_gre";
our $gFldDummy = "task_dummy";
our $gFldAllData = "all_data_files";

###############################################################################
##                               The Program                                 ##
###############################################################################

print "\n";
my %hFilesFound = check_files();
print "Found : ", %hFilesFound, "\n";

print "\n";

my @arrTn = (1853..1932);
my @arrBad1MD = check_firstStep_convergence(\@arrTn);
print "Convergence check finished. Inappropriate convergence: \n";
foreach my $tnBad (@arrBad1MD) {print "  $tnBad "}
print "\n";

print "\n";

###############################################################################
##                            sub declarations                               ##
###############################################################################

#sub pattern{
#  ##   <<<   Input parameters   >>>   ##
#  ##   <<<   ----------------   >>>   ##
#}

sub change_tasks_status {
  ##   <<<   Input parameters   >>>   ##
	my %arrTNewSt = %{$_[0]};
  ##   <<<   ----------------   >>>   ##
	my $cl; my $clCp;
	my @arrLine; my @arrOutL;
	my $tState;
	open(FPOOL, "<$::gBasePath/$::gFnPool") or die "Could not open /$::gFnPool : $! \n\n";
	open(FTMPP, "<$::gBasePath/$::gFnTmpPool") or die "Could not open /$::gFnTmpPool : $! \n\n";
	while($cl = <FPOOL>){
	}
	close(FPOOL);
	close(FTMPP);


keys %hash; # reset the internal iterator so a prior each() doesn't affect the loop
while(my($k, $v) = each %hash) { ... }
}

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
			while ($cl = <FOUT>){if ($cl =~ /DAV:/){
				chomp $cl; $cl =~ s/^\s+//;
				@arrOutL = split /\s+/, $cl;
				##print "  $arrOutL[0] $arrOutL[1]  $arrOutL[3]\n";
				$dE = $arrOutL[3];
			}}
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

sub write_pool_str{
  my $fh = $_[0];
  my $pStr = sprintf("  %05d     %25s  %2d  %20s  %8s  %8d  %8d\n",
    $_[1],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7]);
  print $fh $pStr;
}



















##

