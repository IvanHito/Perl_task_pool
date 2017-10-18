#!/usr/bin/perl

use strict;
##use Math::Trig;

###############################################################################
##                            Glogal Parameters                              ##
###############################################################################

our $gSubLevel = 0;
our $gFnLog = "/home/ivan/git/Perl_task_pool/log";
our $gFLog;
our $gDbgStr;

my $basePath = "/home/ivan/git/Perl_task_pool";     ## should be same as in add_task.pl
my $fldData = "data";                                     ## should be same as in add_task.pl
my $fnPool = "pool";                                      ## should be same as in add_task.pl
my $maxProc = 56;         ## maximum processors in use
my $fnTmp = ".tmp_pool";
my $fnTmpR = ".tmp_run";

###############################################################################
##                               The Program                                 ##
###############################################################################
my $fnRun = "p_run";
my $fnDone = "p_done";
my $cl;
my @arrLine;
my $tState;
my $pInUse = 0;
my $foundPool = 0;
my $foundData = 0;
my $foundRun = 0;
my $foundDone = 0;
my $fh; my $fh2;
my $tPath;

my $datestring = localtime();

$::gDbgStr = "### ======= ======= $datestring\n";
print $::gDbgStr;
write_log_line($::gDbgStr, "No Arrow");

print "Manage pool \n";
print "\n";

##---##
## 0 ## initial check
##---##
my $dir;
my $cfile;
opendir $dir, $basePath or die "Cannot open directory $basePath: $!";
my @files = readdir $dir;
closedir $dir;
foreach $cfile (@files){
  ##print "$cfile \n";
  if ($cfile eq $fldData){$foundData = 1}
  if ($cfile eq $fnPool){$foundPool = 1}
  if ($cfile eq $fnRun){$foundRun = 1}
  if ($cfile eq $fnDone){$foundDone = 1}
}
if ($foundPool == 0){die "No POOL FILE! $fnPool is not in $basePath";}
if ($foundData == 0){die "No DATA FOLDER! $fldData is not in $basePath";}
if ($foundRun == 0){
  print "Create Run\n";
  open(FRUN, ">$basePath/$fnRun") or die "Could not open $fnRun: $!";
  print FRUN ("## \$gBasePath = $basePath\n## All paths are in Base Path. \n");
  print FRUN ("## nTask              name     nc          path           status     dNBeg     sNEnd\n");
  close FRUN;
  open(FRUN, ">$basePath/$fnTmpR") or die "Could not open $fnTmpR: $!";
  print FRUN ("## \$gBasePath = $basePath\n## All paths are in Base Path. \n");
  print FRUN ("## nTask              name     nc          path           status     dNBeg     sNEnd\n");
  close FRUN;
}
if ($foundDone == 0){
  open(FDONE, ">$basePath/$fnDone") or die "Could not open $fnDone: $!";
  print FDONE ("## \$gBasePath = $basePath\n## All paths are in Base Path. \n");
  print FDONE ("## nTask              name     nc          path           status     dNBeg     sNEnd\n");
  close FDONE;
}

##---##
## 1 ## check run
##---##
my @arrDone = ();
my $clCp;
if ($foundRun == 1){
  print "Found Run\n";
  open(FRUN, "<$basePath/$fnRun") or die "Could not open $fnRun: $!";
  open($fh, ">$basePath/$fnTmpR") or die "Could not open $fnTmpR: $!";
  while($cl = <FRUN>){
    if ($cl =~ /## Processors/) {next;}
    if ($cl =~ /#/) {print $fh $cl; next;}
    $clCp = $cl;
    chomp $cl; $cl =~ s/^\s+//;    ## get rid of space in the beginning and new line in the end
    @arrLine = split /\s+/, $cl;
    $tPath = $arrLine[3];
    opendir $dir, "$basePath/$tPath" or die "Cannot open directory $basePath/$tPath: $!";
    my @files = readdir $dir;
    closedir $dir;
    $foundDone = 0;
    foreach $cfile (@files){if ($cfile eq "DONE"){
      $foundDone = 1;
      push @arrDone, $arrLine[0];
      $arrLine[4] = "done";
      open($fh2, ">>$basePath/$fnDone") or die "Could not open $fnDone: $!";
      write_pool_str($fh2,$arrLine[0],$arrLine[1],$arrLine[2],$arrLine[3],$arrLine[4],$arrLine[5],$arrLine[6]);
      close $fh2;
      $::gDbgStr = "Task " . $arrLine[0] . " is done \n";
      print $::gDbgStr;
      write_log_line($::gDbgStr, "No Arrow");
      last;
    }}
    if ($foundDone == 0){
      print "Task ", $arrLine[0], " goes back to tmp_run\n";
      print $fh $clCp;
      $pInUse += $arrLine[2];
    }
  }
  close $fh;
  close FRUN;
}
if (@arrDone == 0) {
  $::gDbgStr = "Nothing is done yet. \n";
  print $::gDbgStr;
  write_log_line($::gDbgStr, "No Arrow");
}

##---##
## 2 ## put new tasks from pool to run
##---##
if ($pInUse >= $maxProc){
  $::gDbgStr = "\n                 All procs are buisy. \n\n";
  print $::gDbgStr;
  write_log_line($::gDbgStr, "No Arrow");
  if ($pInUse > $maxProc){
    $::gDbgStr = "\n                 HEY! WTF?! We are using more procs than we want to! \n\n";
    print $::gDbgStr;
    write_log_line($::gDbgStr, "No Arrow");
  }
  if (@arrDone == 0) {print "Loop file will not be affected. \n"; exit;}
}
my $npDemand;
my $ftmp;
my @args;
my $nFinishedRun = 0;
my $nParams = 0;
open(FPOOL, "<$basePath/$fnPool") or die "Could not open $fnPool: $!";
open($ftmp, ">$basePath/$fnTmp") or die "Could not open $fnTmp: $!";
while($cl = <FPOOL>){
  if ($cl =~ /#/) {print $ftmp $cl; next;}
  if (($pInUse == $maxProc) and ($nFinishedRun == @arrDone)){print $ftmp $cl; next;}
  $clCp = $cl;
  chomp $cl; $cl =~ s/^\s+//;    ## get rid of space in the beginning and new line in the end
  @arrLine = split /\s+/, $cl;
  $nParams = @arrLine;
  if ($nParams != 7){print "WARNING!!! Wrong pool string length: $nParams \n";next;}
  $tState = $arrLine[4];
  if ($tState eq "run"){
    # check done array
    #print "Check in "; foreach my $doneN (@arrDone){print $doneN," ";}print "\n";
    $foundDone = 0;
    foreach my $doneN (@arrDone){if ($doneN == $arrLine[0]){
      $foundDone = 1;
      $nFinishedRun++;
      #print "Found ",$doneN," :\n";
      $arrLine[4] = "done";
      write_pool_str($ftmp,$arrLine[0],$arrLine[1],$arrLine[2],$arrLine[3],$arrLine[4],$arrLine[5],$arrLine[6]);
      last;
    }}
    if ($foundDone == 0){print $ftmp $clCp;}
    next;
  }
  if ($tState eq "done"){print $ftmp $clCp; next;}
  $npDemand = $arrLine[2];
  if ($pInUse + $npDemand > $maxProc){print $ftmp $clCp; next;}
  $tPath = $arrLine[3];
  ## do real things
  #print "Run steward in  $basePath/$tPath\n";
  @args = ("$basePath/start_steward.sh","$basePath/$tPath",$arrLine[0],$arrLine[2]);
  system(@args) == 0 or die "system @args failed: $?";
  $::gDbgStr = "\n --> Set to Run ".$arrLine[1]." ".$arrLine[0]." \n\n";
  print $::gDbgStr;
  write_log_line($::gDbgStr, "No Arrow");
  $arrLine[4] = "run";
  write_pool_str($ftmp,$arrLine[0],$arrLine[1],$arrLine[2],$arrLine[3],$arrLine[4],$arrLine[5],$arrLine[6]);
  ## put it to run! (using tmp_run which will be copied to run after)
  open($fh, ">>$basePath/$fnTmpR") or die "Could not open $fnTmpR: $!";
  write_pool_str($fh,$arrLine[0],$arrLine[1],$arrLine[2],$arrLine[3],$arrLine[4],$arrLine[5],$arrLine[6]);
  close $fh;
  $pInUse += $npDemand;
  #if ($pInUse == $maxProc){last;}
}
close $ftmp;
close FPOOL;

print "Copy .tmp_pool > pool\n";
## copy from pool tmp to pool
open(FPOOL, ">$basePath/$fnPool") or die "Could not open $fnPool: $!";
open($ftmp, "<$basePath/$fnTmp") or die "Could not open $fnTmp: $!";
while($cl = <$ftmp>){print FPOOL $cl;}
close $ftmp;
close FPOOL;

print "Copy .tmp_run > run\n";
## copy tmp_run to run file
open(FRUN, ">$basePath/$fnRun") or die "Could not open $fnRun: $!";
open($fh, "<$basePath/$fnTmpR") or die "Could not open $fnTmpR: $!";
print FRUN "## Processors in use : $pInUse / $maxProc\n";
while($cl = <$fh>){print FRUN $cl;}
close $fh;
close FRUN;

$::gDbgStr = "Processors in use : $pInUse / $maxProc\n";
print $::gDbgStr;
write_log_line($::gDbgStr, "No Arrow");

print "Manage pool finished\n";
print "\n";


###############################################################################
##                            sub declarations                               ##
###############################################################################

#sub pattern{
#  ##   <<<   Input parameters   >>>   ##
#  ##   <<<   ----------------   >>>   ##
#}

sub write_pool_str{
  my $fh = $_[0];
  my $pStr = sprintf("  %05d     %25s  %2d  %20s  %8s  %8d  %8d\n",
    $_[1],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7]);
  print $fh $pStr;
}

sub write_log_line{
  ##   <<<   Input parameters   >>>   ##
  my $toPrint = $_[0];
  #print "From write_log_line: ", $toPrint;
  my $ifArr = "Ok"; if (@_ == 2) {$ifArr = $_[1];}
  ##   <<<   ----------------   >>>   ##
  my $logtab = "   ";
  my $logarr = ($ifArr eq "Ok") ? " --> " : "";
  my $i;
  open($::gFLog, ">>$::gFnLog") or die "Could not open $::gFnLog: $!";
  for ($i = 0; $i < $gSubLevel; $i++){ print $::gFLog ($logtab); }
  print $::gFLog ("manage_pool.pl > ",$logarr, $toPrint);
  close $::gFLog;
}

































##
