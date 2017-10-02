#!/usr/bin/perl

use strict;
use Math::Trig;

###############################################################################
##                            Glogal Parameters                              ##
###############################################################################

our $gSubLevel = 0;
our $gFnLog = "log";
our $gFLog;
our $gDbgStr;
our $gBasePath = "/home/ivan/VASP_data/Perl_task_pool";      ## should be same as in manage_pool.pl
our $gFldData = "data";                                      ## should be same as in manage_pool.pl
our $gFnPool = "pool";                                       ## should be same as in manage_pool.pl
our $gFldDummy = "task_dummy_gre";

###############################################################################
##                               The Program                                 ##
###############################################################################
my $datestring = localtime();
$::gDbgStr = "### -------- -------- $datestring\n";
print $::gDbgStr;
write_log_line($::gDbgStr, "No Arrow");

my $rho = 1.412;
my $nCellsX = 8;
my $nCellsY = 12;
my $lStr = 1.10;
my $fnPscr;
my $tInfo;

print "Ok, let's start \n";
print "\n";

my $tPath;

for ($lStr = 0.98; $lStr <= 1.05; $lStr+=0.01)
{
  $tInfo = sprintf("gre_armCh_ribbon_str%4.2f",$lStr);
  $tPath = new_task($tInfo,56,1);
  $fnPscr = "$::gBasePath/$::gFldData/$tPath/VASP/POSCAR";
  pscr_gre_stripe_armCh($fnPscr,$rho,$nCellsX,$nCellsY,$lStr);
  print "$tInfo\n";
}

print "\n";


###############################################################################
##                            sub declarations                               ##
###############################################################################

#sub pattern{
#  ##   <<<   Input parameters   >>>   ##
#  ##   <<<   ----------------   >>>   ##
#}

sub test_1{
  ##   <<<   Input parameters   >>>   ##
  ##   <<<   ----------------   >>>   ##
  my $datestring = localtime();
  $::gDbgStr = "### -------- -------- $datestring\n";
  print $::gDbgStr;
  write_log_line($::gDbgStr, "No Arrow");
  print "Ok, let's start \n";
  print "\n";
  my $tPath;
  $tPath = new_task("new_task",1,1);
  create_INCAR("$::gFldData/$tPath/VASP/INCAR","Test system");
  print "\n";
}

sub write_pool_str{
  my $fh = $_[0];
  my $pStr = sprintf("  %05d     %15s  %2d  %20s  %8s  %8d  %8d\n",
    $_[1],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7]);
  print $fh $pStr;
}

sub new_task_dummy_files{
  ##   <<<   Input parameters   >>>   ##
  my $taskPath = $_[0];      ## folder of the tsk inside data folder
  ##   <<<   ----------------   >>>   ##
  my $destDir = "$::gBasePath/$::gFldData/$taskPath";
  $::gDbgStr = "Copying task_dummy to $destDir \n";
  print $::gDbgStr;
  write_log_line($::gDbgStr, "No Arrow");
  my @filesToCp = glob "$::gBasePath/$::$gFldDummy/*";
  my @args;
  @args = ("touch", "$destDir/WAIT");
  system(@args) == 0 or die "system @args failed: $?";
  foreach my $file (@filesToCp){
    @args = ("cp", "-R", "$file", "$destDir");
    system(@args) == 0 or die "system @args failed: $?";
  }
}

sub new_task{
  ##   <<<   Input parameters   >>>   ##
  my $tName = $_[0];       ## task name
  my $nc = $_[1];          ## number of cores required
  my $nDtRes = $_[2];      ## number of output results
  ##   <<<   ----------------   >>>   ##
  my $foundPool = 0;
  my $foundData = 0;
  my $lastTaskN = 0;
  my $lastDataN = 0;
  opendir my $dir, $::gBasePath or die "Cannot open directory: $!";
  my @files = readdir $dir;
  closedir $dir;
  foreach my $cfile (@files){
    ##print "$cfile \n";
    if ($cfile eq $::gFldData){$foundData = 1}
    if ($cfile eq $::gFnPool){$foundPool = 1}
  }
  if (($foundData == 1) && ($foundPool == 0)){
    $::gDbgStr = "Data folder is here but no Pool file! \n Now EXIT\n\n";
    print $::gDbgStr;
    write_log_line($::gDbgStr, "No Arrow");
    die;
  }
  if ($foundPool == 0){
    $::gDbgStr = "Creating pool file $::gBasePath/$::gFnPool \n";
    print $::gDbgStr;
    write_log_line($::gDbgStr, "No Arrow");
    open(FPOOL, ">>$::gBasePath/$::gFnPool") or die "Could not open $::gFnPool: $!";
    print FPOOL ("## \$gBasePath = $::gBasePath\n## All paths are in Base Path. \n");
    print FPOOL ("## nTask              name     nc          path           status     dNBeg     sNEnd\n");
    close FPOOL;
  }
  if ($foundData == 0){
    $::gDbgStr = "Creating data folder $::gBasePath/$::gFldData \n";
    print $::gDbgStr;
    write_log_line($::gDbgStr, "No Arrow");
    my @args = ("mkdir", "$::gBasePath/$::gFldData");
    system(@args) == 0 or die "system @args failed: $?";
  }
  ## === check existing tasks
  if ($foundPool == 1){
    my $cl;
    my @arrLine;
    open(FPOOL, "<$::gBasePath/$::gFnPool") or die "Could not open $::gFnPool: $!";
    while($cl = <FPOOL>) {
      if ($cl =~ /#/) {next;}
      #print $cl;
      chomp $cl; $cl =~ s/^\s+//;    ## get rid of space in the beginning and new line in the end
      @arrLine = split /\s+/, $cl;
      $lastTaskN = $arrLine[0];
      $lastDataN = $arrLine[6];
    }
    close FPOOL;
  }
  ## === add new task finaly!
  $lastTaskN += 1;
  my $dBeg = $lastDataN+1;
  my $dEnd = $dBeg+$nDtRes-1;
  my $taskPath = sprintf("task_%05d",$lastTaskN);
  my $fPool;
  open($fPool, ">>$::gBasePath/$::gFnPool") or die "Could not open $::gFnPool: $!";
  write_pool_str($fPool, $lastTaskN, $tName, $nc, "$::gFldData/$taskPath","waiting",$dBeg,$dEnd);
  close $fPool;
  ## === creating task folder
  my @args = ("mkdir", "$::gBasePath/$::gFldData/$taskPath");
  system(@args) == 0 or die "system @args failed: $?";
  new_task_dummy_files($taskPath);
  return $taskPath;
}

sub write_log_line{
  ##   <<<   Input parameters   >>>   ##
  my $toPrint = $_[0];
  my $ifArr = "Ok"; if (@_ == 2) {$ifArr = $_[1];}
  ##   <<<   ----------------   >>>   ##
  my $logtab = "   ";
  my $logarr = ($ifArr eq "Ok") ? " --> " : "";
  my $i;
  open($::gFLog, ">>$::gFnLog") or die "Could not open $::gFnLog: $!";
  for ($i = 0; $i < $gSubLevel; $i++){ print $::gFLog ($logtab); }
  print $::gFLog ("add_task.pl > ",$logarr, $toPrint);
  close $::gFLog;
}

###############################################################################
##                             VASP subroutines                              ##
###############################################################################

sub create_INCAR{
  ## Default parameter values
  my $fn = "INCAR";
  my $sn = "no_specific_name";
  my $ibrion = 2;
  my $isif = 2;
  my $nsw = 500;
  my $encut = 410;
  my $npar = 12;
  my $ismear = 0;
  my $sigma = 0.1;
  ##   <<<   Input parameters   >>>   ##
  my $npars = @_;
  if ($npars > 0){ $fn      = $_[0]; }         ## INCAR file name
  if ($npars > 1){ $sn      = $_[1]; }         ## SYSTEM name
  if ($npars > 2){ $ibrion  = $_[2]; }
  if ($npars > 3){ $isif    = $_[3]; }
  if ($npars > 4){ $nsw     = $_[4]; }
  if ($npars > 5){ $encut   = $_[5]; }
  if ($npars > 6){ $npar    = $_[6]; }
  if ($npars > 7){ $ismear  = $_[7]; }
  if ($npars > 8){ $sigma   = $_[8]; }
  ##   <<<   ----------------   >>>   ##
  print "File ",$fn, " will be created.\n";
  open(my $ff, ">", $fn) or die "Can't open $fn: $!";
  print $ff "SYSTEM = $sn\n";
  print $ff "\n";
  print $ff "PREC = HIGH\n";
  print $ff "\n";
  print $ff "# # Smearing\n";
  print $ff "ISMEAR = $ismear\n";
  print $ff "SIGMA  = $sigma\n";
  print $ff "\n";
  print $ff "# # Algorithm\n";
  print $ff "IBRION = $ibrion\n";
  print $ff "NFREE  = 16\n";
  print $ff "ISIF   = $isif\n";
  print $ff "POTIM  = 0.4\n";
  print $ff "EDIFF  = 1.0E-7\n";
  print $ff "EDIFFG = -1.0E-3\n";
  print $ff "NSW    = $nsw\n";
  print $ff "LAECHG = .TRUE.\n";
  print $ff "\n";
  print $ff "# # Parallel\n";
  print $ff "LPLANE = .FALSE.\n";
  print $ff "NPAR   = $npar\n";
  print $ff "NSIM   = 1\n";
  print $ff "\n";
  print $ff "# # Others\n";
  print $ff "AMIX  = 0.2\n";
  print $ff "BMIX  = 0.0001\n";
  print $ff "ENCUT = $encut\n";
  print $ff "\n";
  close($ff) or die "$ff: $!";
} ## end create_INCAR

## --------------------------------------------------- ##
##                      structures                     ##
## --------------------------------------------------- ##

sub pscr_hcp_Mg{
  ##   <<<   Input parameters   >>>   ##
  my $hcpA = $_[0];         ## hcp parameter a
  my $nCells1 = $_[1];      ## number of cells in direction of vector 1
  my $nCells2 = $_[2];      ## number of cells in direction of vector 2
  my $nCells3 = $_[3];      ## number of cells in direction of vector 3
  my $fnPscr = $_[4];       ## file name for POSCAR file
  ##   <<<   ----------------   >>>   ##
  my $hcpC = $hcpA*1.633;   ## hcp parameter c
  my $hcpG = 120*pi/180;    ## hcp angle (initially in degrees)
  my $tp = "Mg";            ## atom type
  pscr_hcp($hcpA,$hcpC,$hcpG,$nCells1,$nCells2,$nCells3,$tp,$fnPscr);
}


sub pscr_hcp{
  ##   <<<   Input parameters   >>>   ##
  my $hcpA = $_[0];         ## hcp parameter a
  my $hcpC = $_[1];         ## hcp parameter c
  my $hcpG = $_[2]*pi/180;  ## hcp angle (initially in degrees)
  my $nCells1 = $_[3];      ## number of cells in direction of vector 1
  my $nCells2 = $_[4];      ## number of cells in direction of vector 2
  my $nCells3 = $_[5];      ## number of cells in direction of vector 3
  my $tp = $_[6];           ## atom type
  my $fnPscr = $_[7];       ## file name for POSCAR file
  ##   <<<   ----------------   >>>   ##
  my $i = 0; my $j = 0; my $k = 0; my $ii = 0;
  my $fixStyle = " T  T  T";
  my $cg = cos($hcpG);
  my $sg = sin($hcpG);
  my @baseVectors = ([$hcpA,0,0], [$hcpA*$cg,$hcpA*$sg,0], [0,0,$hcpC]);
  my @bvAtoms = ([0,0,0], [1/3, 2/3, 0.5]);
  my @cvA = ([0,0,0], [0,0,0]);
  my $nAtoms = $nCells1*$nCells2*$nCells3*2;
  open(FPSCR, ">$fnPscr") or die "Could not open $fnPscr: $!";
  print FPSCR "Pure hcp structure \n";
  print FPSCR " 1.0000000000000000 \n";
  my @vNcs = ($nCells1,$nCells2,$nCells3);
  for ($i = 0; $i<3; $i++){
    printf FPSCR '   ';
    for ($j = 0; $j<3; $j++){printf FPSCR '%20.15f  ', $baseVectors[$i][$j]*$vNcs[$i];}
    print FPSCR "\n";
  }
  print FPSCR '   ', $tp, "\n";
  print FPSCR '   ', $nAtoms, "\n";
  print FPSCR "Selective dynamics \n";
  print FPSCR "Direct \n";
  for ($i = 0; $i<$nCells1; $i++){
    for ($j = 0; $j<$nCells2; $j++){
      for ($k = 0; $k<$nCells3; $k++){
        for ($ii = 0; $ii<2; $ii++){
          $cvA[$ii][0] = $bvAtoms[$ii][0]/$nCells1 + $i/$nCells1;
          $cvA[$ii][1] = $bvAtoms[$ii][1]/$nCells2 + $j/$nCells2;
          $cvA[$ii][2] = $bvAtoms[$ii][2]/$nCells3 + $k/$nCells3;
          printf FPSCR ' %15.12f %15.12f %15.12f %s', $cvA[$ii][0], $cvA[$ii][1], $cvA[$ii][2], $fixStyle;
          print FPSCR "\n";
        }
      }
    }
  }
  close FPSCR;
}## end

sub pscr_gre_stripe_armCh{
  ##   <<<   Input parameters   >>>   ##
  my $fnPscr = $_[0];       ## file name for POSCAR file
  my $rho = $_[1];
  my $nCells1 = $_[2];
  my $nCells2 = $_[3];
  my $lStr = 1; if (@_ >= 5) {$lStr = $_[4]}
  ##   <<<   ----------------   >>>   ##
  my $ax = $rho*sqrt(3.0)/2.0;
  my $ay = $lStr*$rho*3.0/2.0;
  ##print $ay, "\n";
  my $i = 0; my $j = 0; my $ix = 0; my $iy = 0; my $ii = 0;
  my $fixStyle = " T  T  T";
  my $padd = 15.0;
  my @baseVectors = ([$nCells1*$ax+$padd,0,0], [0,$nCells2*$ay+$padd,0], [0,0,15]);
  my $paddX = $padd/$baseVectors[0][0]*0.5;
  my $paddY = $padd/$baseVectors[1][1]*0.5;
  my $nAtoms = $nCells1*$nCells2;
  my @cvA = ([0,0,0], [0,0,0]);
  ##print "$nAtoms \n";
  open(FPSCR, ">$fnPscr") or die "Could not open $fnPscr: $!";
  print FPSCR "Graphene 2D square structure \n";
  print FPSCR " 1.0000000000000000 \n";
  for ($i = 0; $i<3; $i++){
    printf FPSCR '   ';
    for ($j = 0; $j<3; $j++){printf FPSCR '%20.15f  ', $baseVectors[$i][$j];}
    print FPSCR "\n";
  }
  print FPSCR "  C  \n";
  print FPSCR '   ', $nAtoms, "\n";
  print FPSCR "Selective dynamics \n";
  print FPSCR "Direct \n";
  $ii = 0;
  for ($iy = 0; $iy<$nCells2; $iy++){
    $fixStyle = " T  T  T";
    if (($iy < 2) or ($iy >= $nCells2-2)){$fixStyle = " F  F  F";}
    for ($ix = 0; $ix<$nCells1; $ix++){
      $ii++;
      $cvA[$ii][0] = (($ix+1)*$ax - 0.5*$ax)/$baseVectors[0][0] + $paddX;
      $cvA[$ii][1] = ((-1)**($ix+$iy)*$rho/4 + ($iy+1)*$ay - 0.5*$ay)/$baseVectors[1][1] + $paddY;
      $cvA[$ii][2] = 0.5;
      printf FPSCR ' %15.12f %15.12f %15.12f %s', $cvA[$ii][0], $cvA[$ii][1], $cvA[$ii][2], $fixStyle;
      print FPSCR "\n";
    }
  }
  ##print "$ii \n";
} ## end pscr_gre_stripe_armCh































##
