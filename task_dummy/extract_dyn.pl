#!/Users/ivanl/perl5/perlbrew/perls/perl-5.16.0/bin/perl

use strict;
use autodie;
use open ":utf8";

###############################################################################
##                               The Program                                 ##
###############################################################################
my $dataN = 1;
my $reqIter = "last";
my $num_args = $#ARGV + 1;
if ($num_args >= 1) {
    $dataN = $ARGV[0];
}
if ($num_args >= 2) {
    $reqIter = $ARGV[1];
}

## test extraction
print "\n";
extract_crystal_from_dyn("VASP/OUTCAR", $reqIter, "EXTRACTED_".$dataN, \&get_type_allMg);
print "\n";
###############################################################################
##                            sub declarations                               ##
###############################################################################

sub extract_crystal_from_dyn{
  if (@_ < 4) {die "ERROR! Not enough parameters for extract_crystal_from_dyn";}
  ##   <<<   Input parameters   >>>   ##
  my $outFn = $_[0];         ## Name of OUTCAR file
  my $iterN = $_[1];         ## number of iteration (it may be "last" to extract last iteration)
  my $fnOut = $_[2];         ## file name to write results of extraction
  my $funcGetType = $_[3];   ## function for atom Type recognition
  ##   <<<   ----------------   >>>   ##
  my $itFound = 0;
  my $cl = "";
  my $nAtoms = 0;
  my @atomsInfo = ();
  my @arrLine = ();
  my $i = 0;
  my $j = 0;
  my @baseVectors = ([0,0,0], [0,0,0], [0,0,0]);
  my $E;
  my $unitsE;
  my $fo;
  my $type;
  my $reopen = 0;
  #my $dbgStr;

  $::gSubLevel += 1;

  open($fo, $outFn) or die "Could not open $outFn: $!";

  print "Start extraction \n";

  while($cl = <$fo>) {
    #print $cl;
    #last if $. == 2700;
    if ($cl =~ /Iteration/) {
      $cl =~ /(\d+).+\(\ +(\d+)/; ## $1 - Itertation Number, $2 - scf step
      $itFound = $1;
      if (($itFound == $iterN) and ($iterN ne "last")) {last;}
    }
  }
  if ($itFound == 0){die "ERROR!!! THE FILE CONTAINS NOTHING!\n\n";}
  if ($iterN eq "last"){
    close $fo;
    $iterN = $itFound;
    $reopen = 1;
    print "Last iteration is $itFound \n";
  } else {
    if ($itFound == $iterN) {
      print "$. lines were read\n";
      print "Found $itFound \n";
    } else {
      print "Warning! Iteration $iterN is not found! \n";
      print "Instead, iteration $itFound (which is the last one) will be read. \n";
      close $fo;
      $iterN = $itFound;
      $reopen = 1;
    }
  }

  if ($reopen > 0){
    open($fo, $outFn) or die "Could not open $outFn: $!";
    while($cl = <$fo>) {
      if ($cl =~ /Iteration/) {
        $cl =~ /(\d+).+\(\ +(\d+)/; ## $1 - Itertation Number, $2 - scf step
        $itFound = $1;
        if ($itFound == $iterN) {last;}
      }
    }
    print "$. lines were read\n";
    if ($itFound) {print "Found $itFound \n";}
  }

  ## Find vectors inside current iteration
  while($cl = <$fo>) {if ($cl =~ /lattice vectors/) {last;}}
  for ($i=0; $i<3; $i++){
    $cl = <$fo>;
    chomp $cl; $cl =~ s/^\s+//;
    @arrLine = split /\s+/, $cl;
    for ($j=0; $j<3; $j++){ $baseVectors[$i][$j] = $arrLine[$j] }
  }

  ## Find position inside current iteration
  while($cl = <$fo>) {if ($cl =~ /POSITION/) {last;}}
  $cl = <$fo>;
  $i = 0;
  while($cl = <$fo>) {
    if ($cl =~ /-----------------/) {last;}
    $nAtoms += 1;
    chomp $cl;
    $cl =~ s/^\s+//;
    @arrLine = split /\s+/, $cl;
    #print $nAtoms, "  ", join(", ",  @arrLine), "\n";
    foreach my $num (@arrLine){ push @{$atomsInfo[$i]}, $num; }
    $i++;
  }

  ## Find Energy for current iteration
  while($cl = <$fo>) {if ($cl =~ /FREE ENERGIE/) {last;}}
  $cl = <$fo>;
  $cl = <$fo>;
  $cl =~ /TOTEN\s+=\s+([-+\d\.]+)\s+([a-zA-Z]+)/;
  $E = $1;
  $unitsE = $2;

  close $fo;

  print "\n finished atoms POSITION reading\n\nWriting it to $fnOut\n";

  open(FOUT, ">$fnOut") or die "Could not open $fnOut: $!";
  print FOUT "# total energy = $E $unitsE\n";
  print FOUT "CRYSTAL\n";
  print FOUT "PRIMVEC\n";
  for ($i = 0; $i<3; $i++){
    for ($j = 0; $j<3; $j++){printf FOUT '   %9.4f ', $baseVectors[$i][$j];}
    print FOUT "\n";
  }
  print FOUT "PRIMCOORD\n";
  print FOUT "$nAtoms 1 \n";
  for ($i = 0; $i<$nAtoms; $i++){
    $type = $funcGetType->($i+1);
    print FOUT "  $type    ";
    for ($j = 0; $j<3; $j++){printf FOUT '%9.4f  ', $atomsInfo[$i][$j];}
    printf FOUT '  ';
    for ($j = 3; $j<6; $j++){printf FOUT '%9.4f  ', $atomsInfo[$i][$j];}
    print FOUT "\n";
  }
  close FOUT;

  $::gSubLevel -= 1;
} ## end extract_crystal_from_dyn

sub get_type_cz96{
  ##   <<<   Input parameters   >>>   ##
  my $cn = $_[0];  ## current atom number starting from 1
  ##   <<<   ----------------   >>>   ##
  my $tp = "Unknown";
  $tp = ($cn <= 48) ? "Cu" : "Zr";
  return $tp;
}

sub get_type_allMg{
  ##   <<<   Input parameters   >>>   ##
  my $cn = $_[0];  ## current atom number starting from 1
  ##   <<<   ----------------   >>>   ##
  my $tp = "Mg";
  return $tp;
}

































##
