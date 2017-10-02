#!/usr/bin/perl

package vaspSTRUCT;
use Class::Accessor 'antlers';
use strict;
use Math::Trig;

#### ================================= ####
####             Attributes            ####
#### ================================= ####

has isDefined        => ( is => "ro", isa => "Str" ); ## "OK" means it is defined
has sourceFile       => ( is => "rw", isa => "Str" );
has infoLine         => ( is => "rw", isa => "Str" );
has factor           => ( is => "rw" );
has baseVs           => ( is => "rw" );
has baseAvs          => ( is => "rw" );
has atomTypeNames    => ( is => "rw" );
has atomTypeNums     => ( is => "rw" );
has nAtoms           => ( is => "rw" );
has dynType          => ( is => "rw" );
has aPosType         => ( is => "rw" );
has fixedNums        => ( is => "rw" );
has nCellsArr        => ( is => "rw" );
has ainfD            => ( is => "rw" );
has ainfC            => ( is => "rw" );


#### ================================= ####
####          Inner variables          ####
#### ================================= ####

my $VALUNDEF  = "undefined";
my $VALOK     = "OK";
my $FNPOSCAR  = "POSCAR";
my $FNXSF     = "COORDS.xsf";
my $PRINTDBG  = 0;
my @aRefNames = ( "baseVs",         # array of names which can be recopied
                  "baseAvs",
                  "atomTypeNames",
                  "atomTypeNums",
                  "ainfD",
                  "ainfC" );


#### ================================= ####
####            Subroutines            ####
#### ================================= ####

#sub dummy {
#  my $self = shift;
  ##   <<<   Input parameters   >>>   ##
  ##   <<<   ----------------   >>>   ##
#}

sub my_default_vals {
  my $self = shift;
  $self->set("isDefined",$VALUNDEF);
  ##print "Me: $self. State: $VALUNDEF \n";
  ##print "\n";
  $self->sourceFile($VALUNDEF);
  $self->infoLine($VALUNDEF);
  $self->factor(1.0);
  my @arrArr1 = ([0,0,0], [0,0,0], [0,0,0]);
  $self->baseVs(@arrArr1);
  my @arr1;
  $self->atomTypeNames(@arr1);
  my @arr2;
  $self->atomTypeNums(@arr2);
  $self->nAtoms(0);
  $self->dynType("S");
  $self->aPosType("Direct");
  my @arr0 = (0,0,0);
  $self->nCellsArr(@arr0);
  my @arr3;
  $self->ainfD(@arr3);
  my @arr4;
  $self->ainfC(@arr4);
  my @arr5;
  $self->fixedNums(@arr5);
  my @arr6;
  $self->baseAvs(@arr6);
}

sub recopy_arr {
  my $self = shift;
  my $nm = shift // $VALUNDEF;
  if ( $nm ~~ @aRefNames ){
    ##print "Recopy $nm \n";
    my @newArray;
    map { push(@newArray,$_); } @{$self->SUPER::get($nm)};
    $self->SUPER::set($nm, \@newArray);
  } else { print "WARNING!!! Illigal recopy $nm \n"; }
}

##sub new_

sub test_shift {
  my $self = shift;
}

sub new_struct_hcp_Mg {
  my $self = shift;
  ##   <<<   Input parameters   >>>   ##
  my $hcpA = $_[0];         ## hcp parameter a
  my $nCells1 = $_[3];      ## number of cells in direction of vector 1
  my $nCells2 = $_[4];      ## number of cells in direction of vector 2
  my $nCells3 = $_[5];      ## number of cells in direction of vector 3
  ##   <<<   ----------------   >>>   ##
  my $hcpC = $hcpA*1.633;   ## hcp parameter c
  my $hcpG = 120*pi/180;    ## hcp angle (initially in degrees)
  my $tp = "Mg";            ## atom type
  $self->new_struct_hcp($hcpA,$hcpC,$hcpG,$nCells1,$nCells2,$nCells3,$tp);
}

sub new_struct_hcp {
  my $self = shift;
  ##   <<<   Input parameters   >>>   ##
  my $hcpA = $_[0];         ## hcp parameter a
  my $hcpC = $_[1];         ## hcp parameter c
  my $hcpG = $_[2]*pi/180;  ## hcp angle (initially in degrees)
  my $nCells1 = $_[3];      ## number of cells in direction of vector 1
  my $nCells2 = $_[4];      ## number of cells in direction of vector 2
  my $nCells3 = $_[5];      ## number of cells in direction of vector 3
  my $tp = $_[6];           ## atom type
  ##   <<<   ----------------   >>>   ##
  my $cg = cos($hcpG);
  my $sg = sin($hcpG);
  my @baseVectors = ([$hcpA,0,0], [$hcpA*$cg,$hcpA*$sg,0], [0,0,$hcpC]);
  my @bvAtoms = ([0,0,0], [1/3, 2/3, 0.5]);
  $self->set("isDefined",$VALOK);
  $self->infoLine("Pure $tp hcp structure");
  $self->baseVs(@baseVectors);
  $self->baseAvs(@bvAtoms);
  $self->atomTypeNames([$tp]);
  $self->nCellsArr([$nCells1,$nCells2,$nCells3]);
  $self->struct_base();
}

sub new_struct_fcc {
  my $self = shift;
  ##   <<<   Input parameters   >>>   ##
  my $hcpA = $_[0];         ## hcp parameter a
  my $hcpB = $_[1];         ## hcp parameter b
  my $hcpC = $_[2];         ## hcp parameter c
  my $nCells1 = $_[3];      ## number of cells in direction of vector 1
  my $nCells2 = $_[4];      ## number of cells in direction of vector 2
  my $nCells3 = $_[5];      ## number of cells in direction of vector 3
  my $tp = $_[6];           ## atom type
  ##   <<<   ----------------   >>>   ##
  my @baseVectors = ([$hcpA,0,0], [0,$hcpB,0], [0,0,$hcpC]);
  my @bvAtoms = ([0,0,0], [1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2]);
  $self->set("isDefined",$VALOK);
  $self->infoLine("Pure $tp fcc structure");
  $self->baseVs(@baseVectors);
  $self->baseAvs(@bvAtoms);
  $self->atomTypeNames([$tp]);
  $self->nCellsArr([$nCells1,$nCells2,$nCells3]);
  $self->struct_base();
}

sub new_struct_bcc {
  my $self = shift;
  ##   <<<   Input parameters   >>>   ##
  my $hcpA = $_[0];         ## hcp parameter a
  my $hcpB = $_[1];         ## hcp parameter b
  my $hcpC = $_[2];         ## hcp parameter c
  my $nCells1 = $_[3];      ## number of cells in direction of vector 1
  my $nCells2 = $_[4];      ## number of cells in direction of vector 2
  my $nCells3 = $_[5];      ## number of cells in direction of vector 3
  my $tp = $_[6];           ## atom type
  ##   <<<   ----------------   >>>   ##
  my @baseVectors = ([$hcpA,0,0], [0,$hcpB,0], [0,0,$hcpC]);
  my @bvAtoms = ([0,0,0], [1/2, 1/2, 1/2]);
  $self->set("isDefined",$VALOK);
  $self->infoLine("Pure $tp fcc structure");
  $self->baseVs(@baseVectors);
  $self->baseAvs(@bvAtoms);
  $self->atomTypeNames([$tp]);
  $self->nCellsArr([$nCells1,$nCells2,$nCells3]);
  $self->struct_base();
}

sub new_struct_diamond {
  my $self = shift;
  ##   <<<   Input parameters   >>>   ##
  my $diaA = $_[0];         ## dia parameter a
  my $diaB = $_[1];         ## dia parameter b
  my $diaC = $_[2];         ## dia parameter c
  my $nCells1 = $_[3];      ## number of cells in direction of vector 1
  my $nCells2 = $_[4];      ## number of cells in direction of vector 2
  my $nCells3 = $_[5];      ## number of cells in direction of vector 3
  my $tp = $_[6];           ## atom type
  ##   <<<   ----------------   >>>   ##
  my @baseVectors = ([$diaA,0,0], [0,$diaB,0], [0,0,$diaC]);
  my @bvAtoms = ([0,0,0], [1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2],
      [1/4, 1/4, 1/4], [3/4, 3/4, 1/4], [3/4, 1/4, 3/4], [1/4, 3/4, 3/4]);
  $self->set("isDefined",$VALOK);
  $self->infoLine("Pure $tp diamond structure");
  $self->baseVs(@baseVectors);
  $self->baseAvs(@bvAtoms);
  $self->atomTypeNames([$tp]);
  $self->nCellsArr([$nCells1,$nCells2,$nCells3]);
  $self->struct_base();
}

sub new_struct_single {
  my $self = shift;
  ##   <<<   Input parameters   >>>   ##
  my $diaA = $_[0];         ## dia parameter a
  my $diaB = $_[1];         ## dia parameter b
  my $diaC = $_[2];         ## dia parameter c
  my $nCells1 = $_[3];      ## number of cells in direction of vector 1
  my $nCells2 = $_[4];      ## number of cells in direction of vector 2
  my $nCells3 = $_[5];      ## number of cells in direction of vector 3
  my $tp = $_[6];           ## atom type
  ##   <<<   ----------------   >>>   ##
  my @baseVectors = ([$diaA,0,0], [0,$diaB,0], [0,0,$diaC]);
  my @bvAtoms = ([1/2, 1/2, 1/2]);
  $self->set("isDefined",$VALOK);
  $self->infoLine("Single $tp atom");
  $self->baseVs(@baseVectors);
  $self->baseAvs(\@bvAtoms);
  $self->atomTypeNames([$tp]);
  $self->nCellsArr([$nCells1,$nCells2,$nCells3]);
  $self->struct_base();
}

sub struct_base {
  my $self = shift;
  if ($self->isDefined ne $VALOK){
    print "ERROR!!! Illigal call for <vaspSTRUCT.struct_base>\n";
    return -1;
  }
  $self->factor(1.0);
  my @baseVectors = @{$self->baseVs()};
  my @bvAtoms = @{$self->baseAvs()};
  my @nCellsArr = @{$self->nCellsArr()};
  my @vNcs = @{$self->nCellsArr};
  my $nCells1 = $nCellsArr[0];
  my $nCells2 = $nCellsArr[1];
  my $nCells3 = $nCellsArr[2];
  my $i = 0; my $j = 0; my $k = 0; my $ii = 0; my $ia; my $ix; my $iy;
  my $nbaseAtoms = @bvAtoms;
  ##print "Number of atoms in the unit cell $nbaseAtoms \n";
  my $nAtoms = $nCells1*$nCells2*$nCells3*$nbaseAtoms;
  my @newAinfC;
  my @newAinfD;
  $ia = 0;
  for ($i = 0; $i<$nCells1; $i++){
    for ($j = 0; $j<$nCells2; $j++){
      for ($k = 0; $k<$nCells3; $k++){
        my @cvA = map {[0,0,0]} 1..$nbaseAtoms;
        my @cvAC = map {[0,0,0]} 1..$nbaseAtoms;
        for ($ii = 0; $ii<$nbaseAtoms; $ii++){
          $ia++;
          $cvA[$ii][0] = $bvAtoms[$ii][0]/$nCells1 + $i/$nCells1;
          $cvA[$ii][1] = $bvAtoms[$ii][1]/$nCells2 + $j/$nCells2;
          $cvA[$ii][2] = $bvAtoms[$ii][2]/$nCells3 + $k/$nCells3;
          ##print $cvA[$ii][0], "\n";
          push(@newAinfD, @cvA[$ii]);
          for ($ix = 0; $ix<3; $ix++){
            for ($iy = 0; $iy<3; $iy++){
              $cvAC[$ii][$ix] += $cvA[$ii][$iy]*$baseVectors[$iy][$ix]*$vNcs[$iy];
            }
          }
          push(@newAinfC, @cvAC[$ii]);
        }
      }
    }
  } ## end of unit cells cicles
  ## setting object values
  $self->ainfC(\@newAinfC);
  $self->ainfD(\@newAinfD);
  $self->atomTypeNums([$nAtoms]);
  $self->nAtoms($nAtoms);
  my $strInfo = $self->infoLine;
  if ($PRINTDBG) {
    print "Structure info = <$strInfo>\n";
    for ($i = 0; $i<$ia; $i++){
      print "$i  [ ";
      for ($j = 0; $j<3; $j++){printf "%9.5f, ", $newAinfC[$i][$j];}
      print " ]    [ ";
      for ($j = 0; $j<3; $j++){printf "%9.5f, ", $newAinfD[$i][$j];}
      print " ] \n";
    }
  }
}## end struct_base

sub write_poscar {
  my $self = shift;
  if ($self->isDefined ne $VALOK){
    print "ERROR!!! Illigal call for <vaspSTRUCT.write_poscar>\n";
    return -1;
  }
  ##   <<<   Input parameters   >>>   ##
  my $fn = $_[0] // $FNPOSCAR;
  ##   <<<   ----------------   >>>   ##
  my $i = 0; my $j = 0;
  my @baseVectors = @{$self->baseVs};
  my @vNcs = @{$self->nCellsArr};
  my @atpNames = @{$self->atomTypeNames};
  my @atpNums = @{$self->atomTypeNums};
  my $nAtomTypes = @atpNames;
  my $fixStyle = " T  T  T";
  my @ainfD = @{$self->ainfD};
  open(FPSCR, ">$fn") or die "Could not open $fn: $!";
  print FPSCR $self->infoLine, "\n";
  print FPSCR $self->factor ," \n";
  for ($i = 0; $i<3; $i++){
    printf FPSCR '   ';
    for ($j = 0; $j<3; $j++){printf FPSCR '%20.14f  ', $baseVectors[$i][$j]*$vNcs[$i];}
    print FPSCR "\n";
  }
  print FPSCR '   ';
  for ($i=0; $i<=$nAtomTypes; $i++){print FPSCR $atpNames[$i],'   ';}
  print FPSCR '   ', "\n";
  print FPSCR '   ';
  for ($i=0; $i<=$nAtomTypes; $i++){print FPSCR $atpNums[$i],'   ';}
  print FPSCR '   ', "\n";
  print FPSCR "Selective dynamics \n";
  print FPSCR "Direct \n";
  for ($i = 0; $i<$self->nAtoms; $i++){
    printf FPSCR ' %15.12f %15.12f %15.12f %s', $ainfD[$i][0], $ainfD[$i][1], $ainfD[$i][2], $fixStyle;
    print FPSCR "\n";
  }
  close FPSCR;
} ## write_poscar

sub write_xsf {
  my $self = shift;
  if ($self->isDefined ne $VALOK){
    print "ERROR!!! Illigal call for <vaspSTRUCT.write_xsf>\n";
    return -1;
  }
  ##   <<<   Input parameters   >>>   ##
  my $fn = $_[0] // $FNXSF;
  ##   <<<   ----------------   >>>   ##
  my $i = 0; my $j = 0;
  my @baseVectors = @{$self->baseVs};
  my @vNcs = @{$self->nCellsArr};
  my @atpNames = @{$self->atomTypeNames};
  my @atpNums = @{$self->atomTypeNums};
  my $nAtomTypes = @atpNames;
  my @ainfC = @{$self->ainfC};
  my $type;
  open(FOUT, ">$fn") or die "Could not open $fn: $!";
  print FOUT "## ", $self->infoLine, "\n";
  print FOUT "CRYSTAL\n";
  print FOUT "PRIMVEC\n";
  for ($i = 0; $i<3; $i++){
    for ($j = 0; $j<3; $j++){printf FOUT '   %20.14f ', $baseVectors[$i][$j]*$vNcs[$i];}
    print FOUT "\n";
  }
  print FOUT "PRIMCOORD\n";
  print FOUT $self->nAtoms," 1 \n";
  for ($i = 0; $i<$self->nAtoms; $i++){
    $type = $self->get_atom_type($i+1);
    print FOUT "  $type    ";
    for ($j = 0; $j<3; $j++){printf FOUT '%20.14f  ', $ainfC[$i][$j];}
    print FOUT "\n";
  }
  close FOUT;
} ## end write_xsf

## na starts from 1 !!!
sub get_atom_type {
  my $self = shift;
  if ($self->isDefined ne $VALOK){
    print "ERROR!!! Illigal call for vaspSTRUCT.get_type \n";
    return 0;
  }
  ##   <<<   Input parameters   >>>   ##
  my $na = $_[0] // 1;
  ##   <<<   ----------------   >>>   ##
  my $ct = 0;
  my $atBord = 0;
  my @atomTypeNums = @{$self->atomTypeNums};
  my @atomTypeNames = @{$self->atomTypeNames};
  while($na > $atBord){
    $atBord += $atomTypeNums[$ct];
    $ct += 1;
  }
  my $tp = $atomTypeNames[$ct-1];
  return $tp;
}

sub new_from_poscar {
  my $self = shift;
  if ($self->isDefined eq $VALOK){
    print "Warning!!! Previous structure will be destroyed\n";
  }
  ##   <<<   Input parameters   >>>   ##
  my $fn = $_[0] // $FNPOSCAR;
  ##   <<<   ----------------   >>>   ##
  my $i = 0; my $j = 0; my $ix = 0; my $iy = 0;
  my $cl;
  my @arrLine;
  my @baseVectors = map {[0,0,0]} 1..3;
  my $nAtoms;
  my $nTypes;
  my @atpNames;
  my @atpNums;
  my @ainfD;
  my @ainfC;
  $self->my_default_vals();
  open(FPSCR, "<$fn") or die "Could not open $fn: $!";
  $cl = <FPSCR>;
  chomp $cl; $cl =~ s/^\s+//;
  $self->infoLine($cl);
  $cl = <FPSCR>;
  chomp $cl; $cl =~ s/^\s+//;
  $self->factor($cl);
  for ($i=0; $i<3; $i++){
    $cl = <FPSCR>;
    chomp $cl; $cl =~ s/^\s+//;
    @arrLine = split /\s+/, $cl;
    for ($j=0; $j<3; $j++){$baseVectors[$i][$j] = $arrLine[$j]}
  }
  $self->baseVs(@baseVectors);
  $self->nCellsArr([1,1,1]);
  ## Atoms types
  $cl = <FPSCR>;
  chomp $cl; $cl =~ s/^\s+//;
  @arrLine = split /\s+/, $cl;
  $nTypes = @arrLine;
  $self->atomTypeNames(\@arrLine);
  $self->recopy_arr("atomTypeNames");
  ## Atoms number
  $cl = <FPSCR>;
  chomp $cl; $cl =~ s/^\s+//;
  @arrLine = split /\s+/, $cl;
  $self->atomTypeNums(\@arrLine);
  $self->recopy_arr("atomTypeNums");
  $nAtoms = 0;
  foreach $i (@arrLine) {$nAtoms += $i}
  $self->nAtoms($nAtoms);
  ## Dynamic type
  $cl = <FPSCR>;
  chomp $cl; $cl =~ s/^\s+//;
  $self->dynType($cl);
  ## Atoms coords type
  $cl = <FPSCR>;
  chomp $cl; $cl =~ s/^\s+//; $cl =~ s/\s+$//;
  $self->aPosType($cl);
  ##print "<$cl>\n";
  if ($cl ne "Direct"){print "ERROR!!! I AM DESIGNED TO READ ONLY Direct COORDINATES\n"}
  $i = 0;
  while ($cl = <FPSCR>){
    if ($i == $nAtoms){last;}
    $i++;
    chomp $cl; $cl =~ s/^\s+//;
    @arrLine = split /\s+/, $cl;
    push(@ainfD,[$arrLine[0],$arrLine[1],$arrLine[2]]);
    my @cvAC = (0,0,0);
    for ($ix = 0; $ix<3; $ix++){
      for ($iy = 0; $iy<3; $iy++){$cvAC[$ix] += $ainfD[$i-1][$iy]*$baseVectors[$iy][$ix]}
    }
    push(@ainfC, \@cvAC);
  }
  close FPSCR;
  $self->ainfC(@ainfC);
  $self->ainfD(@ainfD);
  $self->set("isDefined",$VALOK);

  ## check
  #print "\n";
  #print "atoms Types: \n";
  #map { print $_,"   " } @{$self->atomTypeNames};
  #print "\n";
  #print "atoms Numbers: \n";
  #map { print $_,"   " } @{$self->atomTypeNums};
  #print "\n";
  #print "Atom 2    [ ", $ainfD[1][0]," ",$ainfD[1][1]," ",$ainfD[1][2], " ]\n";
  #print "Atom 2 d  [ ", $ainfC[1][0]," ",$ainfC[1][1]," ",$ainfC[1][2], " ]\n";
  #print "\n";
} ## new_from_poscar

sub distort_bvs {
  my $self = shift;
  if ($self->isDefined ne $VALOK){
    print "ERROR!!! Illigal call for vaspSTRUCT.distort_bvs \n";
    return 0;
  }
  ##   <<<   Input parameters   >>>   ##
  my $dType = $_[0] // "C11";
  my $dAmp  = $_[1] // 0.0;   ## in base vectors coordinates
  my $dAx   = $_[2] // 0;     ## 0 - x,  1 - y,  2 - z
  ##   <<<   ----------------   >>>   ##
  my @T = ([1,0,0],[0,1,0],[0,0,1]);
  my @HD;
  my @bVs = @{$self->baseVs};
  my $i = 0; my $j = 0; my $ix; my $iy;
  my @ainfD = @{$self->ainfD};
  my @ainfC = @{$self->ainfC};
  my @vNcs = @{$self->nCellsArr};
  my @bVMods = (0,0,0);
  for($j=0;$j<3;$j++) { $bVMods[$j] = sqrt($bVs[$j][0]**2+$bVs[$j][1]**2+$bVs[$j][2]**2); }
  $dAmp *= $bVMods[$dAx];    ## !!! CAREFULL HERE !!!
  if ($dType eq "C11"){
    @T = map {[0,0,0]} 1..3;
    $T[$dAx][$dAx] = 1;
  }
  if ($dType eq "C44"){
    @T = map {[0.5,0.5,0.5]} 1..3;
    $T[0][0] = 0;
    $T[1][1] = 0;
    $T[2][2] = 0;
  }
  if ($dType eq "CP"){
    $T[1][1] = -0.5;
    $T[2][2] = -0.5;
  }
  ## V0 is default
  @HD = $self->mat3_mult($self->baseVs,\@T);
  for($i=0; $i<3; $i++){
    for($j=0; $j<3; $j++){$HD[$i][$j] = $bVs[$i][$j] + $HD[$i][$j]*$dAmp}
  }
  $self->baseVs(@HD);
  for ($i = 0; $i<$self->nAtoms; $i++){
    for ($ix = 0; $ix<3; $ix++){
      $ainfC[$i][$ix] = 0.0;
      for ($iy = 0; $iy<3; $iy++){
        $ainfC[$i][$ix] += $ainfD[$i][$iy]*$HD[$iy][$ix]*$vNcs[$iy];
      }
    }
  }
} ## distort_bvs

sub distort_pos {
  my $self = shift;
  if ($self->isDefined ne $VALOK){
    print "ERROR!!! Illigal call for vaspSTRUCT.distort_pos \n";
    return 0;
  }
  ##   <<<   Input parameters   >>>   ##
  my $dAmp  = $_[0] // 0.0;  ## in Angstroms
  ##   <<<   ----------------   >>>   ##
  my @bVs = @{$self->baseVs};
  my $j; my $ix; my $iy;
  my @bVMods = (0,0,0);
  my @vNcs = @{$self->nCellsArr};
  for($j=0;$j<3;$j++) { $bVMods[$j] = sqrt($bVs[$j][0]**2+$bVs[$j][1]**2+$bVs[$j][2]**2); }
  #print "Base vectors mods : @bVMods \n";
  my @ainfD = @{$self->ainfD};
  my @ainfC = @{$self->ainfC};
  for (my $i=0; $i<$self->nAtoms; $i++){
    for ($j=0; $j<3; $j++)
    {
      $ainfD[$i][$j] += 2*(rand() - 0.5)*$dAmp/$bVMods[$j];
      if ($ainfD[$i][$j] > 1.0) { $ainfD[$i][$j] -= 1.0; }
      if ($ainfD[$i][$j] < 0.0) { $ainfD[$i][$j] += 1.0; }
      $ainfD[$i][$j] = sprintf("%20.17f",$ainfD[$i][$j]);
    }
    for ($ix = 0; $ix<3; $ix++){
      $ainfC[$i][$ix] = 0.0;
      for ($iy = 0; $iy<3; $iy++){
        $ainfC[$i][$ix] += $ainfD[$i][$iy]*$bVs[$iy][$ix]*$vNcs[$iy];
      }
    }
    ##print "$i Pos vector : ",$ainfD[$i][0]," ",$ainfD[$i][1]," ",$ainfD[$i][2],"\n";
  }
}

## multiplication
sub mat3_mult {
  my $self = shift;
  ##   <<<   Input parameters   >>>   ##
  my @m1 = @{$_[0]};
  my @m2 = @{$_[1]};
  ##   <<<   ----------------   >>>   ##
  my $i = 0; my $j = 0; my $k = 0;
  my @m = map {[0,0,0]} 1..3;
  for ($i=0; $i < 3; $i++){
    for ($j=0; $j < 3; $j++){
      for ($k=0; $k < 3; $k++){$m[$i][$j] += $m1[$i][$k]*$m2[$k][$j]}
    }
  }
  return @m;
}
























##
