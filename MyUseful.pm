#!/usr/bin/perl

package MyUseful;
use Exporter;
@ISA = ('Exporter');
@EXPORT = ('print_v', 'print_m',
  'v3_mult', 'v3_mod', 'v3_m3_mult',
  'float_eq');

our $FLOATACC = 1.0e-10;

#sub dummy {
  ##   <<<   Input parameters   >>>   ##
  ##   <<<   ----------------   >>>   ##
#}

sub print_v {
  print "[  ";
  foreach my $val (@_){
    printf " $val ";}
  print "]\n";
}

sub print_m {
  my @m;
  my $n = @_;
  if ($n == 1){@m = @{$_[0]}}
  else {@m = @_}
  foreach my $v (@m){
    foreach my $val (@{$v}){
      printf " $val ";}
    print "\n";
  }
}

sub v3_mult {
  ##   <<<   Input parameters   >>>   ##
  my @v1 = @{$_[0]};
  my @v2 = @{$_[1]};
  ##   <<<   ----------------   >>>   ##
  my @v = (0.0, 0.0, 0.0);
  $v[0] = $v1[1]*$v2[2] - $v2[1]*$v1[2];
  $v[1] = $v2[0]*$v1[2] - $v1[0]*$v2[2];
  $v[2] = $v1[0]*$v2[1] - $v2[0]*$v1[1];
  return @v;
}

sub v3_mod {
  my @v = @{$_[0]};
  return sqrt($v[0]*$v[0] + $v[1]*$v[1] + $v[2]*$v[2]);
}

sub v3_m3_mult {
  ##   <<<   Input parameters   >>>   ##
  my @v = @{$_[0]};
  my @m = @{$_[1]};
  ##   <<<   ----------------   >>>   ##
  my @nv;
  for(my $i = 0; $i < 3; $i++){
    for(my $i = 0; $i < 3; $i++){ $nv[$i] = $v[$j]*$m[$j][$i]; }
  }
  return @nv;
}

sub m3alaV_Invert {
  my @m = @{$_[0]};
  my @mInv;
  my $mDet = ($m[0]*$m[4]*$m[8] + $m[1]*$m[5]*$m[6] + $m[2]*$m[3]*$m[7]
    - $m[2]*$m[4]*$m[6] - $m[1]*$m[3]*$m[8] - $m[0]*$m[5]*$m[7]);
  if ( float_eq($mDet,0)){
    print("PAY ATTANTION determinant is 0");
    return @mInv;
  }
  $mInv[0] =  ($m[4]*$m[8]-$m[5]*$m[7])/$mDet;
  $mInv[1] = -($m[1]*$m[8]-$m[2]*$m[7])/$mDet;
  $mInv[2] =  ($m[1]*$m[5]-$m[2]*$m[4])/$mDet;
  $mInv[3] = -($m[3]*$m[8]-$m[5]*$m[6])/$mDet;
  $mInv[4] =  ($m[0]*$m[8]-$m[2]*$m[6])/$mDet;
  $mInv[5] = -($m[0]*$m[5]-$m[2]*$m[3])/$mDet;
  $mInv[6] =  ($m[3]*$m[7]-$m[4]*$m[6])/$mDet;
  $mInv[7] = -($m[0]*$m[7]-$m[1]*$m[6])/$mDet;
  $mInv[8] =  ($m[0]*$m[4]-$m[1]*$m[3])/$mDet;
  return @mInv;
}

sub m3_Invert {
  my @m = @{$_[0]};
  my @mInv;
  my $mDet = ($m[0][0]*$m[1][1]*$m[2][2] + $m[0][1]*$m[1][2]*$m[2][0] + $m[0][2]*$m[1][0]*$m[2][1]
    - $m[0][2]*$m[1][1]*$m[2][0] - $m[0][1]*$m[1][0]*$m[2][2] - $m[0][0]*$m[1][2]*$m[2][1]);
  if ( float_eq($mDet,0)){
    print("PAY ATTANTION determinant is 0");
    return @mInv;
  }
  $mInv[0][0] =  ($m[1][1]*$m[2][2]-$m[1][2]*$m[2][1])/$mDet;
  $mInv[0][1] = -($m[0][1]*$m[2][2]-$m[0][2]*$m[2][1])/$mDet;
  $mInv[0][2] =  ($m[0][1]*$m[1][2]-$m[0][2]*$m[1][1])/$mDet;
  $mInv[1][0] = -($m[1][0]*$m[2][2]-$m[1][2]*$m[2][0])/$mDet;
  $mInv[1][1] =  ($m[0][0]*$m[2][2]-$m[0][2]*$m[2][0])/$mDet;
  $mInv[1][2] = -($m[0][0]*$m[1][2]-$m[0][2]*$m[1][0])/$mDet;
  $mInv[2][0] =  ($m[1][0]*$m[2][1]-$m[1][1]*$m[2][0])/$mDet;
  $mInv[2][1] = -($m[0][0]*$m[2][1]-$m[0][1]*$m[2][0])/$mDet;
  $mInv[2][2] =  ($m[0][0]*$m[1][1]-$m[0][1]*$m[1][0])/$mDet;
  return @mInv;
}

sub float_eq {
  my $val1 = $_[0];
  my $val2 = $_[1];
  my $acc = $FLOATACC;
  if ( @_ > 2 ) { $acc = $_[2]; }
  if ( abs($val1 - $val2) < $acc ){ return 1; }
  else { return 0; }
}

## ========================================================================== ##
##                          Subroutines not for export                        ##
## ========================================================================== ##























##
