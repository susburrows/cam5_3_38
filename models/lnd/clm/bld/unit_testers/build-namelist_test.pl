#!/usr/bin/env perl

# Test command line options of the build-namelist script.
# Try to test that all the different options at least work.
# Test that inconsistentcies are appropriately caught.

#########################

use Test::More;
use xFail::expectedFail;

#########################

use strict;
use Getopt::Long;
use NMLTest::CompFiles;
use English;

sub usage {
    die <<EOF;
SYNOPSIS
     build-namelist_test.pl [options]

     Test the the CLM build-namelist 
OPTIONS
     -help [or -h]                 Print usage to STDOUT.                               
     -compare <directory>          Compare namelists for this version to namelists
                                   created by another version.
     -generate                     Leave the namelists in place to do a later compare.
     -test                         Use the -test option to make sure datasets exist.
     -csmdata "dir"                Root directory of CESM input data.

EOF
}


#
# Process command-line options.
#
my %opts = ( help     => 0,
             generate => 0,
             test     => 0,
             compare  => undef,
             csmdata  => undef,
            );

GetOptions(
    "h|help"     => \$opts{'help'},
    "compare=s"  => \$opts{'compare'},
    "generate"   => \$opts{'generate'},
    "test"       => \$opts{'test'},
    "csmdata=s"  => \$opts{'csmdata'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

# Check that the CESM inputdata root directory has been specified.  This must be
# a local or nfs mounted directory.
my $inputdata_rootdir = undef;
if (defined($opts{'csmdata'})) {
    $inputdata_rootdir = $opts{'csmdata'};
} elsif (defined $ENV{'CSMDATA'} ) { 
    $inputdata_rootdir = $ENV{'CSMDATA'};
} else {
   # use yellowstone location as default
   $inputdata_rootdir="/glade/p/cesm/cseg/inputdata";
   print("WARNING:  -csmdata nor CSMDATA are set, using default yellowstone location: $inputdata_rootdir\n");
}

###################################
#_# read in expected fail test list
###################################
my $compGen;
if ( $opts{'generate'} eq 1 && !(defined($opts{'compare'}) )) {
   $compGen='generate';
} elsif ( defined($opts{'compare'}) ) {
   $compGen='compare';
} elsif ( defined($opts{'compare'} && ($opts{'generate'} eq 1 ))) {
   #_# if compare and generate are both given, use compare
   $compGen='compare'; 
}

my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!;
my $testType="namelistTest";

#
# Figure out number of tests that will run
#
my $ntests = 315;
if ( defined($opts{'compare'}) ) {
   $ntests += 175;
}
plan( tests=>$ntests );

#_# ============================================================
#_# setup for xFail module
#_# ============================================================
my $xFail = xFail::expectedFail->new($ProgName,$compGen,$ntests);
my $captOut="";  #_# variable to capture Test::More output
Test::More->builder->output(\$captOut);
#_# ============================================================
#_# 
#_# ============================================================

# Check for unparsed arguments
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}
my $mode = "standard";
system( "../configure -s" );

my $DOMFILE = "$inputdata_rootdir/atm/datm7/domain.lnd.T31_gx3v7.090928.nc";
my $bldnml = "../build-namelist -verbose -csmdata $inputdata_rootdir -lnd_frac $DOMFILE -no-note";
if ( $opts{'test'} ) {
   $bldnml .= " -test";
}

my $tempfile = "temp_file.txt";
if ( -f $tempfile ) {
  system( "/bin/rm $tempfile" );
}

my @files = ( "lnd_in", $tempfile );
my $cwd = `pwd`;
chomp( $cwd );
my $cfiles = NMLTest::CompFiles->new( $cwd, @files );

print "\n==================================================\n";
print "Run simple tests \n";
print "==================================================\n";

# Simple test -- just run build-namelist with -help option
eval{ system( "$bldnml -help > $tempfile 2>&1 " ); };
   is( $@, '', "help" );
   &cleanup();
# Simple test -- just run build-namelist with -version option
eval{ system( "$bldnml -version > $tempfile 2>&1 " ); };
   is( $@, '', "version" );
   system( "/bin/cat $tempfile" );
   &cleanup();
# Simple test -- just run build-namelist
eval{ system( "$bldnml > $tempfile 2>&1 " ); };
   is( $@, '', "plain build-namelist" );
   $cfiles->checkfilesexist( "default", $mode ); 
   # Compare to baseline
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "default", $mode );
      $cfiles->comparefiles( "default", $mode, $opts{'compare'} );
   }

print "\n==================================================\n";
print "Run simple tests with all list options \n";
print "==================================================\n";

$cfiles->copyfiles( "default", $mode );
&cleanup();
# Simple test -- run all the list options
foreach my $options ( "clm_demand", "rcp",      "res", 
                      "sim_year",   "use_case" ) {
   eval{ system( "$bldnml -${options} list > $tempfile 2>&1 " ); };
   my $result = `cat $tempfile`;
   my $expect;
   if ( $options =~ /use_case/ ) {
      $expect = "use cases :";
   } else {
      $expect = "valid values for $options";
   }
   $expect    = "/CLM build-namelist : $expect/";
   like( $result, $expect, "$options list" );
   is( (-f "lnd_in"), undef, "Check that lnd_in file does NOT exist" );
   &cleanup();
}

print "\n==================================================\n";
print "Run simple tests with additional options \n";
print "==================================================\n";

# Exercise a bunch of options
my $options = "-co2_ppmv 250 -glc_nec 10 -glc_grid gland5 -glc_smb .false.";
   $options .= " -res 0.9x1.25 -rcp 2.6";

   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
      $cfiles->checkfilesexist( "default", $mode );
      $cfiles->copyfiles( "most_options", $mode );
   # Compare to default
      $cfiles->doNOTdodiffonfile( "lnd_in",    "default", $mode );
      $cfiles->doNOTdodiffonfile( "$tempfile", "default", $mode );
      $cfiles->comparefiles( "default", $mode );
   # Compare to baseline
   if ( defined($opts{'compare'}) ) {
      $cfiles->dodiffonfile(      "lnd_in",    "most_options", $mode );
      $cfiles->doNOTdodiffonfile( "$tempfile", "most_options", $mode );
      $cfiles->comparefiles( "most_options", $mode, $opts{'compare'} );
   }
   &cleanup();

print "\n==================================================\n";
print "Test drydep and megan namelists  \n";
print "==================================================\n";

# drydep and megan namelists
my @mfiles = ( "lnd_in", "drv_flds_in", $tempfile );
my $mfiles = NMLTest::CompFiles->new( $cwd, @mfiles );
foreach my $options ( "-drydep", "-megan", "-drydep -megan" ) {
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
   $mfiles->checkfilesexist( "$options", $mode);
   if ( $options ne "-drydep" ) {
     $mfiles->shownmldiff( "-drydep", $mode );
   }
   if ( defined($opts{'compare'}) ) {
      $mfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $mfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $mfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}

print "\n==================================================\n";
print "Test irrig, verbose, clm_demand, rcp, test, sim_year, use_case, l_ncpl\n";
print "==================================================\n";

# irrig, verbose, clm_demand, rcp, test, sim_year, use_case, l_ncpl
my $startfile = "clmrun.clm2.r.1964-05-27-00000.nc";
foreach my $options ( "-irrig .true. ", "-verbose", "-rcp 2.6", "-test", "-sim_year 1850",
                      "-use_case 1850_control", "-l_ncpl 1", 
                      "-clm_start_type startup", 
                      "-envxml_dir . -infile myuser_nl_clm", 
                     ) {
   my $file = $startfile;
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
   $cfiles->checkfilesexist( "$options", $mode );
   system( "diff lnd_in lnd_in.default" );
   $cfiles->shownmldiff( "default", $mode );
   my $finidat = `grep finidat lnd_in`;
   if (      $options eq "-l_ncpl 1" ) {
      my $dtime = `grep dtime lnd_in`;
      like( $dtime, "/ 86400\$/", "$options" );
   } elsif ( $options =~ /myuser_nl_clm/ ) {
      my $fsurdat =  `grep fsurdat lnd_in`;
      like( $fsurdat, "/MYDINLOCROOT/lnd/clm2/PTCLMmydatafiles/1x1pt_US-UMB/surfdata_1x1pt_US-UMB_simyr2000_clm4_5_c131122.nc/", "$options" );
   }
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}


print "\n==================================================\n";
print "Start Failure testing.  These should fail \n";
print "==================================================\n";

# Failure testing, do things that SHOULD fail
my $finidat  = "thing.nc";
system( "touch $finidat" );

my %failtest = ( 
     "coldstart but with IC file"=>{ options=>"-clm_start_type cold",
                                     namelst=>"finidat='$finidat'",
                                     conopts=>"",
                                   },
     "l_ncpl is zero"            =>{ options=>"-l_ncpl 0",
                                     namelst=>"",
                                     conopts=>"",
                                   },
     "l_ncpl not integer"        =>{ options=>"-l_ncpl 1.0",
                                     namelst=>"",
                                     conopts=>"",
                                   },
     "both l_ncpl and dtime"     =>{ options=>"-l_ncpl 24",
                                     namelst=>"dtime=1800",
                                     conopts=>"",
                                   },
     "use_crop without -crop"    =>{ options=>"",
                                     namelst=>"use_crop=.true.",
                                     conopts=>"-phys clm4_5",
                                   },
     "CNDV with fpftdyn"         =>{ options=>"",
                                     namelst=>"fpftdyn='myfpftdynfile.nc'",
                                     conopts=>"-bgc cndv",
                                   },
     "CNDV with fpftdyn - clm4_5"=>{ options=>"-bgc bgc -dynamic_vegetation",
                                     namelst=>"fpftdyn='myfpftdynfile.nc'",
                                     conopts=>"-phys clm4_5",
                                   },
     "use_cndv=T without bldnml op"=>{ options=>"-bgc cn",
                                     namelst=>"use_cndv=.true.",
                                     conopts=>"-phys clm4_5",
                                   },
     "use_cndv=F with dyn_veg op"=>{ options=>"-bgc cn -dynamic_vegetation",
                                     namelst=>"use_cndv=.false.",
                                     conopts=>"-phys clm4_5",
                                   },
     "crop with use_crop false"  =>{ options=>"-crop -bgc bgc",
                                     namelst=>"use_crop=.false.",
                                     conopts=>"-phys clm4_5",
                                   },
     "crop without CN"           =>{ options=>"-crop -bgc sp",
                                     namelst=>"",
                                     conopts=>"-phys clm4_5",
                                   },
     "irrigate=T without -irr op"=>{ options=>"-crop -bgc cn",
                                     namelst=>"irrigate=.true.",
                                     conopts=>"-phys clm4_5",
                                   },
     "irrigate=F with -irrg op"  =>{ options=>"-crop -bgc cn -irrig .true.",
                                     namelst=>"irrigate=.false.",
                                     conopts=>"-phys clm4_5",
                                   },
     "-irrig without -crop"      =>{ options=>"-bgc cn -irrig .true.",
                                     namelst=>"",
                                     conopts=>"-phys clm4_5",
                                   },
     "sp and c13"                =>{ options=>"-bgc sp",
                                     namelst=>"use_c13=.true.",
                                     conopts=>"-phys clm4_5",
                                   },
     "sp and c14"                =>{ options=>"-bgc sp",
                                     namelst=>"use_c14=.true.",
                                     conopts=>"-phys clm4_5",
                                   },
     "crop and c13"              =>{ options=>"-crop -bgc bgc",
                                     namelst=>"use_c13=.true.",
                                     conopts=>"-phys clm4_5",
                                   },
     "crop and c14"              =>{ options=>"-crop -bgc cn",
                                     namelst=>"use_c14=.true.",
                                     conopts=>"-phys clm4_5",
                                   },
     "bgc=cn and bgc settings"   =>{ options=>"-bgc cn",
                                     namelst=>"use_lch4=.true.,use_nitrif_denitrif=.true.,use_vertsoilc=.true.,use_century_decomp=.true.",
                                     conopts=>"-phys clm4_5",
                                   },
     "bgc=bgc and cn-only set"   =>{ options=>"-bgc bgc",
                                     namelst=>"use_lch4=.false.,use_nitrif_denitrif=.false.,use_vertsoilc=.false.,use_century_decomp=.false.",
                                     conopts=>"-phys clm4_5",
                                   },
     "use_cn=true bgc=sp"        =>{ options=>"-bgc sp",
                                     namelst=>"use_cn=.true.",
                                     conopts=>"-phys clm4_5",
                                   },
     "use_cn=false bgc=cn"       =>{ options=>"-bgc cn",
                                     namelst=>"use_cn=.false.",
                                     conopts=>"-phys clm4_5",
                                   },
     "vichydro without clm4_5"   =>{ options=>"-vichydro",
                                     namelst=>"",
                                     conopts=>"-phys clm4_0",
                                   },
     "use_vic=F with -vic op"    =>{ options=>"-vichydro",
                                     namelst=>"use_vichydro=.false.",
                                     conopts=>"-phys clm4_5",
                                   },
     "bgc without clm4_5"        =>{ options=>"-bgc sp",
                                     namelst=>"",
                                     conopts=>"-phys clm4_0",
                                   },
     "bgc_spinup without clm4_5" =>{ options=>"-bgc_spinup on",
                                     namelst=>"",
                                     conopts=>"-phys clm4_0",
                                   },
     "DV without clm4_5"         =>{ options=>"-dynamic_vegetation",
                                     namelst=>"",
                                     conopts=>"-phys clm4_0",
                                   },
     "bgc_spinup without cn"     =>{ options=>"-bgc_spinup on -bgc sp",
                                     namelst=>"",
                                     conopts=>"-phys clm4_5",
                                   },
     "spinup=1 without bldnml op"=>{ options=>"-bgc bgc",
                                     namelst=>"spinup_state=1",,
                                     conopts=>"-phys clm4_5",
                                   },
     "spinup=0 with bldnml op"   =>{ options=>"-bgc bgc -bgc_spinup on",
                                     namelst=>"spinup_state=0",
                                     conopts=>"-phys clm4_5",
                                   },
     "both co2_type and on nml"  =>{ options=>"-co2_type constant",
                                     namelst=>"co2_type='prognostic'",
                                     conopts=>"",
                                   },
     "both lnd_frac and on nml"  =>{ options=>"-lnd_frac domain.nc",
                                     namelst=>"fatmlndfrc='frac.nc'",
                                     conopts=>"",
                                   },
     "branch but NO nrevsn"      =>{ options=>"-clm_start_type branch",
                                     namelst=>"",
                                     conopts=>"",
                                   },
     "glc_nec inconsistent"      =>{ options=>"-glc_nec 10",
                                     namelst=>"maxpatch_glcmec=5",
                                     conopts=>"",
                                   },
     "glc_smb inconsistent"      =>{ options=>"-glc_nec 10 -glc_smb .true.",
                                     namelst=>"glc_smb=.false.",
                                     conopts=>"",
                                   },
     "glc_grid inconsistent"     =>{ options=>"-glc_nec 10 -glc_grid gland10",
                                     namelst=>"glc_grid='gland5'",
                                     conopts=>"",
                                   },
     "envxml_not_dir"            =>{ options=>"-envxml_dir myuser_nl_clm",
                                     namelst=>"",
                                     conopts=>"",
                                   },
     "envxml_emptydir"           =>{ options=>"-envxml_dir xFail",
                                     namelst=>"",
                                     conopts=>"",
                                   },
               );
foreach my $key ( keys(%failtest) ) {
   system( "../configure -s ".$failtest{$key}{"conopts"});
   my $options  = $failtest{$key}{"options"};
   my $namelist = $failtest{$key}{"namelst"};
   eval{ system( "$bldnml $options -namelist \"&clmexp $namelist /\" > $tempfile 2>&1 " ); };
   isnt( $?, 0, $key );
   system( "cat $tempfile" );
}

print "\n==================================================\n";
print "Test ALL resolutions with CN \n";
print "==================================================\n";

# Check for ALL resolutions with CN
my $mode = "CN";
system( "../configure -s -bgc cn" );
my $reslist = `../queryDefaultNamelist.pl -res list -s`;
my @resolutions = split( / /, $reslist );
$options = "";
my @regional;
foreach my $res ( @resolutions ) {
   chomp($res);
   print "=== Test $res === \n";
   my $options  = "-res $res";

   if ( $res eq "512x1024" ) { 
      $options .= " -sim_year 1850"; 
   } elsif ( $res =~ /^([0-9]+x[0-9]+_[a-zA-Z]+)$/ ) {
      push( @regional, $res );
      next;
   } elsif ( $res eq "0.5x0.5"     ||
             $res eq "0.1x0.1"     ||
             $res eq "3x3min"      ||
             $res eq "5x5min"      ||
             $res eq "10x10min"    ||
             $res eq "0.33x0.33"   ||
             $res eq "1km-merge-10min" ) {
      next;
   }

   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );

   $cfiles->checkfilesexist( "$options", $mode );
   system( "diff lnd_in lnd_in.default.standard" );

   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }

   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup(); print "\n";
}

print "\n==================================================\n";
print " Test all use-cases \n";
print "==================================================\n";

# Run over all use-cases...
my $list = `$bldnml -use_case list 2>&1 | grep "use case"`;
my @usecases;
if ( $list =~ /build-namelist : use cases : (.+)$/ ) {
  my @usecases  = split( / /, $list );
} else {
  die "ERROR:: Trouble getting list of use-cases\n";
}
foreach my $usecase ( @usecases ) {
   $options = "-use_case $usecase ";
   eval{ system( "$bldnml $options  > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
   $cfiles->checkfilesexist( "$options", $mode );
   system( "diff lnd_in lnd_in.default.standard" );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}

print "\n==================================================\n";
print "Test single-point regional cases \n";
print "==================================================\n";

# Run over single-point regional cases
foreach my $res ( @regional ) {
   $mode = "$res";
   system( "../configure -s -sitespf_pt $res" );
   eval{ system( "$bldnml > $tempfile 2>&1 " ); };
   is( $@, '', "$res" );
   $cfiles->checkfilesexist( "$res", $mode );
   system( "diff lnd_in lnd_in.default.standard" );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$res", $mode );
      $cfiles->comparefiles( "$res", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$res", $mode );
   }
   &cleanup();
}

print "\n==================================================\n";
print "Test crop resolutions \n";
print "==================================================\n";

# Check for crop resolutions
my $mode = "crop";
system( "../configure -s -crop on -bgc cn" );
my @crop_res = ( "10x15", "1.9x2.5" );
foreach my $res ( @crop_res ) {
   $options = "-res $res";
   eval{ system( "$bldnml $options  > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   system( "diff lnd_in lnd_in.default.standard" );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}
print "\n==================================================\n";
print " Test glc_mec resolutions \n";
print "==================================================\n";

# Check for glc_mec resolutions
my $mode = "standard";
system( "../configure -s" );
my @glc_res = ( "48x96", "0.9x1.25", "1.9x2.5" );
my @use_cases = ( "1850-2100_rcp2.6_glacierMEC_transient",
                  "1850-2100_rcp4.5_glacierMEC_transient",
                  "1850-2100_rcp6_glacierMEC_transient",
                  "1850-2100_rcp8.5_glacierMEC_transient",
                  "1850_glacierMEC_control",
                  "2000_glacierMEC_control",
                  "20thC_glacierMEC_transient",
                 );
my $GLC_NEC         = 10;
foreach my $res ( @glc_res ) {
   foreach my $usecase ( @usecases ) {
      $options = "-glc_nec $GLC_NEC -res $res -use_case $usecase ";
      eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
      is( $@, '', "$options" );
      $cfiles->checkfilesexist( "$options", $mode );
      system( "diff lnd_in lnd_in.default.standard" );
      $cfiles->shownmldiff( "default", "standard" );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
      }
      if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$options", $mode );
      }
      &cleanup();
   }
}
# Transient 20th Century simulations
my $mode = "standard";
system( "../configure -s" );
my @tran_res = ( "48x96", "0.9x1.25", "1.9x2.5", "ne30np4", "ne60np4", "ne120np4", "10x15", "1x1_tropicAtl" );
my $usecase  = "20thC_transient";
my $GLC_NEC         = 0;
foreach my $res ( @tran_res ) {
   $options = "-res $res -use_case $usecase ";
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   system( "diff lnd_in lnd_in.default.standard" );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}
# Transient rcp scenarios
my $mode = "standard";
system( "../configure -s" );
my @tran_res = ( "48x96", "0.9x1.25", "1.9x2.5", "ne30np4", "10x15" );
foreach my $usecase ( "1850-2100_rcp2.6_transient", "1850-2100_rcp4.5_transient", "1850-2100_rcp6_transient", "1850-2100_rcp8.5_transient" ) {
   foreach my $res ( @tran_res ) {
      $options = "-res $res -use_case $usecase ";
      eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
      is( $@, '', "$options" );
      $cfiles->checkfilesexist( "$options", $mode );
      system( "diff lnd_in lnd_in.default.standard" );
      $cfiles->shownmldiff( "default", "standard" );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
      }
      if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$options", $mode );
      }
      &cleanup();
   }
}

print "\n==================================================\n";
print "Test clm4.5/clm5.0 resolutions \n";
print "==================================================\n";

foreach my $phys ( "clm4_5", 'clm5_0' ) {
  my $mode;
  if ( $phys eq "clm4_5" ) {
     $mode = "phys45";
  } else {
     $mode = "phys50";
  }
  system( "../configure -s -phys ".$phys );
  my $clmoptions = "-bgc bgc";
  my @clmres = ( "10x15", "48x96", "0.9x1.25", "1.9x2.5", "360x720cru" );
  foreach my $res ( @clmres ) {
     $options = "-res $res";
     eval{ system( "$bldnml $options $clmoptions  > $tempfile 2>&1 " ); };
     is( $@, '', "$options" );
     $cfiles->checkfilesexist( "$options", $mode );
     system( "diff lnd_in lnd_in.default.standard" );
     $cfiles->shownmldiff( "default", "standard" );
     if ( defined($opts{'compare'}) ) {
        $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
        $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
     }
     if ( defined($opts{'generate'}) ) {
        $cfiles->copyfiles( "$options", $mode );
     }
     &cleanup();
  }
  my $mode = "phys-crop";
  system( "../configure -s -phys ".$phys );
  my $clmoptions = "-bgc cn -crop";
  my $res = "1.9x2.5";
  $options = "-res $res -irrig .true. -crop -bgc cn";
  eval{ system( "$bldnml $options $clmoptions  > $tempfile 2>&1 " ); };
  is( $@, '', "$options" );
  $cfiles->checkfilesexist( "$options", $mode );
  system( "diff lnd_in lnd_in.default.standard" );
  $cfiles->shownmldiff( "default", "standard" );
  if ( defined($opts{'compare'}) ) {
     $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
     $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
  }
  if ( defined($opts{'generate'}) ) {
     $cfiles->copyfiles( "$options", $mode );
  }
  &cleanup();
}

system( "/bin/rm $finidat" );

print "\n==================================================\n";
print " Dumping output  \n";
print "==================================================\n";

$xFail->parseOutput($captOut);

print "Successfully ran all testing for build-namelist\n\n";

&cleanup( "config" );
system( "/bin/rm lnd_in.default" );
system( "/bin/rm $tempfile" );

sub cleanup {
#
# Cleanup files created
#
  my $type = shift;

  print "Cleanup files created\n";
  if ( defined($type) ) {
     if ( $type eq "config" ) {
        system( "/bin/rm Filepath config_cache.xml CESM_cppdefs" );
     }
  } else {
     system( "/bin/rm $tempfile *_in" );
  }
}

