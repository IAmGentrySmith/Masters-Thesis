#!/usr/bin/perl -w

# Converts an xyz file to a series of povray files

use Getopt::Std;
use constant M_PI => 3.1415926536;
use POSIX;
# use Math::Trig;


# default variables
#$bondTol = 1.5;		# bonds formed to atoms less than this distance apart
#$vdwFrac = 0.2;			# spheres at 20% of their Van der Waals radius
$buffer = 4.0;			# image atoms within this distance from the walls for clipping
#$bondSize = 0.125;  # default bond size
#$shrinkFrac = 1.0; 	# shrink factor for water molecules
#$nonWordChar = 0;   # a counter for unknown atom types
$bWalls = 0;        # default bWalls to off
#$outlineWidth = 0.075;
#$outlining = 0; 
$tolerance = 0.001;
$sep_dist = 1;
$n_radius = 2.0;
$p_radius = 2.0;
$r_radius = 0;
$dist_tol2 = 25;
$cull_wat = 0;
$single_ion = 0;
$make_neutral = 0;
$solute_count = 0;
$x_shift = 0; 
$y_shift = 0;
$channel_halfwidth = 5.5;
$scale_factor = 20; 
push(@color_list, "1.0,0.98,0.8");
push(@color_list, "1.0,0.94,0.4");
push(@color_list, "1.0,0.90,0.2");
push(@color_list, "1.0,0.85,0.0");
push(@color_list, "1.0,0.75,0.0");
$color_list_factor = 2;
$findHBonds = 0;
$ion_pairs_count = -1;
$total_ions = $ion_pairs_count * 2;
$total_ions_count = -1;
$hbond_len_tol = 3.625;
$hbond_ang_tol = 25;
$hbond_ang_tol2 = 120 - $hbond_ang_tol;



# get our options
#getopts('acdhob:s:t:v:w:');
getopts('bchiqwa:d:f:l:m:n:p:r:s:t:x:y:z:');

# if we don't have a filename, either write a ATOMPARAM or drop to -h
if ($#ARGV != 0){
#	if ($opt_a){
## write out a default ATOMPARAM file if -a flag is used
#  		writeAtomParams();
#		die "\nFor more options, use the -h command.\n\n";
#	} else {
    $opt_h = 'true';
#	}
}

if ($opt_h) {
    print "\n$0: converts a 2D-OOPSE .dump coordinate file (.dump) to a povray input file (.pov)\n\n";
    print "usage: $0 -[options] [input dump file]\n\n";
    print "Options:\n";
    print "  -b : turn on hbond and nucleation identification\n";
    print "  -c : center the configurations about the origin\n";
    print "  -h : show this message\n";
    print "  -i : treat the first molecule as a single ion \n";
    print "  -q : treat ions as neutrals \n";
    print "  -w : include bounding walls\n\n";

    print "  -a (real) : alpha transparency (default: 0)\n";
    print "  -d (real) : distance tolerance in Å (default: 5)\n";
    print "  -f (real) : scale factor for render resolution (default: 20)\n";
    print "  -l (int)  : oil color value 0-4 (default: 3)\n";
    print "  -m (int)  : MB styling as 1 - black, 2 - dipole colored, 3 - colored spokes, 4 - dipole dots, 5 - no dipole, 6 - black spokes\n";
    print "  -n (real) : negative ion radius size, when used with option -s (default: 2.0 Å)\n";
    print "  -p (real) : positive ion radius size, when used with option -s (default: 2.0 Å)\n";
    print "  -r (real) : neutral particle radius size (default: 2.0 Å)\n";
    print "  -s (real) : treat the first molecule as an ion pair with this separation distance\n";
    print "  -t (int)  : rose type - 3, 4, or 5 (default: 3)\n";
    print "  -x (real) : shift x-coordinates by this amount (default: 0 Å)\n";
    print "  -y (real) : shift y-coordinates by this amount (default: 0 Å)\n";
    print "  -z (int)  : number of ion pairs for salt imaging\n";
    die "\n";
}

$fileName = $ARGV[0];

# some input checking
if ($opt_w) {
    $bWalls = $opt_w;
} else {
    $bWalls = 0;
}

if ($opt_b){
    $findHBonds = $opt_b;
} else {
    $findHBonds = 0;
}

if ($opt_c){
    $centerConfig = $opt_c;
} else {
    $centerConfig = 0;
}

if ($opt_i){
    $single_ion = $opt_i;
    $single_ion = 3;
} else {
    $single_ion = 0;
}

if ($opt_q){
    $make_neutral = $opt_q;
    $make_neutral = 1;
} else {
    $make_neutral = 0;
}

if (defined($opt_a)) {
    if ($opt_a >= 0 && $opt_a <= 1) {
        $alpha_t = $opt_a;
        $alpha_t2 = 0.5*(1.0 - $alpha_t) + $alpha_t;
    } else {
        die "\tThe alpha transparency (-a) value ($opt_a) is not a valid real number.\n\tPlease choose a decimal number between 0 and 1.\n";
    }
} else {
    $alpha_t = 0;
    $alpha_t2 = 0;
}

if (defined($opt_d)) {
    if ($opt_d >= 0) {
        $dist_tol2 = $opt_d*$opt_d;
        $cull_wat = 1;
    } else {
        die "\tThe distance tolerance (-d) value ($opt_d) is not a valid real number.\n\tPlease choose a decimal number greater than 0.\n";
    }
}

if (defined($opt_f)) {
    if ($opt_f >= 0) {
        $scale_factor = $opt_f;
    } else {
        die "\tThe resolution scale factor (-f) value ($opt_f) is not a valid real number.\n\tPlease choose a decimal number greater than 0.\n";
    }
}

if (defined($opt_l)) {
    if ($opt_l == 0 || $opt_l == 1 || $opt_l == 2 || $opt_l == 3 || $opt_l == 4) {
        $color_list_factor = $opt_l;
    } else {
        die "\tThe oil color option (-l) value ($opt_l) is not a valid integer option.\n\tPlease choose 0, 1, 2, 3, or 4.\n";
    }
} 

if (defined($opt_m)) {
    if ($opt_m == 1 || $opt_m == 2 || $opt_m == 3 || $opt_m == 4 || $opt_m == 5 || $opt_m == 6) {
        $mb_style = $opt_m;
    } else {
        die "\tThe MB style (-m) value ($opt_m) is not a valid integer option.\n\tPlease choose 1, 2, 3, 4, 5, or 6.\n";
    }
} else {
    $mb_style = 0;
}

if (defined($opt_n)) {
    if ($opt_n > 0) {
        $n_radius = $opt_n;
    } else {
        die "\tThe negative ion radius (-n) value ($opt_n) is not a valid number.\n\tPlease choose a real number greater than 0.\n";
    }
}

if (defined($opt_p)) {
    if ($opt_p > 0) {
        $p_radius = $opt_p;
    } else {
        die "\tThe positive ion radius (-p) value ($opt_p) is not a valid number.\n\tPlease choose a real number greater than 0.\n";
    }
}

if (defined($opt_r)) {
    if ($opt_r > 0) {
        $r_radius = $opt_r;
    } else {
        die "\tThe neutral particle radius (-r) value ($opt_r) is not a valid number.\n\tPlease choose a real number greater than 0.\n";
    }
}

if (defined($opt_s)) {
    if ($opt_s > 0) {
        $sep_dist = 0.5*$opt_s;
        $channel_halfwidth = $sep_dist;
    } else {
        die "\tThe ion separation distance (-s) value ($opt_s) is not a valid number.\n\tPlease choose a real number greater than 0.\n";
    }
}

if (defined($opt_t)) {
    if ($opt_t == 3 || $opt_t == 4 || $opt_t == 5) {
        $rose_type = $opt_t;
    } else {
        die "\tThe rose type (-t) value ($opt_t) is not a valid integer option.\n\tPlease choose 3, 4, or 5.\n";
    }
} else {
    $rose_type = 3;
}

if (defined($opt_x)) {
    if ($opt_x > 0 || $opt_x < 0) {
        $x_shift = $opt_x;
    } else {
        die "\tThe x shift (-x) value ($opt_x) is not a valid number.\n\tPlease choose a real number greater than or less than 0.\n";
    }
}

if (defined($opt_y)) {
    if ($opt_y > 0 || $opt_y < 0) {
        $y_shift = $opt_y;
    } else {
        die "\tThe y shift (-y) value ($opt_y) is not a valid number.\n\tPlease choose a real number greater than or less than 0.\n";
    }
}

if (defined($opt_z)) {
    if ($opt_z > 0){
        $ion_pairs_count = $opt_z;
        $total_ions = $ion_pairs_count * 2;
        $total_ions_count = 0;
    } else {
        die "\tThe (-z) # of ion pairs value ($opt_z) is not a valid number.\n\tPlease choose an integer number greater than 0.\n";
    }
}


open(DUMPFILE, "./$fileName") || die "Error: Unable to open the file $fileName for reading.\n";
#loadAtomProps();

mkdir "pov", 0777 || die "Error: can't mkdir pov: $!";

if ($single_ion > 0){
    $single_ion = 1 if defined($opt_p);
    $single_ion = 2 if defined($opt_n);
}

#$count = 0;
$frameCount = 0;
$wrapFlag = 0;

$framedat_region = 0;
$stuntdoubles_region = 0;

while (<DUMPFILE>) {
    @line = split;

    if (defined($line[0])){
        $framedat_region = 1 if $line[0] eq '<FrameData>'; 
        $framedat_region = 0 if $line[0] eq '</FrameData>'; 
        $stuntdoubles_region = 1 if $line[0] eq '<StuntDoubles>'; 
        $stuntdoubles_region = 0 if $line[0] eq '</StuntDoubles>'; 

        if ($framedat_region == 1 && $line[0] ne '<FrameData>'){
            $frameComment = $_ if $line[0] eq 'Time:';
            if ($line[0] eq 'Hmat:'){
                chop($line[2]);
                chop($line[8]);
                $hxx = $line[2];
                $hyy = $line[8];
                $hzz = 100;
                $halfx = 0.5*$hxx;
                $halfy = 0.5*$hyy;
                $halfz = 0.5*$hzz;
            }
        }

        if ($stuntdoubles_region == 1 && $line[0] ne '<StuntDoubles>'){
            shiftCoordinates();
            wrapCoordinates();# if $wrapFlag == 0;
            centerCoordinates();
            saveCoordinates();
        }

        if ($line[0] eq '</StuntDoubles>') {
            # write a povray file
            writePOVRAY();
            wipeArrays();
        }

    }
}

# the following are our subroutines

sub writePOVRAY {
    open(POVOUT, ">pov/$frameCount.pov") || die "Error: Unable to open pov/$frameCount.pov for writing.\n";
    povHeader();
    povOptions();

    if ($bWalls == 1){
        if ($wrapFlag == 0){ 
            povBoundWalls();
        } else {
            print STDERR "Warning: The boundry box cannot be included due to the lack of box dimensions.\n";
        }
    }
    povAtoms();

    $frameCount++;
}

sub povHeader {
    print POVOUT "
    //**************************************
    // generated by rose2pov.pl
    // $frameComment
    //**************************************\n";
}

sub povOptions {
    $ratio = $hxx/$hyy;
    $camera_pos = 5*$hxx;
    #printf STDOUT "povray -w%-4i -h%-4i +a0.1 +UA -D $frameCount.pov\n", $scale_factor*$hxx, $scale_factor*$hyy;
    #printf STDOUT "../../../bin/povray -w%-4i -h%-4i +a0.1 -D $frameCount.pov\n", $scale_factor*$hxx, $scale_factor*$hyy;
    printf STDOUT "povray -w%-4i -h%-4i +a0.1 -D $frameCount.pov\n", $scale_factor*$hxx, $scale_factor*$hyy;
    print POVOUT "
    //**************************************
    // Lights, camera, resolution!
    //**************************************

#declare Ratio = $ratio ;
#declare zoom = $camera_pos; 
#declare RAD = off;

global_settings {
#if(RAD)
radiosity {
pretrace_start 0.08
pretrace_end   0.01
count 500

nearest_count 10
error_bound 0.02
recursion_limit 1

low_error_factor 0.2
gray_threshold 0.0
minimum_reuse 0.015
brightness 1.4

adc_bailout 0.01/2
}
#end
}


camera{
orthographic
location < 0, 0, zoom >
direction < 0, 0, 2 >
up < 0, $hyy, 0 >
// Ratio is negative to switch povray to a right hand coordinate system.
right < -$hxx, 0, 0 >
look_at < 0, 0, 0 >
}

background { color rgb 1 }

global_settings { ambient_light rgb 1 }
light_source { < 0, 0, zoom > rgb 1.45 }
//light_source { < -zoom/2, zoom/2, zoom >  rgb < 0.5, 0.5, 0.5 > }

#macro make_3_lobe_rose (center_x, center_y, center_z, angle_theta)
difference{
sphere{
< center_x, center_y, center_z >,
1.0
texture{
pigment{ rgb .95 }
}
}
plane{ 
<0, 0, -1>, -center_z+0.1
texture{
pigment { rgb 0.95 }
}               
}
}
difference {
blob {
threshold .65
sphere{<0.000, 0.000, 0.000>,1.2, 1 scale .8}   
sphere{<0.000, 1.33, 0.000>,1.1, 1 scale .8}   
sphere{<1.152, -0.665, 0.000>,1.1, 1 scale .8}   
sphere{<-1.152,-0.665, 0.000>,1.1, 1 scale .8}

texture{
pigment{ rgb < 1.0, 0.5, 0.8 > }
finish{
ambient .2
diffuse .6
specular 1
roughness .001
}
}
rotate <0, 0, angle_theta> 
translate<center_x, center_y, center_z>
}

plane{ 
<0, 0, -1>, -center_z-0.1
texture{
pigment{ rgb < 1.0, 0.5, 0.9 > }
}               
}
}
difference{
blob {
threshold .55
sphere{<0.000, 0.000, 0.000>,1.2, 1 scale .8}   
sphere{<0.000, 1.33, 0.000>,1.1, 1 scale .8}   
sphere{<1.152, -0.665, 0.000>,1.1, 1 scale .8}   
sphere{<-1.152,-0.665, 0.000>,1.1, 1 scale .8}
texture{
pigment { rgb 0 }
}               
rotate <0, 0, angle_theta> 
translate<center_x, center_y, center_z>
}
plane{ 
<0, 0, -1>, -center_z
texture{
pigment { rgb 0 }
}               
}
}
#end

#macro make_4_lobe_rose (center_x, center_y, center_z, angle_theta)
difference{
sphere{
< center_x, center_y, center_z >,
1.0
texture{
pigment{ rgb .95 }
}
}
plane{ 
<0, 0, -1>, -center_z+0.1
texture{
pigment { rgb 0.95 }
}               
}
}
difference {
blob {
threshold .65
sphere{<0.0, 0.0, 0.0>,1.2, 1 scale 0.8}   
sphere{<0.0+0.50897, 0.0+1.22876,0.0>,1.1, 1 scale 0.8}   
sphere{<0.0-1.22876,0.0+0.50897, 0.0>,1.1, 1 scale 0.8}   
sphere{<0.0-0.50897,0.0-1.22876, 0.0>,1.1, 1 scale 0.8}
sphere{<0.0+1.22876, 0.0-0.50897,0.0>,1.1, 1 scale 0.8}

texture{
pigment{ rgb < 0.5, 0.5, 1.0 > }
finish{
ambient .2
diffuse .6
specular 1
roughness .001
}
}
rotate <0, 0, angle_theta> 
translate<center_x, center_y, center_z>
}

plane{ 
<0, 0, -1>, -center_z-0.1
texture{
pigment{ rgb < 0.5, 0.75, 1.0 > }
}               
}
}
difference{
blob {
threshold .55
sphere{<0.0, 0.0, 0.0>,1.2, 1 scale 0.8}   
sphere{<0.0+0.50897, 0.0+1.22876,0.0>,1.1, 1 scale 0.8}   
sphere{<0.0-1.22876,0.0+0.50897, 0.0>,1.1, 1 scale 0.8}   
sphere{<0.0-0.50897,0.0-1.22876, 0.0>,1.1, 1 scale 0.8}
sphere{<0.0+1.22876, 0.0-0.50897,0.0>,1.1, 1 scale 0.8}
texture{
pigment { rgb 0 }
}               
rotate <0, 0, angle_theta> 
translate<center_x, center_y, center_z>
}
plane{ 
<0, 0, -1>, -center_z
texture{
pigment { rgb 0 }
}               
}
}

#end

#macro make_5_lobe_rose (center_x, center_y, center_z, angle_theta)
difference{
sphere{
< center_x, center_y, center_z >,
1.0
texture{
pigment{ rgb .95 }
}
}
plane{ 
<0, 0, -1>, -center_z+0.1
texture{
pigment { rgb 0.95 }
}               
}
}
difference {
blob {
threshold .65
sphere{<0.0, 0.0, 0.0>,1.2, 1 scale 0.8}   
sphere{<0.0-0.78175, 0.0+1.07599,0.0>,1.1, 1 scale 0.8}   
sphere{<0.0-1.26491,0.0-0.41099, 0.0>,1.1, 1 scale 0.8}   
sphere{<0.0,0.0-1.33, 0.0>,1.1, 1 scale 0.8}
sphere{<0.0+1.26491, 0.0-0.41099,0.0>,1.1, 1 scale 0.8}
sphere{<0.0+0.78175, 0.0+1.07599,0.0>,1.1, 1 scale 0.8}

texture{
pigment{ rgb < 1.0, 0.75, 0.5 > }
finish{
ambient .2
diffuse .6
specular 1
roughness .001
}
}
rotate <0, 0, angle_theta> 
translate<center_x, center_y, center_z>
}

plane{ 
<0, 0, -1>, -center_z-0.1
texture{
pigment{ rgb < 1.0, 0.75, 0.5 > }
}               
}
}
difference{
blob {
threshold .55
sphere{<0.0, 0.0, 0.0>,1.13, 1 scale 0.8}   
sphere{<0.0-0.78175, 0.0+1.07599,0.0>,1.1, 1 scale 0.8}   
sphere{<0.0-1.26491,0.0-0.41099, 0.0>,1.1, 1 scale 0.8}   
sphere{<0.0,0.0-1.33, 0.0>,1.1, 1 scale 0.8}
sphere{<0.0+1.26491, 0.0-0.41099,0.0>,1.1, 1 scale 0.8}
sphere{<0.0+0.78175, 0.0+1.07599,0.0>,1.1, 1 scale 0.8}
texture{
pigment { rgb 0 }
}               
rotate <0, 0, angle_theta> 
translate<center_x, center_y, center_z>
}
plane{ 
<0, 0, -1>, -center_z
texture{
pigment { rgb 0 }
}               
}
}

#end

#macro make_mbnd_black (center_x, center_y, center_z, angle_theta)
union{
difference{
cylinder{
< 0,0,-0.19 >,
< 0,0,0.19 >,
1.3
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
cylinder{
< 0,0,-0.29 >,
< 0,0,+0.29 >,
1.1
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
}
cylinder{
< 0,0,-0.19 >,
< 0,0,+0.19 >,
0.10
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
box{
< -0.1, 0, +0.1 >,
< +0.1, +1.2, -0.1 >
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
box{
< -0.1, 0, +0.1 >,
< +0.1, +1.2, -0.1 >
rotate <0, 0, 240>
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
box{
< -0.1, 0, +0.1 >,
< +0.1, +1.2, -0.1 >
rotate <0, 0, 120>
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
rotate <0, 0, angle_theta> 
translate<center_x, center_y, center_z>
}

#end

#macro make_mb_black (center_x, center_y, center_z, angle_theta)
union{
difference{
cylinder{
< 0,0,-0.19 >,
< 0,0,0.19 >,
1.3
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
cylinder{
< 0,0,-0.29 >,
< 0,0,+0.29 >,
1.1
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
}
cylinder{
< 0,0,-0.19 >,
< 0,0,+0.19 >,
0.10
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
box{
< -0.1, 0, +0.1 >,
< +0.1, +1.2, -0.1 >
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
intersection{
difference{
cone {
<0, 0.7, 0>, 0.1    
<0, 0.4, 0>, 0.25
texture{
pigment{ rgbt <0,0,0,$alpha_t> }
}
}
cone {
    <0, 0.5, 0>, 0.1    
    <0, 0.3, 0>, 0.40
    texture{
        pigment{ rgbt <0,0,0,$alpha_t> }
    }
}
}
box{
    < -0.3, 0, +0.01 >,
    < +0.3, +1.2, -0.01 >
    texture{
        pigment{ rgbt <0,0,0,$alpha_t> }
    }
}
}

box{
    < -0.1, 0, +0.1 >,
    < +0.1, +1.2, -0.1 >
    rotate <0, 0, 240>
    texture{
        pigment{ rgbt <0,0,0,$alpha_t> }
    }
}
box{
    < -0.1, 0, +0.1 >,
    < +0.1, +1.2, -0.1 >
    rotate <0, 0, 120>
    texture{
        pigment{ rgbt <0,0,0,$alpha_t> }
    }
}
rotate <0, 0, angle_theta> 
translate<center_x, center_y, center_z>
}

#end

#macro make_mb_gold (center_x, center_y, center_z, angle_theta)
union{
    difference{
        cylinder{
            < 0,0,-0.19 >,
            < 0,0,0.19 >,
            1.3
            texture{
                pigment{ rgbt <0,0,0,$alpha_t> }
            }
        }
        cylinder{
            < 0,0,-0.29 >,
            < 0,0,+0.29 >,
            1.1
            texture{
                pigment{ rgbt <0,0,0,$alpha_t> }
            }
        }
    }
    difference{
        cylinder{
            < 0,0,-0.195 >,
            < 0,0,0.195 >,
            1.25
            texture{
                pigment{ rgbt <1,0.8,0,$alpha_t> }
            }
        }
        cylinder{
            < 0,0,-0.295 >,
            < 0,0,+0.295 >,
            1.15
            texture{
                pigment{ rgbt <1,0.8,0,$alpha_t> }
            }
        }
    }
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        texture{
            pigment{ rgbt <0,0,0,$alpha_t> }
        }
    }
    box{
        < -0.05, 0, +0.2>,
        < +0.05, +1.2, -0.2>
        texture{
            pigment{ rgbt <1,0.8,0,$alpha_t> }
        }
    }
    intersection{
        difference{
            cone {
                <0, 0.7, 0>, 0.1    
                <0, 0.4, 0>, 0.25
                texture{
                    pigment{ rgbt <0,0,0,$alpha_t> }
                }
            }
            cone {
                <0, 0.5, 0>, 0.1    
                <0, 0.3, 0>, 0.40
                texture{
                    pigment{ rgbt <0,0,0,$alpha_t> }
                }
            }
        }
        box{
            < -0.3, 0, +0.01 >,
            < +0.3, +1.2, -0.01 >
            texture{
                pigment{ rgbt <0,0,0,$alpha_t> }
            }
        }
    }
    intersection{
        difference{
            cone {
                <0, 0.7, 0>, 0.090  
                <0, 0.475, 0>, 0.2
                texture{
                    pigment{ rgbt <1,0.8,0,$alpha_t> }
                }
            }
            cone {
                <0, 0.55, 0>, 0.1   
                <0, 0.35, 0>, 0.40
                texture{
                    pigment{ rgbt <1,0.8,0,$alpha_t> }
                }
            }
        }
        box{
            < -0.3, 0, +0.11>,
            < +0.3, +1.2, -0.11>
            texture{
                pigment{ rgbt <1,0.8,0,$alpha_t> }
            }
        }
    }

    box{
        < -0.1, 0, +0.191 >,
        < +0.1, +1.2, -0.191 >
        rotate <0, 0, 240>
        texture{
            pigment{ rgbt <0,0,0,$alpha_t> }
        }
    }
    box{
        < -0.05, 0, +0.2>,
        < +0.05, +1.2, -0.2>
        rotate <0, 0, 240>
        texture{
            pigment{ rgbt <1,0.85,0,$alpha_t> }
        }
    }
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        rotate <0, 0, 120>
        texture{
            pigment{ rgbt <0,0,0,$alpha_t> }
        }
    }
    box{
        < -0.05, 0, +0.2>,
        < +0.05, +1.2, -0.2>
        rotate <0, 0, 120>
        texture{
            pigment{ rgbt <1,0.8,0,$alpha_t> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_mb_dipole (center_x, center_y, center_z, angle_theta)
union{
    difference{
        cylinder{
            < 0,0,-0.19 >,
            < 0,0,0.19 >,
            1.3
            texture{
                pigment{ rgbt <0,0,0,$alpha_t> }
            }
        }
        cylinder{
            < 0,0,-0.29 >,
            < 0,0,+0.29 >,
            1.1
            texture{
                pigment{ rgbt <0,0,0,$alpha_t> }
            }
        }
    }
    cylinder{
        < 0,0,-0.19 >,
        < 0,0,+0.19 >,
        0.10
        texture{
            pigment{ rgbt <0,0,0,$alpha_t> }
        }
    }
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        texture{
            pigment{ rgbt <0,0.68,1,$alpha_t> }
        }
    }
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        rotate <0, 0, 240>
        texture{
            pigment{ rgbt <1,0.25,0.25,$alpha_t2> }
        }
    }
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        rotate <0, 0, 120>
        texture{
            pigment{ rgbt <1,0.25,0.25,$alpha_t2> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_mb_spokes (center_x, center_y, center_z, angle_theta)
union{
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        texture{
            pigment{ rgbt <0,0.68,1,$alpha_t> }
        }
    }
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        rotate <0, 0, 240>
        texture{
            pigment{ rgbt <1,0.25,0.25,$alpha_t2> }
        }
    }
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        rotate <0, 0, 120>
        texture{
            pigment{ rgbt <1,0.25,0.25,$alpha_t2> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_mb_black_spokes (center_x, center_y, center_z, angle_theta)
union{
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        // < -0.025, 0, +0.025 >,
        // < +0.025, +1.2, -0.025 >
        texture{
            pigment{ rgbt <0,0,0,(0.5*(1-$alpha_t)+$alpha_t)> }
        }
    }
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        // < -0.025, 0, +0.025 >,
        // < +0.025, +1.2, -0.025 >
        rotate <0, 0, 240>
        texture{
            pigment{ rgbt <0,0,0,(0.5*(1-$alpha_t)+$alpha_t)> }
        }
    }
    box{
        < -0.1, 0, +0.1 >,
        < +0.1, +1.2, -0.1 >
        // < -0.025, 0, +0.025 >,
        // < +0.025, +1.2, -0.025 >
        rotate <0, 0, 120>
        texture{
            pigment{ rgbt <0,0,0,(0.5*(1-$alpha_t)+$alpha_t)> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_mb_orange_spokes(center_x, center_y, center_z, angle_theta)
union{
    box{
        < -0.025, 0, +0.025 >,
        < +0.025, +1.2, -0.025 >
        rotate <0, 0, 240>
        texture{
            pigment{ rgbt <1,0.55,0,$alpha_t> }
        }
    }
    box{
        < -0.025, 0, +0.025 >,
        < +0.025, +1.2, -0.025 >
        rotate <0, 0, 120>
        texture{
            pigment{ rgbt <1,0.55,0,$alpha_t> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_mb_orange_dots(center_x, center_y, center_z, angle_theta)
union{
    cylinder{
        < 0,1.3,-0.1 >,
        < 0,1.3,+0.1 >,
        0.15
        rotate <0, 0, 240>
        texture{
            pigment{ rgbt <1,0.575,0,$alpha_t> }
        }
    }
    cylinder{
        < 0,1.3,-0.1 >,
        < 0,1.3,+0.1 >,
        0.15
        rotate <0, 0, 120>
        texture{
            pigment{ rgbt <1,0.575,0,$alpha_t> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_mb_orange_tabs(center_x, center_y, center_z, angle_theta)
union{
    box{
        < -0.05, 0.8, +0.025 >,
        < +0.05, +1.4, -0.025 >
        rotate <0, 0, 240>
        texture{
            pigment{ rgbt <1,0.55,0,$alpha_t> }
        }
    }
    box{
        < -0.05, 0.8, +0.025 >,
        < +0.05, +1.4, -0.025 >
        rotate <0, 0, 120>
        texture{
            pigment{ rgbt <1,0.55,0,$alpha_t> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_mb_dots (center_x, center_y, center_z, angle_theta)
union{
    cylinder{
        < 0,0,-0.1 >,
        < 0,0,+0.1 >,
        0.25
        texture{
            pigment{ rgbt <1,0.25,0.25,$alpha_t> }
        }
    }
    cylinder{
        < 0,0.495,-0.1 >,
        < 0,0.495,+0.1 >,
        0.25
        texture{
            pigment{ rgbt <0,0.68,1,$alpha_t> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_neut_dots (center_x, center_y, center_z)
union{
    cylinder{
        < 0,0,-0.1 >,
        < 0,0,+0.1 >,
        0.5
        texture{
            pigment{ rgbt <1,0.55,0.0,$alpha_t> }
        }
    }
    translate<center_x, center_y, center_z>
}

#end

#macro make_ion_pair (center_x, center_y, center_z, angle_theta)
cylinder{
    < -$sep_dist,0,0 >,
    < -$sep_dist,0,+0.21 >,
    $p_radius
    texture{
        pigment{ rgbf <0,0.68,1.0,$alpha_t> }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}
union{
    difference{
        cylinder{
            < -$sep_dist,0,0 >,
            < -$sep_dist,0,+0.205 >,
            $p_radius+0.2
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
        cylinder{
            < -$sep_dist,0,0 >,
            < -$sep_dist,0,+0.21 >,
            $p_radius
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
    }
    union{
        box{
            < -$sep_dist-0.5,0.075,0>
            < -$sep_dist+0.5,-0.075,0.23>
            texture{
                pigment{ rgbf <1,1,1,$alpha_t> }
            }
        }
        box{
            < -$sep_dist-0.075,-0.5,0>
            < -$sep_dist+0.075,0.5,0.23>
            texture{
                pigment{ rgbf <1,1,1,$alpha_t> }
            }
        }
    }

    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}
cylinder{
    < $sep_dist,0,0 >,
    < $sep_dist,0,+0.21 >,
    $n_radius
    texture{
        pigment{ rgbf <1.0,0.25,0.25,$alpha_t> }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}
union{
    difference{
        cylinder{
            < $sep_dist,0,0 >,
            < $sep_dist,0,+0.205 >,
            $n_radius+0.2
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
        cylinder{
            < $sep_dist,0,0 >,
            < $sep_dist,0,+0.21 >,
            $n_radius
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
    }
    box{
        < $sep_dist-0.5,0.075,0>
        < $sep_dist+0.5,-0.075,0.23>
        texture{
            pigment{ rgbf <1,1,1,$alpha_t> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_neut_pair (center_x, center_y, center_z, angle_theta)
cylinder{
    < -$sep_dist,0,0 >,
    < -$sep_dist,0,+0.21 >,
    $p_radius
    texture{
        pigment{ rgbf <1.0,0.85,0.0,$alpha_t> }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}
union{
    difference{
        cylinder{
            < -$sep_dist,0,0 >,
            < -$sep_dist,0,+0.205 >,
            $p_radius+0.2
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
        cylinder{
            < -$sep_dist,0,0 >,
            < -$sep_dist,0,+0.21 >,
            $p_radius
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}
cylinder{
    < $sep_dist,0,0 >,
    < $sep_dist,0,+0.21 >,
    $n_radius
    texture{
        pigment{ rgbf <1.0,0.85,0.0,$alpha_t> }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}
union{
    difference{
        cylinder{
            < $sep_dist,0,0 >,
            < $sep_dist,0,+0.205 >,
            $n_radius+0.2
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
        cylinder{
            < $sep_dist,0,0 >,
            < $sep_dist,0,+0.21 >,
            $n_radius
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_neut_ion (center_x, center_y, center_z)
cylinder{
    < 0,0,0 >,
    < 0,0,+0.21 >,
    $r_radius
    texture{
        pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
    }
    translate<center_x, center_y, center_z>
}
union{
    difference{
        cylinder{
            < 0,0,0 >,
            < 0,0,+0.205 >,
            $r_radius+0.2
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
        cylinder{
            < 0,0,0 >,
            < 0,0,+0.21 >,
            $r_radius
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
    }
    translate<center_x, center_y, center_z>
}
#end

#macro make_pos_ion (center_x, center_y, center_z)
cylinder{
    < 0,0,0 >,
    < 0,0,+0.21 >,
    $p_radius
    texture{
        pigment{ rgbf <0,0.68,1.0,$alpha_t> }
    }
    translate<center_x, center_y, center_z>
}
union{
    difference{
        cylinder{
            < 0,0,0 >,
            < 0,0,+0.205 >,
            $p_radius+0.2
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
        cylinder{
            < 0,0,0 >,
            < 0,0,+0.21 >,
            $p_radius
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
    }
    union{
        box{
            < -0.5,0.075,0>
            < +0.5,-0.075,0.23>
            texture{
                pigment{ rgbf <1,1,1,$alpha_t> }
            }
        }
        box{
            < -0.075,-0.5,0>
            < +0.075,0.5,0.23>
            texture{
                pigment{ rgbf <1,1,1,$alpha_t> }
            }
        }
    }

    translate<center_x, center_y, center_z>
}
#end

#macro make_neg_ion (center_x, center_y, center_z)
cylinder{
    < 0,0,0 >,
    < 0,0,+0.21 >,
    $n_radius
    texture{
        pigment{ rgbf <1.0,0.25,0.25,$alpha_t> }
    }
    translate<center_x, center_y, center_z>
}
union{
    difference{
        cylinder{
            < 0,0,0 >,
            < 0,0,+0.205 >,
            $n_radius+0.2
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
        cylinder{
            < 0,0,0 >,
            < 0,0,+0.21 >,
            $n_radius
            texture{
                pigment{ rgbf <0,0,0,$alpha_t> }
            }
        }
    }
    box{
        < -0.5,0.075,0>
        < +0.5,-0.075,0.23>
        texture{
            pigment{ rgbf <1,1,1,$alpha_t> }
        }
    }
    translate<center_x, center_y, center_z>
}

#end

#macro make_wall (center_x, center_y, center_z, angle_theta)
union{
    cylinder{
        < -10,0,0 >,
        < -10,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -8,0,0 >,
        < -8,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -6,0,0 >,
        < -6,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -4,0,0 >,
        < -4,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -2,0,0 >,
        < -2,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -0,0,0 >,
        < -0,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 2,0,0 >,
        < 2,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 4,0,0 >,
        < 4,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 6,0,0 >,
        < 6,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 8,0,0 >,
        < 8,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 10,0,0 >,
        < 10,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}
union{
    difference{
        union{
            cylinder{
                < -10,0,0 >,
                < -10,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -8,0,0 >,
                < -8,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -6,0,0 >,
                < -6,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -4,0,0 >,
                < -4,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -2,0,0 >,
                < -2,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -0,0,0 >,
                < -0,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 2,0,0 >,
                < 2,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 4,0,0 >,
                < 4,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 6,0,0 >,
                < 6,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 8,0,0 >,
                < 8,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 10,0,0 >,
                < 10,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
        }
        union{
            cylinder{
                < -10,0,0 >,
                < -10,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -8,0,0 >,
                < -8,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -6,0,0 >,
                < -6,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -4,0,0 >,
                < -4,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -2,0,0 >,
                < -2,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -0,0,0 >,
                < -0,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 2,0,0 >,
                < 2,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 4,0,0 >,
                < 4,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 6,0,0 >,
                < 6,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 8,0,0 >,
                < 8,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 10,0,0 >,
                < 10,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
        }
    }

    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

#macro make_channel (center_x, center_y, center_z, angle_theta)
union{
    cylinder{
        < -10,$channel_halfwidth,0 >,
        < -10,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -8,$channel_halfwidth,0 >,
        < -8,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -6,$channel_halfwidth,0 >,
        < -6,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -4,$channel_halfwidth,0 >,
        < -4,$channel_halfwidth,0 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -2,$channel_halfwidth,0 >,
        < -2,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -0,$channel_halfwidth,0 >,
        < -0,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 2,$channel_halfwidth,0 >,
        < 2,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 4,$channel_halfwidth,0 >,
        < 4,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 6,$channel_halfwidth,0 >,
        < 6,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 8,$channel_halfwidth,0 >,
        < 8,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 10,$channel_halfwidth,0 >,
        < 10,$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -10,-$channel_halfwidth,0 >,
        < -10,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -8,-$channel_halfwidth,0 >,
        < -8,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -6,-$channel_halfwidth,0 >,
        < -6,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -4,-$channel_halfwidth,0 >,
        < -4,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -2,-$channel_halfwidth,0 >,
        < -2,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -0,-$channel_halfwidth,0 >,
        < -0,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 2,-$channel_halfwidth,0 >,
        < 2,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 4,-$channel_halfwidth,0 >,
        < 4,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 6,-$channel_halfwidth,0 >,
        < 6,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 8,-$channel_halfwidth,0 >,
        < 8,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 10,-$channel_halfwidth,0 >,
        < 10,-$channel_halfwidth,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}
union{
    difference{
        union{
            cylinder{
                < -10,$channel_halfwidth,0 >,
                < -10,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -8,$channel_halfwidth,0 >,
                < -8,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -6,$channel_halfwidth,0 >,
                < -6,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -4,$channel_halfwidth,0 >,
                < -4,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -2,$channel_halfwidth,0 >,
                < -2,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -0,$channel_halfwidth,0 >,
                < -0,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 2,$channel_halfwidth,0 >,
                < 2,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 4,$channel_halfwidth,0 >,
                < 4,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 6,$channel_halfwidth,0 >,
                < 6,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 8,$channel_halfwidth,0 >,
                < 8,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 10,$channel_halfwidth,0 >,
                < 10,$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -10,-$channel_halfwidth,0 >,
                < -10,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -8,-$channel_halfwidth,0 >,
                < -8,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -6,-$channel_halfwidth,0 >,
                < -6,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -4,-$channel_halfwidth,0 >,
                < -4,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -2,-$channel_halfwidth,0 >,
                < -2,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -0,-$channel_halfwidth,0 >,
                < -0,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 2,-$channel_halfwidth,0 >,
                < 2,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 4,-$channel_halfwidth,0 >,
                < 4,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 6,-$channel_halfwidth,0 >,
                < 6,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 8,-$channel_halfwidth,0 >,
                < 8,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 10,-$channel_halfwidth,0 >,
                < 10,-$channel_halfwidth,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
        }
        union{
            cylinder{
                < -10,$channel_halfwidth,0 >,
                < -10,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -8,$channel_halfwidth,0 >,
                < -8,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -6,$channel_halfwidth,0 >,
                < -6,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -4,$channel_halfwidth,0 >,
                < -4,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -2,$channel_halfwidth,0 >,
                < -2,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -0,$channel_halfwidth,0 >,
                < -0,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 2,$channel_halfwidth,0 >,
                < 2,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 4,$channel_halfwidth,0 >,
                < 4,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 6,$channel_halfwidth,0 >,
                < 6,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 8,$channel_halfwidth,0 >,
                < 8,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 10,$channel_halfwidth,0 >,
                < 10,$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -10,-$channel_halfwidth,0 >,
                < -10,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -8,-$channel_halfwidth,0 >,
                < -8,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -6,-$channel_halfwidth,0 >,
                < -6,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -4,-$channel_halfwidth,0 >,
                < -4,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -2,-$channel_halfwidth,0 >,
                < -2,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -0,-$channel_halfwidth,0 >,
                < -0,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 2,-$channel_halfwidth,0 >,
                < 2,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 4,-$channel_halfwidth,0 >,
                < 4,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 6,-$channel_halfwidth,0 >,
                < 6,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 8,-$channel_halfwidth,0 >,
                < 8,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 10,-$channel_halfwidth,0 >,
                < 10,-$channel_halfwidth,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
        }
    }

    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end
#macro make_surface (center_x, center_y, center_z, angle_theta)
union{
    cylinder{
        < -10,0,0 >,
        < -10,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -8,0,0 >,
        < -8,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -6,0,0 >,
        < -6,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -4,0,0 >,
        < -4,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -2,0,0 >,
        < -2,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < -0,0,0 >,
        < -0,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 2,0,0 >,
        < 2,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 4,0,0 >,
        < 4,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 6,0,0 >,
        < 6,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 8,0,0 >,
        < 8,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    cylinder{
        < 10,0,0 >,
        < 10,0,+0.21 >,
        $p_radius
        texture{
            pigment{ rgbf <$color_list[$color_list_factor],$alpha_t> }
        }
    }
    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}
union{
    difference{
        union{
            cylinder{
                < -10,0,0 >,
                < -10,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -8,0,0 >,
                < -8,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -6,0,0 >,
                < -6,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -4,0,0 >,
                < -4,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -2,0,0 >,
                < -2,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -0,0,0 >,
                < -0,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 2,0,0 >,
                < 2,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 4,0,0 >,
                < 4,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 6,0,0 >,
                < 6,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 8,0,0 >,
                < 8,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 10,0,0 >,
                < 10,0,+0.205 >,
                $p_radius+0.2
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
        }
        union{
            cylinder{
                < -10,0,0 >,
                < -10,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -8,0,0 >,
                < -8,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -6,0,0 >,
                < -6,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -4,0,0 >,
                < -4,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -2,0,0 >,
                < -2,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < -0,0,0 >,
                < -0,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 2,0,0 >,
                < 2,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 4,0,0 >,
                < 4,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 6,0,0 >,
                < 6,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 8,0,0 >,
                < 8,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
            cylinder{
                < 10,0,0 >,
                < 10,0,+0.21 >,
                $p_radius
                texture{
                    pigment{ rgbf <0,0,0,$alpha_t> }
                }
            }
        }
    }

    rotate <0, 0, angle_theta> 
    translate<center_x, center_y, center_z>
}

#end

\n";
}

sub povBoundWalls {
print POVOUT "
//**************************************
// make boundry wall
//**************************************\n\n

box{
    <-$halfx,-$halfy,-$hzz> <$halfx,$halfy,$hzz>
    texture{
        pigment{ rgbt < 0.9, 0.9, 0.9, 0.8 > }
        finish{
            ambient .2
            diffuse .6
            specular 1
            roughness .001
        }
    }
    clipped_by{ box{ <-$hxx,-$hyy,-$halfz> <$hxx,$hyy,$halfz> }}
}\n";

}

sub povAtoms {
    if ($findHBonds == 1){
        buildHBondMap();
    }
    print POVOUT "
    //**************************************
    // List of all the atoms
    //**************************************\n\n";

    for ($k=0; $k<=$#angle_thetas; $k++) {
        if ($include_mol[$k] == 1){
            if ($centerConfig == 1){
                if ($sep_dist > 0 && $k==0){
                    if ($make_neutral == 1 && $bwalls == 0){
                        printf POVOUT ("make_neut_pair( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                    } elsif ($make_neutral == 1 && $bwalls == 1){
                    #    printf POVOUT ("make_channel( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                        printf POVOUT ("make_surface( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                    } else {
                        printf POVOUT ("make_ion_pair( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                    }
                } elsif ($single_ion == 1 && $k==0){
                    printf POVOUT ("make_pos_ion( %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0);
                } elsif ($single_ion == 2 && $k==0){
                    printf POVOUT ("make_neg_ion( %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0);
                } elsif ($r_radius > 0 && $k>=0 && $is_solute[$k] == 1){
                    if ($alpha_t > 0){
                        printf POVOUT ("make_neut_dots( %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0);
                    } else {
                        printf POVOUT ("make_neut_ion( %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0);
                    }
                } elsif ($is_solute[$k] > 0){
                    if ($is_solute[$k] == 2){
                        if ($make_neutral == 1){
                            printf POVOUT ("make_neut_ion( %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0);
                        } else {
                            printf POVOUT ("make_pos_ion( %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0);
                        }
                    } elsif ($is_solute[$k] == 3) {
                        if ($make_neutral == 1){
                            printf POVOUT ("make_neut_ion( %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0);
                        } else {
                            printf POVOUT ("make_neg_ion( %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0);
                        }
                    }
                } else {
                    if ($mb_style == 1 || $mb_style == 5){
                        if ($include_gold[$k] == 0 || $single_ion > 0){
                            if ($mb_style == 1){
                                printf POVOUT ("make_mb_black( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                            } else {
                                printf POVOUT ("make_mbnd_black( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                            }
                        } else {
                            if ($make_neutral == 1){
                                printf POVOUT ("make_mbnd_black( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                            } else {
                                printf POVOUT ("make_mb_gold( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                            }
                        }
                    } elsif ($mb_style == 2){
                        printf POVOUT ("make_mb_dipole( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                    } elsif ($mb_style == 3){
                        printf POVOUT ("make_mb_spokes( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                    } elsif ($mb_style == 6){
                        printf POVOUT ("make_mb_black_spokes( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                    } elsif ($mb_style == 4){
                        if ($make_neutral == 1){
                            printf POVOUT ("make_neut_dots( %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0);
                        } elsif ($r_radius == 0) {
                            printf POVOUT ("make_mb_dots( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                            printf POVOUT ("make_mb_orange_tabs( %f, %f, %f, %f)\n", $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                        }
                    } else {
                        printf POVOUT ("make_%i\_lobe_rose( %f, %f, %f, %f)\n", $rose_type, $configs_x[$k]-$halfx, $configs_y[$k]-$halfy, 0, $angle_thetas[$k]);
                    }
                }
            } else {
                if ($sep_dist > 0 && $k==0){
                    if ($make_neutral == 1 && $bWalls == 0){
                        printf POVOUT ("make_neut_pair( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                    } elsif ($make_neutral == 1 && $bWalls == 1){
                    # printf POVOUT ("make_channel( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                        printf POVOUT ("make_surface( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                    } else {
                        printf POVOUT ("make_ion_pair( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                    }
                } elsif ($single_ion == 1 && $k==0){
                    printf POVOUT ("make_pos_ion( %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0);
                } elsif ($single_ion == 2 && $k==0){
                    printf POVOUT ("make_neg_ion( %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0);
                } elsif ($r_radius > 0 && $k>=0 && $is_solute[$k] == 1){
                    if ($alpha_t > 0){
                        printf POVOUT ("make_neut_dots( %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0);
                    } else {
                        printf POVOUT ("make_neut_ion( %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0);
                    }
                } elsif ($is_solute[$k] > 0){
                    if ($is_solute[$k] == 2){
                        if ($make_neutral == 1){
                            printf POVOUT ("make_neut_ion( %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0);
                        } else {
                            printf POVOUT ("make_pos_ion( %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0);
                        }
                    } elsif ($is_solute[$k] == 3) {
                        if ($make_neutral == 1){
                            printf POVOUT ("make_neut_ion( %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0);
                        } else {
                            printf POVOUT ("make_neg_ion( %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0);
                        }
                    }
                
                } else {
                    if ($mb_style == 1 || $mb_style == 5){
                        if ($include_gold[$k] == 0 || $single_ion > 0){
                            if ($mb_style == 1){
                                printf POVOUT ("make_mb_black( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                            } else {
                                printf POVOUT ("make_mbnd_black( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                            }
                        } else {
                            if ($make_neutral == 1){
                                printf POVOUT ("make_mbnd_black( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                            } else {
                                printf POVOUT ("make_mb_gold( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                            }
                        }
                    } elsif ($mb_style == 2){
                        printf POVOUT ("make_mb_dipole( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                    } elsif ($mb_style == 3){
                        printf POVOUT ("make_mb_spokes( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                    } elsif ($mb_style == 6){
                        printf POVOUT ("make_mb_black_spokes( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                    } elsif ($mb_style == 4){
                        if ($make_neutral == 1){
                            printf POVOUT ("make_neut_dots( %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0);
                        } elsif ($r_radius == 0) {
                            printf POVOUT ("make_mb_dots( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                            printf POVOUT ("make_mb_orange_tabs( %f, %f, %f, %f)\n", $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                        }
                    } else {
                        printf POVOUT ("make_%i\_lobe_rose( %f, %f, %f, %f)\n", $rose_type, $configs_x[$k], $configs_y[$k], 0, $angle_thetas[$k]);
                    }
                }
            }
        }
        print POVOUT "\n";
    }
}

sub wrapCoordinates {
    if ($line[2] != 0){
        $sign = round($line[2]/abs($line[2]));    
        $box_num = floor((abs($line[2]-$halfx)/$hxx) + 0.5);
    } else {
        $sign = 1;
        $box_num = 0;
    }
    $wrap = $line[2] - ($sign*$hxx*$box_num);
    $line[2] = $wrap;
    if ($line[3] != 0){
        $sign = round($line[3]/abs($line[3]));    
        $box_num = floor((abs($line[3]-$halfy)/$hyy) + 0.5);
    } else {
        $sign = 1;
        $box_num = 0;
    }
    $wrap = $line[3] - ($sign*$hyy*$box_num);
    $line[3] = $wrap;

}

sub shiftCoordinates{
    $line[2] += $x_shift;
    $line[3] += $y_shift;
}

sub centerCoordinates{
    $line[2] -= $halfx;
    $line[3] -= $halfy;
}

sub saveCoordinates {
    $solute_flag = 0; 
    push(@configs_x, $line[2]);
    push(@configs_y, $line[3]);
    # determine the rotation angle from the quaternions
    if ($line[1] eq 'pv'){
        $q0 = 1; #$line[5];
        $q3 = 0; #$line[8];
        if ($line[0] < $ion_pairs_count){
            $solute_flag = 2;
        } elsif ($line[0] < $total_ions){
            $solute_flag = 3;
        } else {
            $solute_flag = 1; # note: we consider this a solute particle if it has no orientation
        }

    } elsif ($line[1] eq 'pvqj'){
        $q0 = $line[8];
        $q3 = $line[11];
    } else {
        #die "\nError: Can't determine rose angle for a line of type $line[1]\n\n";
        $q0 = 1;
        $q3 = 0;
    }
    $angle_val = atan2(2*$q0*$q3, 1-2*$q3*$q3);
    $angle_val = $angle_val*180/M_PI;
    push(@angle_thetas,$angle_val);
    if ($solute_flag >= 1) {
        push(@is_solute, $solute_flag); 
    } else {
        push(@is_solute, 0); 
    }

    if ($cull_wat == 1){
        if ($#configs_x > 0){
            checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);

            if ($temp_dist2 <= $dist_tol2){
                push(@include_mol, 1);
                checkIncludeGold();
            } else {
                push(@include_mol, 0);
                checkIncludeGold();
            }
        } else {
            push(@include_mol, 1);
            push(@include_gold, 0);
        }
    } else {
        push(@include_mol, 1);
        push(@include_gold, 0);
    }

    if ($wrapFlag == 0) {
# now do some imaging as needed for the clipping rendering
        if ($line[2] >= ($halfx - $buffer)) {
            push(@configs_x, $line[2]-$hxx);
            push(@configs_y, $line[3]);
            push(@angle_thetas,$angle_val);
            if ($solute_flag >= 1) {
                push(@is_solute, $solute_flag); 
            } else {
                push(@is_solute, 0); 
            }
            if ($cull_wat == 1){
                if ($#configs_x > 0){
                    checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                    if ($temp_dist2 <= $dist_tol2){
                        push(@include_mol, 1);
                    } else {
                        push(@include_mol, 0);
                    }
                    checkIncludeGold();
                    if ($line[3] >= ($halfy - $buffer)){
                        push(@configs_x, $line[2]-$hxx);
                        push(@configs_y, $line[3]-$hyy);
                        push(@angle_thetas,$angle_val);
                        if ($solute_flag >= 1) {
                            push(@is_solute, $solute_flag); 
                        } else {
                            push(@is_solute, 0); 
                        }
                        checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                        if ($temp_dist2 <= $dist_tol2){
                            push(@include_mol, 1);
                        } else {
                            push(@include_mol, 0);
                        }
                        checkIncludeGold();
                    } elsif ($line[3] <= (-$halfy + $buffer)){
                        push(@configs_x, $line[2]-$hxx);
                        push(@configs_y, $line[3]+$hyy);
                        push(@angle_thetas,$angle_val);
                        if ($solute_flag >= 1) {
                            push(@is_solute, $solute_flag); 
                        } else {
                            push(@is_solute, 0); 
                        }
                        checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                        if ($temp_dist2 <= $dist_tol2){
                            push(@include_mol, 1);
                        } else {
                            push(@include_mol, 0);
                        }
                        checkIncludeGold();
                    }
                }
            } else {
                push(@include_mol, 1);
                push(@include_gold, 0);
                if ($line[3] >= ($halfy - $buffer)){
                    push(@configs_x, $line[2]-$hxx);
                    push(@configs_y, $line[3]-$hyy);
                    push(@angle_thetas,$angle_val);
                    if ($solute_flag >= 1) {
                        push(@is_solute, $solute_flag); 
                    } else {
                        push(@is_solute, 0); 
                    }
                    push(@include_mol, 1);
                    push(@include_gold, 0);
                } elsif ($line[3] <= (-$halfy + $buffer)){
                    push(@configs_x, $line[2]-$hxx);
                    push(@configs_y, $line[3]+$hyy);
                    push(@angle_thetas,$angle_val);
                    if ($solute_flag >= 1) {
                        push(@is_solute, $solute_flag); 
                    } else {
                        push(@is_solute, 0); 
                    }
                    push(@include_mol, 1);
                    push(@include_gold, 0);
                }
            }
        }

        if ($line[2] <= (-$halfx + $buffer)) {
            push(@configs_x, $line[2]+$hxx);
            push(@configs_y, $line[3]);
            push(@angle_thetas,$angle_val);
            if ($solute_flag >= 1) {
                push(@is_solute, $solute_flag); 
            } else {
                push(@is_solute, 0); 
            }
            if ($cull_wat == 1){
                if ($#configs_x > 0){
                    checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                    if ($temp_dist2 <= $dist_tol2){
                        push(@include_mol, 1);
                    } else {
                        push(@include_mol, 0);
                    }
                    checkIncludeGold();
                    if ($line[3] <= (-$halfy + $buffer)){
                        push(@configs_x, $line[2]+$hxx);
                        push(@configs_y, $line[3]+$hyy);
                        push(@angle_thetas,$angle_val);
                        if ($solute_flag >= 1) {
                            push(@is_solute, $solute_flag); 
                        } else {
                            push(@is_solute, 0); 
                        }
                        checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                        if ($temp_dist2 <= $dist_tol2){
                            push(@include_mol, 1);
                        } else {
                            push(@include_mol, 0);
                        }
                        checkIncludeGold();
                    } elsif ($line[3] >= ($halfy - $buffer)){
                        push(@configs_x, $line[2]+$hxx);
                        push(@configs_y, $line[3]-$hyy);
                        push(@angle_thetas,$angle_val);
                        if ($solute_flag >= 1) {
                            push(@is_solute, $solute_flag); 
                        } else {
                            push(@is_solute, 0); 
                        }
                        checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                        if ($temp_dist2 <= $dist_tol2){
                            push(@include_mol, 1);
                        } else {
                            push(@include_mol, 0);
                        }
                        checkIncludeGold();
                    }
                }
            } else {
                push(@include_mol, 1);
                push(@include_gold, 0);
                if ($line[3] <= (-$halfy + $buffer)){
                    push(@configs_x, $line[2]+$hxx);
                    push(@configs_y, $line[3]+$hyy);
                    push(@angle_thetas,$angle_val);
                    if ($solute_flag >= 1) {
                        push(@is_solute, $solute_flag); 
                    } else {
                        push(@is_solute, 0); 
                    }
                    push(@include_mol, 1);
                    push(@include_gold, 0);
                } elsif ($line[3] >= ($halfy - $buffer)){
                    push(@configs_x, $line[2]+$hxx);
                    push(@configs_y, $line[3]-$hyy);
                    push(@angle_thetas,$angle_val);
                    if ($solute_flag >= 1) {
                        push(@is_solute, $solute_flag); 
                    } else {
                        push(@is_solute, 0); 
                    }
                    push(@include_mol, 1);
                    push(@include_gold, 0);
                }

            }
        }

        if ($line[3] >= ($halfy - $buffer)) {
            push(@configs_x, $line[2]);
            push(@configs_y, $line[3]-$hyy);
            push(@angle_thetas,$angle_val);
            if ($solute_flag >= 1) {
                push(@is_solute, $solute_flag); 
            } else {
                push(@is_solute, 0); 
            }
            # push(@configs_x, $line[2]);
            # push(@configs_y, $line[3]-2*$hyy);
            # push(@angle_thetas,$angle_val);
            # if ($solute_flag == 1) {
            #     push(@is_solute, 1); 
            # } else {
            #     push(@is_solute, 0); 
            # }
            if ($cull_wat == 1){
                if ($#configs_x > 0){
                    checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                    if ($temp_dist2 <= $dist_tol2){
                        push(@include_mol, 1);
                    } else {
                        push(@include_mol, 0);
                    }
                    checkIncludeGold();
                    if ($line[2] >= ($halfx - $buffer)){
                        push(@configs_x, $line[2]-$hxx);
                        push(@configs_y, $line[3]-$hyy);
                        push(@angle_thetas,$angle_val);
                        if ($solute_flag >= 1) {
                            push(@is_solute, $solute_flag); 
                        } else {
                            push(@is_solute, 0); 
                        }
                        checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                        if ($temp_dist2 <= $dist_tol2){
                            push(@include_mol, 1);
                        } else {
                            push(@include_mol, 0);
                        }
                        checkIncludeGold();
                    } elsif ($line[2] <= (-$halfx + $buffer)){
                        push(@configs_x, $line[2]+$hxx);
                        push(@configs_y, $line[3]-$hyy);
                        push(@angle_thetas,$angle_val);
                        if ($solute_flag >= 1) {
                            push(@is_solute, $solute_flag); 
                        } else {
                            push(@is_solute, 0); 
                        }
                        checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                        if ($temp_dist2 <= $dist_tol2){
                            push(@include_mol, 1);
                        } else {
                            push(@include_mol, 0);
                        }
                        checkIncludeGold();
                    }
                }
            } else {
                push(@include_mol, 1);
                push(@include_gold, 0);
                if ($line[2] >= ($halfx - $buffer)){
                    push(@configs_x, $line[2]-$hxx);
                    push(@configs_y, $line[3]-$hyy);
                    push(@angle_thetas,$angle_val);
                    if ($solute_flag >= 1) {
                        push(@is_solute, $solute_flag); 
                    } else {
                        push(@is_solute, 0); 
                    }
                    push(@include_mol, 1);
                    push(@include_gold, 0);
                } elsif ($line[2] <= (-$halfx + $buffer)){
                    push(@configs_x, $line[2]+$hxx);
                    push(@configs_y, $line[3]-$hyy);
                    push(@angle_thetas,$angle_val);
                    if ($solute_flag >= 1) {
                        push(@is_solute, $solute_flag); 
                    } else {
                        push(@is_solute, 0); 
                    }
                    push(@include_mol, 1);
                    push(@include_gold, 0);
                }
            }
        }

        if ($line[3] <= (-$halfy + $buffer)) {
            push(@configs_x, $line[2]);
            push(@configs_y, $line[3]+$hyy);
            push(@angle_thetas,$angle_val);
            if ($solute_flag >= 1) {
                push(@is_solute, $solute_flag); 
            } else {
                push(@is_solute, 0); 
            }
            if ($cull_wat == 1){
                if ($#configs_x > 0){
                    checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                    if ($temp_dist2 <= $dist_tol2){
                        push(@include_mol, 1);
                    } else {
                        push(@include_mol, 0);
                    }
                    checkIncludeGold();
                    if ($line[2] <= (-$halfx + $buffer)){
                        push(@configs_x, $line[2]+$hxx);
                        push(@configs_y, $line[3]+$hyy);
                        push(@angle_thetas,$angle_val);
                        if ($solute_flag >= 1) {
                            push(@is_solute, $solute_flag); 
                        } else {
                            push(@is_solute, 0); 
                        }
                        checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                        if ($temp_dist2 <= $dist_tol2){
                            push(@include_mol, 1);
                        } else {
                            push(@include_mol, 0);
                        }
                        checkIncludeGold();
                    } elsif ($line[2] >= ($halfx - $buffer)){
                        push(@configs_x, $line[2]-$hxx);
                        push(@configs_y, $line[3]+$hyy);
                        push(@angle_thetas,$angle_val);
                        if ($solute_flag >= 1) {
                            push(@is_solute, $solute_flag); 
                        } else {
                            push(@is_solute, 0); 
                        }
                        checkCull($configs_x[$#configs_x], $configs_y[$#configs_y]);
                        if ($temp_dist2 <= $dist_tol2){
                            push(@include_mol, 1);
                        } else {
                            push(@include_mol, 0);
                        }
                        checkIncludeGold();
                    }
                }
            } else {
                push(@include_mol, 1);
                push(@include_gold, 0);

                if ($line[2] <= (-$halfx + $buffer)){
                    push(@configs_x, $line[2]+$hxx);
                    push(@configs_y, $line[3]+$hyy);
                    push(@angle_thetas,$angle_val);
                    if ($solute_flag >= 1) {
                        push(@is_solute, $solute_flag); 
                    } else {
                        push(@is_solute, 0); 
                    }
                    push(@include_mol, 1);
                    push(@include_gold, 0);
                } elsif ($line[2] >= ($halfx - $buffer)){
                    push(@configs_x, $line[2]-$hxx);
                    push(@configs_y, $line[3]+$hyy);
                    push(@angle_thetas,$angle_val);
                    if ($solute_flag >= 1) {
                        push(@is_solute, $solute_flag); 
                    } else {
                        push(@is_solute, 0); 
                    }
                    push(@include_mol, 1);
                    push(@include_gold, 0);
                }
            }
        }
    }
}

sub checkCull {
    $temp_x = $_[0]-$configs_x[0]-$sep_dist;
    $temp_x2 = $_[0]-$configs_x[0]+$sep_dist;
    $temp_y = $_[1]-$configs_y[0];
    $temp_dist1 = $temp_x*$temp_x + $temp_y*$temp_y;
    $temp_dist2 = $temp_dist1;
    $temp_dist12 = $temp_x2*$temp_x2 + $temp_y*$temp_y;
    $temp_dist2 = $temp_dist12 if $temp_dist12 < $temp_dist2;
#    $avg_rad = 0.5*($n_radius+$p_radius);
#    $avg_rad2 = $avg_rad*$avg_rad;
    $small_rad = $p_radius;
    $small_rad = $n_radius if $small_rad > $n_radius;
    $avg_rad3 = ($small_rad+1.5)**2;
}

sub checkIncludeGold {
    if ((abs($temp_y) <= $avg_rad3) && ($temp_dist1 + $temp_dist12) < 2*($sep_dist*$sep_dist + $avg_rad3)){
        push(@include_gold, 1);
    } else {
        push(@include_gold, 0);
    }
}

sub wipeArrays {
    $#configs_x = -1;
    $#configs_y = -1;
    $#angle_thetas = -1;
    $#include_mol = -1;
    $#include_gold = -1;
    $#is_solute = -1;
    $solute_count = 0;
    $#hbonds = -1;
    $#five_member_rings = -1;
    $#six_member_rings = -1;
    $#seven_member_rings = -1;
}

sub round {
    return int( $_[0] + 0.5 * ($_[0] <=> 0) );
}

sub buildHBondMap {
    $hbond_len2_tol = $hbond_len_tol*$hbond_len_tol;
    #print STDERR "angle thetas: $#angle_thetas\n";

    $start = 0;
    $start = 1 if defined($single_ions) && $single_ions > 0;
    $start = 1 if defined($sep_dist) && $sep_dist > 0;
    #print STDERR "start #: $start\n";

    for ($k=0; $k<=$#angle_thetas; $k++) {
        push @{ $hbonds[$k] }, 0;
    }
    # loop over all waters and identify HBonded neighbors
#    $count = 0;
    for ($k=$start; $k<=$#angle_thetas; $k++) {
#        for ($l=$start; $l<=$#angle_thetas; $l++){
        for ($l=$k+1; $l<=$#angle_thetas; $l++){
            if ($l != $k){
                # test if distance and angle tolerances pass    
                $diff_x = $configs_x[$l]-$configs_x[$k];
                $diff_y = $configs_y[$l]-$configs_y[$k];
                $dist2 = $diff_x*$diff_x + $diff_y*$diff_y;
#print STDERR "dist2 = $dist2 \n";
                if ($dist2 <= $hbond_len2_tol && $dist2 > $tolerance){
#print STDERR "$count\n"; $count++;
                    $ang_2_j = atan2($diff_y,$diff_x)*180/M_PI + 30.0;
#print STDERR "ang_2j = $ang_2_j \n";
                    $ang_2_i = $ang_2_j + 180;
#print STDERR "ang_2i = $ang_2_i \n";
                    $ang_2_j = $ang_2_j + 360 if $ang_2_j < 0;
#print STDERR "ang_2j = $ang_2_j \n";
                    $mod_val1 = int($ang_2_j / 120); 
                    $mod_val2 = int($ang_2_i / 120); 

                    $ang_i = $angle_thetas[$k];
                    $ang_i += 360.0 if $ang_i < 0;
                    $ang_j = $angle_thetas[$l];
                    $ang_j += 360.0 if $ang_j < 0;
                    $mod_val3 = int($ang_i / 120.0);
                    $mod_val4 = int($ang_j / 120.0);
#print STDERR "theta_i: $angle_thetas[$k] ($ang_i)  theta_j: $angle_thetas[$l]($ang_j)\n";
#print STDERR "mods 1 2 3 4: $mod_val1 $mod_val2 $mod_val3 $mod_val4\n";

                    $ang_val1 = $ang_2_j - $mod_val1*120.0;
                    $ang_val2 = $ang_2_i - $mod_val2*120.0;
                    $ang_val3 = $ang_i - $mod_val3*120.0;
                    $ang_val4 = $ang_j - $mod_val4*120.0;

                    #$ang_val1 = $ang_val1 - 120.0 if $ang_val1 > 60.0;
                    #$ang_val2 = $ang_val2 - 120.0 if $ang_val2 > 60.0;
                    #$ang_val3 = $ang_val3 - 120.0 if $ang_val3 > 60.0;
                    #$ang_val4 = $ang_val4 - 120.0 if $ang_val4 > 60.0;
                    
#print STDERR "angs 1 2 3 4: $ang_val1 $ang_val2 $ang_val3 $ang_val4\n";
                    $ang_diff1 = abs(($ang_val1 - $ang_val3));
                    $ang_diff2 = abs(($ang_val2 - $ang_val4));
                    #$ang_diff3 = abs(($ang_val3 - $ang_val1));
                    #$ang_diff4 = abs(($ang_val4 - $ang_val2));
                    #$ang_diff1 -= 60.0 if $ang_diff1 > 60;
                    #$ang_diff2 -= 60.0 if $ang_diff2 > 60;
                    #$ang_diff1 += 60.0 if $ang_diff1 < -60;
                    #$ang_diff2 += 60.0 if $ang_diff2 < -60;
#print STDERR "av1: $ang_val1, av2: $ang_val2, av3: $ang_val3, av4: $ang_val4\n";
#print STDERR "angle diff 1: $ang_diff1, angle diff 2: $ang_diff2\n";
#print STDERR "angle diff 3: $ang_diff3, angle diff 4: $ang_diff4\n";
                    if (($ang_diff1 <= $hbond_ang_tol || $ang_diff1 > $hbond_ang_tol2) && ($ang_diff2 <= $hbond_ang_tol || $ang_diff2 > $hbond_ang_tol2)){
                        #if (($ang_diff1 <= $hbond_ang_tol || $ang_diff3 <= $hbond_ang_tol) && ($ang_diff2 <= $hbond_ang_tol || $ang_diff4 <= $hbond_ang_tol)){
#$diff_diff = abs($ang_diff1 - $ang_diff2);
#if ($diff_diff <= $hbond_ang_tol){ 
#print STDOUT "passed!: $k $l $ang_diff1 $dist2 $configs_x[$l] $configs_x[$k] $configs_y[$l] $configs_y[$k] \n";
                        push @{ $hbonds[$k] }, $l;
                        push @{ $hbonds[$l] }, $k;
                        # uncomment below to show hbonds
                        print POVOUT "// k: $k, l: $l\n";
                        print POVOUT "intersection{union{cylinder{<$configs_x[$l],$configs_y[$l],-0.5>, <$configs_x[$k],$configs_y[$k],-0.5>, 0.35\ntexture{pigment{ rgb <1.0,0.75,0.2> }}} sphere{<$configs_x[$k],$configs_y[$k],-0.5>,0.35\ntexture{pigment{ rgb <1.0,0.75,0.2> }}} sphere{<$configs_x[$l],$configs_y[$l],-0.5>,0.35\ntexture{pigment{ rgb <1.0,0.75,0.2> }}}} box{<-1000,-1000,-1.0>,<1000,1000,-0.5>\ntexture{pigment{ rgb <1.0,0.75,0.2> }}}}\n";
                    }
                }
            }
        }
    }
    # for ($k=0; $k<=$#angle_thetas; $k++){
    #     $array_len1 = $#{$hbonds[$k]};
    #     print "$k ($array_len1): ";
    #     for ($l=0; $l<=$array_len1; $l++){
    #         print "$hbonds[$k][$l] ";
    #     }
    #     print "\n";
    # }
    #identifyNucleationRings();

    identifyNucleationRings();
}

sub identifyNucleationRings {
    # loop over water and identify HBond rings
    # This could be slow.  Worst case scenario: if each had 5 neighbors, there would be ~80K "tests" per particle.
    $start = 0;
    $start = 1 if defined($single_ions) && $single_ions > 0;
    $start = 1 if defined($sep_dist) && $sep_dist > 0;

    for ($k=$start; $k<=$#angle_thetas; $k++){
        $array_len1 = $#{$hbonds[$k]};
        if ($array_len1 > 0){
            for ($l=1; $l<=$array_len1; $l++){  
                $curr_l_mol = $hbonds[$k][$l];
                if ($curr_l_mol == $k){
                    next;
                } else {
                    $array_len2 = $#{$hbonds[$curr_l_mol]};
                    # print "array_len2: $array_len2 ; l_mol: $curr_l_mol\n" if $k == 0;
                }
                if ($array_len2 > 0){
                    for ($m=1; $m<=$array_len2; $m++){  
                        $curr_m_mol = $hbonds[$curr_l_mol][$m];
                        if ($curr_m_mol == $k || $curr_m_mol == $curr_l_mol){
                            next;
                        } else {
                            $array_len3 = $#{$hbonds[$curr_m_mol]};
                            # print "array_len3: $array_len3 ; m_mol: $curr_m_mol\n" if $k == 0;
                        }
                        if ($array_len3 > 0){
                            for ($n=1; $n<=$array_len3; $n++){  
                                $curr_n_mol = $hbonds[$curr_m_mol][$n];
                                if ($curr_n_mol == $k || $curr_n_mol == $curr_l_mol || $curr_n_mol == $curr_m_mol){
                                    next;
                                } else {
                                    $array_len4 = $#{$hbonds[$curr_n_mol]};
                                }
                                if ($array_len4 > 0){
                                    for ($o=1; $o<=$array_len4; $o++){  
                                        $curr_o_mol = $hbonds[$curr_n_mol][$o];
                                        if ($curr_o_mol == $k || $curr_o_mol == $curr_l_mol || $curr_o_mol == $curr_m_mol || $curr_o_mol == $curr_n_mol){
                                            next;
                                        } else {
                                            $array_len5 = $#{$hbonds[$curr_o_mol]};
                                        }
                                        if ($array_len5 > 0){
                                            for ($p=1; $p<=$array_len5; $p++){
# test if we have a 5 member ring   
                                                #print "test 6: $curr_l_mol, $curr_m_mol, $curr_n_mol, $curr_o_mol, $hbonds[$o][$p]\n";
                                                if ($hbonds[$curr_o_mol][$p] == $k){
                                                    $val0 = $k;
                                                    $val1 = $curr_l_mol;
                                                    $val2 = $curr_m_mol;
                                                    $val3 = $curr_n_mol;
                                                    $val4 = $curr_o_mol;

                                                    #if ($val0 < $val1 && $val0 != $val2 && $val0 != $val3 && $val0 != $val4 && $val1 != $val2 && $val1 != $val3 && $val1 != $val4 && $val2 != $val3 && $val2 != $val4 && $val3 != $val4){
                                                        @ring_array = ($val0, $val1, $val2, $val3, $val4);
                                                        push(@five_member_rings, [@ring_array]);
                                                        #                                
#                print STDOUT "@ring_array\n";
                                                        #}
                                                }
# continue the nested search
                                                $curr_p_mol = $hbonds[$curr_o_mol][$p];
                                                if ($curr_p_mol == $k || $curr_p_mol == $curr_l_mol || $curr_p_mol == $curr_m_mol || $curr_p_mol == $curr_n_mol || $curr_p_mol == $curr_o_mol){
                                                    next;
                                                } else {
                                                    $array_len6 = $#{$hbonds[$curr_p_mol]};
                                                }
                                                if ($array_len6 > 0){
                                                    for ($q=1; $q<=$array_len6; $q++){
# test if we have a 6 member ring   
                                                        if ($hbonds[$curr_p_mol][$q] == $k){
                                                            $val0 = $k;
                                                            $val1 = $curr_l_mol;
                                                            $val2 = $curr_m_mol;
                                                            $val3 = $curr_n_mol;
                                                            $val4 = $curr_o_mol;
                                                            $val5 = $curr_p_mol;
                                                            #    if ($val0 < $val1 && $val0 != $val2 && $val0 != $val3 && $val0 != $val4 && $val0 != $val5 && $val1 != $val2 && $val1 != $val3 && $val1 != $val4 && $val1 != $val5 && $val2 != $val3 && $val2 != $val4 && $val2 != $val5 && $val3 != $val4 && $val3 != $val5 && $val4 != $val5){
                                                                @ring_array = ($val0, $val1, $val2, $val3, $val4, $val5);
                                                                push(@six_member_rings, [@ring_array]);
                                                                # print "@ring_array\n";
                                                                #}
                                                        }
# continue the nested search
                                                        $curr_q_mol = $hbonds[$curr_p_mol][$q];
                                                        if ($curr_q_mol == $k || $curr_q_mol == $curr_l_mol || $curr_q_mol == $curr_m_mol || $curr_q_mol == $curr_n_mol || $curr_q_mol == $curr_o_mol || $curr_q_mol == $curr_p_mol){
                                                            next;
                                                        } else {
                                                            $array_len7 = $#{$hbonds[$curr_q_mol]};
                                                        }
                                                        if ($array_len7 > 0){
                                                            for ($r=1; $r<=$array_len7; $r++){
# test if we have a 7 member ring   
                                                                if ($hbonds[$curr_q_mol][$r] == $k){
                                                                    $val0 = $k;
                                                                    $val1 = $curr_l_mol;
                                                                    $val2 = $curr_m_mol;
                                                                    $val3 = $curr_n_mol;
                                                                    $val4 = $curr_o_mol;
                                                                    $val5 = $curr_p_mol;
                                                                    $val6 = $curr_q_mol;
                                                                    #    if ($val0 < $val1 && $val0 != $val2 && $val0 != $val3 && $val0 != $val4 && $val0 != $val5 && $val0 != $val6 && $val1 != $val2 && $val1 != $val3 && $val1 != $val4 && $val1 != $val5 && $val1 != $val6 && $val2 != $val3 && $val2 != $val4 && $val2 != $val5 && $val2 != $val6 && $val3 != $val4 && $val3 != $val5 && $val3 != $val6 && $val4 != $val5 && $val4 != $val6 && $val5 != $val6){
                                                                        @ring_array = ($val0, $val1, $val2, $val3, $val4, $val5, $val6);
                                                                        push(@seven_member_rings, [@ring_array]);
                                                                        #}
                                                                }
# end of the nested search
                                                            }
                                                        }
                                                    }
                                                }
                                            }   
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

# make 5 member polygons
    $five_ring_num = @five_member_rings;
    # print STDOUT "number of 5s: $five_ring_num\n";
    for ($k=0; $k<$five_ring_num; $k++) {
        print POVOUT "polygon{ 5, <$configs_x[$five_member_rings[$k][0]],$configs_y[$five_member_rings[$k][0]]>,<$configs_x[$five_member_rings[$k][1]],$configs_y[$five_member_rings[$k][1]]>,<$configs_x[$five_member_rings[$k][2]],$configs_y[$five_member_rings[$k][2]]>,<$configs_x[$five_member_rings[$k][3]],$configs_y[$five_member_rings[$k][3]]>,<$configs_x[$five_member_rings[$k][4]],$configs_y[$five_member_rings[$k][4]]> texture{pigment{ rgb <0.6,0.6,1.0> }} translate<0,0,-0.505>}\n";
    }
# make 6 member polygons
    $six_ring_num = @six_member_rings;
    # print STDOUT "number of 6s: $six_ring_num\n";
    for ($k=0; $k<$six_ring_num; $k++) {
        print POVOUT "polygon{ 6, <$configs_x[$six_member_rings[$k][0]],$configs_y[$six_member_rings[$k][0]]>,<$configs_x[$six_member_rings[$k][1]],$configs_y[$six_member_rings[$k][1]]>,<$configs_x[$six_member_rings[$k][2]],$configs_y[$six_member_rings[$k][2]]>,<$configs_x[$six_member_rings[$k][3]],$configs_y[$six_member_rings[$k][3]]>,<$configs_x[$six_member_rings[$k][4]],$configs_y[$six_member_rings[$k][4]]>,<$configs_x[$six_member_rings[$k][5]],$configs_y[$six_member_rings[$k][5]]> texture{pigment{ rgb <0.6,0.75,1.0> }} translate<0,0,-0.51>}\n";
    }
# make 7 member polygons
    $seven_ring_num = @seven_member_rings;
    # print STDOUT "number of 7s: $seven_ring_num\n";
    for ($k=0; $k<$seven_ring_num; $k++) {
        print POVOUT "polygon{ 7, <$configs_x[$seven_member_rings[$k][0]],$configs_y[$seven_member_rings[$k][0]]>,<$configs_x[$seven_member_rings[$k][1]],$configs_y[$seven_member_rings[$k][1]]>,<$configs_x[$seven_member_rings[$k][2]],$configs_y[$seven_member_rings[$k][2]]>,<$configs_x[$seven_member_rings[$k][3]],$configs_y[$seven_member_rings[$k][3]]>,<$configs_x[$seven_member_rings[$k][4]],$configs_y[$seven_member_rings[$k][4]]>,<$configs_x[$seven_member_rings[$k][5]],$configs_y[$seven_member_rings[$k][5]]>,<$configs_x[$seven_member_rings[$k][6]],$configs_y[$seven_member_rings[$k][6]]>  texture{pigment{ rgb <0.6,0.9,1.0> }} translate<0,0,-0.515>}\n";
    }

}

