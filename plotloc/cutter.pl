#!/usr/bin/perl

# REQUIREMENTS:
# You'll need the Perl module Image::Magick to run this script.
# If you don't have it, but have adminstrator access on your
# system, you can install this from the command line using:
# 	perl -MCPAN -e shell
# Once you answered any CPAN set up questions, at the prompt
# you should type:
# 	install Image::Magick
#
#
# USAGE:
# perl cutter.pl file="sky.png" minzoom=3 maxzoom=4
#
#
# OPTIONS:
# tilesize    - default = 256
# minzoom     - default = 1
# maxzoom     - default is calculated from image size
# ext         - default is "jpg"
# compression - default is'JPEG'
# quality     - default is 85

use Image::Magick;


# Set some options
getArgs();
$path = $data{'file'} ? $data{'file'} : "sky.png";


# Build the output directory name
$dir = $path;
$dir =~ s/.[^.]*$//g;
$dir .= '-tiles';
print "Saving tiles to $dir.\n";
mkdir "$dir", 0755 unless -d "$dir";



# Read in the input file
my($image, $x);
$image = Image::Magick->new;
$x = $image->Read($path);
warn "$x" if "$x";

$wide = $image->Get('width');
$tall = $image->Get('height');
$maxdim = ($wide < $tall) ? $tall : $wide;


# Now that we know the size we can set some more arguments
$tileSize = $data{'tilesize'} ? $data{'tilesize'} : 256;
$minZoom = $data{'minzoom'} ? $data{'minzoom'} : 1;
$maxZoom = $data{'maxzoom'} ? $data{'maxzoom'} : int(log($wide/$tileSize)/log(2));
$fileExt = $data{'ext'} ? $data{'ext'} : "jpg";
$compres = $data{'compression'} ? $data{'compression'} : 'JPEG';
$quality = $data{'quality'} ? $data{'quality'} : 85;
$verbose = $data{'verbose'} ? $data{'verbose'} : 0;


# Here we assume the image is the entire sky. Could be adapted
$map = Image::Magick->new(size=>$maxdim.'x'.$maxdim);
$map->ReadImage('xc:black');
$offx = 0;
$offy = 0;
$offx = ($wide > $tall) ? 0 : ($maxdim-$wide)/2.;
$offy = ($wide > $tall) ? ($maxdim-$tall)/2.: 0;
$map->Composite(image=>$image,compose=>'over',x=>$offx,y=>$offy);

if($verbose){print "The original image is ".$wide."x".$tall."px. Expanding to ".$maxdim."x".$maxdim."\n";}

for ($z = $minZoom; $z <= $maxZoom ; $z++){

	if($verbose){int "Creating zoom level $z:\n";}
	$xpix = 2**$z;
	$ypix = 2**$z;

	# Scale the original to the appropriate size for the tiles
	$zoommap = $map->Clone();
	if($verbose){print "Resizing to ".($xpix*$tileSize)."x".($ypix*$tileSize)."px\n";}
	$zoommap->Resize(width=>$xpix*$tileSize,height=>$ypix*$tileSize);

	# Use the Transform() command to create an array of tiles of the appropriate size
	$tiles = $zoommap->Transform(crop=>($tileSize).'x'.($tileSize));

	# Loop over the tiles
	for ($i = 0; $tiles->[$i]; $i++){
		$x = ($i % $xpix);
		$y = int($i / $xpix);
		# Get the appropriate filename
		$file = getTile($x,$y,$z).".".$fileExt;
		# print "$i => $dir/$file\n";
		$tiles->[$i]->Write(filename=>$dir.'/'.$file, compression=>$compres,quality=>$quality);
	}

}




sub getTile() {
	my ($x,$y,$z) = @_;
	my $pixels = (2**$z);
	my $x=($x+$pixels)%($pixels);
	my $y=($y+$pixels)%($pixels);
	#print "$pixels\n";
	my $f="t";
	my $g;
	for($g=0 ; $g < $z ; $g++){
		$pixels=$pixels/2;
		if($y<$pixels){
			if($x<$pixels){$f.="q";}
			else{$f.="r";$x-=$pixels;}
		}else{
			if($x<$pixels){$f.="t";$y-=$pixels;}
			else{$f.="s";$x-=$pixels;$y-=$pixels;}
		}
	}
	return $f;
}


sub getArgs {

	local($temp,@pairs,$field,$value,@fields,$idx,%numbers,$my);
	@fields = @_;

	foreach $my (@ARGV){
		$temp .= $my."\&";
	}

	# no point going further
	if(!$temp){ return 1; }

	@pairs = split("&", $temp);

	# For each name-value pair:
	foreach $pair (@pairs) {
		# Split the pair up into individual variables.
		($field, $value) = split("=", $pair);

		#$value = &parseQuery($value);

		$data{$field} = $value;
	}
	return 0;
}
