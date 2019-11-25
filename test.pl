use strict;
use warnings;
use diagnostics;

use feature 'say';

use feature "switch";

use v5.14.4; # Make sure we have the correct version



my $new_sample="123";
my $new_lane="5";
my $RGID=join("_",$new_sample,$new_lane);

print "$RGID\n"
