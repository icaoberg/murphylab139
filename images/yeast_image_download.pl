#! /usr/bin/perl -s
###############################
# This script is to download the yeast image archives to the currect directory
###############################

use URI::URL;
use LWP::Simple qw(getstore);

open(FILE, "yeast_image_download.txt");
@lines = <FILE>;
foreach(@lines)
{
    $url=new URI::URL($_);
    $name2 = substr($_,33,80);
    print ("$name2\n");
    getstore($url,$name2);
}
close(FILE);
