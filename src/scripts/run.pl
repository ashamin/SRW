use strict;
use warnings;

chdir('..');
system('cmake CMakeLists.txt -B../build');
chdir('../build');
system('make');
system('./src.exe');