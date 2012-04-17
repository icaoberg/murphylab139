#!/usr/bin/perl -w 
for ($i=1; $i<=6; $i++)
{
    $command = ('matlab5 -display 0 < svm_new_rbf(\'SLF' . $i . '\') > temp ');
    print STDOUT $command;
    system($command);
}
