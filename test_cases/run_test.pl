#!/usr/bin/perl


#################################################################################################
###                                 James T. MacDonald 2013                                   ###
#################################################################################################


#
# Runs quick test examples for PD2
#

$user = getlogin();
$user_home = $ENV{HOME};
$host=`uname -n`;
chomp $host;
$hg_id = `hg id`;
chomp $hg_id;
$date = `date`;
chomp $date;

$debug_mode = 1;              #debug mode 1=yes 0=no

if ($#ARGV != 0) {
    select STDERR;
    print "ERROR: Wrong number of arguments\n";
    print "\nUsage: <test directory>\n\n";
    print "Runs tests\n";
    exit;
}


$dir = $ARGV[0];




$cmd = `cat $dir/cmd`;
chomp $cmd;
$log = "$dir/run.log";
$cmd = $cmd." > $log 2>&1";

print "running command : $cmd\n";
print `$cmd`;

$elapsed = `grep elapsed $log | cut -f 5 -d" "`;
chomp $elapsed;

$error = `grep -i error $log `;
chomp $error;

if ($error ne ""){
    $error = "ERROR";
}
else {
    $error = "RUN_OK";
}

$outline = "$dir\t$elapsed\t$host\t$user\t$hg_id\t$date\t$error\n";
print $outline;

$resultfile = ">>$dir/results.txt";

open (OUTFILE, $resultfile);
print OUTFILE $outline;


