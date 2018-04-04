

`echo > src/prodart_env/prodart_env_var.h`;

open VARFILE, ">src/prodart_env/prodart_env_var.h"
    or die "ERROR: Can't open file\n\n";


$svn_rev = "";
chomp $svn_rev;

$hg_rev = "";
chomp $hg_rev;

$comp_date = `date`;
chomp $comp_date;

$pd2_ver = `cat src/Version.txt`;
chomp $pd2_ver;

print VARFILE "//DO NOT INCLUDE THIS FILE OR ADD THIS TO ANY REPOSITORY. IT IS DESIGNED FOR VARIABLES ONLY AND IS CREATED BY THE SCRIPT: src/setup.pl\n";
print VARFILE "#ifdef USING_SCONS_BUILDER\n";
print VARFILE "#define SVN_REV \"".$svn_rev."\"\n";
print VARFILE "#define HG_REV \"".$hg_rev."\"\n";
print VARFILE "#define _COMPILE_DATE_ \"".$comp_date."\"\n";
print VARFILE "#define _PRODART_VERSION_ \"".$pd2_ver."\"\n";
print VARFILE "#endif\n";

