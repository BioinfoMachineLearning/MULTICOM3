#!/usr/bin/perl -w
 use FileHandle; # use FileHandles instead of open(),close()
 use Cwd;
 use Cwd 'abs_path';

######################## !!! customize settings here !!! ############################
#																					#
# Set directory of DeepRank databases and tools								        #


######################## !!! End of customize settings !!! ##########################

######################## !!! Don't Change the code below##############


$install_dir = getcwd;
$install_dir=abs_path($install_dir);


if(!-s $install_dir)
{
	die "The installation directory ($install_dir) is not existing, please revise the customize settings part inside the configure.pl.\n";
}

if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
        $install_dir .= "/";
}


print "checking whether the configuration file run in the installation folder ...";
$cur_dir = `pwd`;
chomp $cur_dir;
$configure_file = "$cur_dir/configure.pl";
if (! -f $configure_file || $install_dir ne "$cur_dir/")
{
        die "\nPlease check the installation directory setting and run the configure program under the installation directory.\n";
}
print " OK!\n";



if (! -d $install_dir)
{
	die "can't find installation directory.\n";
}
if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
	$install_dir .= "/"; 
}

$option_list = "$install_dir/configure_scripts_list";

if (! -f $option_list)
{
        die "\nOption file $option_list not exists.\n";
}

configure_file($option_list);

print "#########  Configuring scripts, done\n\n";

sub configure_file 
{
	my ($option_list) = @_;
	open(IN, $option_list) || die "Failed to open file $option_list\n";
	$file_indx = 0;
	while (<IN>) 
	{
		$file = $_;
		chomp $file;
		@tmparr = split('/', $file);
		$filename = pop @tmparr;
		chomp $filename;
		$filepath = join('/', @tmparr);
		$option_default = $install_dir . $filepath . '/.' . $filename . '.default';
		$option_new = $install_dir . $file;
		$file_indx++;
		print "$file_indx: Configuring $option_new\n";
		if (!-f $option_default) 
		{
			die "\nOption file $option_default not exists.\n";
		}

		open(IN1, $option_default) || die "Failed to open file $option_default\n";
		open(OUT1, ">$option_new") || die "Failed to open file $option_new\n";
		while (<IN1>) 
		{
			$line = $_;
			chomp $line;
			if (index($line, 'SOFTWARE_PATH') >= 0) 
			{
				$line =~ s/SOFTWARE_PATH/$install_dir/g;
				$line =~ s/\/\//\//g;
				#print OUT1 $line."\n";
			}
			print OUT1 $line . "\n";
		}
		close IN1;
		close OUT1;
	}
	close IN;
}
