#!/usr/bin/perl
#
# Usage: sfmakedepend [-s] [-e ext] [-f file] [-I incdir]
#                     [-m mod_ext] [-c] [-p] [-g] [-h] [-o obj_ext]
#                     [-a add_ext] file ...
#
# This is a makedepend script for Fortran, including Fortran 90.
# It searches for Fortran style includes, C preprocessor includes,
# and module dependencies to the extent that I understand them.
#
# Your files must have an extension listed in the @suffixes list
# below.  You might also want to modify $compile_string if you use
# the -c or -p option.  I call the compiler $(CFT) for historical
# reasons.
#
# The switch [-s] is for stupid Fortran compilers that don't know
# how to automatically send things through the C preprocessor.
# It is possible to force 'make' to invoke 'cpp' to create a .f
# file from a .F file (which has cpp directives), but make won't
# know that a .f file will depend on the files that the .F file
# included.  This option will provide those dependencies.
#
# The [-e ext] switch is used with the [-s] switch for compilers
# which expect an extension other than .f on source files.  For
# instance, for the Connection Machine one would use "-e fcm".
#
# The [-f file] switch is used to change the name of the current
# Makefile.
#
# The [-I incdir] option tells sfmakedepend to look in alternate
# directories for the include files.  There can be several "-I dir"
# options used at once.  The current directory is still searched
# first.
#
# The [-m mod_ext] option tells sfmakedepend what extension to
# use for Fortran 90 module files.  The default for "use My_Mod"
# is to list the dependency as "my_mod.mod" since this is what
# NAG f90 and IBM xlf both use.  Let me know if other compilers
# use a different filename for the module information.
#
# The [-c] option specifies that the Cray compiler is being used.
# This compiler requires a "-p file.o" option for each object file
# that contains a module used by your current source file.
#
# The [-p] option specifies that the Parasoft compiler is being used.
# This compiler requires a "-module file.o" option for each object file
# that contains a module used by your current source file.
#
# The [-g] option specifies that the SGI compiler is being used.
# This compiler names the module file in uppercase with the extension
# .kmo.
#
# The [-h] option specifies that the HP compiler is being used.
# This compiler names the module file in uppercase with the extension
# .mod (added by Patrick Jessee who also fixed an include bug).
#
# The [-o obj_ext] option tells sfmakedepend what extension to use for
# object files.  The default is "o", but "obj", for instance, is
# appropriate on MS-DOG etc.  This option was added by Dave Love.
#
# The [-a add_ext] option (also added by Dave Love) tells sfmakedepend
# to add targets with extension add_ext to the rules for object files.
# For instance, to operate with (f77) ftnchek .prj files, you could use
# `-a prj' to get rules like:
# foo.prj foo.o: ...
#
# The final arguments contain the list of source files to be
# searched for dependencies.
#
# EXAMPLE
#       sfmakedepend -I /usr/local/include *.F
#
# NOTES
#	This makedepend script is my first attempt at using perl 5
#	objects.  Therefore, it may not be the best example of how
#	to do this.  Also, it requires perl 5 and will die if you
#	to use it with an older perl.  The latest version is
#	available from:
#
#		http://marine.rutgers.edu/po/perl.html
#		ftp://ahab.rutgers.edu/pub/perl/sfmakedepend
#
#	Fortran 90 introduces some interesting dependencies.  Two
#	compilers I have access to (NAG f90 and IBM xlf) produce a
#	private "mod_name.mod" file if you define "module mod_name"
#	in your code.  This file is used by the compiler when you
#	use the module as a consistency check (type-safe).  On the
#	other hand, the Cray and Parasoft compilers store the module
#	information in the object file and then files which use the
#	modules need to be compiled with extra flags pointing to the
#	module object files.
#
#	This script assumes that all the files using and defining
#	modules are in the same directory and are all in the list of
#	files to be searched.  It seems that the industry has not
#	settled on a practical way to deal with a separate modules
#	directory, anyway.
#
#	I sometimes include non-existent files as a compile time
#	consistency check:
#
#	   #ifndef PLOTS
#	   #include "must_define_PLOTS"       /* bogus include */
#	   #endif
#
#	This program warns about include files it can't find, but
#	not if there is a "bogus" on the same line.
#
#     *	The f90 module dependencies can confuse some versions of
#	make, especially of the System V variety.  We use gnu
#	make because it has no problems with these dependencies.
#
# BUGS
#	It can sometimes produce duplicate dependencies.
#
#	It treats C preprocessor includes the same as Fortran
#	includes.  This can add unnecessary dependencies if you
#	use the -s flag and both kinds of includes.
#
#       Please let me know if you find any others.
#	Kate Hedstrom
#	kate@ahab.rutgers.edu
#

package source_file;

# hashes containing names of included files, modules in files
%inc_files = ();
%main::mod_files = ();

# Constructor
sub new {
    my $type = shift;
    my $filename = shift;
    my $path = shift;
    my $self = {};
    $self->{'source_file'} = $filename;
    $self->{'filepath'} = $path;
    $self->{'includes'} = {};
    $self->{'uses'} = {};
    $self->{'modules'} = {};
    bless $self;
}

sub find_includes {
    my $self = shift;
    my $file = $self->{'filepath'};
    my($after, $filepath, $ref, $included, $use, $modname);
    local(*FILE);
    local($_);

    if (-f $file) {
        open(FILE, $file) || warn "Can't open $file: $!\n";
    } elsif (-f "RCS/$file,v" || -f "$file,v" ) {
        system("co $file");
        open(FILE, $file) || warn "Can't open $file: $!\n";
        $main::rcs{$file} = 1;
    } else {
	return;
    }
    while (<FILE>) {
	$included = "";
	$use = "";
        # look for Fortran style includes
        if (/^\s*include\s*['"]([^"']*)["']/i) {
            $included = $1;
	    $after = $';
	# C preprocessor style includes
        } elsif (/^#\s*include\s*["<]([^">]*)[">]/) {
            $included = $1;
	    $after = $';
        # Fortran 90 "use"
        } elsif (/^\s*use\s+(\w+)/i) {   
	    $use = $1;
# Make the module name lowercase except for SGI & HP.
# May be compiler dependent!!
	    if ($main::sgi || $main::hp) {
		$use = uc($use);
	    } else {
		$use = lc($use);
	    }
	    $self->{'uses'}{$use} = 1;
	# Fortran 90 module
	} elsif (/^\s*module\s+(\w+)/i) {
	    $modname = $1;
	    if ($main::sgi || $main::hp) {
		$modname = uc($modname);
	    } else {
		$modname = lc($modname);
	    }
	    unless (lc($modname) eq "procedure") {
		$main::mod_files{$modname} = $file;
		$self->{'modules'}{$modname} = 1;
	    }
	}
	if ($included) {
	    if ( $inc_files{$included} ) {
		$filepath = $inc_files{$included}{'filepath'};
	    } else {
                $filepath = &main::findfile($included);
		$ref = new source_file($included, $filepath);
                $inc_files{$included} = $ref;
# Search included file for includes
		$ref->find_includes();
	    }
            if ( $filepath ) {
		$self->{'includes'}{$included} = 1;
	    } else {
                if ($after !~ /bogus/i) {
		    warn "Can't find file: $included\n";
		}
	    }
	}
    }
    close FILE;
}

sub print_includes {
    my $self = shift;
    my $target = shift;
    my $len_sum = shift;
    my($file, $ref);
    my %printed = ();

    foreach $file (keys %{$self->{'includes'}}) {
	next if $printed{$file};
	$ref = $inc_files{$file};
	my $len = length($ref->{'filepath'}) + 1;
	if (($len_sum + $len > 80) &&
	       (length($target) + 1 < $len_sum)) {
	    print "\n$target:";
	    $len_sum = length($target) + 1;
	}
	print " " . $ref->{'filepath'};
	$printed{$file} = 1;
	$len_sum += $len;
	$len_sum = $ref->print_includes($target, $len_sum);
    }
    $len_sum;
}

# return list of modules used by included files
sub inc_mods {
    my $self = shift;
    my($file, $ref, $mod, @sub_list);
    my @list = ();
    my %printed = ();

    foreach $mod (keys %{$self->{'uses'}}) {
	push(@list, $mod);
    }

    foreach $file (keys %{$self->{'includes'}}) {
	next if $printed{$file};
	$ref = $inc_files{$file};
	$printed{$file} = 1;
	@sub_list = $ref->inc_mods();
	@list = (@list, @sub_list);
    }
    @list;
}

# filenames containing the list of modules used by file and all its includes
sub find_mods {
    my $self = shift;
    my($ref, $modname, $file, @list, $base);
    my @module_files = ();
    my @mod_list = ();
    my @tmp_list = ();

# find modules used by include files
    if (%{$self->{'includes'}}) {
	foreach $file (keys %{$self->{'includes'}}) {
	    $ref = $inc_files{$file};
	    @list = $ref->inc_mods();
            @tmp_list = @mod_list;
	    @mod_list = (@tmp_list, @list);
	}
    }

# add them to the uses list (hash ensures uniqueness)
    foreach $modname (@mod_list) {
	$self->{'uses'}{$modname} = 1;
    }
    
# now find the filename that contains the module information
    foreach $modname (keys %{$self->{'uses'}}) {
	if ($main::cray || $main::parasoft) {
	    if ($file = $main::mod_files{$modname}) {
                $base = &main::basename($file, @main::suffixes);
		$file = $base . "." . $main::obj_ext;
		push(@module_files, $file);
	    } else {
		warn "Don't know where module $modname lives.\n";
	    }
	} else {
	    $modname .= "." . $main::mod_ext;
	    push(@module_files, $modname);
	}
    }
    sort(@module_files);
}

sub print {
    my $self = shift;
    my $source = $self->{'source_file'};
    my $compile_string = "\t" . '$(CFT) $(FFLAGS) -c';
    my($base, $object, $modname, $flag, $target, $ftarget);

    $base = &main::basename($source, @main::suffixes);
    $target = $base . "." . $main::obj_ext;
    if ($main::stupid) {
	$ftarget = $base . "." . $main::ext;
    }

    if ($main::cray) {
	$flag = " -p ";
    } elsif ($main::parasoft) {
	$flag = " -module ";
    }

# print out "include" dependencies
    if (%{$self->{'includes'}}) {
	my $len_sum = length($target) + 1;
	if ($main::add_ext) {
	    $target .= " $base.$main::add_ext";
	    $len_sum += length($base) + length($main::add_ext) + 2;
	}
        print "$target:";
        $self->print_includes($target, $len_sum);
        print "\n";
        if ($main::stupid) {
	    $len_sum = length($ftarget) + 1;
            print "$ftarget:";
            $self->print_includes($ftarget, $len_sum);
            print "\n";
        }
    }

# clean out "use" of modules in own file
    foreach $mod ( keys %{$self->{'uses'}} ) {
	if ( ${$self->{'modules'}}{$mod} ) {
	    delete ${$self->{'uses'}}{$mod};
	}
    }

# print out "use" dependencies
    if (%{$self->{'uses'}} || %{$self->{'includes'}}) {
	@module_files = $self->find_mods();
	my $len_sum = length($target) + 1;
	print "$target:";
	foreach $file (@module_files) {
	    my $len = length($file) + 1;
	    if (($len_sum + $len > 80) &&
		     (length($target) + 1 < $len_sum)) {
		print "\n$target:";
		$len_sum = length($target) + 1;
	    }
	    $len_sum += $len;
	    print " " . $file;
	}
	if ($main::need_f) {
	    my $len = length($ftarget) + 1;
	    if (($len_sum + $len > 80) &&
		     (length($target) + 1 < $len_sum)) {
		print "\n$target:";
		$len_sum = length($target) + 1;
	    }
	    print " " . $ftarget;
	}
	print "\n";
# extra Cray / Parasoft stuff
        if ($main::cray || $main::parasoft) {
	    print $compile_string;
	    foreach $file (@module_files) {
		print $flag . $file;
	    }
	    if ($main::stupid) {
	        print " " . $ftarget . "\n";
	    } else {
	        print " " . $source . "\n";
	    }
	}
    }
}


# Start of main program
package main;

if ($] < 5.000) { die "Need perl 5.000 or newer\n"; }
use File::Basename;
use Getopt::Long;
@suffixes = qw( .c .C .cc .cxx .cpp .f .F .fcm .FCM .f90 .F90 .for);

GetOptions("s", "e=s", "f=s", "I=s@", "m=s", "c", "p", "g", "h", "o=s", "a=s")
		       || die "problem in GetOptions";

# For compilers that don't invoke cpp for you
if ($opt_s) {
	$stupid = 1;
}
if ($opt_e) {
	$ext = $opt_e;
} else {
	$ext = "f";
}

# list of directories to search, starting with current directory
if (@opt_I) {
    @incdirs = @opt_I;
} elsif (@opt_i) {
    @incdirs = @opt_i;
}

if ($opt_f) {
    $mf = $opt_f;
} elsif (-f "makefile") {
    $mf = 'makefile';
} else {
    $mf = 'Makefile';
}
if ( !(-f $mf)) {
    system "touch $mf";
}

# extension used for compiler's private module information
if ($opt_m) {
    $mod_ext = $opt_m;
} else {
    $mod_ext = 'mod';
}

if ($opt_c) {
    $cray = 1;
}

if ($opt_p) {
    $parasoft = 1;
}

if ($opt_g) {
    $sgi = 1;
    $mod_ext = 'kmo';
}

if ($opt_h) {
    $hp = 1;
}

# need to add some more dependencies so the .f file gets created
if ($stupid && ($cray || $parasoft)) {
    $need_f = 1;
}

if ($opt_c && $opt_p) {
    die "Doesn't make sense to have both Cray and Parasoft options!";
}

# object file extension
if ($opt_o) {
    $obj_ext = $opt_o;
} else {
    $obj_ext = 'o';
}

# extension for additional targets (like .prj)
if ($opt_a) {
    $add_ext = $opt_a;
}

$mystring = '# DO NOT DELETE THIS LINE - used by make depend';


# Search for the includes in all the files
foreach $file (@ARGV) {
    $sources{$file} = new source_file($file, $file);
    $sources{$file}->find_includes();
}

# Create new Makefile with new dependencies.

open(MFILE, $mf) || die "can't read Makefile $mf: $!\n";
open(NMFILE, "> Makefile.new") || die "can't write Makefile.new: $!\n";
select(NMFILE);

while (<MFILE>) {
    if (!/$mystring/) {
        print;
    } else {
        last;
    }
}

print $mystring, "\n";

# Now print out include and use dependencies in sorted order.
foreach $target (sort keys(%sources)) {
    $sources{$target}->print();
}

# print out module dependencies
if ( !( $cray || $parasoft) ) {
    foreach $modname (sort keys(%mod_files)) {
        ($name, $path, $suffix) =
		&fileparse($sources{$mod_files{$modname}}->{'filepath'}, @suffixes);
	$object = $name . "." . $obj_ext;
##	$object = $path . "/" . $name . "." . $obj_ext;
	 print "$modname.$mod_ext: $object\n";
    }
}

# Sort out the Makefiles

rename($mf, "$mf.old") || warn "can't overwrite $mf.old: $!\n";
rename('Makefile.new', $mf) ||
     warn "can't move Makefile.new to $mf: $!\n";

# Delete those RCS files we checked out
foreach $file (keys %rcs) {
    unlink($file);
}

#
# End of main
#

sub findfile {
# Let's see if we can find the included file.  Look in current
# directory first, then in directories from -I arguments.  Finally,
# look for RCS files in current directory.
    my $file = shift;
    my($found, $i, $filepath);

    $found = 0;

    if ( -f $file ) {
        $found = 1;
        $file =~ s#^\./##;          # convert ./foo.h to foo.h
        return $file;
    }
    foreach $i (0 .. $#incdirs) {
        $filepath = $incdirs[$i]."/".$file;
        if ( -f $filepath ) {
	    $found = 1;
            $filepath =~ s#^\./##;          # convert ./foo.h to foo.h
            return $filepath;
        }
    }
#see if it is a checked-in RCS file
    if (-f "RCS/$file,v" || -f "$file,v" ) {
        $found = 1;
        system("co $file");
	$filepath = $file;
        $rcs{$file} = 1;
    }
    if ( ! $found ) {
	$filepath = "";
    }
    $filepath;
}
