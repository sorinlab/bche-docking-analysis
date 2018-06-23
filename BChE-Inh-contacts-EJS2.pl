#!/usr/bin/perl


### Subroutines and Modules ###
use Cwd;
use File::Copy;
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


########## USAGE and global variables ####################
$usage="\nUsage: \.\/kmeans_structural_analysis\.pl  \[Project\#\]  \[Inhibitor 3-Letter Code\]  \[contact cutoff (A)\]  \[Keyfile\]\n\n";
$proj     = $ARGV[0] || die "$usage\n";
$code     = $ARGV[1] || die "$usage\n";
$cutoff   = $ARGV[2] || die "$usage\n";
$keyfile  = $ARGV[3] || die "$usage\n";
chomp $proj;
chomp $code;
chomp $keyfile;
$dir = getcwd;
$numatoms = 8309;

#### Read in Inhibitor Key File ####
# and idetnify each chemical group in the inhibitor
open(KI,"<$keyfile")|| die "\n\tError reading from $keyfile\n\n";
while (<KI>){
# for every line in the file KI
    if (/Proj $proj/i) {
        $spool = 1;  # track spool counter
        next;
    }
    # is this the proj you're looking for ...	
    elsif ($_ =~ /Proj\s\d/ && $_ !~ /Proj $proj/) {
        $spool = 0;
    }
    elsif ($spool) {
	    	# initialize @singroup array
      	undef @singroup;
      	$groupline = $_;
      	# remove white space from line ...
		for($groupline) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
      	# assign array elements to new array
		@logline = split(/ /,$groupline);
      	$group = "@logline[0]";
      	chomp $group;
		# define which atoms are in the group $group	
      	$atoms = "@logline[1..$#logline]";
          	@singroup = split(/ /,$atoms);
		# assign those atom numbers to an array for that group
           if ($group eq "PheCOP"){
                push (@PheCOP,@singroup);
                undef @singroup;
           }
           elsif ($group eq "OCH1"){
                push (@OCH1, @singroup);
                undef @singroup;
           }
           elsif ($group eq "OCH2"){
                push (@OCH2, @singroup);
                undef @singroup;
           }
           elsif ($group eq "Phos"){
                push (@Phos, @singroup);
           	     undef @singroup;
           }
    }
}


### identify all pdb's 
my @pdbfiles = `ls *.pdb`;
$count = scalar(@pdbfiles);
$pdbtotal = $count;
	


### Find contacts in each pdb ###
foreach $pdb(@pdbfiles){
	chomp $pdb;
	
	### Open/Read pdb line and save protein coordinates ###
	open(PDB,"<$pdb") || die "\n\tError reading from $pdb\n\n";
	$i = 0;
	while(defined($origline=<PDB>)) {
		$line = $origline;
		chomp $line;
		for($line) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
		@input = split(/ /,$line);

		### Save protein coordinates ###
		if(@input[0] eq 'ATOM'){
			$atomnum = @input[1];
			$atomname = @input[2];
			$resname[$atomnum] = @input[3];
			$resnum[$atomnum]  = @input[5];
			$x[$atomnum]= @input[6];
			$y[$atomnum]= @input[7];
			$z[$atomnum]= @input[8];
			if((@input[3] ne 'SOL')&&(@input[3] ne '$code')&&(@input[3] ne 'UNK')&&(@input[3] ne 'CL')&&(@input[3] ne 'NA')){ 
				$protx[$atomnum] = $x[$atomnum];
				$proty[$atomnum] = $y[$atomnum];
				$protz[$atomnum] = $z[$atomnum];
				$i++;
			}
		}
	}
	close(PDB);


	### Find the inhibitor group's atom distance from protein ###
	### let's get rid of these blocks/condense them into a single loop ###
	### add "." delimiters between amino acid name, res number, and add on inhibitor chemical group name
	### also add if loop to detect BB versus SC interaction

	foreach $num (@PheCOP){
		for ( $i=1; $i<=$numatoms; $i++ ){
			$X = $protx[$i]-$x[$num];
			$Y = $proty[$i]-$y[$num];
			$Z = $protz[$i]-$z[$num];
			$D2 = ($X**2) + ($Y**2) + ($Z**2);
			$D = sqrt($D2);
			if ( $D < $cutoff ){
				$phe = "$resname[$i]$resnum[$i]";
				for($phe) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
				push @pheline, $phe;
				@pheorg = uniq(@pheline);
			}
		}
	}


	foreach $num (@Phos){
		#print $num;
		for ( $i=1; $i<=$numatoms; $i++){
			$X = $protx[$i]-$x[$num];
			$Y = $proty[$i]-$y[$num];
			$Z = $protz[$i]-$z[$num];
			$D2 = $X**2 + $Y**2 + $Z**2;
			$D = sqrt($D2);
			if ($D < $cutoff ){
				@phoorg=();
				$jj = $j - 1;
				$pho = "$resname[$i]$resnum[$i]";
				for($pho) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
				push @pholine, $pho;
				@phoorg = uniq(@pholine);
			}
		}
	}
	
	foreach $num (@OCH1){
		#print $num;
		for ( $i=1; $i<=$numatoms; $i++){
			$X = $protx[$i]-$x[$num];
			$Y = $proty[$i]-$y[$num];
			$Z = $protz[$i]-$z[$num];
			$D2 = $X**2 + $Y**2 + $Z**2;
			$D = sqrt($D2);
			if ($D < $cutoff ){
				@oc1org=();
				$jj = $j - 1;
				$oc1 = "$resname[$i]$resnum[$i]";
				for($oc1) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
				push @oc1line, $oc1;
				@oc1org = uniq(@oc1line);
			}
		}
	}

	foreach $num (@OCH2){
		for ( $i=1; $i<=$numatoms; $i++){
			$X = $protx[$i]-$x[$num];
			$Y = $proty[$i]-$y[$num];
			$Z = $protz[$i]-$z[$num];
			$D2 = $X**2 + $Y**2 + $Z**2;
			$D = sqrt($D2);
			if ($D < $cutoff ){
				@oc2org=();
				$jj = $j - 1;
				$oc2 = "$resname[$i]$resnum[$i]";
				for($oc2) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
				push @oc2line, $oc2;
				@oc2org = uniq(@oc2line);
			}
		}
	}


	### Add these array values to master-arrays and re-initialize ###
	push(@phecluster,@pheorg);
	@pheorg=();
	push(@phocluster,@phoorg);
	@phoorg=();
	push(@oc1cluster,@oc1org);
	@oc1org=();
	push(@oc2cluster,@oc2org);
	@oc2org=();
}



####Finding total interaction time###
print "CLUSTER GROUP: $cluster\n";
print "----------------------------------\n";
print "----------------------------------\n";


print "Interactions with Phenyl Group:\n";
$phecount{$_}++ foreach @phecluster;
foreach my $name (sort { $phecount{$b} <=> $phecount{$a} or $a cmp $b } keys %phecount){
	printf"%-8s%d%s\n", $name, $phecount{$name}*100/$pdbtotal,"%";
}
@phecluster = ();
print "----------------------------------\n";


print "Interactions with Phosphate Group:\n";
my %phocount;
$phocount{$_}++ foreach @phocluster;
foreach my $name (sort { $phocount{$b} <=> $phocount{$a} or $a cmp $b } keys %phocount){
	printf"%-8s%d%s\n", $name, $phocount{$name}*100/$pdbtotal,"%";
}
@phocluster =();
print "----------------------------------\n";


print "Interactions with Alkyl (Left) Group:\n";
my %oc1count;
$oc1count{$_}++ foreach @oc1cluster;
foreach my $name (sort { $oc1count{$b} <=> $oc1count{$a} or $a cmp $b } keys %oc1count){
	printf"%-8s%d%s\n", $name, $oc1count{$name}*100/$pdbtotal,"%";
}
@oc1cluster = ();
print "----------------------------------\n";


print "Interactions with Alkyl (Right) Group:\n";
my %oc2count;
$oc2count{$_}++ foreach @oc2cluster;
foreach my $name (sort { $oc2count{$b} <=> $oc2count{$a} or $a cmp $b } keys %oc2count){
	printf"%-8s%d%s\n", $name, $oc2count{$name}*100/$pdbtotal,"%";
}
@oc2cluster=();
print "----------------------------------\n";


print"\n";
