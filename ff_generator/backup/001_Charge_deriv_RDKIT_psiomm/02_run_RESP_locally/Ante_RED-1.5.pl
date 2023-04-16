#!/usr/bin/perl -w
# use warnings;
use strict;
use FileHandle;

###############################################################################
##                            Ante_R.E.D. 1.5                                ##
##                 http://q4md-forcefieldtools.org/RED/                      ##
##                                                                           ##
##   Ante_R.E.D. developments were initiated in Prof. D.A. Case's lab. by    ##
##         R. Lelong,(1,2) P. Cieplak,(3) & F.-Y. Dupradeau (1,2)            ##
##                                                                           ##
##             Ante_R.E.D. 1.x were developed in Amiens by                   ##
##                  F.-Y. Dupradeau(4) & P. Cieplak(3)                       ##
##                                                                           ##
##          Ante_R.E.D. 2.x is currently developed in Amiens by              ##
##           G. Klimerak,(4) P. Cieplak(3) & F.-Y. Dupradeau(4)              ##
##                                                                           ##
##  (1) DMAG EA-3901 & Faculte de Pharmacie, UPJV, Amiens, France            ##
##  (2) The Scripps Research Institute, La Jolla, CA, USA                    ##
##  (3) Sanford|Burnham Institute for Medical Research, La Jolla, CA, USA    ##
##  (4) UMR CNRS 6219 & UFR de Pharmacie & UPJV, Amiens, France              ##
##                                                                           ##
##            Distributed under the GNU General Public License               ##
###############################################################################   

# Information:
# One can modify the "$sort" & "$notlig" variables in the main section of the program (line 591)
#
######## Chemical elements recognized by Ante_R.E.D.: up to Bromine (Z = 35)...
# 'Terminal' chemical elements that may generate vdW bumps in a bad initial structure
my %Terminal = (H =>"1",F =>"9",CL =>"17",BR =>"35");

# Elements: Warning: Calcium is XX !!!
my %Elements = (H =>"1",LI =>"3",BE =>"4",B =>"5",C =>"6",N =>"7",O =>"8",F =>"9",NA =>"11",MG =>"12",AL =>"13",SI =>"14",P =>"15",S =>"16",CL =>"17",K =>"19",XX =>"20",SC=>"21",TI => "22",V=>"23",CR =>"24",MN=>"25",FE =>"26",CO =>"27",NI =>"28",CU =>"29",ZN =>"30",GA =>"31","GE" =>"32",AS =>"33",SE =>"34",BR =>"35");

# Maximal number of bonds for each chemical element
my %Valence = (H =>"1",LI =>"1",BE =>"2",B =>"3",C =>"4",N =>"4",O =>"2",F =>"1",NA =>"1",MG =>"2",AL =>"6",SI =>"6",P =>"5",S =>"6",CL =>"1",K =>"1",XX =>"2",SC=>"6",TI => "6",V=>"6",CR =>"6",MN=>"8",FE =>"6",CO =>"6",NI =>"6",CU =>"6",ZN =>"6",GA =>"3","GE" =>"4",AS =>"3",SE =>"2",BR =>"1");

# Radii used by Ante_R.E.D.
my %Radius = (H =>"0.230",LI =>"0.680",BE =>"0.350",B =>"0.830",C =>"0.680",N =>"0.680",O =>"0.680",F =>"0.640",NA =>"0.970",MG =>"1.100",AL =>"1.350",SI =>"1.200",P =>"0.750",S =>"1.020",CL =>"0.990",K =>"1.330",CA =>"0.990",SC=>"1.440",TI => "1.470",V=>"1.330",CR =>"1.350",MN=>"1.350",FE =>"1.340",CO =>"1.330",NI =>"1.500",CU =>"1.520",ZN =>"1.450",GA =>"1.220","GE" =>"1.170",AS =>"1.210",SE =>"1.220",BR =>"1.210");

# Empirical radii: You might try the "Empirical radii" instead of the "Radii used by Ante_R.E.D."...
#my %Radius = (H =>"0.250",LI =>"1.450",BE =>"1.050",B =>"0.850",C =>"0.700",N =>"0.650",O =>"0.600",F =>"0.500",NA =>"1.800",MG =>"1.500",AL =>"1.250",SI =>"1.100",P =>"1.000",S =>"1.000",CL =>"1.000",K =>"2.200",CA =>"1.800",SC=>"1.600",TI => "1.400",V=>"1.350",CR =>"1.400",MN=>"1.400",FE =>"1.400",CO =>"1.350",NI =>"1.350",CU =>"1.350",ZN =>"1.350",GA =>"1.300","GE" =>"1.250",AS =>"1.150",SE =>"1.150",BR =>"1.150");

# Calculed radii: You might try the "Calculated radii" instead of the "Radii used by Ante_R.E.D."...
# my %Radius = (H =>"0.530",LI =>"1.670",BE =>"1.120",B =>"0.870",C =>"0.670",N =>"0.560",O =>"0.480",F =>"0.420",NA =>"1.900",MG =>"1.450",AL =>"1.180",SI =>"1.110",P =>"0.980",S =>"0.880",CL =>"0.790",K =>"2.430",CA =>"1.940",SC=>"1.840",TI => "1.760",V=>"1.710",CR =>"1.660",MN=>"1.610",FE =>"1.560",CO =>"1.520",NI =>"1.490",CU =>"1.450",ZN =>"1.420",GA =>"1.360","GE" =>"1.250",AS =>"1.140",SE =>"1.030",BR =>"0.94");

sub SEARCH_ATOM
{
	my ($brut,$nbatom) = @_;
	my ($error,$carac,$ret,$len,$at) = 0;
	$brut = uc($brut);
	($carac,$ret) = ($brut =~ m/[a-zA-Z]+/g); # Number of letters for each chemical element
	if (defined($carac)) {
		$len = length($carac);
		if ($len >= 2) { # Take out the two first letters for the chemical element
			($carac,$ret) = ($brut =~ m/[a-zA-Z][a-zA-Z]/g);
			if (exists $Elements{$carac}) { # Test if the chemical element exists (with 2 letters)
				if ($carac eq "XX")	{ $at = "CA"; }
				else	{ $at = $carac; }
			}
			else {	# Take out the first letter if it does not exist
				($carac,$ret) = ($brut =~ m/[a-zA-Z]/g);
				# Test if the chemical element exists (with 1 letter)
				if (exists $Elements{$carac})	{ $at = $carac; }
				else {
					print "\n  ERROR: The atom number ",$nbatom+1," is not handled by Ante_R.E.D.\n\n"; $error = 1;
				}
			}
		}
		elsif ($len == 1) {
			($carac,$ret) = ($brut =~ m/[a-zA-Z]/g);
			if (exists $Elements{$carac})	{ $at = $carac; }
			else {
				print "\n  ERROR: The atom number ",$nbatom+1," is not handled by Ante_R.E.D.\n\n";
				$error = 1;
			}
		}
		else {
			print "\n  ERROR: The atom number ",$nbatom+1," is not handled by Ante_R.E.D.\n\n";
			$error = 1;
		}
	}
	else {
		print "\n  ERROR: There is no letter in the atom name number: ",$nbatom+1,"\n\n";
		$error = 1;
	}
	if ($error == 1) { exit(1); }
	return $at;
}

############## Calculate distances between two points in 3D ###############
sub DISTANCE
{
	my ($x,$y,$z,$xb,$yb,$zb) = @_;
	my $distance = sqrt(($xb-$x)**2+($yb-$y)**2+($zb-$z)**2);
	return $distance;
}

############### Test if the two atoms are connected ###############
sub CONNECTED
{
	my ($radiusA,$radiusB,$x,$y,$z,$xb,$yb,$zb) = @_;
	my ($distance,$res) = 0;
	my $CONSTANT = 0.45;
	$distance = DISTANCE($x,$y,$z,$xb,$yb,$zb);
	$res = $radiusA + $radiusB + $CONSTANT;
	if ($distance < $res)	{ return 1; }
	else	{ return 0; }
}

############### Establish the matrix of connections ###############
sub MATRIX_CONNECTION
{
	my ($nbatom,$ref_tab,$ref_matrix,$ref_H_bond_H,$nb_H_bond_H) = @_;
	my ($i,$j,$atom,$atom2) = 0;
	$nb_H_bond_H = 0;
	while ($i < $nbatom) {
		$j = 0;
		while ($j < $nbatom) {
			if ($i != $j) {
				$atom = $$ref_tab [$i][0];
				$atom = SEARCH_ATOM ($atom,$nbatom);
				$atom2 = $$ref_tab [$j][0];
				$atom2 = SEARCH_ATOM ($atom2,$nbatom);
				if (($atom eq "H") && ($atom2 eq "H"))
				{
					$$ref_matrix [$i][$j] = 0;
					if (CONNECTED ($$ref_tab [$i][1],$$ref_tab [$j][1],$$ref_tab [$i][2],$$ref_tab [$i][3],$$ref_tab [$i][4],$$ref_tab [$j][2],$$ref_tab [$j][3],$$ref_tab [$j][4]) == 1) { # Store H info to print a warning later
						$$ref_H_bond_H[0][$nb_H_bond_H] = $i;
						$$ref_H_bond_H[1][$nb_H_bond_H] = $j;
						$nb_H_bond_H++;
					}
				}
				else {
					$$ref_matrix [$i][$j] = CONNECTED ($$ref_tab [$i][1],$$ref_tab [$j][1],$$ref_tab [$i][2],$$ref_tab [$i][3],$$ref_tab [$i][4],$$ref_tab [$j][2],$$ref_tab [$j][3],$$ref_tab [$j][4]);
				}
			}
			else {	# An atom cannot be conected to itself
				$$ref_matrix [$i][$j] = 0;
			}
			$j++;
		}
		$i++;
	}
	return $nb_H_bond_H;
}

############### For each H atom, add the same number as the heavy atom it is bound to  ###############
sub PUT_NUMBER_TO_H
{
	my ($nbatom,$nbH,$ref_tab,$ref_matrix,$ref_tabH,$i,$compteur) = @_;
	my $j = 0;
	while ($j<$nbH) {
		if($$ref_matrix[$i][$$ref_tabH[0][$j]] == 1) {
			$$ref_tab[$$ref_tabH[0][$j]][0] = $$ref_tab[$$ref_tabH[0][$j]][0].$compteur;
		}
		$j++;
	}
}

############### Add a number for each chemical element  ###############
sub PUT_NUMBER
{
	my ($nbatom,$nbH,$ref_tab,$ref_matrix,$ref_tabH) = @_;
	my ($i,$bool,$j) = 0;
	my $compteur = 1;
	while($i<$nbatom) {
		$bool = 0; $j = 0;
		while($j<$nbH) {
			if($$ref_matrix[$i][$$ref_tabH[0][$j]] == 1)	{ $bool = 1; }
			$j++;
		}
		if ($bool == 1) {
			if($$ref_tab[$i][0] =~ /^[^H]/) {
				$$ref_tab[$i][0] = $$ref_tab[$i][0].$compteur;
				PUT_NUMBER_TO_H($nbatom,$nbH,$ref_tab,$ref_matrix,$ref_tabH,$i,$compteur);
				$compteur ++;
			}
		}
		elsif ($$ref_tab[$i][0] =~ /^[^H]/) {
			$$ref_tab[$i][0] = $$ref_tab[$i][0].$compteur;
			$compteur ++;
		}
		$i++;
	}
	print "  Add a number to each new atom name...\t\t\t\t\t[Done]\n";
}

############### Locate the C atoms with two or more H connected atoms ###############
sub PUT_CT
{
	my ($nbatom,$ref_tab,$ref_matrix) = @_;
	my ($i,$j,$nbcolC,$nbH) = 0;
	my (@tabC,@tabH);
	$nbcolC = 0; $nbH = 0;
	while($i<$nbatom) {	# Store the line number of my table for C & H atoms
		if ($$ref_tab[$i][0] eq "C") {
			$tabC[0][$nbcolC] = $i;
			$tabC[1][$nbcolC] = 0;
			$nbcolC++;
		}
		elsif ($$ref_tab[$i][0] eq "H") {
			$tabH[0][$nbH] = $i;
			$nbH++;
		}
		$i++;
	}
	$i = 0;		# Count the number of connected H atom for each C atom
	while($i<$nbcolC) {
		$j = 0;	
		while($j<$nbH) {
			if($$ref_matrix[$tabC[0][$i]][$tabH[0][$j]] == 1) {
				$tabC[1][$i]++;
			}
			$j++;
		}
		$i++;
	}
	$i = 0;		# Add a 'T' for each C atom connected to 2 or more H atoms
	while($i<$nbcolC) {
		if($tabC[1][$i]>=2) {
			$$ref_tab[$tabC[0][$i]][0] = "CT";
		}
	$i++;
	}
	print "  Use 'CT' for C atoms connected to 2 or 3 H atoms...\t\t\t[Done]\n";
	PUT_NUMBER($nbatom,$nbH,$ref_tab,$ref_matrix,\@tabH);
}

############### Sort by number ################
sub SORT_BY_NUMBER
{
	my ($nbatom,$ref_tab) = @_;
	my ($i,$j,$num,$ret,$num2) = 1;
	my @temp;
	while ($i<$nbatom-1) {
		$j = 0;
		while($j<($nbatom - $i)) {
			($num,$ret) = ($$ref_tab[$j][0] =~ m/[0-9]+/g);
			($num2,$ret) = ($$ref_tab[$j+1][0] =~ m/[0-9]+/g);
			if (!(defined($num))) { $num=0; }
			if (!(defined($num2))) { $num2=0; }
			if ($num > $num2) { # Debug
				$temp[0] = $$ref_tab[$j];
				$$ref_tab[$j] = $$ref_tab[$j+1];
				$$ref_tab[$j+1] = $temp[0];
			}
			$j++;
		}
		$i++;
	}
	print "  Sort atoms by numbers...\t\t\t\t\t\t[Done]\n";
}

############### Sort by letter & number ###############
sub SORT_LETTER_BY_NUMBER
{
	my($nbatom,$ref_tab) = @_;
	my ($i,$num,$ret,$num2,$first) = 0;
	my @temp;
	while ($i<$nbatom-1) {
		($num,$ret) = ($$ref_tab[$i][0] =~ m/[0-9]+/g);
		($num2,$ret) = ($$ref_tab[$i+1][0] =~ m/[0-9]+/g);
		if (!(defined($num))) { $num=0; }
		if (!(defined($num2))) { $num2=0; }
		if ($num != $num2) { $i++; }	# Debug
		else {
			$first = $i;
			while (($$ref_tab[$i+1][0] =~ /^H/)&&($num == $num2)&&($i<($nbatom-2))) { # Debug
				$i++;
				($num,$ret) = ($$ref_tab[$i][0] =~ m/[0-9]+/g);
				($num2,$ret) = ($$ref_tab[$i+1][0] =~ m/[0-9]+/g);
				if (!(defined($num))) { $num=0; }
				if (!(defined($num2))) { $num2=0; }
			}
			if (!(defined($num))) { $num=0; }
			if (!(defined($num2))) { $num2=0; }
			if (($$ref_tab[$i+1][0] =~ /^[^H]/)&&($num == $num2)) { # Debug
				$temp[0] = $$ref_tab[$first];
				$$ref_tab[$first] = $$ref_tab[$i+1];
				$$ref_tab[$i+1] = $temp[0];
			}
			$i++;
		}
	}
	print "  Sort atoms by atom names...\t\t\t\t\t\t[Done]\n";
}

################ Write lines with 'CONECT' in the PDB output #################
sub WRITE_CONNECT
{
	my ($nameout,$nameout2,$nbatom,$ref_matrix) = @_;
	my ($i,$k,$j) = 0;
	open(FILE,">> $nameout") or die "can't open this file";
	open(FILE2,">> $nameout2") or die "can't open this file";
	while($i<$nbatom) {
		print FILE "CONECT  ";
		print FILE2 "CONECT  ";
		$k = $i+1;
		printf FILE "%3d  ",$k;
		printf FILE2 "%3d  ",$k;
		$j = 0;
		while($j<$nbatom) {
			if($$ref_matrix[$i][$j] == 1) {
			$k = $j+1;
			printf FILE "%3d  ",$k;
			printf FILE2 "%3d  ",$k;
			}
			$j++;
		}
		print FILE "\n";
		print FILE2 "\n";
		$i++;
	}
	print FILE "END \n";
	print FILE2 "END \n";
	close(FILE);
	close(FILE2);
}
################ Write results in the P2N & PDB outputs #################
sub WRITE_FILE_OUT
{
	my ($name,$nameout,$nameout2,$ref_tab,$nbatom) = @_;
	my $number = 1;
	my ($atom,$x,$y,$z,$rest,$nb,$element) = 0;
# Format of the P2N File
format ATOMP2N =
ATOM   @### @||| @<<< @>>>    @##.### @##.### @##.###                    @<<<<<
$number,$atom,$rest,$nb,$x,$y,$z,$element
.

# Format of the PDB File
format ATOMPDB =
ATOM   @### @||| @<<< @>>>    @##.### @##.### @##.###
$number,$element,$rest,$nb,$x,$y,$z
.

	# Format for the PDB output file
	format_name FILEOUT "ATOMP2N";
	format_name FILEOUT2 "ATOMPDB";
	open(FILE,$name) or die "can't open this file";
	open(FILEOUT,"> $nameout") or die "can't write in this file";
	open(FILEOUT2,"> $nameout2") or die "can't write in this file";
	print FILEOUT  "REMARK\n";
	print FILEOUT  "REMARK TITLE MOLECULE\n";
	print FILEOUT  "REMARK CHARGE-VALUE 0\n";
	print FILEOUT  "REMARK MULTIPLICITY-VALUE 1\n";
	print FILEOUT  "REMARK\n";
	print FILEOUT2 "REMARK\n";
	while(<FILE>) {
		if ((/^ATOM/ig) || (/^HETATM/ig)) {
			$atom = $$ref_tab[$number-1][0];
			$x = $$ref_tab[$number-1][2];
			$y = $$ref_tab[$number-1][3];
			$z = $$ref_tab[$number-1][4];
			$rest = $$ref_tab[$number-1][5];;
			$nb = $$ref_tab[$number-1][6];;
			$element = $$ref_tab[$number-1][7];;
			write FILEOUT; # Write in -out.p2n output file with the right mask
			write FILEOUT2; # Write in -out.pdb output file with the right mask
			$number ++;
		}
		elsif (/^REMARK/ig) {
			print FILEOUT $_;
			print FILEOUT2 $_;
		}
	}
	close (FILE);
	close (FILEOUT);
	close (FILEOUT2);
	print "  Create the P2N output file:\t\t\t$nameout...\t[Done]\n";
	print "  Create a new PDB output file:\t\t\t$nameout2...\t[Done]\n";
}

############### Create the summary for bond information and warnings ###############
sub INFORMATION
{
	my ($nameinfo,$nbatom,$ref_tab,$ref_matrix,$ref_H_bond_H,$nb_H_bond_H) = @_;
	my ($ret,$i,$j,$bond,$at,$k,$nb,$atom) = 0;
	$nameinfo = (split /\./,$nameinfo)[0];
	$nameinfo = $nameinfo."-info.txt";
	open(FILEINFO,"> $nameinfo") or die "can't write in this file";
	print FILEINFO "  NUM  ATOM_NAME  NUM_BOND  CONECT TO XX,XX,X.... ATOMS \n\n";
	$i = 0;
	while($i<$nbatom) {
		$j = 0; $bond = 0;
		while($j<$nbatom) {
			if($$ref_matrix[$i][$j] == 1) {
				$bond++;
			}
			$j++;
		}
		$atom = $$ref_tab[$i][0];
		printf FILEINFO "  %3d       %4s         %1d  CONECT TO ",$i+1,$atom,$bond;
		$j = 0;
		while($j<$nbatom) {
			if($$ref_matrix[$i][$j] == 1) {
				$at = $$ref_tab[$j][0];
				printf FILEINFO "%4s ",$at;
			}
			$j++;
		}
		$atom = SEARCH_ATOM ($atom,$i);
		$k = 0; $nb = 0;
		while($k<$nb_H_bond_H) {
			if ($$ref_H_bond_H[0][$k] == $i) {
				$nb++;
			}
			if (($nb != 0) && ($k == $nb_H_bond_H-1)) {
				print FILEINFO "  WARNING !!! This atom is close to ",$nb," H atom(s). Check your initial PDB file, & you might re-built your initial structure. Atom connectivities between two hydrogens are never considered by Ante_R.E.D.";
			}
			$k++;
		}
		if ($bond>$Valence{$atom}) {
			print FILEINFO "  WARNING !!! Number of bond(s) too high for this atom. Check your initial PDB file, & re-built your initial structure. ";
		}
		if ($bond == 0) { # Debug
			print FILEINFO "  WARNING !!! Isolated atom detected. Check your initial PDB file, re-build your initial structure, or try the \"Empirical radii\" or \"Calculed radii\" instead of the \"Radii used by Ante_R.E.D.\" (see code lines 23, 26 & 29).";
		}
		print FILEINFO "\n";
		$i++;
	}
	print FILEINFO "END";
	close(FILEINFO);
	print "  Create the information output file:\t\t$nameinfo...\t[Done]\n\n";
}

############### Avoid atom connectivities between two 'terminal' atoms ###############
sub TERMINAL_ATOM
{
	my ($nbatom,$ref_tab,$ref_matrix) = @_;
	my ($i,$j,$atom,$atom2) = 0;
	$i = 0;
	while ($i < $nbatom) {
		$j = 0;
		while ($j < $nbatom) {
			if ($i != $j) {
				$atom = $$ref_tab [$i][0];
				$atom = SEARCH_ATOM ($atom,$nbatom);
				$atom2 = $$ref_tab [$j][0];
				$atom2 = SEARCH_ATOM ($atom2,$nbatom);
				if ((exists $Terminal{$atom}) && (exists $Terminal{$atom2})) {
					$$ref_matrix[$i][$j] = 0;
				}
			}
			$j++;
		}
		$i++;
	}
}

############### Create the input files for Gaussian & GAMESS-US ###############
sub WRITE_INPUT_GG
{
	my ($nameout3,$nameout4,$nameout5,$ref_tab,$nbatom) = @_;
	my ($atom,$residu,$x,$y,$z,$i) = 0;

# Format of the GAMESS-US + PC-GAMESS .inp File
format ATOMINP =
 @<<  @#.#    @##.### @##.### @##.###
$atom,$residu,$x,$y,$z
.

# Format of the Gaussian .com File
format ATOMCOM =
 @<<    @##.### @##.### @##.###
$atom,$x,$y,$z
.

	format_name FILEOUT3 "ATOMINP"; # Keywords for the input files (.inp & .com)
	format_name FILEOUT4 "ATOMINP"; 
	format_name FILEOUT5 "ATOMCOM";
	open(FILEOUT3,"> $nameout3") or die "can't write in this file";
	open(FILEOUT4,"> $nameout4") or die "can't write in this file";
	open(FILEOUT5,"> $nameout5") or die "can't write in this file";
	print FILEOUT3 " ! Keywords below are only useful for providing a starting input...\n";
	print FILEOUT3 " ! See the GAMESS-US documentation...\n";
	print FILEOUT3 " \$CONTRL  ICHARG=0 MPLEVL=0\n";
	print FILEOUT3 "          RUNTYP=OPTIMIZE SCFTYP=RHF EXETYP=RUN\n";
	print FILEOUT3 "          MULT=1 UNITS=ANGS MAXIT=200\n";
	print FILEOUT3 " !        INTTYP=HONDO QMTTOL=1.0E-08 ITOL=30 ICUT=20\n";
	print FILEOUT3 "          COORD=CART                           \$END\n";
	print FILEOUT3 " \$SCF     DIRSCF=.T. CONV=1.0E-08              \$END\n";
	print FILEOUT3 " \$SYSTEM  TIMLIM=50000 MWORDS=32 MEMDDI=0      \$END\n";
	print FILEOUT3 " \$STATPT  NSTEP=200 OPTTOL=1.0E-06\n";
        print FILEOUT3 "          HESS=CALC IHREP=10                   \$END\n";
	print FILEOUT3 " \$FORCE   METHOD=ANALYTIC VIBANL=.F.           \$END\n";
	print FILEOUT3 " \$BASIS   GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=0\n";
	print FILEOUT3 "          DIFFSP=.F.                           \$END\n";
	print FILEOUT3 " \$GUESS   GUESS=HUCKEL                         \$END\n";
	print FILEOUT3 " \$DATA\n";
	print FILEOUT3 " GAMESS-US optimization output to be used by R.E.D.\n";
	print FILEOUT3 " C1 \n";

	print FILEOUT4 " ! Keywords below are only useful for providing a starting input...\n";
	print FILEOUT4 " ! See the PC-GAMESS documentation...\n";
	print FILEOUT4 " \$CONTRL  ICHARG=0 MPLEVL=0\n";
	print FILEOUT4 "          RUNTYP=OPTIMIZE SCFTYP=RHF EXETYP=RUN\n";
	print FILEOUT4 "          MULT=1 UNITS=ANGS MAXIT=200\n";
	print FILEOUT4 " !        INTTYP=HONDO ITOL=30 ICUT=20\n";
	print FILEOUT4 "          COORD=CART                           \$END\n";
	print FILEOUT4 " \$SCF     DIRSCF=.T. NCONV=1.0E-08             \$END\n";
	print FILEOUT4 " ! \$P2P   P2P=.T. DLB=.T.                      \$END\n";
	print FILEOUT4 " \$SYSTEM  TIMLIM=50000 MWORDS=32               \$END\n";
	print FILEOUT4 " \$STATPT  NSTEP=200 OPTTOL=1.0E-06\n";
	print FILEOUT4 "          HESS=CALC IHREP=10                   \$END\n";
	print FILEOUT4 " \$FORCE   METHOD=ANALYTIC VIBANL=.F.           \$END\n";
	print FILEOUT4 " \$BASIS   GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=0\n";
	print FILEOUT4 "          DIFFSP=.F.                           \$END\n";
	print FILEOUT4 " \$GUESS   GUESS=HUCKEL                         \$END\n";
	print FILEOUT4 " \$DATA\n";
	print FILEOUT4 " PC-GAMESS optimization output to be used by R.E.D.\n";
	print FILEOUT4 " C1 \n";

	print FILEOUT5 "\%Chk=Your-chkfile.chk\n";
	print FILEOUT5 "\%Mem=256MB\n";
	print FILEOUT5 "\%NProc=1\n";
	print FILEOUT5 "\n";
	print FILEOUT5 "#P hf/6-31G* Opt=(Tight,CalcFC) Freq SCF(Conver=8) Test\n";
	print FILEOUT5 "\n";
	print FILEOUT5 "Gaussian optimization output to be used by R.E.D.\n";
	print FILEOUT5 "\n";
	print FILEOUT5 " 0 1\n";
	$i = 0;
	while($i<$nbatom) {
		$atom = $$ref_tab[$i][0];
		$x = $$ref_tab[$i][2];
		$y = $$ref_tab[$i][3];
		$z = $$ref_tab[$i][4];
		$atom = SEARCH_ATOM ($atom,$nbatom);
		$residu = $Elements{$atom};
		write FILEOUT3;
		write FILEOUT4;
		write FILEOUT5;
		$i++;
	}
	print FILEOUT3 " \$END\n\n\n";
	print FILEOUT4 " \$END\n\n\n";
	print FILEOUT5 "\n\n\n";
	if ((-e $nameout3) && (-e $nameout4) && (-e $nameout5)) {
		print "  Create a GAMESS-US input file:\t\t$nameout3...   \t[Done]\n";
		print "  Create a PC-GAMESS input file:\t\t$nameout4...   \t[Done]\n";
		print "  Create a GAUSSIAN input file:\t\t\t$nameout5...   \t[Done]\n";
	}
}
################## Print matrix #####################
##   This function is used to debug Ante_R.E.D.    ##
##   Print the matrix which is given in parameter  ##
sub Printab
{
	my ($nbatom,$ref_tab,$ref_matrix) = @_;
	my ($i,$j) = 0;
	my ($k,$l) = 0;
	while ($i < $nbatom) {
		$j = 0;
		while ($j < 8) {
			print "$$ref_tab[$i][$j],\t";
			$j++;
		}
		print "\n";
		$i++;
	}
	print "\n";
	while ($k < $nbatom) {
		$l = 0;
		while ($l < $nbatom) {
			print "$$ref_matrix[$k][$l],\t";
			$l++;
		}
		print "\n";
		$k++;
	}
	print "\n";
	return 0;
}

############### MAIN PROGRAM SECTION ###############
	my $name = $ARGV[0];
	my ($nameout,$nameout2,$nameout3,$nameout4,$nameout5,$ret,$nbatom,$element,$rest,$nb,$x,$y,$z,$nb_H_bond_H,$rad,$atom,$sort,$notlig,$debug);
	my (@tab,@matrix,@H_bond_H);
	$nbatom = 0;
	system("clear");
	print " ****************************************************************************\n";
	print " **                   Ante_R.E.D. program (version 1.5)                    **\n";
	print " **                 http://q4md-forcefieldtools.org/RED/                   **\n";
	print " **                                                                        **\n";
	print " **             Distributed under the GNU General Public License           **\n";
	print " **                                                                        **\n";
	print " **  Ante_R.E.D. developments were initiated in Prof. D.A. Case's lab. by  **\n";
	print " **          R. Lelong,(1,2) P. Cieplak,(3) & F.-Y. Dupradeau (1,2)        **\n";
	print " **                                                                        **\n";
	print " **             Ante_R.E.D. 1.x were developed in Amiens by                **\n";
	print " **                   F.-Y. Dupradeau(4) & P. Cieplak,(3)                  **\n";
	print " **                                                                        **\n";
	print " **          Ante_R.E.D. 2.0 is currently developed in Amiens by           **\n";
	print " **           G. Klimerak,(4) P. Cieplak(3) & F.-Y. Dupradeau(4)           **\n";
	print " ****************************************************************************\n";
	print " ** (1) DMAG EA-3901 & Faculte de Pharmacie, UPJV, Amiens, France          **\n";
	print " ** (2) The Scripps Research Institute, La Jolla, CA, USA                  **\n";
	print " ** (3) Sanford|Burnham Institute for Medical Research, La Jolla, CA, USA  **\n";
	print " ** (4) UMR CNRS 6219 & UFR de Pharmacie, UPJV, Amiens, France             **\n";
	print " ****************************************************************************\n\n";
## WARNING: If one wishes to sort the atoms of a molecule choose $sort="ON"; otherwise use "OFF" ##
## WARNING: If one wishes to avoid the creation of an atom connectivity between two 'terminal' atoms (in a bad initial geometry) choose $nolig="ON"; otherwise use "OFF" ##
	$sort = "ON";	$notlig = "ON";  $debug = "OFF";     # $sort=$notlig="ON", $debug=off: default

	$sort=~ s/^\s*(.*?)\s*$/$1/;    $sort=uc($sort);
	$notlig=~ s/^\s*(.*?)\s*$/$1/;  $notlig=uc($notlig);
	$debug=~ s/^\s*(.*?)\s*$/$1/;   $debug=uc($debug);
	if (defined($name)) {
		chomp($name);
		open(FILE,$name) or die "can't open this file";
		$nameout = (split /\./,$name)[0];
		$nameout = $nameout.".p2n";
		$nameout2 = (split /\./,$name)[0];
		$nameout2 = $nameout2."-out.pdb";
		$nameout3 = (split /\./,$name)[0];
		$nameout3 = $nameout3."-gam.inp";
		$nameout4 = (split /\./,$name)[0];
		$nameout4 = $nameout4."-pcg.inp";
		$nameout5 = (split /\./,$name)[0];
		$nameout5 = $nameout5."-gau.com";
		while(<FILE>) {
			if ((/^ATOM/ig) || (/^HETATM/ig)) {	# Take only the columns one is interested in
				($element,$rest,$nb,$x,$y,$z) = unpack("x12 A4 x1 A3 x2 A4 x4 A8 A8 A8",$_);
				$atom = SEARCH_ATOM ($element,$nbatom);
				$rad = $Radius{$atom};
				$tab [$nbatom][0] = $atom;
				$tab [$nbatom][1] = $rad;
				$tab [$nbatom][2] = $x;
				$tab [$nbatom][3] = $y;
				$tab [$nbatom][4] = $z;
				$tab [$nbatom][5] = $rest;
				$tab [$nbatom][6] = $nb;
				$tab [$nbatom][7] = $element;
				$nbatom++;
			}
		}
		close(FILE);
		print "  Replace atom names by names used by R.E.D.-III.x\/IV... \t\t[Done]\n";
		$nb_H_bond_H = MATRIX_CONNECTION($nbatom,\@tab,\@matrix,\@H_bond_H,$nb_H_bond_H);
		PUT_CT($nbatom,\@tab,\@matrix);
		if ($sort eq "ON") {
			SORT_BY_NUMBER($nbatom,\@tab);
			SORT_LETTER_BY_NUMBER($nbatom,\@tab);
		}
		$nb_H_bond_H = MATRIX_CONNECTION($nbatom,\@tab,\@matrix,\@H_bond_H,$nb_H_bond_H);
		if ($notlig eq "ON") { TERMINAL_ATOM ($nbatom,\@tab,\@matrix); }
		print "  Establish new connectivities between atoms...\t\t\t\t[Done]\n";
		WRITE_FILE_OUT($name,$nameout,$nameout2,\@tab,$nbatom);
		WRITE_INPUT_GG($nameout3,$nameout4,$nameout5,\@tab,$nbatom);
		WRITE_CONNECT($nameout,$nameout2,$nbatom,\@matrix);
		INFORMATION($name,$nbatom,\@tab,\@matrix,\@H_bond_H,$nb_H_bond_H);
		if ($debug eq "ON") { Printab($nbatom,\@tab,\@matrix) };
		print "  QM program inputs generated by Ante_R.E.D. only provides starting points: \n";
		print "  They should be manually adapted to consider the specificity of each case... \n\n";

		print " ****************************************************************************\n";
		print " ** Do you need a new feature which is not yet implemented in Ante_R.E.D.? **\n";
		print " **               contact the q4md force field tools team \@                **\n";
		print " **                 contact\@q4md-forcefieldtools.org                       **\n";
		print " **          Regularly look for bug fixes at the R.E.D. Home Page          **\n";
		print " **                                 ----                                   **\n";
		print " **  Visit R.E.D. Server \@ http://q4md-forcefieldtools.org/REDS/ as well   **\n";
		print " **                                 ----                                   **\n";
		print " **      Please, submit your force field library(ies) to R.E.DD.B. \@       **\n";
		print " **                 http://q4md-forcefieldtools.org/REDDB/                 **\n";
		print " **      to freely share your results with the scientific community        **\n";
		print " **                                 ----                                   **\n";
		print " **           Do you need help about the q4md force field tools?           **\n";
		print " **        Please, use the q4md-forcefieldtools.org mailing list \@         **\n";
		print " **                http://lists.q4md-forcefieldtools.org/                  **\n";
		print " ****************************************************************************\n";
		print " **                   We thank you for using Ante_R.E.D.                   **\n";
		print " ****************************************************************************\n\n";
	}
	else { print "\t  The PDB input file requested does not exist or is not provided! \n\n"; exit(0); }

