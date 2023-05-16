#!/usr/local/bin/perl -w
### generate reference .fa for STR repeat and flanking sequences # Bookmarks: 0,0 0,0 0,0 0,192

sub ratioN { ### calculate ratio of N in $_[0]
    my $i=0;
    my $count=0;
    my $letter="";
    for ($i=0;$i<=length($_[0])-1;$i++) {
          $letter=substr ($_[0], $i, 1);
          if ($letter eq "N") {
               $count++;
          }
    }
    $count=$count/(length($_[0]));
    return $count;
}

sub complementary { ### return complementary sequence of input
    my $temp=reverse($_[0]);
    $temp=~tr/ATCGatcg/TAGCtagc/;
    return $temp;
}

sub AddUnit {
    ### check if $_[0] (with/without complementary) can be transformed to any exist index of %UniUnit
    ### if not, add $_[0] to %UniUnit index list
    my $i=0;
    my $temp="";
    for ($i=length($_[0])-1;$i>=0;$i--) {
        $temp=join "",(substr($_[0],$i,length($_[0])-$i),substr($_[0],0,$i));
        if (defined($UniUnit{$temp})) {
            return;
        }
    }
    my $unit=&complementary($_[0]);
    for ($i=length($unit)-1;$i>=0;$i--) {
        $temp=join "",(substr($unit,$i,length($unit)-$i),substr($unit,0,$i));
        if (defined($UniUnit{$temp})) {
            return;
        }
    }
    $UniUnit{$_[0]}=1;
}

### read parameters ###
my $i=0;
@Para=();
### Para[0]=main result file from Finder
### Para[1] minimum length of _Repeat & _Fake, suggestion = expected reads length + ~50
### generate .fa with flank sequences, repeat sequences & fake repeat sequences
### output .fa 60nt per line
### Para[2] optional to specific output file name, otherwise output file will be Para[0].fa
### Para[3] N ratio to discard sequences, suggestion = 0.3
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-h')or($ARGV[$i] eq '--help')) {
        print "Instruction:\n";
        print "\tThe RefCreater module processes the repeat and flanking sequences of STRs to generate both real and artificial references to remap raw reads to.\n";
        print "\tIt outputs a fasta format file with 60nt lines. Sequences with too many Ns can be filtered according to user setting.\n";
        print "Usage Sample:\n";
        print "\tperl LUSTR_RefCreater.pl -i <STRsequences.txt> -o <STRreferences.fa>\n";
        print "Options:\n";
        print "\t--help/-h\t\tPrint usage instructions\n";
        print "\t--input/-i <file>\tThe STR sequences obtained by the Finder module. See format details in README.txt Step (1.2)\n";
        print "\t--output/-o <file>\tOutput file in fasta format. If not appointed, LUSTR will generate a file with \".fa\" extension after the input file name\n";
        print "\t--length/-l <value>\tMinimum length when generating the references for repeats. A setting with 50nt longer than read length is recommended\n";
        print "\t\t\t\t(Default = 200)\n";
        print "\t--ratio/-r <value>\tA value between 0 and 1 as the threshold for the ratio of Ns. Sequences with Ns beyond this ratio will be discarded\n";
        print "\t\t\t\t(Default = 0.3)\n";
        exit;
    }
}
$Para[0]='';
$Para[1]=200;
$Para[2]='';
$Para[3]=0.3;
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-r')or($ARGV[$i] eq '--ratio')) {
        $i++;
        $Para[3]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-i')or($ARGV[$i] eq '--input')) {
        $i++;
        $Para[0]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-l')or($ARGV[$i] eq '--length')) {
        $i++;
        $Para[1]=$ARGV[$i];
    }
}
$Para[2]=join '',($Para[0],'.fa');
for ($i=0;$i<=$#ARGV;$i++) {
    if ($ARGV[$i] eq '--output') {
        $i++;
        $Para[2]=$ARGV[$i];
    }
}
if (($Para[0] eq '')or($Para[2] eq '')) {
    print ("Input/Output not appointed. Stop running.\n");
    print "Please use --help or -h to see instructions.\n";
    exit;
}

### main ###
my $filename=$Para[2];
open (OUTPUT_RESULT, ">$filename") or die "Couldn't open: $!";

my $limit=$Para[3];
my $line="";
my $id="";
my $temp="";
my $count=0;
%UniUnit=();
open (INPUT_DATA, "<$Para[0]") or die "Couldn't open: $!";
if (-s INPUT_DATA) {
    chomp ($line=<INPUT_DATA>);
    while (1) {
        chomp ($line=<INPUT_DATA>);
        $line=~/\s+/;
        $id=join '',('>',$`);
        $temp=$';
        $temp=~/\s+/;
        &AddUnit($`);
        chomp ($line=<INPUT_DATA>);
        chomp ($line=<INPUT_DATA>);
        if (&ratioN($line)<=$limit) {
            $count++;
            $temp=$line;
            while ((length($temp))<$Para[1]) {
                $temp.=$line;
            }
            $line=$temp;
            $temp=join '',($id,'_Repeat');
            print OUTPUT_RESULT ("$temp\n");
            while ((length($line))>60) {
                $temp=substr($line,0,60);
                print OUTPUT_RESULT ("$temp\n");
                $line=substr($line,60,length($line)-60);
            }
            print OUTPUT_RESULT ("$line\n");
        }
        chomp ($line=<INPUT_DATA>);
        chomp ($line=<INPUT_DATA>);
        if ($line ne "") {
            if (&ratioN($line)<=$limit) {
                $count++;
                $temp=join '',($id,'_UpStream');
                print OUTPUT_RESULT ("$temp\n");
                while ((length($line))>60) {
                    $temp=substr($line,0,60);
                    print OUTPUT_RESULT ("$temp\n");
                    $line=substr($line,60,length($line)-60);
                }
                print OUTPUT_RESULT ("$line\n");
            }
        }
        chomp ($line=<INPUT_DATA>);
        chomp ($line=<INPUT_DATA>);
        if ($line ne '') {
            if (&ratioN($line)<=$limit) {
                $count++;
                $temp=join '',($id,'_DownStream');
                print OUTPUT_RESULT ("$temp\n");
                while ((length($line))>60) {
                    $temp=substr($line,0,60);
                    print OUTPUT_RESULT ("$temp\n");
                    $line=substr($line,60,length($line)-60);
                }
                print OUTPUT_RESULT ("$line\n");
            }
        }
        chomp ($line=<INPUT_DATA>);
        if (eof) {
            last;
        }
    }
}
close (INPUT_DATA);
foreach (keys(%UniUnit)) {
    $line='';
    for ($i=int($Para[1]/(length($_)))+1;$i>=1;$i--) {
        $line.=$_;
    }
    $count++;
    $temp=join '',('>',$_,'_Fake');
    print OUTPUT_RESULT ("$temp\n");
    while ((length($line))>60) {
        $temp=substr($line,0,60);
        print OUTPUT_RESULT ("$temp\n");
        $line=substr($line,60,length($line)-60);
    }
    print OUTPUT_RESULT ("$line\n");
}
close (OUTPUT_RESULT);
print ("Done, generated $filename with $count sequences\n");
############

exit;

