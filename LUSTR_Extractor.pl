#!/usr/local/bin/perl -w
### extract reads from pair-end .sam/.bam to generate .fastq # Bookmarks: 0,0 0,0 0,0 0,384

sub complementary { ### return complementary sequence of input
    my $temp=reverse($_[0]);
    $temp=~tr/ATCGatcg/TAGCtagc/;
    return $temp;
}

sub ProcessFlag { ### explain sam format FLAG (column 2) $_[0], return "priority(p/s)_pair(1/2/0)_strand(+/-)_map(y/n)"
    my $num=$_[0];
    my @tempdata=();
    while ($num>0) {
        $tempdata[$#tempdata+1]=$num%2;
        $num=int($num/2);
    }

    my $priority=""; ### "p"=primary, "s"=secondary or supplementary
    if ($#tempdata<8) {
        $priority="p";
    } elsif ($tempdata[8]==1) {
        $priority="s";
    } elsif ($#tempdata<11) {
        $priority="p";
    } elsif ($tempdata[11]==1) {
        $priority="s";
    } else {
        $priority="p";
    }

    my $pair=""; ### 0=unpaired
    if ($#tempdata==-1) {
        $pair="0";
    } elsif ($tempdata[0]==0) {
        $pair="0";
    } elsif ($#tempdata<6) {
        print ("Warning: flag $_[0] unexplainable\n");
        $pair="0";
    } elsif ($#tempdata==6) {
        $pair="1";
    } elsif (($tempdata[6]==1)and($tempdata[7]==0)) {
        $pair="1";
    } elsif (($tempdata[6]==0)and($tempdata[7]==1)) {
        $pair="2";
    } else {
        print ("Warning: flag $_[0] unexplainable\n");
        $pair="0";
    }

    my $strand=""; ### no matter mapped or unmapped
    if ($#tempdata<4) {
        $strand="+";
    } elsif ($tempdata[4]==1) {
        $strand="-";
    } else {
        $strand="+";
    }

    my $map=""; ### "y"=mapped, "n"=unmapped
    if ($#tempdata<2) {
        $map="y";
    } elsif ($tempdata[2]==1) {
        $map="n";
    } else {
        $map="y";
    }
    my $temp=join "_",($priority,$pair,$strand,$map);
    return ($temp);
}

sub sortarray {
    my $mid=int(($_[0]+$_[1])/2);
    if (($_[1]-$_[0])<=1) {
        if ($ID[$Order[$_[0]]] gt $ID[$Order[$_[1]]]) {
            ($Order[$_[0]],$Order[$_[1]])=($Order[$_[1]],$Order[$_[0]]);
        }
        return;
    }
    &sortarray($_[0],$mid);
    if (($_[1]-$_[0])>2) {
        &sortarray($mid+1,$_[1]);
    }
    my $i=$_[0];
    my $j=$mid+1;
    my @tempOrder=();
    while (($i<=$mid)and($j<=$_[1])) {
        if ($ID[$Order[$i]] gt $ID[$Order[$j]]) {
            push @tempOrder, $Order[$j];
            $j++;
        } else {
            push @tempOrder, $Order[$i];
            $i++;
        }
    }
    if ($j>($_[1])) {
        while ($i<=$mid) {
            push @tempOrder, $Order[$i];
            $i++;
        }
    } else {
        while ($j<=$_[1]) {
            push @tempOrder, $Order[$j];
            $j++;
        }
    }
    for ($i=$_[0];$i<=$_[1];$i++) {
        $Order[$i]=$tempOrder[$i-$_[0]];
    }
}

### read parameters ###
my $i=0;
@Para=();
### Para[0] .sam or .bam no matter sorted or not when raw .fastq not available
### Para[1] output .fastq
### Para[2] max reads num of each intermediate file (10M reads ~ 5G file)
### will generate a set of intermediate files, total size = final output Para[1]
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-h')or($ARGV[$i] eq '--help')) {
        print "Instruction:\n";
        print "\tThe Extractor module extracts raw read pairs from given .sam/.bam files.\n";
        print "\tIt can process both sorted or unsorted .sam/.bam files.\n";
        print "\tIt outputs a single .fastq file written by read pairs.\n";
        print "\tIn case the .sam/.bam is incomplete due to filteration of unmapped reads, a warning message will be given.\n";
        print "Usage Sample:\n";
        print "\tperl LUSTR_Extractor.pl -i <file_name.bam> -o <file_name.fastq>\n";
        print "Options:\n";
        print "\t--help/-h\t\tPrint usage instructions\n";
        print "\t--input/-i <file>\tInput .sam/.bam file\n";
        print "\t--output/-o <file>\tOutput .fastq file of raw read pairs\n";
        print "\t--max/-m <value>\tMaximum pair number of intermediate files, which will be automatically deleted when finish\n";
        print "\t\t\t\tThis setting affects memory usage. Setting of 10,000,000 usually costs 5~10G memory\n";
        print "\t\t\t\t(Default = 10000000)\n";
        exit;
    }
}
$Para[0]='';
$Para[1]='';
$Para[2]=10000000;
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-i')or($ARGV[$i] eq '--input')) {
        $i++;
        $Para[0]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-o')or($ARGV[$i] eq '--output')) {
        $i++;
        $Para[1]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-m')or($ARGV[$i] eq '--max')) {
        $i++;
        $Para[2]=$ARGV[$i];
    }
}
if (($Para[0] eq '')or($Para[1] eq '')) {
    print ("Input/Output not appointed. Stop running.\n");
    print "Please use --help or -h to see instructions.\n";
    exit;
}

### main ###
my $line="";
my @tempdata=();
my $temp="";
@ID=();
@Seq=();
@Qual=();
@Order=();
my $count=0;
my $unmap=0;
my $filenum=0;
my $filename="";
my $index=0;
if ($Para[0]=~/\.sam$/) { ### extract reads and store reads into intermediate files, each file is sorted by reads name
    open (INPUT_DATA, "<$Para[0]") or die "Couldn't open: $!";
    if (-s INPUT_DATA) {
        while (1) {
            chomp ($line=<INPUT_DATA>);
            if (!($line=~/^\@/)) {
                @tempdata=();
                @tempdata=split /\s+/, $line;
                $temp=&ProcessFlag($tempdata[1]);
                if ($temp=~/_n$/) {
                    $unmap=1;
                }
                if ($temp=~/^p_/) {
                    if ($temp=~/_0_/) {
                        print ("Warning: $tempdata[0] unpaired, possibly because alignment was not done in pair-end mode\n");
                    } else {
                        $count++;
                        $Order[$count]=$count;
                        if ($temp=~/_-_/) {
                            $Seq[$count]=&complementary($tempdata[9]);
                            $Qual[$count]=reverse($tempdata[10]);
                        } else {
                            $Seq[$count]=$tempdata[9];
                            $Qual[$count]=$tempdata[10];
                        }
                        if ($temp=~/_1_/) {
                            $ID[$count]=join "",('@',$tempdata[0],'/1');
                        } elsif ($temp=~/_2_/) {
                            $ID[$count]=join "",('@',$tempdata[0],'/2');
                        }
                        if ($count==$Para[2]) { ### sort and export
                            &sortarray(1,$count);
                            $filenum++;
                            $filename=join "",($Para[1],".im",$filenum);
                            open (OUTPUT_INTER, ">$filename") or die "Couldn't open: $!";
                            for ($i=1;$i<=$count;$i++) {
                                $index=$Order[$i];
                                print OUTPUT_INTER ("$ID[$index]\n");
                                print OUTPUT_INTER ("$Seq[$index]\n");
                                print OUTPUT_INTER ("$Qual[$index]\n");
                            }
                            close (OUTPUT_INTER);
                            print ("$filenum group of $count reads sorted...\n");
                            $count=0;
                            @ID=();
                            @Seq=();
                            @Qual=();
                            @Order=();
                        }
                    }
                }
            }
            if (eof INPUT_DATA) {
                last;
            }
        }
    }
    close (INPUT_DATA);
} elsif ($Para[0]=~/\.bam$/) { ### extract reads and store reads into intermediate files, each file is sorted by reads name
    open INPUT_DATA, "samtools view -h $ARGV[0] |";
    while (<INPUT_DATA>) {
        chomp;
        @tempdata=();
        @tempdata=split /\s+/;
        if (!($tempdata[0]=~/^\@/)) {
            $temp=&ProcessFlag($tempdata[1]);
            if ($temp=~/_n$/) {
                $unmap=1;
            }
            if ($temp=~/^p_/) {
                if ($temp=~/_0_/) {
                    print ("Warning: $tempdata[0] unpaired, possibly because alignment was not done in pair-end mode\n");
                } else {
                    $count++;
                    $Order[$count]=$count;
                    if ($temp=~/_-_/) {
                        $Seq[$count]=&complementary($tempdata[9]);
                        $Qual[$count]=reverse($tempdata[10]);
                    } else {
                        $Seq[$count]=$tempdata[9];
                        $Qual[$count]=$tempdata[10];
                    }
                    if ($temp=~/_1_/) {
                        $ID[$count]=join "",('@',$tempdata[0],'/1');
                    } elsif ($temp=~/_2_/) {
                        $ID[$count]=join "",('@',$tempdata[0],'/2');
                    }
                    if ($count==$Para[2]) { ### sort and export
                        &sortarray(1,$count);
                        $filenum++;
                        $filename=join "",($Para[1],".temp",$filenum);
                        open (OUTPUT_INTER, ">$filename") or die "Couldn't open: $!";
                        for ($i=1;$i<=$count;$i++) {
                            $index=$Order[$i];
                            print OUTPUT_INTER ("$ID[$index]\n");
                            print OUTPUT_INTER ("$Seq[$index]\n");
                            print OUTPUT_INTER ("$Qual[$index]\n");
                        }
                        close (OUTPUT_INTER);
                        print ("Group $filenum $count reads sorted...\n");
                        $count=0;
                        @ID=();
                        @Seq=();
                        @Qual=();
                        @Order=();
                    }
                }
            }
        }
    }
    close (INPUT_DATA);
} else {
    print ("Error: $Para[0] is not .sam or .bam file\n");
}

if ((($Para[0]=~/\.sam$/)or($Para[0]=~/\.bam$/))and($count>0)) { ### export last group of reads
    &sortarray(1,$count);
    $filenum++;
    $filename=join "",($Para[1],".temp",$filenum);
    open (OUTPUT_INTER, ">$filename") or die "Couldn't open: $!";
    for ($i=1;$i<=$count;$i++) {
        $index=$Order[$i];
        print OUTPUT_INTER ("$ID[$index]\n");
        print OUTPUT_INTER ("$Seq[$index]\n");
        print OUTPUT_INTER ("$Qual[$index]\n");
    }
    close (OUTPUT_INTER);
    print ("Group $filenum $count reads sorted...\n");
}

$count=0;
@ID=();
@Seq=();
@Qual=();
@Order=();
my $command="";
my $idpre="";
my $idnow="";
my $id1="";
my $id1o="";
my $seq1="";
my $qual1="";
if (0==$filenum) {
    print ("No reads extracted, stop\n");
} else {
    print ("Merging all $filenum groups...\n");
    open (OUTPUT_RESULT, ">$Para[1]") or die "Couldn't open: $!";
    for ($i=1;$i<=$filenum;$i++) { ### read first reads from each file
        $filename=join "",($Para[1],".temp",$i);
        $temp=join "",("FH",$i);
        open ($temp, "<$filename") or die "Couldn't open: $!";
        chomp ($ID[$i]=<$temp>);
        chomp ($Seq[$i]=<$temp>);
        chomp ($Qual[$i]=<$temp>);
        $Order[$i]=$i;
    }
    &sortarray(1,$filenum); ### sort all first reads
    while ($ID[$Order[1]] ne "X") { ### process the reads with the smallest ID until all files reach end
        if ($ID[$Order[1]] eq $idpre) {
            print ("Warning: $ID[$Order[1]] has multiple primary alignment\n");
        } elsif ($ID[$Order[1]]=~/1$/) {
            $idpre=$ID[$Order[1]];
            $id1o=$ID[$Order[1]];
            $id1=substr($ID[$Order[1]],0,length($ID[$Order[1]])-1);
            $seq1=$Seq[$Order[1]];
            $qual1=$Qual[$Order[1]];
        } elsif ($ID[$Order[1]]=~/2$/) {
            $idpre=$ID[$Order[1]];
            $idnow=substr($ID[$Order[1]],0,length($ID[$Order[1]])-1);
            if ($idnow eq $id1) {
                $count++;
                print OUTPUT_RESULT ("$id1o\n");
                print OUTPUT_RESULT ("$seq1\n");
                print OUTPUT_RESULT ("+\n");
                print OUTPUT_RESULT ("$qual1\n");
                print OUTPUT_RESULT ("$ID[$Order[1]]\n");
                print OUTPUT_RESULT ("$Seq[$Order[1]]\n");
                print OUTPUT_RESULT ("+\n");
                print OUTPUT_RESULT ("$Qual[$Order[1]]\n");
            }
        }
        $temp=join "",("FH",$Order[1]); ### replace the processed reads with the next from the same file
        if (eof $temp) {
            $ID[$Order[1]]="X"; ### "X" is defined greater than any reads ID
        } else {
            chomp ($ID[$Order[1]]=<$temp>);
            chomp ($Seq[$Order[1]]=<$temp>);
            chomp ($Qual[$Order[1]]=<$temp>);
        }
        for ($i=1;$i<$filenum;$i++) { ### insert the new reads
            if ($ID[$Order[$i+1]] eq "X") {
                last;
            } elsif ($ID[$Order[$i]] eq "X") {
                ($Order[$i],$Order[$i+1])=($Order[$i+1],$Order[$i]);
            } elsif ($ID[$Order[$i]] gt $ID[$Order[$i+1]]) {
                ($Order[$i],$Order[$i+1])=($Order[$i+1],$Order[$i]);
            } else {
                last;
            }
        }
    }
    for ($i=1;$i<=$filenum;$i++) { ### close and delete intermediate files
        $filename=join "",($Para[1],".temp",$i);
        $temp=join "",("FH",$i);
        close ($temp);
        $command=join " ",("rm",$filename);
        system($command);
    }
    close (OUTPUT_RESULT);
    print ("Job done. Total $count pairs extracted into $Para[1]\n");
    if (0==$unmap) {
        print ("Warning: NO unmapped reads found in $Para[0]\n");
    }
}
############

exit;

