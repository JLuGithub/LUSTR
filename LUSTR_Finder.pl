#!/usr/local/bin/perl -w
### process basic information of STRs to retrieve standarized sequences # Bookmarks: 0,0 0,0 0,0 0,1937

sub read_fai { ### output @FaiInfo, %ChrIndex
    ### [0]=.fai
    my $line="";
    my $i=0;
    my @tempdata=();
    @FaiInfo=();
    %ChrIndex=();
    open (INPUT_DATA, "<$_[0]") or die "Couldn't open: $!";
    if (-s INPUT_DATA) {
        while (1) {
            chomp ($line=<INPUT_DATA>);
            $FaiInfo[$i]=$line;
            @tempdata=();
            @tempdata=split /\s+/, $line;
            $ChrIndex{$tempdata[0]}=$i;
            $i++;
            if (eof) {
                last;
            }
        }
    }
    close (INPUT_DATA);
}

sub complementary { ### return complementary sequence of input
    my $temp=reverse($_[0]);
    $temp=~tr/ATCGatcg/TAGCtagc/;
    return $temp;
}

sub standardizeDNA { ### make sequence ATCGN only
    my $temp=uc($_[0]);
    $temp=~s/[^ATCGN]/N/gi;
    return $temp;
}

sub innerrepeat { ### if $_[0] is composed by multiple shortest string X, return X; otherwise return $_[0]
    my $i=0;
    my $j=0;
    my $result=$_[0];
    my $tempresult="";
    my $temp="";
    my $len=length($_[0]);
    for ($i=1;$i<=int($len/2);$i++) {
        if ($len%$i==0) {
            $temp=substr($_[0],0,$i);
            $tempresult="";
            for ($j=1;$j<=($len/$i);$j++) {
                $tempresult.=$temp;
            }
            if ($tempresult eq $_[0]) {
                $result=$temp;
                last;
            }
        } else {
            next;
        }
    }
    if ($result ne $_[0]) {
        print ("Warning: $StrID repeat unit $_[0] replaced by $result\n");
        print OUTPUT_WARN ("Warning: $StrID repeat unit $_[0] replaced by $result\n");
    }
    return ($result);
}

sub adjust_matchscore { ### return adjusted match score
    ### [0]=length of STR unit [1]=original match score [2]=mismatch score [3]=gap score
    my $max=$_[2];
    if (!(defined($dyscore[$_[0]]))) {
        if ($dynamic eq "no") {
            $dyscore[$_[0]]=$_[1];
        } else {
            if ($max<$_[3]) {
                $max=$_[3];
            }
            if ((($_[0]-1)*$_[1]+$max)<=-0.5) {
                $dyscore[$_[0]]=$_[1];
            } else {
                $dyscore[$_[0]]=(-0.5-$max)/($_[0]-1);
            }
        }
    }
    return ($dyscore[$_[0]]);
}

sub get_seed { ### output $totalinput, get $Seq $StrUnit $Chr $Pos $Center $StrID for &get_seq
    ### [0]=.fa [1]=Para[2]
    my $line="";
    my @tempdata=();
    my @fai=();
    my $start=0;
    my $wantlength=0;
    $Seq="";
    $StrUnit=""; ### unit sequence on genome + strand, can be identical/complementary to input
    $Chr="";
    $Pos=0;
    $StrID="";
    my %checkdup=(); ### check duplicated input
    $totalinput=0;
    $oriinput=0;
    $compleinput=0;
    $bothinput=0;
    $missinput=0;
    $goodinput=0;
    my $templength=0;

    my @read=();
    open (INPUT_DATA2, "<$_[1]") or die "Couldn't open: $!";
    if (-s INPUT_DATA2) {
        while (1) {
            chomp ($line=<INPUT_DATA2>);
            push @read, $line;
            if (eof) {
                last;
            }
        }
    }
    close (INPUT_DATA2);
    $totalinput=$#read+1;

    my $i=0;
    my $seqlimit=20;
    open (INPUT_DATA1, "<$_[0]") or die "Couldn't open: $!";
    if (-s INPUT_DATA1) {
        for ($i=0;$i<=$#read;$i++) {
            $line=$read[$i];
            @tempdata=();
            @tempdata=split /\s+/, $line;
            if ($#tempdata<4) {
                print ("Warning: $tempdata[0] missing information, skipped\n");
                print OUTPUT_WARN ("Warning: $tempdata[0] missing information, skipped\n");
                print OUTPUT_NOTRUN ("$tempdata[0]\n");
            } elsif (defined($checkdup{$tempdata[0]})) {
                print ("Warning: $tempdata[0] duplicated\n");
                print OUTPUT_WARN ("Warning: $tempdata[0] duplicated\n");
            } else {
                if (!defined($ChrIndex{$tempdata[1]})) {
                    print ("Warning: $tempdata[0] Chr $tempdata[1] doesn't exist in reference, skipped\n");
                    print OUTPUT_WARN ("Warning: $tempdata[0] Chr $tempdata[1] doesn't exist in reference, skipped\n");
                    print OUTPUT_NOTRUN ("$tempdata[0]\n");
                } else {
                    @fai=();
                    @fai=split /\s+/, $FaiInfo[$ChrIndex{$tempdata[1]}];
                    if ($tempdata[2]>$tempdata[3]) {
                        ($tempdata[2],$tempdata[3])=($tempdata[3],$tempdata[2]);
                    }
                    if ($tempdata[2]<=0) {
                        $tempdata[2]=1;
                        print ("Warning: $tempdata[0] search range $tempdata[2] out of $tempdata[1], changed\n");
                        print OUTPUT_WARN ("Warning: $tempdata[0] search range $tempdata[2] out of $tempdata[1], changed\n");
                    }
                    if ($tempdata[3]>$fai[1]) {
                        $tempdata[3]=$fai[1];
                        print ("Warning: $tempdata[0] search range $tempdata[3] out of $tempdata[1], changed\n");
                        print OUTPUT_WARN ("Warning: $tempdata[0] search range $tempdata[3] out of $tempdata[1], changed\n");
                    }
                    if (($tempdata[3]-$tempdata[2]+1)<$seqlimit) {
                        print ("Warning: $tempdata[0] search range $tempdata[2] $tempdata[3] too short, extended\n");
                        print OUTPUT_WARN ("Warning: $tempdata[0] search range $tempdata[2] $tempdata[3] too short, extended\n");
                        while ((($tempdata[3]-$tempdata[2]+1)<$seqlimit)and(($tempdata[2]>1)or($tempdata[3]<$fai[1]))) {
                            if (1==$tempdata[2]) {
                                $tempdata[3]=$tempdata[2]+$seqlimit-1;
                            } elsif ($tempdata[3]==$fai[1]) {
                                $tempdata[2]=$tempdata[3]-$seqlimit+1;
                            } else {
                                $tempdata[2]=$tempdata[2]-int(($seqlimit-$tempdata[3]+$tempdata[2]-1)/2)-1;
                                $tempdata[3]=$tempdata[3]+int(($seqlimit-$tempdata[3]+$tempdata[2]-1)/2)+1;
                            }
                            if ($tempdata[2]<=0) {
                                $tempdata[2]=1;
                            }
                            if ($tempdata[3]>$fai[1]) {
                                $tempdata[3]=$fai[1];
                            }
                        }
                    }
                    $checkdup{$tempdata[0]}=1;
                    $StrID=$tempdata[0];
                    $StrUnit=&standardizeDNA($tempdata[4]);
                    $StrUnit=&innerrepeat($StrUnit);
                    $Chr=$tempdata[1];
                    $start=$tempdata[2];
                    $Pos=$start;
                    $wantlength=$tempdata[3]-$start+1;
                    $templength=0;
                    $Seq="";
                    seek (INPUT_DATA1,$fai[2]+int(($start-1)/$fai[3])*$fai[4]+($start-1)%$fai[3],0);
                    while ($templength<$wantlength) {
                        chomp ($line=<INPUT_DATA1>);
                        if ($line ne "") {
                            $line=&standardizeDNA($line);
                            if (($templength+length($line))<=$wantlength) {
                                $templength=$templength+length($line);
                                $Seq.=$line;
                            } else {
                                $line=substr($line,0,$wantlength-$templength);
                                $templength=$templength+length($line);
                                $Seq.=$line;
                            }
                        }
                    } ### $Seq will be used in &get_seq to find seed repeats, $Pos turn pos in $Seq to pos in genome
                    &get_seq($flanklength,$match_score,$mis_score,$gap_score,$limit_score);
                }
            }
        }
    }
    close (INPUT_DATA1);
}

sub get_seq { ### output StrSeq, UpFlank, DownFlank, $StrPos for each run into files, count $goodinput $oriinput $compleinput $bothinput for log
    ### [0]=flank size for retrieve flanking sequences
    ### [1]=match score, [2]=mismatch score, [3]=gap score, [4]=limit to stop for STR extendsion
    ### $Seq $StrUnit $Chr $Pos $StrID are from &get_seed

    my $dymatch_score=&adjust_matchscore(length($StrUnit),$_[1],$_[2],$_[3]);
    my $exgap=abs($_[2]);
    if ($exgap<(abs($_[3]))) {
        $exgap=abs($_[3]);
    }
    $exgap=int($exgap/$dymatch_score);
    my $exFL=0;
    my $seed="";
    my $minrepeat=0;
    if ((length($StrUnit))==0) {
        print ("Warning: $StrID missing repeat unit information, skipped\n");
        print OUTPUT_WARN ("Warning: $StrID missing repeat unit information, skipped\n");
        print OUTPUT_NOSTR ("$StrID\n");
        return;
    } elsif ((length($StrUnit))<=1) {
        $seed=join "",($StrUnit,$StrUnit,$StrUnit,$StrUnit,$StrUnit);
        $minrepeat=10;
        $goodinput++;
    } elsif ((length($StrUnit))<=2) {
        $seed=join "",($StrUnit,$StrUnit,$StrUnit);
        $minrepeat=5;
        $goodinput++;
    } elsif ((length($StrUnit))<=4) {
        $seed=join "",($StrUnit,$StrUnit);
        $minrepeat=3;
        $goodinput++;
    } elsif ((length($StrUnit))<=10) {
        $seed=$StrUnit;
        $minrepeat=2;
        $goodinput++;
    } else {
        $seed=$StrUnit;
        $minrepeat=2;
        print ("Warning: This program may NOT run properly when repeat unit is longer than 10nt\n"); ### STR with long unit needs further test
        print OUTPUT_WARN ("Warning: This program may NOT run properly when repeat unit is longer than 10nt\n");
        $goodinput++;
    }
    my $seedcomple=&complementary($seed);

    my $i=0;
    my $temp="";
    my $tempcompare=0;
    @find=();
    for ($i=0;$i<=(length($Seq)-length($seed));$i++) { ### search $seed
        $temp=substr($Seq,$i,length($seed));
        if ($temp eq $seed) {
            push @find, ($Pos+$i);
        }
    }
    @findcomple=();
    if (($seedcomple ne $seed)and($check_comple eq "yes")) { ### search complementary sequence of $seed when they are different
        for ($i=0;$i<=(length($Seq)-length($seedcomple));$i++) {
            $temp=substr($Seq,$i,length($seedcomple));
            if ($temp eq $seedcomple) {
                push @findcomple, ($Pos+$i);
            }
        }
    }

    my @density=();
    my @densitycomple=();
    my $distance=0;
    my $j=0;
    for ($i=0;$i<=$#find;$i++) { ### get @density
        $density[$i]=0;
        for ($j=0;$j<=$#find;$j++) {
            $distance=abs($find[$i]-$find[$j])-length($seed)+1;
            if ($distance<1) {
                $distance=1;
            }
            $density[$i]=$density[$i]+1/$distance;
        }
    }
    for ($i=1;$i<=$#find;$i++) { ### move $find[] with max $density[] to [0]
        if ($density[$i]>$density[0]) {
            ($find[0],$find[$i])=($find[$i],$find[0]);
            ($density[0],$density[$i])=($density[$i],$density[0]);
        }
    }
    for ($i=0;$i<=$#findcomple;$i++) { ### get @densitycomple
        $densitycomple[$i]=0;
        for ($j=0;$j<=$#findcomple;$j++) {
            $distance=abs($findcomple[$i]-$findcomple[$j])-length($seedcomple)+1;
            if ($distance<1) {
                $distance=1;
            }
            $densitycomple[$i]=$densitycomple[$i]+1/$distance;
        }
    }
    for ($i=1;$i<=$#findcomple;$i++) { ### move $findcomple[] with max $densitycomple[] to [0]
        if ($densitycomple[$i]>$densitycomple[0]) {
            ($findcomple[0],$findcomple[$i])=($findcomple[$i],$findcomple[0]);
            ($densitycomple[0],$densitycomple[$i])=($densitycomple[$i],$densitycomple[0]);
        }
    }

    my $index=-1;
    my $indexcomple=-1;
    my $upext=0;
    my $downext=0;
    my $upextcomple=0;
    my $downextcomple=0;
    my $upextmore=0;
    my $downextmore=0;
    my $upextcomplemore=0;
    my $downextcomplemore=0;
    my $maxgap=0;
    my $maxgapcomple=0;
    if ($#find>=0) {
        while ($density[0]>0) { ### try extend best $seed ($find[0])
            $upext=&up_extension($Chr,$find[0],$StrUnit,$dymatch_score,$_[2],$_[3],$_[4]);
            $downext=&down_extension($Chr,$find[0]+length($seed)-1,$StrUnit,$dymatch_score,$_[2],$_[3],$_[4]);
            if ((length($seed)+$upext+$downext)>=($minrepeat*length($StrUnit))) { ### pass test
                $index=0;
                if ($downgap<$upgap) {
                    $maxgap=$upgap;
                } else {
                    $maxgap=$downgap;
                }
                if ($FL_limit eq "unlimited") {
                    $upextmore=&up_extensionmore($Chr,$find[0],$StrUnit,$dymatch_score,$_[2],$_[3],$_[4],$allow_score-3*$dymatch_score);
                    $downextmore=&down_extensionmore($Chr,$find[0]+length($seed)-1,$StrUnit,$dymatch_score,$_[2],$_[3],$_[4],$allow_score-3*$dymatch_score);
                } elsif ($FL_limit<=0) {
                    $upextmore=$upext;
                    $downextmore=$downext;
                } else {
                    $upextmore=&up_extensionmore_wl($Chr,$find[0],$StrUnit,$dymatch_score,$_[2],$_[3],$_[4],$allow_score-3*$dymatch_score,$upext+$FL_limit);
                    $downextmore=&down_extensionmore_wl($Chr,$find[0]+length($seed)-1,$StrUnit,$dymatch_score,$_[2],$_[3],$_[4],$allow_score-3*$dymatch_score,$downext+$FL_limit);
                }
                last;
            } else {
                $density[0]=0;
                for ($i=1;$i<=$#find;$i++) { ### move $find[] with max $density[] to [0]
                    if ($density[$i]>$density[0]) {
                        ($find[0],$find[$i])=($find[$i],$find[0]);
                        ($density[0],$density[$i])=($density[$i],$density[0]);
                    }
                }
            }
        }
    }
    if ($#findcomple>=0) {
        while ($densitycomple[0]>0) { ### try extend best $seedcomple ($findcomple[0])
            $upextcomple=&up_extension($Chr,$findcomple[0],&complementary($StrUnit),$dymatch_score,$_[2],$_[3],$_[4]);
            $downextcomple=&down_extension($Chr,$findcomple[0]+length($seedcomple)-1,&complementary($StrUnit),$dymatch_score,$_[2],$_[3],$_[4]);
            if ((length($seedcomple)+$upextcomple+$downextcomple)>=($minrepeat*length($StrUnit))) { ### pass test
                $indexcomple=0;
                if ($downgap<$upgap) {
                    $maxgapcomple=$upgap;
                } else {
                    $maxgapcomple=$downgap;
                }
                if ($FL_limit eq "unlimited") {
                    $upextcomplemore=&up_extensionmore($Chr,$findcomple[0],&complementary($StrUnit),$dymatch_score,$_[2],$_[3],$_[4],$allow_score-3*$dymatch_score);
                    $downextcomplemore=&down_extensionmore($Chr,$findcomple[0]+length($seedcomple)-1,&complementary($StrUnit),$dymatch_score,$_[2],$_[3],$_[4],$allow_score-3*$dymatch_score);
                } elsif ($FL_limit<=0) {
                    $upextcomplemore=$upextcomple;
                    $downextcomplemore=$downextcomple;
                } else {
                    $upextcomplemore=&up_extensionmore_wl($Chr,$findcomple[0],&complementary($StrUnit),$dymatch_score,$_[2],$_[3],$_[4],$allow_score-3*$dymatch_score,$upextcomple+$FL_limit);
                    $downextcomplemore=&down_extensionmore_wl($Chr,$findcomple[0]+length($seedcomple)-1,&complementary($StrUnit),$dymatch_score,$_[2],$_[3],$_[4],$allow_score-3*$dymatch_score,$downextcomple+$FL_limit);
                }
                last;
            } else {
                $densitycomple[0]=0;
                for ($i=1;$i<=$#findcomple;$i++) { ### move $findcomple[] with max $densitycomple[] to [0]
                    if ($densitycomple[$i]>$densitycomple[0]) {
                        ($findcomple[0],$findcomple[$i])=($findcomple[$i],$findcomple[0]);
                        ($densitycomple[0],$densitycomple[$i])=($densitycomple[$i],$densitycomple[0]);
                    }
                }
            }
        }
    }

    my $judge="";
    if (($index==-1)and($indexcomple==-1)) {
        $judge="no";
        print ("Warning: $StrID can't find STR sequence, skipped\n");
        print OUTPUT_WARN ("Warning: $StrID can't find STR sequence, skipped\n");
        print OUTPUT_NOSTR ("$StrID\n");
        $missinput++;
    } elsif ($indexcomple==-1) { ### retrieve StrSeq $StrPos UpFlank DownFlank using $find[$index] $seed $upext $downext
        $judge="ori";
        $oriinput++;
    } elsif ($index==-1) { ### retrieve StrSeq $StrPos UpFlank DownFlank using $findcomple[$indexcomple] $seedcomple $upextcomple $downextcomple
        $judge="comple";
        $StrUnit=&complementary($StrUnit);
        $compleinput++;
    } elsif ($density[$index]>$densitycomple[$indexcomple]) { ### retrieve StrSeq $StrPos UpFlank DownFlank using $find[$index] $seed $upext $downext
        $judge="ori";
        $oriinput++;
    } elsif ($density[$index]<$densitycomple[$indexcomple]) { ### retrieve StrSeq $StrPos UpFlank DownFlank using $findcomple[$indexcomple] $seedcomple $upextcomple $downextcomple
        $judge="comple";
        $StrUnit=&complementary($StrUnit);
        $compleinput++;
    } else { ### retrieve StrSeq $StrPos UpFlank DownFlank using $find[$index] $seed $upext $downext
        $judge="ori";
        print ("Warning: $StrID can't determine best STR sequence. Using input repeat.\n");
        print OUTPUT_WARN ("Warning: $StrID can't determine best STR sequence. Using input repeat.\n");
        $bothinput++;
    }

    my @fai=();
    my $start=0;
    my $wantlength=0;
    my $templength=0;
    my $flanksize=$_[0];
    my $line="";
    my $StrPos=0;
    my $ratio=0;
    my $tempstr="";
    if ($judge eq "") {
        print ("Error: Unknown No.1\n");
    } elsif ($judge ne "no") {
        if ($judge eq "ori") {
            $StrPos=$find[$index]-$upext;
        } elsif ($judge eq "comple") {
            $StrPos=$findcomple[$indexcomple]-$upextcomple;
        }
        @fai=();
        @fai=split /\s+/, $FaiInfo[$ChrIndex{$Chr}];

        $start=$StrPos;
        if ($judge eq "ori") {
            $wantlength=$upext+length($seed)+$downext;
            $temp=join "",("$StrID $StrUnit $Chr",":",$StrPos,"-",$StrPos+$wantlength-1," MaxGap:",$maxgap+$exgap," 5FL:",$upextmore-$upext+$exFL," 3FL:",$downextmore-$downext+$exFL);
        } elsif ($judge eq "comple") {
            $wantlength=$upextcomple+length($seedcomple)+$downextcomple;
            $temp=join "",("$StrID $StrUnit $Chr",":",$StrPos,"-",$StrPos+$wantlength-1," MaxGap:",$maxgapcomple+$exgap," 5FL:",$upextcomplemore-$upextcomple+$exFL," 3FL:",$downextcomplemore-$downextcomple+$exFL);
        }
        print OUTPUT_RESULT ("$temp\n");
        print OUTPUT_RESULT ("STR Sequence:\n");
        $templength=0;
        seek (INPUT_DATA1,$fai[2]+int(($start-1)/$fai[3])*$fai[4]+($start-1)%$fai[3],0);
        while ($templength<$wantlength) {
            chomp ($line=<INPUT_DATA1>);
            if ($line ne "") {
                $line=&standardizeDNA($line);
                if (($templength+length($line))<=$wantlength) {
                    $templength=$templength+length($line);
                } else {
                    $line=substr($line,0,$wantlength-$templength);
                    $templength=$templength+length($line);
                }
                print OUTPUT_RESULT ($line);
            }
        }
        print OUTPUT_RESULT ("\n");

        $start=$StrPos-$flanksize;
        if ($start<1) {
            $start=1;
        }
        $wantlength=$StrPos-$start;
        print OUTPUT_RESULT ("UpStream:\n");
        $templength=0;
        seek (INPUT_DATA1,$fai[2]+int(($start-1)/$fai[3])*$fai[4]+($start-1)%$fai[3],0);
        while ($templength<$wantlength) {
            chomp ($line=<INPUT_DATA1>);
            if ($line ne "") {
                $line=&standardizeDNA($line);
                if (($templength+length($line))<=$wantlength) {
                    $templength=$templength+length($line);
                } else {
                    $line=substr($line,0,$wantlength-$templength);
                    $templength=$templength+length($line);
                }
                print OUTPUT_RESULT ($line);
            }
        }
        print OUTPUT_RESULT ("\n");

        if ($judge eq "ori") {
            $start=$find[$index]+length($seed)+$downext;
        } elsif ($judge eq "comple") {
            $start=$findcomple[$indexcomple]+length($seedcomple)+$downextcomple;
        }
        $wantlength=$flanksize;
        if (($start+$wantlength-1)>$fai[1]) {
            $wantlength=$flanksize-$start-$wantlength+1+$fai[1];
        }
        print OUTPUT_RESULT ("DownStream:\n");
        $templength=0;
        seek (INPUT_DATA1,$fai[2]+int(($start-1)/$fai[3])*$fai[4]+($start-1)%$fai[3],0);
        while ($templength<$wantlength) {
            chomp ($line=<INPUT_DATA1>);
            if ($line ne "") {
                $line=&standardizeDNA($line);
                if (($templength+length($line))<=$wantlength) {
                    $templength=$templength+length($line);
                } else {
                    $line=substr($line,0,$wantlength-$templength);
                    $templength=$templength+length($line);
                }
                print OUTPUT_RESULT ($line);
            }
        }
        print OUTPUT_RESULT ("\n\n");
    }
    if ($goodinput%10000==0) {
        print ("$goodinput processed...\n");
    }
}

sub up_extension { ### output nt length can be extended to upstream direction
    ### continue to use INPUT_DATA1 as .fa, and information from @FaiInfo, %ChrIndex
    ### [0]=chr [1]=start position to extend (right after extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion
    my @faiext=split /\s+/, $FaiInfo[$ChrIndex{$_[0]}];
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=length($_[2]);$y>=1;$y--) { ### store repeat unit in reversed order
        $repeat[length($_[2])-$y+1]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence in reversed order
    seek (INPUT_DATA1,$faiext[2]+int($_[1]/$faiext[3])*$faiext[4]+$_[1]%$faiext[3],0);
    if (($_[1]%$faiext[3])==0) {
        seek (INPUT_DATA1,-1,1);
    }
    my $posnow=$_[1]-1;
    my $refchar="";

    my $gapscore=$_[5];
    my $misscore=$_[4]; ### N misscore will = $misscore*3/4
    my $matchscore=$_[3];
    my $limit=$_[6];
    my $count=0;
    my $result=0; ### max successful extension nt number
    ### score matrix ###
    my @bottom=();
    $bottom[0]=0;
    my @tempold=();
    my @tempnew=();
    my @rightold=();
    $rightold[0]=0;
    my @rightnew=();
    my @gapbool=();
    ####################

    my $keepext="no";
    my $value=0;
    while ($posnow>=1) {
        seek (INPUT_DATA1,-2,1);
        $refchar=getc(INPUT_DATA1);
        if ($refchar ne "\n") {
            $refchar=uc($refchar);
            if (($refchar ne "A")and($refchar ne "C")and($refchar ne "T")and($refchar ne "G")) {
                $refchar="N";
            }
            $count++; ### $count=X (=$_[1]-$posnow) means now trying to align the Xst nt
            $ref[$count]=$refchar;

            $tempnew[0]=$bottom[0];
            for ($y=1;$y<=length($_[2]);$y++) {
                if ($tempnew[$y-1]==$limit) { ### any calculation using a stop match is considered as a stop
                    $tempnew[$y]=$limit;
                } else {
                    $tempnew[$y]=$tempnew[$y-1]+$gapscore;
                }
                if ($tempnew[$y]<=$limit) { ### any score less than $limit is considered as a stop
                    $tempnew[$y]=$limit;
                } elsif ($tempnew[$y]>=0) { ### any score >0 will be reset to 0, so previous extension doesn't affect downstream
                    $tempnew[$y]=0;
                }
            }
            @tempold=@tempnew;
            @tempnew=();
            $bottom[0]=$tempold[length($_[2])];

            for ($x=1;$x<$count;$x++) {
                $tempnew[0]=$bottom[$x];
                for ($y=1;$y<=length($_[2]);$y++) { ### get max score from 3 routes
                    if ($tempnew[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempnew[$y-1]+$gapscore;
                    }
                    $tempnew[$y]=$value;

                    if ($tempold[$y]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempold[$y]+$gapscore;
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempold[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        if (($ref[$x] eq "N")or($repeat[$y] eq "N")) {
                            $value=$tempold[$y-1]+$misscore*0.75;
                        } elsif ($ref[$x] ne $repeat[$y]) {
                            $value=$tempold[$y-1]+$misscore;
                        } else {
                            $value=$tempold[$y-1]+$matchscore;
                        }
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempnew[$y]<=$limit) {
                        $tempnew[$y]=$limit;
                    } elsif ($tempnew[$y]>=0) {
                        $tempnew[$y]=0;
                        $gapbool[$x]=1;
                        if ($x>$result) {
                            $result=$x; ### if any score =0, consider it as successful extension
                        }
                    }
                }
                @tempold=@tempnew;
                @tempnew=();
                $bottom[$x]=$tempold[length($_[2])];
            }
            for ($y=1;$y<=length($_[2]);$y++) {
                $rightold[$#rightold+1]=$tempold[$y];
            }
            @tempold=();
            if ($#rightold!=($count*length($_[2]))) {
                print ("Error: Unknown No.2\n");
            }

            $keepext="no";
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext="yes"; ### if any score >limit, continue to try next nt
                if ($rightnew[0]>=0) {
                    $rightnew[0]=0;
                    $gapbool[$count]=1;
                    $result=$count; ### if any score =0, consider it as successful extension
                }
            }

            for ($y=1;$y<=$count*length($_[2]);$y++) { ### get max score from 3 routes
                if ($rightnew[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightnew[$y-1]+$gapscore;
                }
                $rightnew[$y]=$value;

                if ($rightold[$y]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightold[$y]+$gapscore;
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightold[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    if (($ref[$count] eq "N")or($repeat[$y%(length($_[2]))] eq "N")) {
                        $value=$rightold[$y-1]+$misscore*0.75;
                    } elsif ($ref[$count] ne $repeat[$y%(length($_[2]))]) {
                        $value=$rightold[$y-1]+$misscore;
                    } else {
                        $value=$rightold[$y-1]+$matchscore;
                    }
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightnew[$y]<=$limit) {
                    $rightnew[$y]=$limit;
                } else {
                    $keepext="yes"; ### if any score >limit, continue to try next nt
                    if ($rightnew[$y]>=0) {
                        $rightnew[$y]=0;
                        $gapbool[$count]=1;
                        $result=$count; ### if any score =0, consider it as successful extension
                    }
                }
            }
            @rightold=@rightnew;
            @rightnew=();
            $bottom[$count]=$rightold[$#rightold];

            if ($keepext eq "no") {
                last;
            } else {
                $posnow--;
            }
        }
    }

    my @gapstart=();
    my @gapend=();
    my $i=0;
    my $record="no";
    $upgap=0;
    for ($i=1;$i<=$result;$i++) {
        if (defined($gapbool[$i])) {
            if ($record eq "yes") {
                $gapend[$#gapend+1]=$i-1;
                $record="no";
            }
        } else {
            if ($record eq "no") {
                $gapstart[$#gapstart+1]=$i;
                $record="yes";
            }
        }
    }
    if ($record eq "yes") { ### should not happen
        $gapend[$#gapend+1]=$result;
    }
    if ($#gapstart!=$#gapend) {
        print ("Error: Unknown No.3\n");
    }
    for ($i=0;$i<=$#gapstart;$i++) {
        if (($gapend[$i]-$gapstart[$i]+1)>$upgap) {
            $upgap=$gapend[$i]-$gapstart[$i]+1;
        }
    }

    return($result);
}

sub down_extension { ### output nt length can be extended to downstream direction
    ### continue to use INPUT_DATA1 as .fa, and information from @FaiInfo, %ChrIndex
    ### [0]=chr [1]=start position to extend (right before extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion
    my @faiext=split /\s+/, $FaiInfo[$ChrIndex{$_[0]}];
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=1;$y<=length($_[2]);$y++) { ### store repeat unit
        $repeat[$y]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence
    seek (INPUT_DATA1,$faiext[2]+int($_[1]/$faiext[3])*$faiext[4]+$_[1]%$faiext[3],0);
    my $posnow=$_[1]+1;
    my $refchar="";

    my $gapscore=$_[5];
    my $misscore=$_[4]; ### N misscore will = $misscore*3/4
    my $matchscore=$_[3];
    my $limit=$_[6];
    my $count=0;
    my $result=0; ### max successful extension nt number
    ### score matrix ###
    my @bottom=();
    $bottom[0]=0;
    my @tempold=();
    my @tempnew=();
    my @rightold=();
    $rightold[0]=0;
    my @rightnew=();
    my @gapbool=();
    ####################

    my $keepext="no";
    my $value=0;
    while ($posnow<=$faiext[1]) {
        $refchar=getc(INPUT_DATA1);
        if ($refchar ne "\n") {
            $refchar=uc($refchar);
            if (($refchar ne "A")and($refchar ne "C")and($refchar ne "T")and($refchar ne "G")) {
                $refchar="N";
            }
            $count++; ### $count=X (=$posnow-$_[1]) means now trying to align the Xst nt
            $ref[$count]=$refchar;

            $tempnew[0]=$bottom[0];
            for ($y=1;$y<=length($_[2]);$y++) {
                if ($tempnew[$y-1]==$limit) { ### any calculation using a stop match is considered as a stop
                    $tempnew[$y]=$limit;
                } else {
                    $tempnew[$y]=$tempnew[$y-1]+$gapscore;
                }
                if ($tempnew[$y]<=$limit) { ### any score less than $limit is considered as a stop
                    $tempnew[$y]=$limit;
                } elsif ($tempnew[$y]>=0) { ### any score >0 will be reset to 0, so previous extension doesn't affect downstream
                    $tempnew[$y]=0;
                }
            }
            @tempold=@tempnew;
            @tempnew=();
            $bottom[0]=$tempold[length($_[2])];

            for ($x=1;$x<$count;$x++) {
                $tempnew[0]=$bottom[$x];
                for ($y=1;$y<=length($_[2]);$y++) { ### get max score from 3 routes
                    if ($tempnew[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempnew[$y-1]+$gapscore;
                    }
                    $tempnew[$y]=$value;

                    if ($tempold[$y]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempold[$y]+$gapscore;
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempold[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        if (($ref[$x] eq "N")or($repeat[$y] eq "N")) {
                            $value=$tempold[$y-1]+$misscore*0.75;
                        } elsif ($ref[$x] ne $repeat[$y]) {
                            $value=$tempold[$y-1]+$misscore;
                        } else {
                            $value=$tempold[$y-1]+$matchscore;
                        }
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempnew[$y]<=$limit) {
                        $tempnew[$y]=$limit;
                    } elsif ($tempnew[$y]>=0) {
                        $tempnew[$y]=0;
                        $gapbool[$x]=1;
                        if ($x>$result) {
                            $result=$x; ### if any score =0, consider it as successful extension
                        }
                    }
                }
                @tempold=@tempnew;
                @tempnew=();
                $bottom[$x]=$tempold[length($_[2])];
            }
            for ($y=1;$y<=length($_[2]);$y++) {
                $rightold[$#rightold+1]=$tempold[$y];
            }
            @tempold=();
            if ($#rightold!=($count*length($_[2]))) {
                print ("Error: Unknown No.2\n");
            }

            $keepext="no";
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext="yes"; ### if any score >limit, continue to try next nt
                if ($rightnew[0]>=0) {
                    $rightnew[0]=0;
                    $gapbool[$count]=1;
                    $result=$count; ### if any score =0, consider it as successful extension
                }
            }

            for ($y=1;$y<=$count*length($_[2]);$y++) { ### get max score from 3 routes
                if ($rightnew[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightnew[$y-1]+$gapscore;
                }
                $rightnew[$y]=$value;

                if ($rightold[$y]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightold[$y]+$gapscore;
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightold[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    if (($ref[$count] eq "N")or($repeat[$y%(length($_[2]))] eq "N")) {
                        $value=$rightold[$y-1]+$misscore*0.75;
                    } elsif ($ref[$count] ne $repeat[$y%(length($_[2]))]) {
                        $value=$rightold[$y-1]+$misscore;
                    } else {
                        $value=$rightold[$y-1]+$matchscore;
                    }
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightnew[$y]<=$limit) {
                    $rightnew[$y]=$limit;
                } else {
                    $keepext="yes"; ### if any score >limit, continue to try next nt
                    if ($rightnew[$y]>=0) {
                        $rightnew[$y]=0;
                        $gapbool[$count]=1;
                        $result=$count; ### if any score =0, consider it as successful extension
                    }
                }
            }
            @rightold=@rightnew;
            @rightnew=();
            $bottom[$count]=$rightold[$#rightold];

            if ($keepext eq "no") {
                last;
            } else {
                $posnow++;
            }
        }
    }

    my @gapstart=();
    my @gapend=();
    my $i=0;
    my $record="no";
    $downgap=0;
    for ($i=1;$i<=$result;$i++) {
        if (defined($gapbool[$i])) {
            if ($record eq "yes") {
                $gapend[$#gapend+1]=$i-1;
                $record="no";
            }
        } else {
            if ($record eq "no") {
                $gapstart[$#gapstart+1]=$i;
                $record="yes";
            }
        }
    }
    if ($record eq "yes") { ### should not happen
        $gapend[$#gapend+1]=$result;
    }
    if ($#gapstart!=$#gapend) {
        print ("Error: Unknown No.3\n");
    }
    for ($i=0;$i<=$#gapstart;$i++) {
        if (($gapend[$i]-$gapstart[$i]+1)>$downgap) {
            $downgap=$gapend[$i]-$gapstart[$i]+1;
        }
    }

    return($result);
}

sub up_extensionmore { ### output nt length can be extended to upstream direction
    ### continue to use INPUT_DATA1 as .fa, and information from @FaiInfo, %ChrIndex
    ### [0]=chr [1]=start position to extend (right after extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion
    my @faiext=split /\s+/, $FaiInfo[$ChrIndex{$_[0]}];
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=length($_[2]);$y>=1;$y--) { ### store repeat unit in reversed order
        $repeat[length($_[2])-$y+1]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence in reversed order
    seek (INPUT_DATA1,$faiext[2]+int($_[1]/$faiext[3])*$faiext[4]+$_[1]%$faiext[3],0);
    if (($_[1]%$faiext[3])==0) {
        seek (INPUT_DATA1,-1,1);
    }
    my $posnow=$_[1]-1;
    my $refchar="";

    my $gapscore=$_[5];
    my $misscore=$_[4]; ### N misscore will = $misscore*3/4
    my $matchscore=$_[3];
    my $allowance=$_[7];
    my $limit=$_[6]+$allowance;
    my $count=0;
    my $result=0; ### max successful extension nt number
    ### score matrix ###
    my @bottom=();
    $bottom[0]=0;
    my @tempold=();
    my @tempnew=();
    my @rightold=();
    $rightold[0]=0;
    my @rightnew=();
    ####################

    my $keepext="no";
    my $value=0;
    while ($posnow>=1) {
        seek (INPUT_DATA1,-2,1);
        $refchar=getc(INPUT_DATA1);
        if ($refchar ne "\n") {
            $refchar=uc($refchar);
            if (($refchar ne "A")and($refchar ne "C")and($refchar ne "T")and($refchar ne "G")) {
                $refchar="N";
            }
            $count++; ### $count=X (=$_[1]-$posnow) means now trying to align the Xst nt
            $ref[$count]=$refchar;

            $tempnew[0]=$bottom[0];
            for ($y=1;$y<=length($_[2]);$y++) {
                if ($tempnew[$y-1]==$limit) { ### any calculation using a stop match is considered as a stop
                    $tempnew[$y]=$limit;
                } else {
                    $tempnew[$y]=$tempnew[$y-1]+$gapscore;
                }
                if ($tempnew[$y]<=$limit) { ### any score less than $limit is considered as a stop
                    $tempnew[$y]=$limit;
                } elsif ($tempnew[$y]>=0) { ### any score >0 will be reset to 0, so previous extension doesn't affect downstream
                    $tempnew[$y]=0;
                }
            }
            @tempold=@tempnew;
            @tempnew=();
            $bottom[0]=$tempold[length($_[2])];

            for ($x=1;$x<$count;$x++) {
                $tempnew[0]=$bottom[$x];
                for ($y=1;$y<=length($_[2]);$y++) { ### get max score from 3 routes
                    if ($tempnew[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempnew[$y-1]+$gapscore;
                    }
                    $tempnew[$y]=$value;

                    if ($tempold[$y]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempold[$y]+$gapscore;
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempold[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        if (($ref[$x] eq "N")or($repeat[$y] eq "N")) {
                            $value=$tempold[$y-1]+$misscore*0.75;
                        } elsif ($ref[$x] ne $repeat[$y]) {
                            $value=$tempold[$y-1]+$misscore;
                        } else {
                            $value=$tempold[$y-1]+$matchscore;
                        }
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempnew[$y]<=$limit) {
                        $tempnew[$y]=$limit;
                    } elsif ($tempnew[$y]>=0) {
                        $tempnew[$y]=0;
                    }
                    if (($tempnew[$y]>=$allowance)and($x>$result)) { ### if any score >=$allowance, consider it as successful extension
                        $result=$x;
                    }
                }
                @tempold=@tempnew;
                @tempnew=();
                $bottom[$x]=$tempold[length($_[2])];
            }
            for ($y=1;$y<=length($_[2]);$y++) {
                $rightold[$#rightold+1]=$tempold[$y];
            }
            @tempold=();
            if ($#rightold!=($count*length($_[2]))) {
                print ("Error: Unknown No.2\n");
            }

            $keepext="no";
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext="yes"; ### if any score >limit, continue to try next nt
                if ($rightnew[0]>=0) {
                    $rightnew[0]=0;
                }
                if ($rightnew[0]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                    $result=$count;
                }
            }

            for ($y=1;$y<=$count*length($_[2]);$y++) { ### get max score from 3 routes
                if ($rightnew[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightnew[$y-1]+$gapscore;
                }
                $rightnew[$y]=$value;

                if ($rightold[$y]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightold[$y]+$gapscore;
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightold[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    if (($ref[$count] eq "N")or($repeat[$y%(length($_[2]))] eq "N")) {
                        $value=$rightold[$y-1]+$misscore*0.75;
                    } elsif ($ref[$count] ne $repeat[$y%(length($_[2]))]) {
                        $value=$rightold[$y-1]+$misscore;
                    } else {
                        $value=$rightold[$y-1]+$matchscore;
                    }
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightnew[$y]<=$limit) {
                    $rightnew[$y]=$limit;
                } else {
                    $keepext="yes"; ### if any score >limit, continue to try next nt
                    if ($rightnew[$y]>=0) {
                        $rightnew[$y]=0;
                    }
                    if ($rightnew[$y]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                        $result=$count;
                    }
                }
            }
            @rightold=@rightnew;
            @rightnew=();
            $bottom[$count]=$rightold[$#rightold];

            if ($keepext eq "no") {
                last;
            } else {
                $posnow--;
            }
        }
    }
    return($result);
}

sub down_extensionmore { ### output nt length can be extended to downstream direction
    ### continue to use INPUT_DATA1 as .fa, and information from @FaiInfo, %ChrIndex
    ### [0]=chr [1]=start position to extend (right before extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion
    my @faiext=split /\s+/, $FaiInfo[$ChrIndex{$_[0]}];
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=1;$y<=length($_[2]);$y++) { ### store repeat unit
        $repeat[$y]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence
    seek (INPUT_DATA1,$faiext[2]+int($_[1]/$faiext[3])*$faiext[4]+$_[1]%$faiext[3],0);
    my $posnow=$_[1]+1;
    my $refchar="";

    my $gapscore=$_[5];
    my $misscore=$_[4]; ### N misscore will = $misscore*3/4
    my $matchscore=$_[3];
    my $allowance=$_[7];
    my $limit=$_[6]+$allowance;
    my $count=0;
    my $result=0; ### max successful extension nt number
    ### score matrix ###
    my @bottom=();
    $bottom[0]=0;
    my @tempold=();
    my @tempnew=();
    my @rightold=();
    $rightold[0]=0;
    my @rightnew=();
    ####################

    my $keepext="no";
    my $value=0;
    while ($posnow<=$faiext[1]) {
        $refchar=getc(INPUT_DATA1);
        if ($refchar ne "\n") {
            $refchar=uc($refchar);
            if (($refchar ne "A")and($refchar ne "C")and($refchar ne "T")and($refchar ne "G")) {
                $refchar="N";
            }
            $count++; ### $count=X (=$posnow-$_[1]) means now trying to align the Xst nt
            $ref[$count]=$refchar;

            $tempnew[0]=$bottom[0];
            for ($y=1;$y<=length($_[2]);$y++) {
                if ($tempnew[$y-1]==$limit) { ### any calculation using a stop match is considered as a stop
                    $tempnew[$y]=$limit;
                } else {
                    $tempnew[$y]=$tempnew[$y-1]+$gapscore;
                }
                if ($tempnew[$y]<=$limit) { ### any score less than $limit is considered as a stop
                    $tempnew[$y]=$limit;
                } elsif ($tempnew[$y]>=0) { ### any score >0 will be reset to 0, so previous extension doesn't affect downstream
                    $tempnew[$y]=0;
                }
            }
            @tempold=@tempnew;
            @tempnew=();
            $bottom[0]=$tempold[length($_[2])];

            for ($x=1;$x<$count;$x++) {
                $tempnew[0]=$bottom[$x];
                for ($y=1;$y<=length($_[2]);$y++) { ### get max score from 3 routes
                    if ($tempnew[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempnew[$y-1]+$gapscore;
                    }
                    $tempnew[$y]=$value;

                    if ($tempold[$y]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempold[$y]+$gapscore;
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempold[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        if (($ref[$x] eq "N")or($repeat[$y] eq "N")) {
                            $value=$tempold[$y-1]+$misscore*0.75;
                        } elsif ($ref[$x] ne $repeat[$y]) {
                            $value=$tempold[$y-1]+$misscore;
                        } else {
                            $value=$tempold[$y-1]+$matchscore;
                        }
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempnew[$y]<=$limit) {
                        $tempnew[$y]=$limit;
                    } elsif ($tempnew[$y]>=0) {
                        $tempnew[$y]=0;
                    }
                    if (($tempnew[$y]>=$allowance)and($x>$result)) { ### if any score >=$allowance, consider it as successful extension
                        $result=$x;
                    }
                }
                @tempold=@tempnew;
                @tempnew=();
                $bottom[$x]=$tempold[length($_[2])];
            }
            for ($y=1;$y<=length($_[2]);$y++) {
                $rightold[$#rightold+1]=$tempold[$y];
            }
            @tempold=();
            if ($#rightold!=($count*length($_[2]))) {
                print ("Error: Unknown No.2\n");
            }

            $keepext="no";
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext="yes"; ### if any score >limit, continue to try next nt
                if ($rightnew[0]>=0) {
                    $rightnew[0]=0;
                }
                if ($rightnew[0]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                    $result=$count;
                }
            }

            for ($y=1;$y<=$count*length($_[2]);$y++) { ### get max score from 3 routes
                if ($rightnew[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightnew[$y-1]+$gapscore;
                }
                $rightnew[$y]=$value;

                if ($rightold[$y]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightold[$y]+$gapscore;
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightold[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    if (($ref[$count] eq "N")or($repeat[$y%(length($_[2]))] eq "N")) {
                        $value=$rightold[$y-1]+$misscore*0.75;
                    } elsif ($ref[$count] ne $repeat[$y%(length($_[2]))]) {
                        $value=$rightold[$y-1]+$misscore;
                    } else {
                        $value=$rightold[$y-1]+$matchscore;
                    }
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightnew[$y]<=$limit) {
                    $rightnew[$y]=$limit;
                } else {
                    $keepext="yes"; ### if any score >limit, continue to try next nt
                    if ($rightnew[$y]>=0) {
                        $rightnew[$y]=0;
                    }
                    if ($rightnew[$y]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                        $result=$count;
                    }
                }
            }
            @rightold=@rightnew;
            @rightnew=();
            $bottom[$count]=$rightold[$#rightold];

            if ($keepext eq "no") {
                last;
            } else {
                $posnow++;
            }
        }
    }
    return($result);
}

sub up_extensionmore_wl { ### output nt length can be extended to upstream direction
    ### continue to use INPUT_DATA1 as .fa, and information from @FaiInfo, %ChrIndex
    ### [0]=chr [1]=start position to extend (right after extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion
    my @faiext=split /\s+/, $FaiInfo[$ChrIndex{$_[0]}];
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=length($_[2]);$y>=1;$y--) { ### store repeat unit in reversed order
        $repeat[length($_[2])-$y+1]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence in reversed order
    seek (INPUT_DATA1,$faiext[2]+int($_[1]/$faiext[3])*$faiext[4]+$_[1]%$faiext[3],0);
    if (($_[1]%$faiext[3])==0) {
        seek (INPUT_DATA1,-1,1);
    }
    my $posnow=$_[1]-1;
    my $refchar="";

    my $gapscore=$_[5];
    my $misscore=$_[4]; ### N misscore will = $misscore*3/4
    my $matchscore=$_[3];
    my $allowance=$_[7];
    my $limit=$_[6]+$allowance;
    my $count=0;
    my $result=0; ### max successful extension nt number
    my $maxresult=$_[8];
    ### score matrix ###
    my @bottom=();
    $bottom[0]=0;
    my @tempold=();
    my @tempnew=();
    my @rightold=();
    $rightold[0]=0;
    my @rightnew=();
    ####################

    my $keepext="no";
    my $value=0;
    while ($posnow>=1) {
        seek (INPUT_DATA1,-2,1);
        $refchar=getc(INPUT_DATA1);
        if ($refchar ne "\n") {
            $refchar=uc($refchar);
            if (($refchar ne "A")and($refchar ne "C")and($refchar ne "T")and($refchar ne "G")) {
                $refchar="N";
            }
            $count++; ### $count=X (=$_[1]-$posnow) means now trying to align the Xst nt
            $ref[$count]=$refchar;

            $tempnew[0]=$bottom[0];
            for ($y=1;$y<=length($_[2]);$y++) {
                if ($tempnew[$y-1]==$limit) { ### any calculation using a stop match is considered as a stop
                    $tempnew[$y]=$limit;
                } else {
                    $tempnew[$y]=$tempnew[$y-1]+$gapscore;
                }
                if ($tempnew[$y]<=$limit) { ### any score less than $limit is considered as a stop
                    $tempnew[$y]=$limit;
                } elsif ($tempnew[$y]>=0) { ### any score >0 will be reset to 0, so previous extension doesn't affect downstream
                    $tempnew[$y]=0;
                }
            }
            @tempold=@tempnew;
            @tempnew=();
            $bottom[0]=$tempold[length($_[2])];

            for ($x=1;$x<$count;$x++) {
                $tempnew[0]=$bottom[$x];
                for ($y=1;$y<=length($_[2]);$y++) { ### get max score from 3 routes
                    if ($tempnew[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempnew[$y-1]+$gapscore;
                    }
                    $tempnew[$y]=$value;

                    if ($tempold[$y]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempold[$y]+$gapscore;
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempold[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        if (($ref[$x] eq "N")or($repeat[$y] eq "N")) {
                            $value=$tempold[$y-1]+$misscore*0.75;
                        } elsif ($ref[$x] ne $repeat[$y]) {
                            $value=$tempold[$y-1]+$misscore;
                        } else {
                            $value=$tempold[$y-1]+$matchscore;
                        }
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempnew[$y]<=$limit) {
                        $tempnew[$y]=$limit;
                    } elsif ($tempnew[$y]>=0) {
                        $tempnew[$y]=0;
                    }
                    if (($tempnew[$y]>=$allowance)and($x>$result)) { ### if any score >=$allowance, consider it as successful extension
                        $result=$x;
                        if ($result>=$maxresult) {
                            return($maxresult);
                        }
                    }
                }
                @tempold=@tempnew;
                @tempnew=();
                $bottom[$x]=$tempold[length($_[2])];
            }
            for ($y=1;$y<=length($_[2]);$y++) {
                $rightold[$#rightold+1]=$tempold[$y];
            }
            @tempold=();
            if ($#rightold!=($count*length($_[2]))) {
                print ("Error: Unknown No.2\n");
            }

            $keepext="no";
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext="yes"; ### if any score >limit, continue to try next nt
                if ($rightnew[0]>=0) {
                    $rightnew[0]=0;
                }
                if ($rightnew[0]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                    $result=$count;
                    if ($result>=$maxresult) {
                        return($maxresult);
                    }
                }
            }

            for ($y=1;$y<=$count*length($_[2]);$y++) { ### get max score from 3 routes
                if ($rightnew[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightnew[$y-1]+$gapscore;
                }
                $rightnew[$y]=$value;

                if ($rightold[$y]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightold[$y]+$gapscore;
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightold[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    if (($ref[$count] eq "N")or($repeat[$y%(length($_[2]))] eq "N")) {
                        $value=$rightold[$y-1]+$misscore*0.75;
                    } elsif ($ref[$count] ne $repeat[$y%(length($_[2]))]) {
                        $value=$rightold[$y-1]+$misscore;
                    } else {
                        $value=$rightold[$y-1]+$matchscore;
                    }
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightnew[$y]<=$limit) {
                    $rightnew[$y]=$limit;
                } else {
                    $keepext="yes"; ### if any score >limit, continue to try next nt
                    if ($rightnew[$y]>=0) {
                        $rightnew[$y]=0;
                    }
                    if ($rightnew[$y]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                        $result=$count;
                        if ($result>=$maxresult) {
                            return($maxresult);
                        }
                    }
                }
            }
            @rightold=@rightnew;
            @rightnew=();
            $bottom[$count]=$rightold[$#rightold];

            if ($keepext eq "no") {
                last;
            } else {
                $posnow--;
            }
        }
    }
    return($result);
}

sub down_extensionmore_wl { ### output nt length can be extended to downstream direction
    ### continue to use INPUT_DATA1 as .fa, and information from @FaiInfo, %ChrIndex
    ### [0]=chr [1]=start position to extend (right before extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion
    my @faiext=split /\s+/, $FaiInfo[$ChrIndex{$_[0]}];
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=1;$y<=length($_[2]);$y++) { ### store repeat unit
        $repeat[$y]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence
    seek (INPUT_DATA1,$faiext[2]+int($_[1]/$faiext[3])*$faiext[4]+$_[1]%$faiext[3],0);
    my $posnow=$_[1]+1;
    my $refchar="";

    my $gapscore=$_[5];
    my $misscore=$_[4]; ### N misscore will = $misscore*3/4
    my $matchscore=$_[3];
    my $allowance=$_[7];
    my $limit=$_[6]+$allowance;
    my $count=0;
    my $result=0; ### max successful extension nt number
    my $maxresult=$_[8];
    ### score matrix ###
    my @bottom=();
    $bottom[0]=0;
    my @tempold=();
    my @tempnew=();
    my @rightold=();
    $rightold[0]=0;
    my @rightnew=();
    ####################

    my $keepext="no";
    my $value=0;
    while ($posnow<=$faiext[1]) {
        $refchar=getc(INPUT_DATA1);
        if ($refchar ne "\n") {
            $refchar=uc($refchar);
            if (($refchar ne "A")and($refchar ne "C")and($refchar ne "T")and($refchar ne "G")) {
                $refchar="N";
            }
            $count++; ### $count=X (=$posnow-$_[1]) means now trying to align the Xst nt
            $ref[$count]=$refchar;

            $tempnew[0]=$bottom[0];
            for ($y=1;$y<=length($_[2]);$y++) {
                if ($tempnew[$y-1]==$limit) { ### any calculation using a stop match is considered as a stop
                    $tempnew[$y]=$limit;
                } else {
                    $tempnew[$y]=$tempnew[$y-1]+$gapscore;
                }
                if ($tempnew[$y]<=$limit) { ### any score less than $limit is considered as a stop
                    $tempnew[$y]=$limit;
                } elsif ($tempnew[$y]>=0) { ### any score >0 will be reset to 0, so previous extension doesn't affect downstream
                    $tempnew[$y]=0;
                }
            }
            @tempold=@tempnew;
            @tempnew=();
            $bottom[0]=$tempold[length($_[2])];

            for ($x=1;$x<$count;$x++) {
                $tempnew[0]=$bottom[$x];
                for ($y=1;$y<=length($_[2]);$y++) { ### get max score from 3 routes
                    if ($tempnew[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempnew[$y-1]+$gapscore;
                    }
                    $tempnew[$y]=$value;

                    if ($tempold[$y]==$limit) {
                        $value=$limit;
                    } else {
                        $value=$tempold[$y]+$gapscore;
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempold[$y-1]==$limit) {
                        $value=$limit;
                    } else {
                        if (($ref[$x] eq "N")or($repeat[$y] eq "N")) {
                            $value=$tempold[$y-1]+$misscore*0.75;
                        } elsif ($ref[$x] ne $repeat[$y]) {
                            $value=$tempold[$y-1]+$misscore;
                        } else {
                            $value=$tempold[$y-1]+$matchscore;
                        }
                    }
                    if ($value>$tempnew[$y]) {
                        $tempnew[$y]=$value;
                    }

                    if ($tempnew[$y]<=$limit) {
                        $tempnew[$y]=$limit;
                    } elsif ($tempnew[$y]>=0) {
                        $tempnew[$y]=0;
                    }
                    if (($tempnew[$y]>=$allowance)and($x>$result)) { ### if any score >=$allowance, consider it as successful extension
                        $result=$x;
                        if ($result>=$maxresult) {
                            return($maxresult);
                        }
                    }
                }
                @tempold=@tempnew;
                @tempnew=();
                $bottom[$x]=$tempold[length($_[2])];
            }
            for ($y=1;$y<=length($_[2]);$y++) {
                $rightold[$#rightold+1]=$tempold[$y];
            }
            @tempold=();
            if ($#rightold!=($count*length($_[2]))) {
                print ("Error: Unknown No.2\n");
            }

            $keepext="no";
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext="yes"; ### if any score >limit, continue to try next nt
                if ($rightnew[0]>=0) {
                    $rightnew[0]=0;
                }
                if ($rightnew[0]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                    $result=$count;
                    if ($result>=$maxresult) {
                        return($maxresult);
                    }
                }
            }

            for ($y=1;$y<=$count*length($_[2]);$y++) { ### get max score from 3 routes
                if ($rightnew[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightnew[$y-1]+$gapscore;
                }
                $rightnew[$y]=$value;

                if ($rightold[$y]==$limit) {
                    $value=$limit;
                } else {
                    $value=$rightold[$y]+$gapscore;
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightold[$y-1]==$limit) {
                    $value=$limit;
                } else {
                    if (($ref[$count] eq "N")or($repeat[$y%(length($_[2]))] eq "N")) {
                        $value=$rightold[$y-1]+$misscore*0.75;
                    } elsif ($ref[$count] ne $repeat[$y%(length($_[2]))]) {
                        $value=$rightold[$y-1]+$misscore;
                    } else {
                        $value=$rightold[$y-1]+$matchscore;
                    }
                }
                if ($value>$rightnew[$y]) {
                    $rightnew[$y]=$value;
                }

                if ($rightnew[$y]<=$limit) {
                    $rightnew[$y]=$limit;
                } else {
                    $keepext="yes"; ### if any score >limit, continue to try next nt
                    if ($rightnew[$y]>=0) {
                        $rightnew[$y]=0;
                    }
                    if ($rightnew[$y]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                        $result=$count;
                        if ($result>=$maxresult) {
                            return($maxresult);
                        }
                    }
                }
            }
            @rightold=@rightnew;
            @rightnew=();
            $bottom[$count]=$rightold[$#rightold];

            if ($keepext eq "no") {
                last;
            } else {
                $posnow++;
            }
        }
    }
    return($result);
}

### read parameters ###
my $i=0;
@Para=();
### Para[0]=.fa [1]=.fai by samtools faidx
### Para[2]=file containing lines: "STRid chr pos1 pos2 STRunit"
### use range pos1~pos2 to find STR seed; chr can be "X" or "chrX", but need to match .fa and .fai
### Para[3]=output file
### Para[4..7] scores of match mismatch gap stop, suggestion 2 -5 -7 -30
### Para[8] use dynamic match score "yes"/"no", suggestion "no"
### Para[9] check possibility of complementary STR unit "yes"/"no", suggestion "no"
### Para[10] 5FL/3FL max output, optional to speed up, undefined=no limit, suggestion 100
### Para[11] output length of flanking sequences, suggestion 1000
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-h')or($ARGV[$i] eq '--help')) {
        print "Instruction:\n";
        print "\tThe Finder module processes the basic information of the STRs provided by user, and retrieves standarized STR sequences by user defined parameters.\n";
        print "\tThe parameters and other information retrieved will be inherited in subsequent steps.\n";
        print "\tSee README.txt Step (1.1) and (1.2) for more details.\n";
        print "Usage Sample:\n";
        print "\tperl LUSTR_Finder.pl -r <reference.fa> --refindex <reference.fa.fai> -i <STRinfo.txt> -o <STRsequences.txt> -e 100\n";
        print "Options:\n";
        print "\t--help/-h\t\tPrint usage instructions\n";
        print "\t--ref/-r <file>\t\tThe reference file to search for STR sequences\n";
        print "\t--refindex <file>\tIndex .fai file of the reference that can be generated by samtools\n";
        print "\t\t\t\tIf not appointed, LUSTR will search for file with \".fai\" extension after the file name of the reference\n";
        print "\t--input/-i <file>\tThe file for basic information of STRs. See format details in README.txt Step (1.1)\n";
        print "\t--output/-o <file>\tOutput file for STR sequences. See format details in README.txt Step (1.2)\n";
        print "\t--match/-m <value>\tMatch score for repeat extension (Default = 2)\n";
        print "\t--mis/-x <value>\tMismatch penalty for repeat extension (Default = -5)\n";
        print "\t--gap/-g <value>\tGap penalty for repeat extension (Default = -7)\n";
        print "\t--stop/-s <value>\tThreshold to stop repeat extension (Default = -30)\n";
        print "\t--dynamic/-d\t\tSwitch on dynamic adjustment of match score for long repetitive units (Default = off)\n";
        print "\t--comple/-c\t\tSwitch on searching for the STR sequence on reference negative strand when positive strand searching fails (Default = off)\n";
        print "\t--extra/-e <value>\tThis setting helps avoid long time running for certain STR cases when searching for extra repeat extensions allowing mismatches\n";
        print "\t\t\t\t(Default = no limit, but a setting = 100 is recommended)\n";
        print "\t--flank/-f <value>\tLength to export flanking sequences\n";
        print "\t\t\t\tThis length determines whether the reads in pair are close enough to be considered as from the same STR locus\n";
        print "\t\t\t\tPlease make the best setting accordingly to the sequencing library preparation procedure\n";
        print "\t\t\t\t(Default = 1000)\n";
        exit;
    }
}
$Para[0]='';
$Para[1]='';
$Para[2]='';
$Para[3]='';
$Para[4]=2;
$Para[5]=-5;
$Para[6]=-7;
$Para[7]=-30;
$Para[8]='no';
$Para[9]='no';
$Para[10]="unlimited";
$Para[11]=1000;
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-r')or($ARGV[$i] eq '--ref')) {
        $i++;
        $Para[0]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-i')or($ARGV[$i] eq '--input')) {
        $i++;
        $Para[2]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-o')or($ARGV[$i] eq '--output')) {
        $i++;
        $Para[3]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-m')or($ARGV[$i] eq '--match')) {
        $i++;
        $Para[4]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-x')or($ARGV[$i] eq '--mis')) {
        $i++;
        $Para[5]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-g')or($ARGV[$i] eq '--gap')) {
        $i++;
        $Para[6]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-s')or($ARGV[$i] eq '--stop')) {
        $i++;
        $Para[7]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-d')or($ARGV[$i] eq '--dynamic')) {
        $Para[8]='yes';
    } elsif (($ARGV[$i] eq '-c')or($ARGV[$i] eq '--comple')) {
        $Para[9]='yes';
    } elsif (($ARGV[$i] eq '-e')or($ARGV[$i] eq '--extra')) {
        $i++;
        $Para[10]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-f')or($ARGV[$i] eq '--flank')) {
        $i++;
        $Para[11]=$ARGV[$i];
    }
}
$Para[1]=join '',($Para[0],'.fai');
for ($i=0;$i<=$#ARGV;$i++) {
    if ($ARGV[$i] eq '--refindex') {
        $i++;
        $Para[1]=$ARGV[$i];
    }
}
if (($Para[0] eq '')or($Para[1] eq '')or($Para[2] eq '')or($Para[3] eq '')) {
    print ("Input/Output not appointed. Stop running.\n");
    print "Please use --help or -h to see instructions.\n";
    exit;
}

### main ###
my $reffile=$Para[0];
while ($reffile=~/\//) {
    $reffile=$';
}
&read_fai($Para[1]);
$match_score=$Para[4];
$mis_score=$Para[5];
$gap_score=$Para[6];
$limit_score=$Para[7]; ### ~[7]/[5] mismatches will lead to stop
$allow_score=2*$mis_score+$gap_score; ### will further -3*$dymatch_score when use
$dynamic=$Para[8]; ### if "yes", match_score will be dynamically adjusted for long unit STRs, to make the score of a unit with 1mis <0
$check_comple=$Para[9];
$FL_limit=$Para[10];
$flanklength=$Para[11];
$upgap=0; ### max gap length found by &up_extension
$downgap=0; ### max gap length found by &down_extension

my $temp=join "",($Para[3],"_notrun");
open (OUTPUT_NOTRUN, ">$temp") or die "Couldn't open: $!";
$temp=join "",($Para[3],"_noSTR");
open (OUTPUT_NOSTR, ">$temp") or die "Couldn't open: $!";
$temp=join "",($Para[3],"_warning");
open (OUTPUT_WARN, ">$temp") or die "Couldn't open: $!";
open (OUTPUT_RESULT, ">$Para[3]") or die "Couldn't open: $!";
print OUTPUT_RESULT ("Parameters (match/mis/gap/stop/dynamic/maxFL): $match_score $mis_score $gap_score $limit_score $dynamic $FL_limit Reference: $reffile\n");

&get_seed($Para[0],$Para[2]);
close OUTPUT_NOTRUN;
close OUTPUT_NOSTR;
close OUTPUT_WARN;
close OUTPUT_RESULT;

$temp=join "",($Para[3],"_log");
open (OUTPUT_LOG, ">$temp") or die "Couldn't open: $!";
print OUTPUT_LOG ("$totalinput inputs read from $Para[2]\n");
print OUTPUT_LOG ("$goodinput inputs are informative and searched for STR sequences (Check discarded ID list in _notrun)\n");
print OUTPUT_LOG ("\t");
print OUTPUT_LOG ("$oriinput STRs are found by original input repeat\n");
print OUTPUT_LOG ("\t");
if ($check_comple eq "yes") {
    print OUTPUT_LOG ("$compleinput STRs are found by complementary input repeat\n");
    print OUTPUT_LOG ("\t");
    print OUTPUT_LOG ("$bothinput STRs can't determine the best sequence. Using original input repeat\n");
} else {
    print OUTPUT_LOG ("Skipped searching by complementary input repeat\n");
}
print OUTPUT_LOG ("\t");
print OUTPUT_LOG ("$missinput STRs were not found (Check ID list in _noSTR)\n\n");
print OUTPUT_LOG ("Parameters used for STR extension (will be automatically transfered to alignment step):\n\t");
print OUTPUT_LOG ("match score = $match_score\n\t");
print OUTPUT_LOG ("mismatch score = $mis_score\n\t");
print OUTPUT_LOG ("gap score = $gap_score\n\t");
print OUTPUT_LOG ("stop score = $limit_score\n\t");
print OUTPUT_LOG ("dynamic match score : $dynamic\n\t");
print OUTPUT_LOG ("max FL : $FL_limit\n");
close OUTPUT_LOG;
print ("All searching done. Please check _log for IMPORTANT information\n");
############

exit;

