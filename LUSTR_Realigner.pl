#!/usr/local/bin/perl -w
### realigns remapped pairs to each STR, generates a plain text file indicating how each pair is realigned # Bookmarks: 0,0 0,0 0,0 0,5930

sub complementary { ### return complementary sequence of input
    my $temp=reverse($_[0]);
    $temp=~tr/ATCGatcg/TAGCtagc/;
    return $temp;
}

sub PrepareFlag { ### explain sam format FLAG (column 2), store 4 FLAG related arrays FlagPriority(p/s),FlagPair(1/2/0),FlagStrand(+/-),FlagMap(y/n))
    @FlagPriority=(); ### "p"=primary, "s"=secondary or supplementary
    @FlagStrand=(); ### "+"/"-" no matter mapped or unmapped
    @FlagMap=(); ### "y"=mapped, "n"=unmapped
    @FlagPair=(); ### 0=unpaired/unexplainable
    my $i=0;
    my $temp='';
    for ($i=0;$i<=4095;$i++) {
        if ((($i%512)>=256)or(($i%4096)>=2048)) {
            $FlagPriority[$i]='s';
        } else {
            $FlagPriority[$i]='p';
        }
        if (($i%32)>=16) {
            $FlagStrand[$i]='-';
        } else {
            $FlagStrand[$i]='+';
        }
        if (($i%8)>=4) {
            $FlagMap[$i]='n';
        } else {
            $FlagMap[$i]='y';
        }
        $FlagPair[$i]=0;
        if (($i%2)>=1) {
            if (($i%128)>=64) {
                $temp='y';
            } else {
                $temp='n';
            }
            if (($i%256)>=128) {
                $temp.='y';
            } else {
                $temp.='n';
            }
            if ($temp eq 'yn') {
                $FlagPair[$i]=1;
            } elsif ($temp eq 'ny') {
                $FlagPair[$i]=2;
            }
        }
    }
}

sub standardizeDNA { ### make sequence ATCGN only
    my $temp=uc($_[0]);
    $temp=~s/[^ATCGN]/N/gi;
    return $temp;
}

sub TransformSTR { ### transform input STR unit $_[0] to find equivalent index in %SameUnit, return the index (or "" if not found)
    my $i=0;
    my $result='';
    my $temp='';
    for ($i=length($_[0])-1;$i>=0;$i--) {
        $temp=join '',(substr($_[0],$i,length($_[0])-$i),substr($_[0],0,$i));
        if (defined($SameUnit{$temp})) {
            $result=$temp;
            last;
        }
    }
    return ($result);
}

sub PrepareDymscore { ### return adjusted match score $dymscore[unitlength]
    @dymscore=();
    my $i=0;
    my $max=$mis_score;
    if ($max<$gap_score) {
        $max=$gap_score;
    }
    for ($i=$min_ul;$i<=$max_ul;$i++) {
        if ($dynamic eq 'no') {
            $dymscore[$i]=$match_score;
        } else {
            if ((($i-1)*$match_score+$max)<=-0.5) {
                $dymscore[$i]=$match_score;
            } else {
                $dymscore[$i]=(-0.5-$max)/($i-1);
            }
        }
    }
}

sub PrepareSTR { ### $_[0] is the result file by _findseq.pl
    @RawResult=(); ### 2D array for raw realignment results by STRid, first index same to @StrInfo
    %RawResultIndex=(); ### store correlated first [] of @RawResult for each STRid
    @RawResultMax=(); ### max second []s of @RawResult
    my $line='';
    my $count=0;
    my @tempdata=();
    $max_ul=0; ### max length of all STR units, used in &PossibleRepeat & &Prepare_dymscore
    $min_ul=0; ### min length of all STR units, used in &PossibleRepeat & &Prepare_dymscore
    %StrUnit=(); ### Unit of each STR, transformed to %SameUnit index
    %StrUnitcomple=(); ### Unit of each STR, complementary & transformed to %SameUnit index
    %MaxGap=(); ### max gap info of each STR
    %UpFL=(); ### 5FL info of each STR, means length for flank searching at STR upstream
    %DownFL=(); ### 3FL info of each STR, means length for flank searching at STR downstream
    %UpStream=(); ### upstream flanking region of each STR, copied from $_[0] then truncated
    %DownStream=(); ### downstream flanking region of each STR, copied from $_[0] then truncated
    %UpStream_comple=(); ### complementary of %UpStream
    %DownStream_comple=(); ### complementary of %DownStream
    %StrUp=(); ### 5' end region of each STR, copied from $_[0] then truncated
    %StrDown_comple=(); ### complementary of 3' end region of each STR, copied from $_[0] then truncated then complementary
    %SameUnit=(); ### category all STRs by units ("STRid1+ STRid2- .."). Equivalent units after transformation will be in the same category. Value=index for @SameUnitCont
    my $SameUnitIndex=-1;
    @SameUnitCont=(); ### 2d array, [index][1..]="STRid+/-"
    @SameUnitNum=(); ### [index]=number of "STRid+/-" with same index unit
    @StrInfo=(); ### store the info lines of all STRs by original order, for final output
    $match_score=0;
    $mis_score=0;
    $gap_score=0;
    $limit_score=0; ### parameters will be read from 1st line of $_[0]
    $dynamic='';
    $reffile='';
    ### below 3 $skip_ are used to skip impossible positions for flank search
    $skip_minfree=0; ### minimum requirement of free scores
    $skip_misgap=0; ### larger one of $mis_score or $gap_score
    $skip_mislimit=0; ### minimum mismatch score
    my $maxFL=0;
    my $dangerlen=0;
    @FlankMaxMis=(); ### used in &PossibleFlank
    @FlankMaxGap=(); ### used in &flank_search and &PossibleFlank
    @FlankMaxMisScore=(); ### used in &flank_search
    my $i=0;
    $openlimit1=$checklimit;
    if ($openlimit1<10) {
        $openlimit1=10;
    }
    $openlimit2=$openlimit1+20; ### range ($openlimit1 $openlimit2] is used when generating %Open in &PrepareResult
    my $temp='';
    my $longer_length=$search_length;
    my $templen=0;
    open (INPUT_DATA, "<$_[0]") or die "Couldn't open: $!";
    if (-s INPUT_DATA) {
        chomp ($line=<INPUT_DATA>);
        @tempdata=();
        @tempdata=split /\s+/, $line;
        $match_score=$tempdata[2];
        $mis_score=$tempdata[3];
        $gap_score=$tempdata[4];
        $limit_score=$tempdata[5];
        $dynamic=$tempdata[6];
        $maxFL=$tempdata[7];
        $reffile=$tempdata[9];
        if ($maxFL ne 'unlimited') {
            if ($maxFL<$search_length) {
                print ("Warning: flank search range ($search_length, from user input) is larger than max FL ($maxFL, from $Para[0]), may lose realignment results\n");
            }
        }
        if ($gap_score>$mis_score) {
            $dangerlen=0-$gap_score;
            $skip_misgap=$gap_score;
            $skip_minfree=0-$mis_score; ### will +$dymatch_score when use
        } else {
            $dangerlen=0-$mis_score;
            $skip_misgap=$mis_score;
            $skip_minfree=0-$gap_score; ### will +$dymatch_score when use
        }
        $skip_mislimit=$limit_score+$skip_minfree; ### will +$dymatch_score when use, score below $skip_mislimit has risk to be beyond $limit_score when a match turns into mis/gap
        if (0==((0-$limit_score)%$dangerlen)) {
            $dangerlen=(0-$limit_score)/$dangerlen;
        } else {
            $dangerlen=int((0-$limit_score)/$dangerlen)+1;
        }
        for ($i=4;$i>=1;$i--) {
            $FlankMaxMisScore[$i]=0;
        }
        for ($i=$dangerlen;$i>=5;$i--) {
            $FlankMaxMisScore[$i]=3;
        }
        for ($i=1;$i<=3;$i++) {
            push @FlankMaxMisScore,(5);
        }
        for ($i=1;$i<=3;$i++) {
            push @FlankMaxMisScore,(8);
        }
        for ($i=1;$i<=3;$i++) {
            push @FlankMaxMisScore,(9);
        }
        for ($i=1;$i<=5;$i++) {
            push @FlankMaxMisScore,(10);
        }
        for ($i=1;$i<=10;$i++) {
            push @FlankMaxMisScore,(12);
        }
        while ($#FlankMaxMisScore<$search_length) {
            push @FlankMaxMisScore,(15);
        }
        for ($i=1;$i<=$#FlankMaxMisScore;$i++) {
            $FlankMaxGap[$i]=int($FlankMaxMisScore[$i]/(0-$flankgapscore));
            $FlankMaxMis[$i]=int($FlankMaxMisScore[$i]/(0-max($flankmisscore,$flankgapscore)));
        }
        $longer_length=$longer_length+$FlankMaxGap[$#FlankMaxGap]+1;
        if (($flankcheck_length*2)>$longer_length) {
            $longer_length=$flankcheck_length*2;
        }
        while (1) {
            chomp ($line=<INPUT_DATA>);
            $StrInfo[$count]=$line;
            @tempdata=();
            @tempdata=split /\s+/, $line;
            $RawResultIndex{$tempdata[0]}=$count;
            $RawResultMax[$count]=-1;
            $tempdata[3]=~/:/;
            $MaxGap{$tempdata[0]}=$';
            $tempdata[4]=~/:/;
            $UpFL{$tempdata[0]}=$';
            $tempdata[5]=~/:/;
            $DownFL{$tempdata[0]}=$';
            if ((length($tempdata[1]))>$max_ul) {
                $max_ul=length($tempdata[1]);
            }
            if ((0==$min_ul)or((length($tempdata[1]))<$min_ul)) {
                $min_ul=length($tempdata[1]);
            }
            $temp=&TransformSTR($tempdata[1]);
            if ($temp ne '') {
                $StrUnit{$tempdata[0]}=$temp;
                $SameUnitNum[$SameUnit{$temp}]++;
                $SameUnitCont[$SameUnit{$temp}][$SameUnitNum[$SameUnit{$temp}]]=join '',($tempdata[0],'+');
            } else {
                $StrUnit{$tempdata[0]}=$tempdata[1];
                $SameUnitIndex++;
                $SameUnit{$tempdata[1]}=$SameUnitIndex;
                $SameUnitNum[$SameUnitIndex]=1;
                $SameUnitCont[$SameUnitIndex][1]=join '',($tempdata[0],'+');
            }
            $tempdata[1]=&complementary($tempdata[1]);
            $temp=&TransformSTR($tempdata[1]);
            if ($temp ne '') {
                $StrUnitcomple{$tempdata[0]}=$temp;
                $SameUnitNum[$SameUnit{$temp}]++;
                $SameUnitCont[$SameUnit{$temp}][$SameUnitNum[$SameUnit{$temp}]]=join '',($tempdata[0],'-');
            } else {
                $StrUnitcomple{$tempdata[0]}=$tempdata[1];
                $SameUnitIndex++;
                $SameUnit{$tempdata[1]}=$SameUnitIndex;
                $SameUnitNum[$SameUnitIndex]=1;
                $SameUnitCont[$SameUnitIndex][1]=join '',($tempdata[0],'-');
            }
            chomp ($line=<INPUT_DATA>);
            chomp ($line=<INPUT_DATA>);
            $StrUp{$tempdata[0]}=substr($line,0,$flankcheck_length*2);
            $templen=length($line);
            if ($templen>=($flankcheck_length*2)) {
                $StrDown_comple{$tempdata[0]}=substr($line,$templen-$flankcheck_length*2,$flankcheck_length*2);
            } else {
                $StrDown_comple{$tempdata[0]}=$line;
            }
            $StrDown_comple{$tempdata[0]}=&complementary($StrDown_comple{$tempdata[0]});
            chomp ($line=<INPUT_DATA>);
            chomp ($line=<INPUT_DATA>);
            $templen=length($line);
            if ($templen>=$longer_length) {
                $UpStream{$tempdata[0]}=substr($line,$templen-$longer_length,$longer_length);
            } else {
                $UpStream{$tempdata[0]}=$line;
            }
            $UpStream_comple{$tempdata[0]}=&complementary($UpStream{$tempdata[0]});
            chomp ($line=<INPUT_DATA>);
            chomp ($line=<INPUT_DATA>);
            $templen=length($line);
            $DownStream{$tempdata[0]}=substr($line,0,$longer_length);
            $DownStream_comple{$tempdata[0]}=&complementary($DownStream{$tempdata[0]});
            chomp ($line=<INPUT_DATA>);
            $count++;
            if (eof) {
                last;
            }
        }
    }
    close (INPUT_DATA);
    print ("Finish reading information of $count STRs\n");
}

sub up_extension { ### output nt length can be extended to upstream direction, by all methods
    ### [0]=input string [1]=start from Xth nt to extend (right after extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion, [7]=allowance for error and SNP
    ### if result is integer, it means the extension was fully finished; if result is integer+0.1, it means the extension was not finished due to string length
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=length($_[2]);$y>=1;$y--) { ### store repeat unit in reversed order
        $repeat[length($_[2])-$y+1]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence in reversed order
    my $posnow=$_[1]-1;
    my $refchar='';

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

    my $keepext='yes';
    my $value=0;
    while ($posnow>=1) {
        $refchar=substr($_[0],$posnow-1,1);
        if ($refchar ne "\n") {
            $count++; ### $count=X (=$_[1]-$posnow) means now trying to align the Xth nt
            $ref[$count]=$refchar;

            $tempnew[0]=$bottom[0];
            for ($y=1;$y<=length($_[2]);$y++) {
                if ($tempnew[$y-1]==$limit) { ### any calculation using a stop match is considered as a stop
                    $tempnew[$y]=$limit;
                } elsif ($y<length($_[2])) { ### make it able to start from any nt in repeat unit
                    $tempnew[$y]=$tempnew[0];
                } else {
                    $tempnew[$y]=$tempnew[0]+$gapscore*length($_[2]);
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
                        if (($ref[$x] eq 'N')or($repeat[$y] eq 'N')) {
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

            $keepext='no';
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext='yes'; ### if any score >limit, continue to try next nt
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
                    if (($ref[$count] eq 'N')or($repeat[$y%(length($_[2]))] eq 'N')) {
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
                    $keepext='yes'; ### if any score >limit, continue to try next nt
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

            if ($keepext eq 'no') {
                last;
            } else {
                $posnow--;
            }
        }
    }
    if ($keepext eq 'yes') {
        $result=$result+0.1;
    }
    return($result);
}

sub down_extension { ### output nt length can be extended to downstream direction, by all methods
    ### [0]=input string [1]=start from Xth nt to extend (right before extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion, [7]=allowance for error and SNP
    ### if result is integer, it means the extension was fully finished; if result is integer+0.1, it means the extension was not finished due to string length
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=1;$y<=length($_[2]);$y++) { ### store repeat unit
        $repeat[$y]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence
    my $posnow=$_[1]+1;
    my $refchar='';

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

    my $keepext='yes';
    my $value=0;
    while ($posnow<=length($_[0])) {
        $refchar=substr($_[0],$posnow-1,1);
        if ($refchar ne "\n") {
            $count++; ### $count=X (=$posnow-$_[1]) means now trying to align the Xth nt
            $ref[$count]=$refchar;

            $tempnew[0]=$bottom[0];
            for ($y=1;$y<=length($_[2]);$y++) {
                if ($tempnew[$y-1]==$limit) { ### any calculation using a stop match is considered as a stop
                    $tempnew[$y]=$limit;
                } elsif ($y<length($_[2])) { ### make it able to start from any nt in repeat unit
                    $tempnew[$y]=$tempnew[0];
                } else {
                    $tempnew[$y]=$tempnew[0]+$gapscore*length($_[2]);
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
                        if (($ref[$x] eq 'N')or($repeat[$y] eq 'N')) {
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

            $keepext='no';
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext='yes'; ### if any score >limit, continue to try next nt
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
                    if (($ref[$count] eq 'N')or($repeat[$y%(length($_[2]))] eq 'N')) {
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
                    $keepext='yes'; ### if any score >limit, continue to try next nt
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

            if ($keepext eq 'no') {
                last;
            } else {
                $posnow++;
            }
        }
    }
    if ($keepext eq 'yes') {
        $result=$result+0.1;
    }
    return($result);
}

sub up_extension_fix { ### output nt length can be extended to upstream direction, fix start nt & with extension limit
    ### [0]=input string [1]=start from Xth nt to extend (right after extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion, [7]=allowance for error and SNP
    ### if result is integer, it means the extension was fully finished; if result is integer+0.1, it means the extension was not finished due to string length
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=length($_[2]);$y>=1;$y--) { ### store repeat unit in reversed order
        $repeat[length($_[2])-$y+1]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence in reversed order
    my $posnow=$_[1]-1;
    my $refchar='';

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

    my $keepext='yes';
    my $value=0;
    while ($posnow>=1) {
        $refchar=substr($_[0],$posnow-1,1);
        if ($refchar ne "\n") {
            $count++; ### $count=X (=$_[1]-$posnow) means now trying to align the Xth nt
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
                        if (($ref[$x] eq 'N')or($repeat[$y] eq 'N')) {
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
                        if ($result>=$upexlimit) {
                            return(int($upexlimit/2));
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

            $keepext='no';
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext='yes'; ### if any score >limit, continue to try next nt
                if ($rightnew[0]>=0) {
                    $rightnew[0]=0;
                }
                if ($rightnew[0]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                    $result=$count;
                    if ($result>=$upexlimit) {
                        return(int($upexlimit/2));
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
                    if (($ref[$count] eq 'N')or($repeat[$y%(length($_[2]))] eq 'N')) {
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
                    $keepext='yes'; ### if any score >limit, continue to try next nt
                    if ($rightnew[$y]>=0) {
                        $rightnew[$y]=0;
                    }
                    if ($rightnew[$y]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                        $result=$count;
                        if ($result>=$upexlimit) {
                            return(int($upexlimit/2));
                        }
                    }
                }
            }
            @rightold=@rightnew;
            @rightnew=();
            $bottom[$count]=$rightold[$#rightold];

            if ($keepext eq 'no') {
                last;
            } else {
                $posnow--;
            }
        }
    }
    if ($keepext eq 'yes') {
        $result=$result+0.1;
    }
    return($result);
}

sub down_extension_fix { ### output nt length can be extended to downstream direction, fix start nt & with extension limit
    ### [0]=input string [1]=start from Xth nt to extend (right before extended nts) [2]=sequence to extend (5'-3')
    ### [3]=match score, [4]=mismatch score, [5]=gap score, [6]=limit to stop for STR extendsion, [7]=allowance for error and SNP
    ### if result is integer, it means the extension was fully finished; if result is integer+0.1, it means the extension was not finished due to string length
    my $x=0;
    my $y=0;
    my @repeat=();
    for ($y=1;$y<=length($_[2]);$y++) { ### store repeat unit
        $repeat[$y]=substr($_[2],$y-1,1);
    }
    $repeat[0]=$repeat[length($_[2])];
    my @ref=(); ### store aligned reference sequence
    my $posnow=$_[1]+1;
    my $refchar='';

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

    my $keepext='yes';
    my $value=0;
    while ($posnow<=length($_[0])) {
        $refchar=substr($_[0],$posnow-1,1);
        if ($refchar ne "\n") {
            $count++; ### $count=X (=$posnow-$_[1]) means now trying to align the Xth nt
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
                        if (($ref[$x] eq 'N')or($repeat[$y] eq 'N')) {
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
                        if ($result>=$downexlimit) {
                            return(int($downexlimit/2));
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

            $keepext='no';
            if (($rightold[0])==$limit) {
                $rightnew[0]=$limit;
            } else {
                $rightnew[0]=$rightold[0]+$gapscore;
            }
            if ($rightnew[0]<=$limit) {
                $rightnew[0]=$limit;
            } else {
                $keepext='yes'; ### if any score >limit, continue to try next nt
                if ($rightnew[0]>=0) {
                    $rightnew[0]=0;
                }
                if ($rightnew[0]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                    $result=$count;
                    if ($result>=$downexlimit) {
                        return(int($downexlimit/2));
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
                    if (($ref[$count] eq 'N')or($repeat[$y%(length($_[2]))] eq 'N')) {
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
                    $keepext='yes'; ### if any score >limit, continue to try next nt
                    if ($rightnew[$y]>=0) {
                        $rightnew[$y]=0;
                    }
                    if ($rightnew[$y]>=$allowance) { ### if any score >=$allowance, consider it as successful extension
                        $result=$count;
                        if ($result>=$downexlimit) {
                            return(int($downexlimit/2));
                        }
                    }
                }
            }
            @rightold=@rightnew;
            @rightnew=();
            $bottom[$count]=$rightold[$#rightold];

            if ($keepext eq 'no') {
                last;
            } else {
                $posnow++;
            }
        }
    }
    if ($keepext eq 'yes') {
        $result=$result+0.1;
    }
    return($result);
}

sub skip_extension { ### align $_[0] with multiple perfect & complete $_[1]
    ### [2]=match score, [3]=mismatch score, [4]=gap score
    ### put scores of each nt to $SingleScore[$_[5] ~ ..]
    my $x=0;
    my $y=0;
    my $xlen=length($_[0]);
    my $len=length($_[1]);
    my @repeat=();
    for ($y=1;$y<=$len;$y++) { ### store repeat unit
        $repeat[$y]=substr($_[1],$y-1,1);
    }
    $repeat[0]=$repeat[$len];
    my $ylen=$len*$xlen;
    for ($y=$len+1;$y<=$ylen;$y++) {
        $repeat[$y]=$repeat[$y%$len];
    }
    my @ref=();
    for ($x=1;$x<=$xlen;$x++) { ### store sequence to be aligned
        $ref[$x]=substr($_[0],$x-1,1);
    }
    my $gapscore=$_[4];
    my $misscore=$_[3]; ### N misscore will = $misscore*3/4
    my $matchscore=$_[2];
    ### score matrix ###
    my @ScoreMatrix=();
    my @PathMatrix=();
    my $value=0;
    $ScoreMatrix[0][0]=0;
    for ($y=1;$y<=$ylen;$y++) {
        $ScoreMatrix[0][$y]=$ScoreMatrix[0][$y-1]+$gapscore;
        $PathMatrix[0][$y]='y';
    }
    for ($x=1;$x<=$xlen;$x++) {
        $ScoreMatrix[$x][0]=$ScoreMatrix[$x-1][0]+$gapscore;
        $PathMatrix[$x][0]='x';
    }
    for ($x=1;$x<=$xlen;$x++) {
        for ($y=1;$y<=$ylen;$y++) {
            $value=$ScoreMatrix[$x][$y-1]+$gapscore;
            $ScoreMatrix[$x][$y]=$value;
            $PathMatrix[$x][$y]='y';
            $value=$ScoreMatrix[$x-1][$y]+$gapscore;
            if ($value>$ScoreMatrix[$x][$y]) {
                $ScoreMatrix[$x][$y]=$value;
                $PathMatrix[$x][$y]='x';
            }
            if (($ref[$x] eq 'N')or($repeat[$y] eq 'N')) {
                $value=$ScoreMatrix[$x-1][$y-1]+$misscore*0.75;
            } elsif ($ref[$x] ne $repeat[$y]) {
                $value=$ScoreMatrix[$x-1][$y-1]+$misscore;
            } else {
                $value=$ScoreMatrix[$x-1][$y-1]+$matchscore;
            }
            if ($value>$ScoreMatrix[$x][$y]) {
                $ScoreMatrix[$x][$y]=$value;
                $PathMatrix[$x][$y]='xy';
            }
        }
    }
    my $maxy=0;
    $value=$ScoreMatrix[$xlen][0];
    for ($y=$len;$y<=$ylen;$y=$y+$len) {
        if ($ScoreMatrix[$xlen][$y]>$value) {
            $value=$ScoreMatrix[$xlen][$y];
            $maxy=$y;
        }
    }
    $x=$xlen;
    $y=$maxy;
    my $tempsinglescore=0;
    while (($x>0)or($y>0)) {
        if ($PathMatrix[$x][$y] eq 'x') {
            $SingleScore[$x-1+$_[5]]=$ScoreMatrix[$x][$y]-$ScoreMatrix[$x-1][$y];
            $x--;
        } elsif ($PathMatrix[$x][$y] eq 'xy') {
            $SingleScore[$x-1+$_[5]]=$ScoreMatrix[$x][$y]-$ScoreMatrix[$x-1][$y-1];
            $x--;
            $y--;
        } elsif ($x==$xlen) { ### gap from repeat unit counted into next nt, unless this is the last nt
            $tempsinglescore=$tempsinglescore+$ScoreMatrix[$x][$y]-$ScoreMatrix[$x][$y-1];
            $y--;
        } else {
            $SingleScore[$x+$_[5]]=$SingleScore[$x+$_[5]]+$ScoreMatrix[$x][$y]-$ScoreMatrix[$x][$y-1];
            $y--;
        }
    }
    $SingleScore[$xlen-1+$_[5]]=$SingleScore[$xlen-1+$_[5]]+$tempsinglescore;
}

sub PossibleFlank { ### return 1 if possible to find flank, 0 impossible
    ### [0]=string from reads [1]=string from flank
    ### [2][3]=length range of end string from [0] possible to match head string from [1]
    ### [4]=length of $_[0]
    ### [5]=length of $_[1]
    my $splitlen=int(($_[3]-$FlankMaxGap[$_[2]])/($FlankMaxMis[$_[2]]+1));
    if ($splitlen<=0) { ### too short to judge, return possible
        return(1);
    }
    my $qcheck=substr($_[0],$_[4]-$_[2],$_[2]);
    my $rcheck='';
    my $splitnum=$FlankMaxMis[$_[2]]+1;
    if ((($splitnum-1)*$splitlen)>=$_[5]) { ### return possible if flank is too short
        return(1);
    }
    my $i=0;
    for ($i=$splitnum;$i>=1;$i--) {
        $rcheck=substr($_[1],($i-1)*$splitlen,$splitlen);
        if ($qcheck=~/\Q$rcheck\E/) {
            return(1);
        }
    }
    return(0);
}

sub NarrowFlank { ### use &PossibleFlank to narrow down range to search for flank, return new range or (-1,-1) if impossible to find flank
    ### [0]=string from reads [1]=string from flank
    ### [2][3]=length range of end string from [0] possible to match head string from [1]
    ### [4]=length of $_[0]
    ### [5]=length of $_[1]
    if (($_[2]-$_[3])<=0) { ### too short to narrow down, return original
        return($_[2],$_[3]);
    }
    my $div=int(($_[2]+$_[3])/2);
    my $judge1=&PossibleFlank($_[0],$_[1],$_[2],$div+1,$_[4],$_[5]);
    my $judge2=&PossibleFlank($_[0],$_[1],$div,$_[3],$_[4],$_[5]);
    my $newstart=0;
    my $newend=0;
    if ((0==$judge1)and(0==$judge2)) {
        return(-1,-1);
    } elsif ((1==$judge1)and(1==$judge2)) {
        return($_[2],$_[3]);
    } elsif (1==$judge1) {
        ($newstart,$newend)=&NarrowFlank($_[0],$_[1],$_[2],$div+1,$_[4],$_[5]);
    } else {
        ($newstart,$newend)=&NarrowFlank($_[0],$_[1],$div,$_[3],$_[4],$_[5]);
    }
    return($newstart,$newend);
}

sub flank_search { ### return length of the best flank search result, 0 means not found
    ### result may not be reliable if too short
    ### [0]=1/2 to use full Reads1/2 as query string [1]=STRid+/- to get flank sequence as reference
    ### [2]=5 (search towards 5' direction) or 3 (search towards 3' direction)
    ### [3][4]=length range to search the alignment ([3]>=[4])
    ### [5]="yes"/"no" whether to apply seed filter to speed up search
    if (((1!=$_[0])and(2!=$_[0]))or((5!=$_[2])and(3!=$_[2]))or($_[3]<$_[4])) {
        print ("Warning: Invalid parameters for &flank_search, skipped\n");
        return(0);
    }
    my $qstring=$Reads1;
    if (2==$_[0]) {
        $qstring=$Reads2;
    }
    my $rstring='';
    my $rlength=0;
    my $rmax=$_[3]+$FlankMaxGap[$_[3]]+1;
    my $ref=$_[1];
    my $temp=substr($ref,0,length($ref)-1);
    my $index='';
    if ($ref=~/-$/) {
        $index=$StrUnitcomple{$temp};
        if (5==$_[2]) { ### generate @flank in reverse direction
            if (!defined($DownStream_comple{$temp})) {
                print ("Warning: Invalid flank sequence for &flank_search, skipped\n");
                return(0);
            } else {
                $rstring=$DownStream_comple{$temp};
                $rlength=length($rstring);
                if ($rmax<$rlength) {
                    $rstring=substr($rstring,$rlength-$rmax,$rmax);
                    $rlength=$rmax;
                }
                $rstring=reverse($rstring);
                $qstring=substr($qstring,0,$_[3]);
                $qstring=reverse($qstring);
            }
        } else { ### generate @flank in forward direction
            if (!defined($UpStream_comple{$temp})) {
                print ("Warning: Invalid flank sequence for &flank_search, skipped\n");
                return(0);
            } else {
                $rstring=$UpStream_comple{$temp};
                $rlength=length($rstring);
                if ($rmax<$rlength) {
                    $rstring=substr($rstring,0,$rmax);
                    $rlength=$rmax;
                }
                $qstring=substr($qstring,-$_[3],$_[3]);
            }
        }
    } else {
        $index=$StrUnit{$temp};
        if (5==$_[2]) { ### generate @flank in reverse direction
            if (!defined($UpStream{$temp})) {
                print ("Warning: Invalid flank sequence for &flank_search, skipped\n");
                return(0);
            } else {
                $rstring=$UpStream{$temp};
                $rlength=length($rstring);
                if ($rmax<$rlength) {
                    $rstring=substr($rstring,$rlength-$rmax,$rmax);
                    $rlength=$rmax;
                }
                $rstring=reverse($rstring);
                $qstring=substr($qstring,0,$_[3]);
                $qstring=reverse($qstring);
            }
        } else { ### generate @flank in forward direction
            if (!defined($DownStream{$temp})) {
                print ("Warning: Invalid flank sequence for &flank_search, skipped\n");
                return(0);
            } else {
                $rstring=$DownStream{$temp};
                $rlength=length($rstring);
                if ($rmax<$rlength) {
                    $rstring=substr($rstring,0,$rmax);
                    $rlength=$rmax;
                }
                $qstring=substr($qstring,-$_[3],$_[3]);
            }
        }
    }

    my $narrowstart=-1;
    my $narrowend=-1;
    my $tempseed='';
    my $tempstart=0;
    my $temptempstart=0;
    my $tempjudge=0;
    my $fulistjudge='X';
    if (($_[5] eq 'no')or(($flankseed-1)>length($rstring))) {
        $tempjudge=&PossibleFlank($qstring,$rstring,$_[3],$_[4],$_[3],$rlength);
        if (1==$tempjudge) {
            ($narrowstart,$narrowend)=&NarrowFlank($qstring,$rstring,$_[3],$_[4],$_[3],$rlength);
        }
    } else {
        if (5==$_[2]) {
            $tempseed=substr($rstring,0,$flankseed);
            if (defined($FlankSeed5{$tempseed})) {
                $tempstart=$FlankSeed5{$tempseed};
                $fulistjudge.=" $FlankSeed5fulist{$tempseed}";
            } else {
                $tempstart=$flankseed-1;
            }
            $tempseed=substr($rstring,0,$flankseed-1);
            if (defined($FlankSeed5{$tempseed})) {
                $temptempstart=$FlankSeed5{$tempseed};
                $fulistjudge.=" $FlankSeed5fulist{$tempseed}";
            } else {
                $temptempstart=$flankseed-1;
            }
            if ($temptempstart>$tempstart) {
                $tempstart=$temptempstart;
            }
            $tempseed=substr($rstring,0,$flankseed+1);
            if (defined($FlankSeed5{$tempseed})) {
                $temptempstart=$FlankSeed5{$tempseed};
                $fulistjudge.=" $FlankSeed5fulist{$tempseed}";
            } else {
                $temptempstart=$flankseed-1;
            }
            if ($temptempstart>$tempstart) {
                $tempstart=$temptempstart;
            }
        } else {
            $tempseed=substr($rstring,0,$flankseed);
            if (defined($FlankSeed3{$tempseed})) {
                $tempstart=$FlankSeed3{$tempseed};
                $fulistjudge.=" $FlankSeed3fulist{$tempseed}";
            } else {
                $tempstart=$flankseed-1;
            }
            $tempseed=substr($rstring,0,$flankseed-1);
            if (defined($FlankSeed3{$tempseed})) {
                $temptempstart=$FlankSeed3{$tempseed};
                $fulistjudge.=" $FlankSeed3fulist{$tempseed}";
            } else {
                $temptempstart=$flankseed-1;
            }
            if ($temptempstart>$tempstart) {
                $tempstart=$temptempstart;
            }
            $tempseed=substr($rstring,0,$flankseed+1);
            if (defined($FlankSeed3{$tempseed})) {
                $temptempstart=$FlankSeed3{$tempseed};
                $fulistjudge.=" $FlankSeed3fulist{$tempseed}";
            } else {
                $temptempstart=$flankseed-1;
            }
            if ($temptempstart>$tempstart) {
                $tempstart=$temptempstart;
            }
        }
        if ($tempstart>$_[3]) { ### may happen due to %Up/DownFL
            $tempstart=$_[3];
        }
        if ($tempstart>=$_[4]) {
            $tempjudge=&PossibleFlank($qstring,$rstring,$tempstart,$_[4],$_[3],$rlength);
            if (1==$tempjudge) {
                ($narrowstart,$narrowend)=&NarrowFlank($qstring,$rstring,$tempstart,$_[4],$_[3],$rlength);
            }
        }
    }
    if (-1==$narrowstart) {
        return(0);
    }
    my @result=();
    my @resultscore=();
    my $score=1; ### score (<=0) from current alignment, =1 means not defined yet
    @flank=(); ### store nt processed from $_[1] and $_[2], will be used in &flank_extension
    @query=(); ### store nt processed from $_[0], will be used in &flank_extension
    $flank[0]='X';
    $query[0]='X';
    push @flank,(split('',$rstring));
    push @query,(split('',$qstring));
    $index=join ' ',($_[0],$index);
    my $indexnum=0;
    if (defined($SkipIndex{$index})) {
        $indexnum=$SkipIndex{$index};
    }
    my $bad=0;
    my $i=0;
    my $j=0;
    my @fljudge=();
    for ($i=$flankseed-1;$i>=1;$i--) {
        $fljudge[$i]=1;
    }
    my @tempdata=split / /, $fulistjudge;
    for ($i=$#tempdata;$i>=1;$i--) {
        $fljudge[$tempdata[$i]]=1;
    }
    for ($i=$narrowstart;$i>=$narrowend;$i--) {
        $bad=-1-$FlankMaxMisScore[$i];
        if (defined($FlankScore{"$_[0] $_[1] $_[2] $i"})) { ### %FlankScore defined in &ReAligner and records existed results to save time
            $score=$FlankScore{"$_[0] $_[1] $_[2] $i"};
        } else {
            if (((5==$_[2])and(defined($Skip5[$indexnum][$i])))or((3==$_[2])and(defined($Skip3[$indexnum][$i])))) {
                $score=$bad;
            } elsif (($_[5] eq 'yes')and(!defined($fljudge[$i]))) {
                $score=$bad;
            } elsif (-1==$bad) {
                if (substr($rstring,0,$i) eq substr($qstring,-$i,$i)) {
                    $score=0;
                } else {
                    $score=-1;
                }
            } else {
                $score=&flank_extension($i,$bad);
            }
            $FlankScore{"$_[0] $_[1] $_[2] $i"}=$score;
        }
        if ($score>$bad) {
            if (0==$#result) {
                if (($score>=($resultscore[0]+5))and($i>=($result[0]-1))) { ### replace the first alignment by discard last nt, -5 is gapscore used in &flank_extension
                    $result[0]=$i;
                    $resultscore[0]=$score;
                } else {
                    push @result, ($i);
                    push @resultscore, ($score);
                }
            } else {
                push @result, ($i);
                push @resultscore, ($score);
            }
            if (($score<0)and($score==($bad+1))) { ### set a mark if short alignment is done by most mis/gap
                $result[$#result]-=0.1;
            }
        }
        if ($#result>-1) {
            if (($result[0]-$i)>3) { ### search range 5nt from first good flank
                push @result, ($i);
                last;
            }
        }
    }
    if ($#result==$#resultscore) {
        push @result, ($i+1);
    }
    for ($i=0;$i<=$#result-2;$i++) {
        for ($j=$i+1;$j<=$#result-1;$j++) {
            if ($resultscore[$i]<$resultscore[$j]) {
                ($result[$i],$result[$j])=($result[$j],$result[$i]);
                ($resultscore[$i],$resultscore[$j])=($resultscore[$j],$resultscore[$i]);
            }
        }
    }
    return (@result);
}

sub flank_extension { ### return best alignment score of last query nt, or =$_[1] (means extension stops early)
    ### [0]=tail length of query string to align
    ### [1]=stop value
    ### @query @flank are generated in &flank_search
    my $gapscore=$flankgapscore;
    my $misscore=$flankmisscore;
    my $matchscore=$flankmatchscore; ### N from @flank will = $matchscore, otherwise N from @query will = $misscore
    my $bad=$_[1];
    my $x=0;
    my $y=0;
    my $i=0;
    my $maxy=0; ### for a new column $x, only need to check till cell ($x, $maxy+1)
    my $best=0; ### best score of last column
    my @alignscore=();

    $x=$#query-$_[0];
    for ($y=int($bad/$gapscore);$y>=0;$y--) {
        $alignscore[$x][$y]=$y*$gapscore;
    }
    $maxy=int($bad/$gapscore);
    for ($x=$#query-$_[0]+1;$x<=$#query;$x++) {
        $alignscore[$x][0]=$alignscore[$x-1][0]+$gapscore;
        for ($y=1;$y<=($maxy+1);$y++) {
            if ($y>$#flank) { ### skip calculating the rows beyond flank sequence
                last;
            }
            if (($query[$x] eq $flank[$y])or($flank[$y] eq 'N')) {
                $alignscore[$x][$y]=$alignscore[$x-1][$y-1]+$matchscore;
            } else {
                $alignscore[$x][$y]=$alignscore[$x-1][$y-1]+$misscore;
            }
            if (($alignscore[$x][$y-1]+$gapscore)>$alignscore[$x][$y]) {
                $alignscore[$x][$y]=$alignscore[$x][$y-1]+$gapscore;
            }
            if ($y<=$maxy) {
                if (($alignscore[$x-1][$y]+$gapscore)>$alignscore[$x][$y]) {
                    $alignscore[$x][$y]=$alignscore[$x-1][$y]+$gapscore;
                }
            }
        }
        $y--;
        while ($alignscore[$x][$y]>$bad) {
            $alignscore[$x][$y+1]=$alignscore[$x][$y]+$gapscore;
            $y++;
        }
        $best=$alignscore[$x][0];
        $maxy=0;
        for ($i=0;$i<=$y;$i++) {
            if ($alignscore[$x][$i]>$best) {
                $best=$alignscore[$x][$i];
            }
            if ($alignscore[$x][$i]>$bad) {
                $maxy=$i;
            }
        }
        if ($best<=$bad) {
            last;
        }
    }

    return($best);
}

sub softclip_realign { ### align $_[0] against $_[1], store lengths with best alignment score of last $_[0] nt into @AL for &getRaw
    @AL=();
    my @queryseq=(); ### store nt from $_[0]
    my @refseq=(); ### store nt processed from $_[1]
    my $rlength=length($_[1]);
    my $qlength=length($_[0]);
    if ($rlength>($qlength*2)) {
        $rlength=$qlength*2;
    }
    my $x=0;
    my $y=0;
    for ($y=1;$y<=$rlength;$y++) { ### generate @refseq in forward direction
        $refseq[$y]=substr($_[1],$y-1,1);
    }
    for ($x=1;$x<=$qlength;$x++) { ### generate @queryseq in forward direction
        $queryseq[$x]=substr($_[0],$x-1,1);
    }

    my $gapscore=-1;
    my $misscore=-1;
    my $matchscore=0; ### N-N will = $matchscore, otherwise = $misscore
    my @alignscore=(); ### scores of last checked query nt
    my @alignnow=(); ### scores of the query nt being checked
    $alignscore[0]=0;
    for ($y=1;$y<=$rlength;$y++) {
        $alignscore[$y]=$alignscore[$y-1]+$gapscore;
    }
    for ($x=1;$x<=$qlength;$x++) {
        $alignnow[0]=$alignscore[0]+$gapscore;
        for ($y=1;$y<=$rlength;$y++) {
            if ($queryseq[$x] eq $refseq[$y]) {
                $alignnow[$y]=$alignscore[$y-1]+$matchscore;
            } else {
                $alignnow[$y]=$alignscore[$y-1]+$misscore;
            }
            if (($alignnow[$y-1]+$gapscore)>$alignnow[$y]) {
                $alignnow[$y]=$alignnow[$y-1]+$gapscore;
            }
            if (($alignscore[$y]+$gapscore)>$alignnow[$y]) {
                $alignnow[$y]=$alignscore[$y]+$gapscore;
            }
        }
        @alignscore=@alignnow;
        @alignnow=();
    }

    my $maxscore=$alignscore[0];
    for ($y=1;$y<=$rlength;$y++) {
        if ($alignscore[$y]>$maxscore) {
            $maxscore=$alignscore[$y];
        }
    }
    for ($y=0;$y<=$rlength;$y++) {
        if ($alignscore[$y]==$maxscore) {
            $AL[$#AL+1]=$y;
        }
    }
}

sub getRaw { ### from @Lines1/2 to @Result1/2 1st round
    ### $_[0..4] +/- refID U/D pos Cigar from @Lines1/2
    ### $_[5] $Reads1/$Reads2, $_[6] reads length
    my $noS=$_[4];
    my $temp='';
    my @result=();
    my $clip5=0;
    my $clip3=0;
    while ($noS=~/^\d*[SHI]/) {
        $noS=$';
        $temp=substr($&,0,length($&)-1);
        if ($temp eq '') {
            $temp=1;
        }
        $clip5=$clip5+$temp;
    }
    while ($noS=~/\d*[SHI]$/) {
        $noS=$`;
        $temp=substr($&,0,length($&)-1);
        if ($temp eq '') {
            $temp=1;
        }
        $clip3=$clip3+$temp;
    }
    if ($noS=~/[SH]/) { ### skip if CIGAR has abnormal clip
        return(@result);
    }
    if ((($_[2] eq 'U')and($clip5>$fakesoftin))or(($_[2] eq 'D')and($clip3>$fakesoftin))) { ### skip if clip towards inner flank is too long
        return(@result);
    }

    my $end=$_[3]-1+$_[6]-$clip5-$clip3;
    my $NDI='';
    my $refname=$_[1];
    my $reflen=0;
    my $i=0;
    if ($_[2] eq 'U') {
        if (0==$clip3) { ### fully mapped to upstream
            @result=($_[6],0,0);
        } else {
            while ($noS=~/\d*[NDI]/) { ### get end position of mapped part
                $noS=$';
                $NDI=substr($&,-1,1);
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                if ($NDI eq 'I') {
                    $end=$end-$temp;
                } else {
                    $end=$end+$temp;
                }
            }
            $refname.='_UpStream';
            $reflen=$RefLength{$refname};
            if (($end<($reflen-$fakesoft))and($clip3>$fakesoft)) { ### skip if clip is real but mapped part doesn't reach reference end
                return(@result);
            } elsif ($end>=$reflen) { ### mapped part reaches reference end, clip is real
                @result=($_[6]-$clip3,-1,-1);
            } elsif (($end+$clip3)<=$reflen) { ### fake clip within reference
                @result=($_[6],0,0);
            } else { ### get exact length of fake clip, the rest is real clip
                if ($_[0] eq '+') {
                    &softclip_realign(&complementary(substr($UpStream_comple{$_[1]},0,$reflen-$end)),substr($_[5],$_[6]-$clip3,$clip3));
                } else {
                    &softclip_realign(&complementary(substr($UpStream_comple{$_[1]},0,$reflen-$end)),&complementary(substr($_[5],0,$clip3)));
                }
                for ($i=0;$i<=$#AL;$i++) {
                    if ($AL[$i]==$clip3) {
                        push @result, ($_[6],0,0);
                    } else {
                        push @result, ($_[6]-$clip3+$AL[$i],-1,-1);
                    }
                }
            }
        }
    } else {
        if (0==$clip5) { ### fully mapped to downstream
            @result=(0,0,$_[6]);
        } else {
            if (($_[3]>($fakesoft+1))and($clip5>$fakesoft)) { ### ### skip if clip is real but mapped part doesn't start from reference head
                return(@result);
            } elsif ($_[3]<=1) { ### mapped part starts from reference head, clip is real
                @result=(-1,-1,$_[6]-$clip5);
            } elsif (($_[3]-$clip5)>=1) { ### fake clip within reference
                @result=(0,0,$_[6]);
            } else { ### get exact length of fake clip, the rest is real clip
                if ($_[0] eq '+') {
                    &softclip_realign(&complementary(substr($DownStream{$_[1]},0,$_[3]-1)),&complementary(substr($_[5],0,$clip5)));
                } else {
                    &softclip_realign(&complementary(substr($DownStream{$_[1]},0,$_[3]-1)),substr($_[5],$_[6]-$clip5,$clip5));
                }
                for ($i=0;$i<=$#AL;$i++) {
                    if ($AL[$i]==$clip5) {
                        push @result, (0,0,$_[6]);
                    } else {
                        push @result, (-1,-1,$_[6]-$clip5+$AL[$i]);
                    }
                }
            }
        }
    }
    return(@result);
}

sub getRawRepeat { ### from @Lines1/2Repeat & @Lines1/2Fake to generate %Repeat1/2, %RepeatSTR1/2, %Repeattogo1/2, %STRtogo1/2, $seedpos1/2
    my $i=0;
    my $clip5=0;
    my $clip3=0;
    my $clip5fix=0;
    my $clip3fix=0;
    my $clip5more=0; ### fragile mapping near reads end will not be trusted to skip flank search
    my $clip3more=0; ### fragile mapping near reads end will not be trusted to skip flank search
    my $noS='';
    my $tempnoS='';
    my $temp='';
    my $upext=0;
    my $downext=0;
    my $uplimit=0;
    my $downlimit=0;
    my $upext2=0;
    my $downext2=0;
    my $uplimit2=0;
    my $downlimit2=0;
    my %extended=(); ### index="1/2 unit CIGAR" to skip duplicated mapping results
    my $unit='';
    my $unitrev='';
    my $dymatch_score=0;
    my %Repeat1pre=();
    my %Repeat2pre=();
    ### index=STR unit (same to the index of %SameUnit), value="start end mapstart+$flanklimit_bymapping mapend-$flanklimit_bymapping ..."
    ### means $Reads1/2 can find _Fake/_Repeat map result from start to end (overlapped [start end]s will be combined into one)
    ### mapstart/end+/-$flanklimit_bymapping are used to speed up flank search, make them = reads length/1 when not supported by _Fake mapping
    %Repeat1=();
    %Repeat2=();
    ### index=STR unit (same to the index of %SameUnit), value="start end mapstart+$flanklimit_bymapping mapend-$flanklimit_bymapping ..."
    ### means $Reads1/2 can find STR unit from start to end (overlapped [start end]s will be combined into one)
    ### start=(X-1).9/end=X.1 means the extension was not done till X, but not fully finished due to reads length
    ### mapstart/end+/-$flanklimit_bymapping are used to speed up flank search, make them = reads length/1 when not supported by _Fake mapping
    %Repeat1STR=();
    %Repeat2STR=();
    ### index=STRid+/-, value="start end mapstart+$flanklimit_bymapping mapend-$flanklimit_bymapping ..."
    ### modify %Repeat1/2 to add in _Repeat mapping results for each STRid+/-, used for generating %STRtogo1/2 & to speed up flank search for specific STR in most cases
    %Repeattogo1=();
    %Repeattogo2=();
    ### index=STR unit (same to the index of %SameUnit), value = the longest "start end mapstart+$flanklimit_bymapping mapend-$flanklimit_bymapping" for this unit
    ### or value = "y" for further processing in &ReAligner
    ### will initiate repeat-guided realignment and check all qualified STR flanks
    %Repeattogo1list=();
    %Repeattogo2list=();
    ### index=STR unit (same to the index of %SameUnit), a pre list for %Repeattogo1/2
    %STRtogo1=();
    %STRtogo2=();
    ### index=STRid+/-, value = the longest "start end mapstart+$flanklimit_bymapping mapend-$flanklimit_bymapping" for related repeat unit
    ### will initiate repeat-guided realignment but check only for STRid+/- flanks

    if (($#Lines1Repeat>0)or($#Lines1Fake>0)) {
        $Reads1rev=reverse($Reads1);
    }
    if (($#Lines2Repeat>0)or($#Lines2Fake>0)) {
        $Reads2rev=reverse($Reads2);
    }

    for ($i=$#Lines1Fake-2;$i>=0;$i=$i-3) { ### process @Lines1Fake to generate %Repeat1pre
        $unit=$Lines1Fake[$i+1];
        $Repeattogo1list{$unit}=1;
        $noS=$Lines1Fake[$i+2];
        $noS=~s/H/S/g;
        if (!(defined($extended{"1 $unit $noS"}))) {
            $extended{"1 $unit $noS"}=1;
            $clip5=0;
            while ($noS=~/^\d*[SI]/) {
                $noS=$';
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip5=$clip5+$temp;
            }
            $clip3=0;
            while ($noS=~/\d*[SI]$/) {
                $noS=$`;
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip3=$clip3+$temp;
            }
            if ($noS=~/S/) { ### skip if CIGAR has abnormal clip
                next;
            }
            $clip5more=$clip5;
            if ($noS=~/^\d*M/) {
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                if ($temp<=$untrustmap) {
                    $clip5more+=$temp;
                    $tempnoS=$';
                    if ($tempnoS=~/^\d*I/) {
                        $temp=substr($&,0,length($&)-1);
                        if ($temp eq '') {
                            $temp=1;
                        }
                        $clip5more+=$temp;
                    }
                }
            }
            $clip3more=$clip3;
            if ($noS=~/\d*M$/) {
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                if ($temp<=$untrustmap) {
                    $clip3more+=$temp;
                    $tempnoS=$`;
                    if ($tempnoS=~/\d*I$/) {
                        $temp=substr($&,0,length($&)-1);
                        if ($temp eq '') {
                            $temp=1;
                        }
                        $clip3more+=$temp;
                    }
                }
            }
            if ($Lines1Fake[$i] eq '-') {
                ($clip5,$clip3)=($clip3,$clip5);
                ($clip5more,$clip3more)=($clip3more,$clip5more);
            }
            if (!(defined($Repeat1pre{$unit}))) { ### add a new [start end startlimit endlimit] to $Repeat1pre{$unit}
                $Repeat1pre{$unit}=join ' ',($clip5+1,$RLen1-$clip3,$clip5more+$flanklimit_bymapping,$RLen1-$clip3more+1-$flanklimit_bymapping);
            } else {
                $Repeat1pre{$unit}=&AddRepeat($clip5+1,$RLen1-$clip3,$clip5more+$flanklimit_bymapping,$RLen1-$clip3more+1-$flanklimit_bymapping,$Repeat1pre{$unit});
            }
        }
    }
    for ($i=$#Lines1Repeat-2;$i>=0;$i=$i-3) { ### process @Lines1Repeat to generate %Repeat1pre
        if ($Lines1Repeat[$i] eq '-') {
            $unit=$StrUnitcomple{$Lines1Repeat[$i+1]};
        } else {
            $unit=$StrUnit{$Lines1Repeat[$i+1]};
        }
        $noS=$Lines1Repeat[$i+2];
        $noS=~s/H/S/g;
        if (!(defined($extended{"1 $unit $noS"}))) {
            $extended{"1 $unit $noS"}=1;
            $clip5=0;
            while ($noS=~/^\d*[SI]/) {
                $noS=$';
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip5=$clip5+$temp;
            }
            $clip3=0;
            while ($noS=~/\d*[SI]$/) {
                $noS=$`;
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip3=$clip3+$temp;
            }
            if ($noS=~/S/) { ### skip if CIGAR has abnormal clip
                next;
            }
            if ($Lines1Repeat[$i] eq '-') {
                ($clip5,$clip3)=($clip3,$clip5);
            }
            if (!(defined($Repeat1pre{$unit}))) { ### add a new [start end startlimit endlimit] to $Repeat1pre{$unit}
                $Repeat1pre{$unit}=join ' ',($clip5+1,$RLen1-$clip3,$RLen1,1);
            } else {
                $Repeat1pre{$unit}=&AddRepeat($clip5+1,$RLen1-$clip3,$RLen1,1,$Repeat1pre{$unit});
            }
        }
    }
    for ($i=$#Lines2Fake-2;$i>=0;$i=$i-3) { ### process @Lines2Fake to generate %Repeat2pre
        $unit=$Lines2Fake[$i+1];
        $Repeattogo2list{$unit}=1;
        $noS=$Lines2Fake[$i+2];
        $noS=~s/H/S/g;
        if (!(defined($extended{"2 $unit $noS"}))) {
            $extended{"2 $unit $noS"}=1;
            $clip5=0;
            while ($noS=~/^\d*[SI]/) {
                $noS=$';
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip5=$clip5+$temp;
            }
            $clip3=0;
            while ($noS=~/\d*[SI]$/) {
                $noS=$`;
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip3=$clip3+$temp;
            }
            if ($noS=~/S/) { ### skip if CIGAR has abnormal clip
                next;
            }
            $clip5more=$clip5;
            if ($noS=~/^\d*M/) {
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                if ($temp<=$untrustmap) {
                    $clip5more+=$temp;
                    $tempnoS=$';
                    if ($tempnoS=~/^\d*I/) {
                        $temp=substr($&,0,length($&)-1);
                        if ($temp eq '') {
                            $temp=1;
                        }
                        $clip5more+=$temp;
                    }
                }
            }
            $clip3more=$clip3;
            if ($noS=~/\d*M$/) {
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                if ($temp<=$untrustmap) {
                    $clip3more+=$temp;
                    $tempnoS=$`;
                    if ($tempnoS=~/\d*I$/) {
                        $temp=substr($&,0,length($&)-1);
                        if ($temp eq '') {
                            $temp=1;
                        }
                        $clip3more+=$temp;
                    }
                }
            }
            if ($Lines2Fake[$i] eq '-') {
                ($clip5,$clip3)=($clip3,$clip5);
                ($clip5more,$clip3more)=($clip3more,$clip5more);
            }
            if (!(defined($Repeat2pre{$unit}))) { ### add a new [start end startlimit endlimit] to $Repeat2pre{$unit}
                $Repeat2pre{$unit}=join ' ',($clip5+1,$RLen2-$clip3,$clip5more+$flanklimit_bymapping,$RLen2-$clip3more+1-$flanklimit_bymapping);
            } else {
                $Repeat2pre{$unit}=&AddRepeat($clip5+1,$RLen2-$clip3,$clip5more+$flanklimit_bymapping,$RLen2-$clip3more+1-$flanklimit_bymapping,$Repeat2pre{$unit});
            }
        }
    }
    for ($i=$#Lines2Repeat-2;$i>=0;$i=$i-3) { ### process @Lines2Repeat to generate %Repeat2pre
        if ($Lines2Repeat[$i] eq '-') {
            $unit=$StrUnitcomple{$Lines2Repeat[$i+1]};
        } else {
            $unit=$StrUnit{$Lines2Repeat[$i+1]};
        }
        $noS=$Lines2Repeat[$i+2];
        $noS=~s/H/S/g;
        if (!(defined($extended{"2 $unit $noS"}))) {
            $extended{"2 $unit $noS"}=1;
            $clip5=0;
            while ($noS=~/^\d*[SI]/) {
                $noS=$';
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip5=$clip5+$temp;
            }
            $clip3=0;
            while ($noS=~/\d*[SI]$/) {
                $noS=$`;
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip3=$clip3+$temp;
            }
            if ($noS=~/S/) { ### skip if CIGAR has abnormal clip
                next;
            }
            if ($Lines2Repeat[$i] eq '-') {
                ($clip5,$clip3)=($clip3,$clip5);
            }
            if (!(defined($Repeat2pre{$unit}))) { ### add a new [start end startlimit endlimit] to $Repeat2pre{$unit}
                $Repeat2pre{$unit}=join ' ',($clip5+1,$RLen2-$clip3,$RLen2,1);
            } else {
                $Repeat2pre{$unit}=&AddRepeat($clip5+1,$RLen2-$clip3,$RLen2,1,$Repeat2pre{$unit});
            }
        }
    }

    my @tempdata=();
    foreach $unit (keys(%Repeat1pre)) { ### process %Repeat1pre to generate %Repeat1
        $dymatch_score=$dymscore[length($unit)];
        $unitrev=reverse($unit);
        @tempdata=split / /, $Repeat1pre{$unit};
        $clip5=$tempdata[0]-1; ### add first [start end startlimit endlimit] to $Repeat1{$unit}
        $clip3=$RLen1-$tempdata[1];
        if (0==$clip5) {
            $clip5fix=0;
        } else {
            $Reads1map=substr($Reads1,$clip5,$RLen1-$clip5-$clip3);
            if ($Reads1map=~/\Q$unit\E/) {
                $clip5fix=length($`)+$clip5;
            } else {
                $clip5fix=$clip5;
            }
        }
        if (0==$clip3) {
            $clip3fix=0;
        } else {
            $Reads1revmap=substr($Reads1rev,$clip3,$RLen1-$clip5-$clip3);
            if ($Reads1revmap=~/\Q$unitrev\E/) {
                $clip3fix=length($`)+$clip3;
            } else {
                $clip3fix=$clip3;
            }
        }
        $upext=$clip5fix+1-&up_extension_fix($Reads1,$clip5fix+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
        $downext=$RLen1-$clip3fix+&down_extension_fix($Reads1,$RLen1-$clip3fix,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
        if ($upext>($clip5+1)) {
            $upext=$clip5+1;
        }
        if ($downext<($RLen1-$clip3)) {
            $downext=$RLen1-$clip3;
        }
        $Repeat1{$unit}=join ' ',($upext,$downext,$tempdata[2],$tempdata[3]);
        for ($i=$#tempdata-3;$i>=4;$i=$i-4) { ### add more new [start end startlimit endlimit] to $Repeat1{$unit}
            $clip5=$tempdata[$i]-1;
            $clip3=$RLen1-$tempdata[$i+1];
            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($clip5+1,$Repeat1{$unit});
            ($upext2,$downext2,$uplimit2,$downlimit2)=&SearchRepeat($RLen1-$clip3,$Repeat1{$unit});
            if (($upext!=$upext2)or($downext!=$downext2)or(0==$upext)) {
                if (0==$upext) {
                    if (0==$clip5) {
                        $clip5fix=0;
                    } else {
                        $Reads1map=substr($Reads1,$clip5,$RLen1-$clip5-$clip3);
                        if ($Reads1map=~/\Q$unit\E/) {
                            $clip5fix=length($`)+$clip5;
                        } else {
                            $clip5fix=$clip5;
                        }
                    }
                    $upext=$clip5fix+1-&up_extension_fix($Reads1,$clip5fix+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                    if ($upext>($clip5+1)) {
                        $upext=$clip5+1;
                    }
                    $uplimit=$tempdata[$i+2];
                } elsif ($uplimit>$tempdata[$i+2]) {
                    $uplimit=$tempdata[$i+2];
                }
                if (0==$upext2) {
                    if (0==$clip3) {
                        $clip3fix=0;
                    } else {
                        $Reads1revmap=substr($Reads1rev,$clip3,$RLen1-$clip5-$clip3);
                        if ($Reads1revmap=~/\Q$unitrev\E/) {
                            $clip3fix=length($`)+$clip3;
                        } else {
                            $clip3fix=$clip3;
                        }
                    }
                    $downext2=$RLen1-$clip3fix+&down_extension_fix($Reads1,$RLen1-$clip3fix,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                    if ($downext2<($RLen1-$clip3)) {
                        $downext2=$RLen1-$clip3;
                    }
                    $downlimit2=$tempdata[$i+3];
                } elsif ($downlimit2<$tempdata[$i+3]) {
                    $downlimit2=$tempdata[$i+3];
                }
                $Repeat1{$unit}=&AddRepeat($upext,$downext2,$uplimit,$downlimit2,$Repeat1{$unit});
            } elsif (($uplimit>$tempdata[$i+2])or($downlimit2<$tempdata[$i+3])) {
                $Repeat1{$unit}=&AddRepeat($upext,$downext2,$tempdata[$i+2],$tempdata[$i+3],$Repeat1{$unit});
            }
        }
    }
    foreach $unit (keys(%Repeat2pre)) { ### process %Repeat2pre to generate %Repeat2
        $dymatch_score=$dymscore[length($unit)];
        $unitrev=reverse($unit);
        @tempdata=split / /, $Repeat2pre{$unit};
        $clip5=$tempdata[0]-1; ### add first [start end startlimit endlimit] to $Repeat2{$unit}
        $clip3=$RLen2-$tempdata[1];
        if (0==$clip5) {
            $clip5fix=0;
        } else {
            $Reads2map=substr($Reads2,$clip5,$RLen2-$clip5-$clip3);
            if ($Reads2map=~/\Q$unit\E/) {
                $clip5fix=length($`)+$clip5;
            } else {
                $clip5fix=$clip5;
            }
        }
        if (0==$clip3) {
            $clip3fix=0;
        } else {
            $Reads2revmap=substr($Reads2rev,$clip3,$RLen2-$clip5-$clip3);
            if ($Reads2revmap=~/\Q$unitrev\E/) {
                $clip3fix=length($`)+$clip3;
            } else {
                $clip3fix=$clip3;
            }
        }
        $upext=$clip5fix+1-&up_extension_fix($Reads2,$clip5fix+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
        $downext=$RLen2-$clip3fix+&down_extension_fix($Reads2,$RLen2-$clip3fix,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
        if ($upext>($clip5+1)) {
            $upext=$clip5+1;
        }
        if ($downext<($RLen2-$clip3)) {
            $downext=$RLen2-$clip3;
        }
        $Repeat2{$unit}=join ' ',($upext,$downext,$tempdata[2],$tempdata[3]);
        for ($i=$#tempdata-3;$i>=4;$i=$i-4) { ### add more new [start end startlimit endlimit] to $Repeat2{$unit}
            $clip5=$tempdata[$i]-1;
            $clip3=$RLen2-$tempdata[$i+1];
            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($clip5+1,$Repeat2{$unit});
            ($upext2,$downext2,$uplimit2,$downlimit2)=&SearchRepeat($RLen2-$clip3,$Repeat2{$unit});
            if (($upext!=$upext2)or($downext!=$downext2)or(0==$upext)) {
                if (0==$upext) {
                    if (0==$clip5) {
                        $clip5fix=0;
                    } else {
                        $Reads2map=substr($Reads2,$clip5,$RLen2-$clip5-$clip3);
                        if ($Reads2map=~/\Q$unit\E/) {
                            $clip5fix=length($`)+$clip5;
                        } else {
                            $clip5fix=$clip5;
                        }
                    }
                    $upext=$clip5fix+1-&up_extension_fix($Reads2,$clip5fix+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                    if ($upext>($clip5+1)) {
                        $upext=$clip5+1;
                    }
                    $uplimit=$tempdata[$i+2];
                } elsif ($uplimit>$tempdata[$i+2]) {
                    $uplimit=$tempdata[$i+2];
                }
                if (0==$upext2) {
                    if (0==$clip3) {
                        $clip3fix=0;
                    } else {
                        $Reads2revmap=substr($Reads2rev,$clip3,$RLen2-$clip5-$clip3);
                        if ($Reads2revmap=~/\Q$unitrev\E/) {
                            $clip3fix=length($`)+$clip3;
                        } else {
                            $clip3fix=$clip3;
                        }
                    }
                    $downext2=$RLen2-$clip3fix+&down_extension_fix($Reads2,$RLen2-$clip3fix,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                    if ($downext2<($RLen2-$clip3)) {
                        $downext2=$RLen2-$clip3;
                    }
                    $downlimit2=$tempdata[$i+3];
                } elsif ($downlimit2<$tempdata[$i+3]) {
                    $downlimit2=$tempdata[$i+3];
                }
                $Repeat2{$unit}=&AddRepeat($upext,$downext2,$uplimit,$downlimit2,$Repeat2{$unit});
            } elsif (($uplimit>$tempdata[$i+2])or($downlimit2<$tempdata[$i+3])) {
                $Repeat2{$unit}=&AddRepeat($upext,$downext2,$tempdata[$i+2],$tempdata[$i+3],$Repeat2{$unit});
            }
        }
    }

    my $strid='';
    for ($i=$#Lines1Repeat-2;$i>=0;$i=$i-3) { ### process @Lines1Repeat to generate %Repeat1STR
        $strid=join '',($Lines1Repeat[$i+1],$Lines1Repeat[$i]);
        $noS=$Lines1Repeat[$i+2];
        $noS=~s/H/S/g;
        if (!(defined($extended{"1 $strid $noS"}))) {
            $extended{"1 $strid $noS"}=1;
            $clip5=0;
            while ($noS=~/^\d*[SI]/) {
                $noS=$';
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip5=$clip5+$temp;
            }
            $clip3=0;
            while ($noS=~/\d*[SI]$/) {
                $noS=$`;
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip3=$clip3+$temp;
            }
            if ($noS=~/S/) { ### skip if CIGAR has abnormal clip
                next;
            }
            $clip5more=$clip5;
            if ($noS=~/^\d*M/) {
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                if ($temp<=$untrustmap) {
                    $clip5more+=$temp;
                    $tempnoS=$';
                    if ($tempnoS=~/^\d*I/) {
                        $temp=substr($&,0,length($&)-1);
                        if ($temp eq '') {
                            $temp=1;
                        }
                        $clip5more+=$temp;
                    }
                }
            }
            $clip3more=$clip3;
            if ($noS=~/\d*M$/) {
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                if ($temp<=$untrustmap) {
                    $clip3more+=$temp;
                    $tempnoS=$`;
                    if ($tempnoS=~/\d*I$/) {
                        $temp=substr($&,0,length($&)-1);
                        if ($temp eq '') {
                            $temp=1;
                        }
                        $clip3more+=$temp;
                    }
                }
            }
            if ($Lines1Repeat[$i] eq '-') {
                ($clip5,$clip3)=($clip3,$clip5);
                ($clip5more,$clip3more)=($clip3more,$clip5more);
            }
            if (!(defined($Repeat1STR{$strid}))) { ### copy from $Repeat1 to start
                if ($Lines1Repeat[$i] eq '-') {
                    $unit=$StrUnitcomple{$Lines1Repeat[$i+1]};
                } else {
                    $unit=$StrUnit{$Lines1Repeat[$i+1]};
                }
                $Repeat1STR{$strid}=$Repeat1{$unit};
            }
            $Repeat1STR{$strid}=&AddRepeat($clip5+1,$RLen1-$clip3,$clip5more+$flanklimit_bymapping,$RLen1-$clip3more+1-$flanklimit_bymapping,$Repeat1STR{$strid});
        }
    }
    for ($i=$#Lines2Repeat-2;$i>=0;$i=$i-3) { ### process @Lines2Repeat to generate %Repeat2STR
        $strid=join '',($Lines2Repeat[$i+1],$Lines2Repeat[$i]);
        $noS=$Lines2Repeat[$i+2];
        $noS=~s/H/S/g;
        if (!(defined($extended{"2 $strid $noS"}))) {
            $extended{"2 $strid $noS"}=1;
            $clip5=0;
            while ($noS=~/^\d*[SI]/) {
                $noS=$';
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip5=$clip5+$temp;
            }
            $clip3=0;
            while ($noS=~/\d*[SI]$/) {
                $noS=$`;
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                $clip3=$clip3+$temp;
            }
            if ($noS=~/S/) { ### skip if CIGAR has abnormal clip
                next;
            }
            $clip5more=$clip5;
            if ($noS=~/^\d*M/) {
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                if ($temp<=$untrustmap) {
                    $clip5more+=$temp;
                    $tempnoS=$';
                    if ($tempnoS=~/^\d*I/) {
                        $temp=substr($&,0,length($&)-1);
                        if ($temp eq '') {
                            $temp=1;
                        }
                        $clip5more+=$temp;
                    }
                }
            }
            $clip3more=$clip3;
            if ($noS=~/\d*M$/) {
                $temp=substr($&,0,length($&)-1);
                if ($temp eq '') {
                    $temp=1;
                }
                if ($temp<=$untrustmap) {
                    $clip3more+=$temp;
                    $tempnoS=$`;
                    if ($tempnoS=~/\d*I$/) {
                        $temp=substr($&,0,length($&)-1);
                        if ($temp eq '') {
                            $temp=1;
                        }
                        $clip3more+=$temp;
                    }
                }
            }
            if ($Lines2Repeat[$i] eq '-') {
                ($clip5,$clip3)=($clip3,$clip5);
                ($clip5more,$clip3more)=($clip3more,$clip5more);
            }
            if (!(defined($Repeat2STR{$strid}))) { ### copy from $Repeat2 to start
                if ($Lines2Repeat[$i] eq '-') {
                    $unit=$StrUnitcomple{$Lines2Repeat[$i+1]};
                } else {
                    $unit=$StrUnit{$Lines2Repeat[$i+1]};
                }
                $Repeat2STR{$strid}=$Repeat2{$unit};
            }
            $Repeat2STR{$strid}=&AddRepeat($clip5+1,$RLen2-$clip3,$clip5more+$flanklimit_bymapping,$RLen2-$clip3more+1-$flanklimit_bymapping,$Repeat2STR{$strid});
        }
    }

    my %longest1=(); ### longest extended $Repeat1{}/$Repeat1STR{} if passes $RepeatPercentLimitMap
    my %longest2=(); ### longest extended $Repeat2{}/$Repeat2STR{} if passes $RepeatPercentLimitMap
    my $repeatlength=0;
    foreach (keys(%Repeat1)) {
        $repeatlength=-1;
        @tempdata=split / /, $Repeat1{$_};
        for ($i=$#tempdata-3;$i>=0;$i=$i-4) {
            if ((($tempdata[$i+1]-$tempdata[$i]+1)>=($RLen1*$RepeatPercentLimitMap))and(($tempdata[$i+3]-$tempdata[$i+2]-1+2*$flanklimit_bymapping)>=($RLen1*($RepeatPercentLimitMap-$Ratiomodify)))) {
                if (($tempdata[$i+1]-$tempdata[$i])>$repeatlength) {
                    $repeatlength=$tempdata[$i+1]-$tempdata[$i];
                    ($upext,$downext,$uplimit,$downlimit)=($tempdata[$i],$tempdata[$i+1],$tempdata[$i+2],$tempdata[$i+3]);
                }
            }
        }
        if ($repeatlength>=0) {
            $longest1{$_}="$upext $downext $uplimit $downlimit";
        }
    }
    foreach (keys(%Repeat1STR)) {
        $repeatlength=-1;
        @tempdata=split / /, $Repeat1STR{$_};
        for ($i=$#tempdata-3;$i>=0;$i=$i-4) {
            if ((($tempdata[$i+1]-$tempdata[$i]+1)>=($RLen1*$RepeatPercentLimitMap))and(($tempdata[$i+3]-$tempdata[$i+2]-1+2*$flanklimit_bymapping)>=($RLen1*($RepeatPercentLimitMap-$Ratiomodify)))) {
                if (($tempdata[$i+1]-$tempdata[$i])>$repeatlength) {
                    $repeatlength=$tempdata[$i+1]-$tempdata[$i];
                    ($upext,$downext,$uplimit,$downlimit)=($tempdata[$i],$tempdata[$i+1],$tempdata[$i+2],$tempdata[$i+3]);
                }
            }
        }
        if ($repeatlength>=0) {
            $longest1{$_}="$upext $downext $uplimit $downlimit";
        }
    }
    foreach (keys(%Repeat2)) {
        $repeatlength=-1;
        @tempdata=split / /, $Repeat2{$_};
        for ($i=$#tempdata-3;$i>=0;$i=$i-4) {
            if ((($tempdata[$i+1]-$tempdata[$i]+1)>=($RLen2*$RepeatPercentLimitMap))and(($tempdata[$i+3]-$tempdata[$i+2]-1+2*$flanklimit_bymapping)>=($RLen2*($RepeatPercentLimitMap-$Ratiomodify)))) {
                if (($tempdata[$i+3]-$tempdata[$i+2])>$repeatlength) {
                    $repeatlength=$tempdata[$i+3]-$tempdata[$i+2];
                    ($upext,$downext,$uplimit,$downlimit)=($tempdata[$i],$tempdata[$i+1],$tempdata[$i+2],$tempdata[$i+3]);
                }
            }
        }
        if ($repeatlength>=0) {
            $longest2{$_}="$upext $downext $uplimit $downlimit";
        }
    }
    foreach (keys(%Repeat2STR)) {
        $repeatlength=-1;
        @tempdata=split / /, $Repeat2STR{$_};
        for ($i=$#tempdata-3;$i>=0;$i=$i-4) {
            if ((($tempdata[$i+1]-$tempdata[$i]+1)>=($RLen2*$RepeatPercentLimitMap))and(($tempdata[$i+3]-$tempdata[$i+2]-1+2*$flanklimit_bymapping)>=($RLen2*($RepeatPercentLimitMap-$Ratiomodify)))) {
                if (($tempdata[$i+3]-$tempdata[$i+2])>$repeatlength) {
                    $repeatlength=$tempdata[$i+3]-$tempdata[$i+2];
                    ($upext,$downext,$uplimit,$downlimit)=($tempdata[$i],$tempdata[$i+1],$tempdata[$i+2],$tempdata[$i+3]);
                }
            }
        }
        if ($repeatlength>=0) {
            $longest2{$_}="$upext $downext $uplimit $downlimit";
        }
    }

    for ($i=$#Lines1Repeat-2;$i>=0;$i=$i-3) { ### scan @Lines1Repeat to generate %STRtogo1
        $strid=join '',($Lines1Repeat[$i+1],$Lines1Repeat[$i]);
        if (defined($longest1{$strid})) {
            $STRtogo1{$strid}=$longest1{$strid};
        }
    }
    for ($i=$#Lines2Repeat-2;$i>=0;$i=$i-3) { ### scan @Lines2Repeat to generate %STRtogo2
        $strid=join '',($Lines2Repeat[$i+1],$Lines2Repeat[$i]);
        if (defined($longest2{$strid})) {
            $STRtogo2{$strid}=$longest2{$strid};
        }
    }

    if ($kmer_switch eq 'yes') {
        $seedpos1=&PossibleRepeat($Reads1,1); ### add more candidates to %Repeattogo1list
        $seedpos2=&PossibleRepeat($Reads2,2); ### add more candidates to %Repeattogo2list
    }
    foreach (keys(%Repeattogo1list)) { ### scan %Repeattogo1list to generate %Repeattogo1
        if (defined($longest1{$_})) {
            $Repeattogo1{$_}=$longest1{$_};
        } elsif (!defined($Repeat1{$_})) {
            $Repeattogo1{$_}='y';
        }
    }
    foreach (keys(%Repeattogo2list)) { ### scan %Repeattogo2list to generate %Repeattogo2
        if (defined($longest2{$_})) {
            $Repeattogo2{$_}=$longest2{$_};
        } elsif (!defined($Repeat2{$_})) {
            $Repeattogo2{$_}='y';
        }
    }
}

sub SearchRepeat { ### search if a position $_[0] already within any [start end] from $_[1], return the found (start,end,startlimit,endlimit), or (0,0,0,0) if not found
    ### $_[0]=position, $_[1]=$Repeat1/2{index} ("start end startlimit endlimit ...")
    my $i=0;
    my @tempdata=split / /, $_[1];
    for ($i=$#tempdata-3;$i>=0;$i=$i-4) {
        if (($_[0]>=int($tempdata[$i]+0.5))and($_[0]<=int($tempdata[$i+1]+0.5))) {
            return ($tempdata[$i],$tempdata[$i+1],$tempdata[$i+2],$tempdata[$i+3]);
        }
    }
    return(0,0,0,0);
}

sub AddRepeat { ### process new [start end startlimit endlimit] ($_[0..3]) to see if overlapping with given %Repeat1/2{index} ($_[4])
    ### output processed "start end startlimit endlimit ..." as a new $Repeat1/2{index}
    my $start=$_[0];
    my $end=$_[1];
    my $startlimit=$_[2];
    my $endlimit=$_[3];
    my $i=0;
    my $result='';
    my @tempdata=split / /, $_[4];
    my $needcheck='yes';
    while ($needcheck eq 'yes') {
        $needcheck='no';
        for ($i=$#tempdata-3;$i>=0;$i=$i-4) {
            if ((int($start+0.5)<=int($tempdata[$i+1]+0.5))and(int($end+0.5)>=int($tempdata[$i]+0.5))) { ### overlapped
                if ($tempdata[$i]<$start) {
                    $start=$tempdata[$i];
                }
                if ($tempdata[$i+1]>$end) {
                    $end=$tempdata[$i+1];
                }
                if ($tempdata[$i+2]<$startlimit) {
                    $startlimit=$tempdata[$i+2];
                }
                if ($tempdata[$i+3]>$endlimit) {
                    $endlimit=$tempdata[$i+3];
                }
                $tempdata[$i]=-1;
                $tempdata[$i+1]=-1;
                $tempdata[$i+2]=-1;
                $tempdata[$i+3]=-1;
                $needcheck='yes';
                last;
            }
        }
    }
    $result=join ' ',($start,$end,$startlimit,$endlimit);
    for ($i=0;$i<=$#tempdata;$i++) {
        if ($tempdata[$i]!=-1) {
            $result.=" $tempdata[$i]";
        }
    }
    return($result);
}

sub getPairtype { ### get the pairtype judged by given alignment, return a combined string by spaced elements from "fr"/"rf"/"ff"/"rr", or "" if not suitably paired
    ### [0~3] for reads1 from sequencing, [4~7] for reads2 from sequencing, matters for "ff" and "rr"
    ### [0][4]="STRid+/-" or just "+/-", [1][2][3]=[X1 X2 X3] from %Map1, [5][6][7]=[X1 X2 X3] from %Map2
    my $result='';
    my $d1=substr($_[0],-1,1);
    my $d2=substr($_[4],-1,1);
    my $p1=0; ### =11/12/13/22/23/33, means reads1 has start within region 1/2/3 and end within region 1/2/3, 1 upflank, 2 repeat, 3 downflank
    my $p2=0; ### =11/12/13/22/23/33, means reads2 has start within region 1/2/3 and end within region 1/2/3, 1 upflank, 2 repeat, 3 downflank
    my $i=0;
    my $improperlist='12+12+ 12+13+ 13+12+ 13+13+ 13+22+ 13+23+ 13+12- 13+22- 22+13+ 22+13- 23+13+ 23+23+ 23+13- 12-13+ 12-12- 12-13- 13-22+ 13-23+ 13-12- 13-13- 13-22- 13-23- 22-13+ 22-13- 23-13- 23-23-';
    my $checklist='12+12- 12+13- 13+13- 13+23- 23+23- 12-12+ 13-12+ 13-13+ 23-13+ 23-23+';
    my $temp='';
    for ($i=1;$i<=3;$i++) {
        if ($_[$i]>0) {
            $p1=$i*10;
            last;
        }
    }
    for ($i=3;$i>=1;$i--) {
        if ($_[$i]>0) {
            $p1=$p1+$i;
            last;
        }
    }
    for ($i=5;$i<=7;$i++) {
        if ($_[$i]>0) {
            $p2=($i-4)*10;
            last;
        }
    }
    for ($i=7;$i>=5;$i--) {
        if ($_[$i]>0) {
            $p2=$p2+$i-4;
            last;
        }
    }

    $temp=join '',($p1,$d1,$p2,$d2);
    if ($improperlist=~/\Q$temp\E/) { ### following situations means improper pairtype (doesn't mean other situations are definitely good)
        return($result);
    } elsif ($checklist=~/\Q$temp\E/) { ### following situations can be further checked
        if (($d1 eq '+')and(($_[1]<$_[5])or($_[3]>$_[7]))) {
            return($result);
        } elsif (($d2 eq '+')and(($_[1]>$_[5])or($_[3]<$_[7]))) {
            return($result);
        }
    }

    if ($p1<$p2) {
        if (($d1 eq '+')and($d2 eq '-')) {
            $result='fr';
        } elsif (($d1 eq '-')and($d2 eq '+')) {
            $result='rf';
        } elsif (($d1 eq '+')and($d2 eq '+')) {
            $result='ff';
        } elsif (($d1 eq '-')and($d2 eq '-')) {
            $result='rr';
        }
    } elsif ($p1>$p2) {
        if (($d1 eq '+')and($d2 eq '-')) {
            $result='rf';
        } elsif (($d1 eq '-')and($d2 eq '+')) {
            $result='fr';
        } elsif (($d1 eq '+')and($d2 eq '+')) {
            $result='rr';
        } elsif (($d1 eq '-')and($d2 eq '-')) {
            $result='ff';
        }
    } else {
        if ($d1 eq $d2) {
            $result='ff rr';
        } elsif (22==$p1) {
            $result='fr rf';
        } else {
            $result='fr';
        }
    }
    return($result);
}

sub getSkip { ### may slow down mapping-guide speed, sofar used mainly for repeat-guide part
    ### $_[0]=1/2 means $Reads1/$Reads2, $_[1] repeat unit to generate @Skip5/3
    ### $_[2] [3] [4]=match/mis/gap score, $_[5] length of $Reads1/$Reads2
    ### $skip_minfree $skip_misgap $skip_mislimit were prepared in &PrepareSTR
    my $temp=$Reads1;
    if (2==$_[0]) {
        $temp=$Reads2;
    }
    $maxskipindex++;
    $SkipIndex{"$_[0] $_[1]"}=$maxskipindex;
    my @SumScore5=();
    my @SumScore3=();
    my $i=0;
    my $donelength=0;
    my $minpos=0;
    my $maxpos=0;
    my $len=length($_[1]);
    my $a=0;
    my $b=0;
    my $tempstr='';
    my $endstr='';
    my $availfree=0;
    if ($temp=~/(\Q$_[1]\E)+/) {
        @SingleScore=();
        $temp=$';
        $donelength=length($`)+length($&);
        $endstr=$`;
        if ($endstr ne '') { ### try extend 5' end
            for ($i=$len-1;$i>=1;$i--) {
                $tempstr=substr($_[1],$len-$i,$i);
                if ($endstr=~/\Q$tempstr\E$/) {
                    last;
                }
            }
        }
        $minpos=length($endstr)+1-$i;
        for ($i=$minpos;$i<=$donelength;$i++) {
            $SingleScore[$i]=$_[2];
        }
        while ($temp=~/(\Q$_[1]\E)+/) {
            $a=$donelength+length($`)+1;
            $b=$a-1+length($&);
            for ($i=$a;$i<=$b;$i++) {
                $SingleScore[$i]=$_[2];
            }
            $a--;
            if (1==$len) {
                for ($i=$donelength+1;$i<=$a;$i++) {
                    $SingleScore[$i]=$skip_misgap;
                }
            } else {
                &skip_extension($`,$_[1],$_[2],$_[3],$_[4],$donelength+1);
            }
            $temp=$';
            $donelength=$donelength+length($`)+length($&);
        }
        if ($temp ne '') { ### try extend 3' end
            for ($i=$len-1;$i>=1;$i--) {
                $tempstr=substr($_[1],0,$i);
                if ($temp=~/^\Q$tempstr\E/) {
                    last;
                }
            }
        } else {
            $i=0;
        }
        $maxpos=$donelength+$i;
        for ($i=$donelength+1;$i<=$maxpos;$i++) {
            $SingleScore[$i]=$_[2];
        }
        $SumScore3[$minpos-1]=0;
        $SumScore5[$maxpos+1]=0;
        for ($i=$minpos;$i<=$maxpos;$i++) {
            $SumScore3[$i]=$SumScore3[$i-1]+$SingleScore[$i];
            if ($SumScore3[$i]>0) {
                $SumScore3[$i]=0;
            }
        }
        for ($i=$maxpos;$i>=$minpos;$i--) {
            $SumScore5[$i]=$SumScore5[$i+1]+$SingleScore[$i];
            if ($SumScore5[$i]>0) {
                $SumScore5[$i]=0;
            }
        }
        $availfree=0;
        $a=$minpos+$len;
        if ($a<($_[5]-$search_length+1)) {
            $a=$_[5]-$search_length+1;
        }
        for ($i=$maxpos;$i>=$a;$i--) {
            if (0==$SumScore3[$i]) {
                $availfree+=$SumScore3[$i-1]+$_[2];
            } elsif ($SumScore3[$i]<=($skip_mislimit+$_[2])) {
                $availfree=0;
            }
            if ($availfree>=($skip_minfree+$_[2])) {
                $Skip3[$maxskipindex][$_[5]-$i+1]=1;
            }
        }
        $availfree=0;
        $a=$maxpos-$len;
        if ($a>$search_length) {
            $a=$search_length;
        }
        for ($i=$minpos;$i<=$a;$i++) {
            if (0==$SumScore5[$i]) {
                $availfree+=$SumScore5[$i+1]+$_[2];
            } elsif ($SumScore5[$i]<=($skip_mislimit+$_[2])) {
                $availfree=0;
            }
            if ($availfree>=($skip_minfree+$_[2])) {
                $Skip5[$maxskipindex][$i]=1;
            }
        }
    }
}

sub getSeed { ### reset and generate %FlankSeed5/3 (value used in &flank_search to replace start)
    ### [0]=$Reads1/2 [1]=length of $Reads1/2
    ### [2][3] start/end to search seed towards 5'
    ### [4][5] start/end to search seed towards 3'
    my %oriseed=();
    my $tempseed='';
    my $i=0;
    my $half1='';
    my $half2='';
    my $seed='';
    my $newseed='';
    %FlankSeed5=();
    %FlankSeed3=();
    %FlankSeed5fulist=();
    %FlankSeed3fulist=();
    my $searchstart=$_[4];
    my $searchend=$_[5];
    if ($searchend<$flankseed) {
        $searchend=$flankseed;
    }
    for ($i=$searchend;$i<=$searchstart;$i++) {
        $tempseed=substr($_[0],$_[1]-$i,$flankseed);
        $oriseed{$tempseed}=$i;
    }
    foreach $seed (keys(%oriseed)) { ### generate %FlankSeed3
        for ($i=1;$i<=$flankseed;$i++) { ### transform original seed by 1nt deletion
            $half1=substr($seed,0,$i-1);
            $half2=substr($seed,$i,$flankseed-$i);
            $newseed=join '',($half1,$half2);
            if (!(defined($FlankSeed3{$newseed}))) {
                $FlankSeed3{$newseed}=$oriseed{$seed};
                $FlankSeed3fulist{$newseed}=$oriseed{$seed};
            } else {
                $FlankSeed3fulist{$newseed}.=" $oriseed{$seed}";
                if ($oriseed{$seed}>$FlankSeed3{$newseed}) {
                    $FlankSeed3{$newseed}=$oriseed{$seed};
                }
            }
        }
        for ($i=1;$i<=$flankseed;$i++) { ### transform original seed by 1nt insertion (no need to process 3' 1-nt insertion)
            $half1=substr($seed,0,$i-1);
            $half2=substr($seed,$i-1,$flankseed-$i+1);
            foreach ('A','T','C','G','N') {
                $newseed=join '',($half1,$_,$half2);
                if (!(defined($FlankSeed3{$newseed}))) {
                    $FlankSeed3{$newseed}=$oriseed{$seed};
                    $FlankSeed3fulist{$newseed}=$oriseed{$seed};
                } else {
                    $FlankSeed3fulist{$newseed}.=" $oriseed{$seed}";
                    if ($oriseed{$seed}>$FlankSeed3{$newseed}) {
                        $FlankSeed3{$newseed}=$oriseed{$seed};
                    }
                }
            }
        }
        for ($i=1;$i<=$flankseed;$i++) { ### transform original seed by 1nt mis (including perfect match)
            $half1=substr($seed,0,$i-1);
            $half2=substr($seed,$i,$flankseed-$i);
            foreach ('A','T','C','G','N') {
                $newseed=join '',($half1,$_,$half2);
                if (!(defined($FlankSeed3{$newseed}))) {
                    $FlankSeed3{$newseed}=$oriseed{$seed};
                    $FlankSeed3fulist{$newseed}=$oriseed{$seed};
                } else {
                    $FlankSeed3fulist{$newseed}.=" $oriseed{$seed}";
                    if ($oriseed{$seed}>$FlankSeed3{$newseed}) {
                        $FlankSeed3{$newseed}=$oriseed{$seed};
                    }
                }
            }
        }
    }
    %oriseed=();
    $searchstart=$_[2];
    $searchend=$_[3];
    if ($searchend<$flankseed) {
        $searchend=$flankseed;
    }
    for ($i=$searchstart;$i>=$searchend;$i--) {
        $tempseed=substr($_[0],$i-$flankseed,$flankseed);
        $tempseed=reverse($tempseed);
        if (!(defined($oriseed{$tempseed}))) {
            $oriseed{$tempseed}=$i;
        }
    }
    foreach $seed (keys(%oriseed)) { ### generate %FlankSeed5
        for ($i=1;$i<=$flankseed;$i++) { ### transform original seed by 1nt deletion
            $half1=substr($seed,0,$i-1);
            $half2=substr($seed,$i,$flankseed-$i);
            $newseed=join '',($half1,$half2);
            if (!(defined($FlankSeed5{$newseed}))) {
                $FlankSeed5{$newseed}=$oriseed{$seed};
                $FlankSeed5fulist{$newseed}=$oriseed{$seed};
            } else {
                $FlankSeed5fulist{$newseed}.=" $oriseed{$seed}";
                if ($oriseed{$seed}>$FlankSeed5{$newseed}) {
                    $FlankSeed5{$newseed}=$oriseed{$seed};
                }
            }
        }
        for ($i=1;$i<=$flankseed;$i++) { ### transform original seed by 1nt insertion (no need to process 3' 1-nt insertion)
            $half1=substr($seed,0,$i-1);
            $half2=substr($seed,$i-1,$flankseed-$i+1);
            foreach ('A','T','C','G','N') {
                $newseed=join '',($half1,$_,$half2);
                if (!(defined($FlankSeed5{$newseed}))) {
                    $FlankSeed5{$newseed}=$oriseed{$seed};
                    $FlankSeed5fulist{$newseed}=$oriseed{$seed};
                } else {
                    $FlankSeed5fulist{$newseed}.=" $oriseed{$seed}";
                    if ($oriseed{$seed}>$FlankSeed5{$newseed}) {
                        $FlankSeed5{$newseed}=$oriseed{$seed};
                    }
                }
            }
        }
        for ($i=1;$i<=$flankseed;$i++) { ### transform original seed by 1nt mis (including perfect match)
            $half1=substr($seed,0,$i-1);
            $half2=substr($seed,$i,$flankseed-$i);
            foreach ('A','T','C','G','N') {
                $newseed=join '',($half1,$_,$half2);
                if (!(defined($FlankSeed5{$newseed}))) {
                    $FlankSeed5{$newseed}=$oriseed{$seed};
                    $FlankSeed5fulist{$newseed}=$oriseed{$seed};
                } else {
                    $FlankSeed5fulist{$newseed}.=" $oriseed{$seed}";
                    if ($oriseed{$seed}>$FlankSeed5{$newseed}) {
                        $FlankSeed5{$newseed}=$oriseed{$seed};
                    }
                }
            }
        }
    }
}

sub ReAligner {
    ### process @Lines1/2, which store mapping of a single pair of reads ($Reads1/2) from .sam
    ### $_[0] ID of the pair
    ### generate raw realignments %Ful1/2, store those contain/bridge STR into @RawResult by paired "+/- $_[0] 1/2" & "realignment"
    $RLen1=length($Reads1);
    $RLen2=length($Reads2);
    &getRawRepeat;

    my %bwa2=(); ### index="STRid+/-", judge if bwa alignments for Reads2 exists, only to filter out Reads1 bwa alignments unnecessary to process
    my $unit='';
    my $dymatch_score=0;
    @Skip5=(); ### 2D array, [index by "1/2 repeat"][i]=1 means flanksearch of Reads"1/2" given "repeat" towards 5' can skip flank length "i"
    @Skip3=(); ### 2D array, [index by "1/2 repeat"][i]=1 means flanksearch of Reads"1/2" given "repeat" towards 3' can skip flank length "i"
    %SkipIndex=(); ### index="1/2 repeat", repeat same to the index of %SameUnit, value (start from 1) is used by @Skip5/3 first []
    $maxskipindex=0; ### max value of %SkipIndex
    foreach $unit (keys(%Repeattogo1)) { ### prepare @Skip5/3 for $Reads1
        $dymatch_score=$dymscore[length($unit)];
        &getSkip(1,$unit,$dymatch_score,$mis_score,$gap_score,$RLen1);
    }
    foreach $unit (keys(%Repeattogo2)) { ### prepare @Skip5/3 for $Reads2
        $dymatch_score=$dymscore[length($unit)];
        &getSkip(2,$unit,$dymatch_score,$mis_score,$gap_score,$RLen2);
    }

    my $i=0;
    for ($i=$#Lines2-4;$i>=0;$i=$i-5) { ### generate %bwa2
        $bwa2{"$Lines2[$i+1]$Lines2[$i]"}=1;
    }

    my @Result1=(); ### 2D array for realignment results for $Reads1
    my @Result2=(); ### 2D array for realignment results for $Reads2
    ### [X][0]=STRid, [X][1]=+/-, [X][2..]="X1 X2 X3 limit X1 X2 X3 limit ..." (X1/2/3 can =0)
    ### means reads mapped to STRid+/-, X1 nt mapped to upstream flank, X2 nt repeats, X3 downstream flank
    ### X1/X2/X3=-1 means not enough information to be sure
    ### X2 can=X.1 if not fully extended due to reads length
    ### X1/3 can =X.9 if with lower priority judged by length $checklimit, but still considered as X+1 when calculating X2
    ### limit is used in 3rd round to speed up flank search
    my %Index1=(); ### for @Result1, $Index1{"STRid +/-"}=X (minX=1)
    my %Index2=(); ### for @Result2, $Index2{"STRid +/-"}=X (minX=1)
    my $IndexNum1=0; ### for @Result1, max X
    my $IndexNum2=0; ### for @Result2, max X
    my @Tail1=(); ### for @Result1, $Tail1[X]=max of [0..]
    my @Tail2=(); ### for @Result2, $Tail2[X]=max of [0..]
    my $pairref='';
    my $pairunit='';
    my $index=0;
    my @Raw=();
    my $j=0;
    my @tempdata=();
    for ($i=$#Lines1-4;$i>=0;$i=$i-5) { ### 1st round, get @Result1 X1/X3 from .sam/.bam
        @tempdata=($Lines1[$i],$Lines1[$i+1],$Lines1[$i+2],$Lines1[$i+3],$Lines1[$i+4]);
        if (($pairtype eq 'ff')or($pairtype eq 'rr')) {
            $pairref=join '',($tempdata[1],$tempdata[0]);
            if ($tempdata[0] eq '-') {
                $pairunit=$StrUnitcomple{$tempdata[1]};
            } else {
                $pairunit=$StrUnit{$tempdata[1]};
            }
        } elsif (($pairtype eq 'fr')or($pairtype eq 'rf')) {
            if ($tempdata[0] eq '+') {
                $pairref=join '',($tempdata[1],'-');
                $pairunit=$StrUnitcomple{$tempdata[1]};
            } else {
                $pairref=join '',($tempdata[1],'+');
                $pairunit=$StrUnit{$tempdata[1]};
            }
        }
        if ((!defined($bwa2{$pairref}))and(!defined($STRtogo2{$pairref}))and(!defined($Repeattogo2{$pairunit}))) { ### skip if impossible to be paired with $Reads2
            next;
        }
        @Raw=&getRaw($tempdata[0],$tempdata[1],$tempdata[2],$tempdata[3],$tempdata[4],$Reads1,$RLen1);
        if ($#Raw>=0) {
            if (!defined($Index1{"$tempdata[1] $tempdata[0]"})) {
                $IndexNum1++;
                $Index1{"$tempdata[1] $tempdata[0]"}=$IndexNum1;
                $index=$IndexNum1;
                $Result1[$index][0]=$tempdata[1];
                $Result1[$index][1]=$tempdata[0];
                $Tail1[$index]=1;
            } else {
                $index=$Index1{"$tempdata[1] $tempdata[0]"};
            }
            for ($j=0;$j<=$#Raw;$j=$j+3) {
                $Result1[$index][$Tail1[$index]+1]=$Raw[$j];
                $Result1[$index][$Tail1[$index]+2]=$Raw[$j+1];
                $Result1[$index][$Tail1[$index]+3]=$Raw[$j+2];
                $Result1[$index][$Tail1[$index]+4]=0;
                $Tail1[$index]+=4;
            }
        }
    }

    my $upext=0;
    my $downext=0;
    my $uplimit=0;
    my $downlimit=0;
    my $pos=0;
    my $tempX1=0;
    my $tempX2=0;
    my $tempX3=0;
    my $tempX4=0;
    my $strand='';
    my $ref='';
    for ($i=$IndexNum1;$i>=1;$i--) { ### 2nd round 1st half, prepare %Repeat1
        $strand=$Result1[$i][1];
        if ($strand eq '-') {
            $unit=$StrUnitcomple{$Result1[$i][0]};
        } else {
            $unit=$StrUnit{$Result1[$i][0]};
        }
        $dymatch_score=$dymscore[length($unit)];
        for ($j=$Tail1[$i]-3;$j>=2;$j=$j-4) {
            if (-1==$Result1[$i][$j+1]) { ### an X2 to be finished
                $tempX1=$Result1[$i][$j];
                $tempX2=$Result1[$i][$j+1];
                $tempX3=$Result1[$i][$j+2];
                if ($strand eq '-') {
                    if ($tempX1>0) {
                        $pos=$RLen1-$tempX1;
                        if (!(defined($Repeat1{$unit}))) { ### add a new [start end] to $Repeat1{$unit}
                            $upext=$pos+1-&up_extension($Reads1,$pos+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            $downext=$pos+&down_extension($Reads1,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            if ($upext<=$downext) {
                                $Repeat1{$unit}=join ' ',($upext,$downext,$RLen1,1);
                            }
                            #&getSkip(1,$unit,$dymatch_score,$mis_score,$gap_score,$RLen1);
                        } else {
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1{$unit});
                            if (0==$upext) { ### add a new [start end] to $Repeat1{$unit}, then check overlap
                                $upext=$pos+1-&up_extension($Reads1,$pos+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                $downext=$pos+&down_extension($Reads1,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                if ($upext<=$downext) {
                                    $Repeat1{$unit}=&AddRepeat($upext,$downext,$RLen1,1,$Repeat1{$unit});
                                }
                            }
                        }
                    } elsif ($tempX3>0) {
                        $pos=$tempX3+1;
                        if (!(defined($Repeat1{$unit}))) { ### add a new [start end] to $Repeat1{$unit}
                            $upext=$pos-&up_extension($Reads1,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            $downext=$pos-1+&down_extension($Reads1,$pos-1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            if ($upext<=$downext) {
                                $Repeat1{$unit}=join ' ',($upext,$downext,$RLen1,1);
                            }
                            #&getSkip(1,$unit,$dymatch_score,$mis_score,$gap_score,$RLen1);
                        } else {
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1{$unit});
                            if (0==$upext) { ### add a new [start end] to $Repeat1{$unit}, then check overlap
                                $upext=$pos-&up_extension($Reads1,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                $downext=$pos-1+&down_extension($Reads1,$pos-1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                if ($upext<=$downext) {
                                    $Repeat1{$unit}=&AddRepeat($upext,$downext,$RLen1,1,$Repeat1{$unit});
                                }
                            }
                        }
                    }
                } else {
                    if ($tempX1>0) {
                        $pos=$tempX1+1;
                        if (!(defined($Repeat1{$unit}))) { ### add a new [start end] to $Repeat1{$unit}
                            $upext=$pos-&up_extension($Reads1,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            $downext=$pos-1+&down_extension($Reads1,$pos-1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            if ($upext<=$downext) {
                                $Repeat1{$unit}=join ' ',($upext,$downext,$RLen1,1);
                            }
                            #&getSkip(1,$unit,$dymatch_score,$mis_score,$gap_score,$RLen1);
                        } else {
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1{$unit});
                            if (0==$upext) { ### add a new [start end] to $Repeat1{$unit}, then check overlap
                                $upext=$pos-&up_extension($Reads1,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                $downext=$pos-1+&down_extension($Reads1,$pos-1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                if ($upext<=$downext) {
                                    $Repeat1{$unit}=&AddRepeat($upext,$downext,$RLen1,1,$Repeat1{$unit});
                                }
                            }
                        }
                    } elsif ($tempX3>0) {
                        $pos=$RLen1-$tempX3;
                        if (!(defined($Repeat1{$unit}))) { ### add a new [start end] to $Repeat1{$unit}
                            $upext=$pos+1-&up_extension($Reads1,$pos+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            $downext=$pos+&down_extension($Reads1,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            if ($upext<=$downext) {
                                $Repeat1{$unit}=join ' ',($upext,$downext,$RLen1,1);
                            }
                            #&getSkip(1,$unit,$dymatch_score,$mis_score,$gap_score,$RLen1);
                        } else {
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1{$unit});
                            if (0==$upext) { ### add a new [start end] to $Repeat1{$unit}, then check overlap
                                $upext=$pos+1-&up_extension($Reads1,$pos+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                $downext=$pos+&down_extension($Reads1,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                if ($upext<=$downext) {
                                    $Repeat1{$unit}=&AddRepeat($upext,$downext,$RLen1,1,$Repeat1{$unit});
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for ($i=$IndexNum1;$i>=1;$i--) { ### 2nd round 2nd half, get @Result1 X2
        $strand=$Result1[$i][1];
        $ref=join '',($Result1[$i][0],$strand);
        if ($strand eq '-') {
            $unit=$StrUnitcomple{$Result1[$i][0]};
        } else {
            $unit=$StrUnit{$Result1[$i][0]};
        }
        for ($j=$Tail1[$i]-3;$j>=2;$j=$j-4) {
            if (-1==$Result1[$i][$j+1]) { ### an X2 to be finished
                $tempX1=$Result1[$i][$j];
                $tempX2=$Result1[$i][$j+1];
                $tempX3=$Result1[$i][$j+2];
                if ($strand eq '-') {
                    if ($tempX1>0) {
                        $pos=$RLen1-$tempX1;
                        ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1{$unit});
                        if (0==$upext) { ### should not happen because of $allowance
                            print ("Warning: cannot find repeats of $unit covering $pos for $Reads1\n");
                            next;
                        }
                        $tempX2=$RLen1-$upext-$tempX1+1;
                        if (($RLen1-$tempX1-int($tempX2))>$MaxGap{$Result1[$i][0]}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                            $tempX2=int($tempX2);
                        }
                        $tempX4=$uplimit;
                        if (defined($Repeat1STR{$ref})) {
                            $pos=$RLen1-$tempX1-int($tempX2)+1;
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1STR{$ref});
                            if (($upext>0)and($uplimit<$tempX4)) {
                                $tempX4=$uplimit;
                            }
                        }
                    } elsif ($tempX3>0) {
                        $pos=$tempX3+1;
                        ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1{$unit});
                        if (0==$downext) { ### should not happen because of $allowance
                            print ("Warning: cannot find repeats of $unit covering $pos for $Reads1\n");
                            next;
                        }
                        $tempX2=$downext-$tempX3;
                        if (($RLen1-$tempX3-int($tempX2))>$MaxGap{$Result1[$i][0]}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                            $tempX2=int($tempX2);
                        }
                        $tempX4=$RLen1-$downlimit+1;
                        if (defined($Repeat1STR{$ref})) {
                            $pos=$tempX3+int($tempX2);
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1STR{$ref});
                            if (($upext>0)and(($RLen1-$downlimit+1)<$tempX4)) {
                                $tempX4=$RLen1-$downlimit+1;
                            }
                        }
                    }
                } else {
                    if ($tempX1>0) {
                        $pos=$tempX1+1;
                        ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1{$unit});
                        if (0==$downext) { ### should not happen because of $allowance
                            print ("Warning: cannot find repeats of $unit covering $pos for $Reads1\n");
                            next;
                        }
                        $tempX2=$downext-$tempX1;
                        if (($RLen1-$tempX1-int($tempX2))>$MaxGap{$Result1[$i][0]}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                            $tempX2=int($tempX2);
                        }
                        $tempX4=$RLen1-$downlimit+1;
                        if (defined($Repeat1STR{$ref})) {
                            $pos=$tempX1+int($tempX2);
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1STR{$ref});
                            if (($upext>0)and(($RLen1-$downlimit+1)<$tempX4)) {
                                $tempX4=$RLen1-$downlimit+1;
                            }
                        }
                    } elsif ($tempX3>0) {
                        $pos=$RLen1-$tempX3;
                        ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1{$unit});
                        if (0==$upext) { ### should not happen because of $allowance
                            print ("Warning: cannot find repeats of $unit covering $pos for $Reads1\n");
                            next;
                        }
                        $tempX2=$RLen1-$upext-$tempX3+1;
                        if (($RLen1-$tempX3-int($tempX2))>$MaxGap{$Result1[$i][0]}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                            $tempX2=int($tempX2);
                        }
                        $tempX4=$uplimit;
                        if (defined($Repeat1STR{$ref})) {
                            $pos=$RLen1-$tempX3-int($tempX2)+1;
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1STR{$ref});
                            if (($upext>0)and($uplimit<$tempX4)) {
                                $tempX4=$uplimit;
                            }
                        }
                    }
                }
                $Result1[$i][$j+1]=$tempX2;
                $Result1[$i][$j+3]=max($tempX4,1);
            }
        }
    }

    my %Ful1=(); ### $Ful1{STRid+/-} will store $Reads1 info processed from @Result1, by finishing all X1/3s
    my %Ful2=(); ### $Ful2{STRid+/-} will store $Reads2 info processed from @Result2, by finishing all X1/3s
    my %Ful1STR=(); ### $Ful1{STRid+/-} contains/bridges STR
    my %Ful2STR=(); ### $Ful2{STRid+/-} contains/bridges STR
    %FlankScore=(); ### index="1/2 STRid+/- 5/3 sublength", value=alignment score by &flank_extension when Reads1/2 mapped to STRid+/-, using sublength at 5'/3' end for correlated flank
    my %FlankAlign=(); ### index="1/2 STRid+/- 5/3", value="sublength1 sublength2 .." sublengths at 5'/3' end for correlated flank aligned by bwa when Reads1/2 mapped to STRid+/-
    my $findflank=0;
    my $findflankraw=0;
    my @tempflank=();
    my $start=0;
    my $end=0;
    my $direct=0;
    my $reliab=0;
    my $k=0;
    my $X2reach=0;
    for ($i=$IndexNum1;$i>=1;$i--) { ### 3rd round 1st half, get %FlankAlign from bwa alignment
        $strand=$Result1[$i][1];
        $ref=join '',($Result1[$i][0],$strand);
        for ($j=$Tail1[$i]-3;$j>=2;$j=$j-4) {
            $tempX1=$Result1[$i][$j];
            $tempX2=$Result1[$i][$j+1];
            $tempX3=$Result1[$i][$j+2];
            if ($strand eq '-') {
                if ($tempX1>0) {
                    if (defined($FlankAlign{"1 $ref 3"})) {
                        $FlankAlign{"1 $ref 3"}.=" $tempX1";
                    } else {
                        $FlankAlign{"1 $ref 3"}=$tempX1;
                    }
                } elsif ($tempX3>0) {
                    if (defined($FlankAlign{"1 $ref 5"})) {
                        $FlankAlign{"1 $ref 5"}.=" $tempX3";
                    } else {
                        $FlankAlign{"1 $ref 5"}=$tempX3;
                    }
                }
            } else {
                if ($tempX1>0) {
                    if (defined($FlankAlign{"1 $ref 5"})) {
                        $FlankAlign{"1 $ref 5"}.=" $tempX1";
                    } else {
                        $FlankAlign{"1 $ref 5"}=$tempX1;
                    }
                } elsif ($tempX3>0) {
                    if (defined($FlankAlign{"1 $ref 3"})) {
                        $FlankAlign{"1 $ref 3"}.=" $tempX3";
                    } else {
                        $FlankAlign{"1 $ref 3"}=$tempX3;
                    }
                }
            }
        }
    }
    for ($i=$IndexNum1;$i>=1;$i--) { ### 3rd round 2nd half, get $Ful1{STRid+/-} by processing @Result1
        $strand=$Result1[$i][1];
        $ref=join '',($Result1[$i][0],$strand);
        for ($j=$Tail1[$i]-3;$j>=2;$j=$j-4) {
            $tempX1=$Result1[$i][$j];
            $tempX2=$Result1[$i][$j+1];
            $tempX3=$Result1[$i][$j+2];
            $findflank=0;
            if (-1==$tempX1) {
                if ($strand eq '-') {
                    $direct=3;
                } else {
                    $direct=5;
                }
                if (defined($FlankAlign{"1 $ref $direct"})) { ### search available X1 in %FlankAlign
                    @tempflank=split / /, $FlankAlign{"1 $ref $direct"};
                    for ($k=$#tempflank;$k>=0;$k--) {
                        if (($tempflank[$k]>$findflank)and($tempflank[$k]<=($RLen1-$tempX3))and($tempflank[$k]>=($RLen1-$tempX3-int($tempX2)))) {
                            $findflank=$tempflank[$k];
                        }
                    }
                }
                if ($findflank!=0) { ### found good X1 from %FlankAlign
                    $tempX1=$findflank;
                    $tempX2=$RLen1-$tempX1-$tempX3;
                } elsif (($RLen1-$tempX3-int($tempX2))<=$search_length) { ### if still possible, continue to obtain X1 by &flank_search
                    $end=$RLen1-$tempX3-int($tempX2);
                    if (0==$end) {
                        $end=1;
                    }
                    $start=min($search_length,$end+$UpFL{$Result1[$i][0]},$RLen1-$tempX3,$Result1[$i][$j+3]);
                    if (((int($tempX2)+$tempX3)==$RLen1)or($tempX2=~/\./)) { ### X2 replacement exists
                        $X2reach=1;
                    } else {
                        $X2reach=0;
                    }
                    while ($start>=$end) {
                        @tempflank=&flank_search(1,$ref,$direct,$start,$end,'no');
                        if (0==$#tempflank) { ### X1 not found, end search
                            last;
                        }
                        for ($k=0;$k<=$#tempflank-1;$k++) {
                            $findflankraw=$tempflank[$k];
                            $findflank=int($findflankraw+0.5);
                            if ($findflank>$flanklimit) { ### found good X1 by &flank_search, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX1=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX1=$findflank-0.1;
                                    }
                                } else {
                                    $tempX1=$findflank;
                                }
                                $tempX2=$RLen1-$findflank-$tempX3;
                                last;
                            } else { ### found not very good X1 by &flank_search
                                $reliab=&reliab_check(1,$ref,$direct,$findflank,$RLen1-$findflank-$tempX3,$tempX3);
                                if ($reliab>=$reliablimit) { ### X1 is trustful, end search
                                    if (1==$X2reach) { ### X2 replacement exists
                                        if ($findflank>$checklimit) {
                                            $tempX1=$findflankraw;
                                        } else { ### set a mark for further check in &PrepareResult
                                            $tempX1=$findflank-0.1;
                                        }
                                    } else { ### X2 cannot reach end, skip marking
                                        $tempX1=$findflank;
                                    }
                                    $tempX2=$RLen1-$findflank-$tempX3;
                                    last;
                                } else { ### X1 is not trustful, continue to search
                                    $findflank=0;
                                }
                            }
                        }
                        if ($findflank!=0) {
                            last;
                        } else {
                            $start=$tempflank[$#tempflank]-1;
                        }
                    }
                    if ((0==$findflank)and(1==$X2reach)) { ### X1 not found, but X2 reaches end or was incomplete due to reads length
                        $tempX1=0;
                        $tempX2=$RLen1-$tempX3;
                    }
                }
            } elsif (-1==$tempX3) {
                if ($strand eq '-') {
                    $direct=5;
                } else {
                    $direct=3;
                }
                if (defined($FlankAlign{"1 $ref $direct"})) { ### search available X3 in %FlankAlign
                    @tempflank=split / /, $FlankAlign{"1 $ref $direct"};
                    for ($k=$#tempflank;$k>=0;$k--) {
                        if (($tempflank[$k]>$findflank)and($tempflank[$k]<=($RLen1-$tempX1))and($tempflank[$k]>=($RLen1-$tempX1-int($tempX2)))) {
                            $findflank=$tempflank[$k];
                        }
                    }
                }
                if ($findflank!=0) { ### found good X3 from %FlankAlign
                    $tempX3=$findflank;
                    $tempX2=$RLen1-$tempX1-$tempX3;
                } elsif (($RLen1-$tempX1-int($tempX2))<=$search_length) { ### if still possible, continue to obtain X3 by &flank_search
                    $end=$RLen1-$tempX1-int($tempX2);
                    if (0==$end) {
                        $end=1;
                    }
                    $start=min($search_length,$end+$DownFL{$Result1[$i][0]},$RLen1-$tempX1,$Result1[$i][$j+3]);
                    if (((int($tempX2)+$tempX1)==$RLen1)or($tempX2=~/\./)) { ### X2 replacement exists
                        $X2reach=1;
                    } else {
                        $X2reach=0;
                    }
                    while ($start>=$end) {
                        @tempflank=&flank_search(1,$ref,$direct,$start,$end,'no');
                        if (0==$#tempflank) { ### X3 not found, end search
                            last;
                        }
                        for ($k=0;$k<=$#tempflank-1;$k++) {
                            $findflankraw=$tempflank[$k];
                            $findflank=int($findflankraw+0.5);
                            if ($findflank>$flanklimit) { ### found good X3 by &flank_search, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX3=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX3=$findflank-0.1;
                                    }
                                } else {
                                    $tempX3=$findflank;
                                }
                                $tempX2=$RLen1-$tempX1-$findflank;
                                last;
                            } else { ### found not very good X3 by &flank_search
                                $reliab=&reliab_check(1,$ref,$direct,$tempX1,$RLen1-$findflank-$tempX1,$findflank);
                                if ($reliab>=$reliablimit) { ### X3 is trustful, end search
                                    if (1==$X2reach) { ### X2 replacement exists
                                        if ($findflank>$checklimit) {
                                            $tempX3=$findflankraw;
                                        } else { ### set a mark for further check in &PrepareResult
                                            $tempX3=$findflank-0.1;
                                        }
                                    } else { ### X2 cannot reach end, skip marking
                                        $tempX3=$findflank;
                                    }
                                    $tempX2=$RLen1-$tempX1-$findflank;
                                    last;
                                } else { ### X3 is not trustful, continue to search
                                    $findflank=0;
                                }
                            }
                        }
                        if ($findflank!=0) {
                            last;
                        } else {
                            $start=$tempflank[$#tempflank]-1;
                        }
                    }
                    if ((0==$findflank)and(1==$X2reach)) { ### X3 not found, but X2 reaches end or was incomplete due to reads length
                        $tempX3=0;
                        $tempX2=$RLen1-$tempX1;
                    }
                }
            }
            if (($tempX1!=-1)and($tempX2!=-1)and($tempX3!=-1)) {
                if (defined($Ful1{$ref})) {
                    $Ful1{$ref}.=" $tempX1 $tempX2 $tempX3";
                } else {
                    $Ful1{$ref}="$tempX1 $tempX2 $tempX3";
                }
                if (($tempX2>0)or(($tempX1>0)and($tempX3>0))) {
                    $Ful1STR{$ref}=1;
                }
            }
        }
    }

    for ($i=$#Lines2-4;$i>=0;$i=$i-5) { ### 1st round, get @Result2 X1/X3 from .sam/.bam
        @tempdata=($Lines2[$i],$Lines2[$i+1],$Lines2[$i+2],$Lines2[$i+3],$Lines2[$i+4]);
        if (($pairtype eq 'ff')or($pairtype eq 'rr')) {
            $pairref=join '',($tempdata[1],$tempdata[0]);
            if ($tempdata[0] eq '-') {
                $pairunit=$StrUnitcomple{$tempdata[1]};
            } else {
                $pairunit=$StrUnit{$tempdata[1]};
            }
        } elsif (($pairtype eq 'fr')or($pairtype eq 'rf')) {
            if ($tempdata[0] eq '+') {
                $pairref=join '',($tempdata[1],'-');
                $pairunit=$StrUnitcomple{$tempdata[1]};
            } else {
                $pairref=join '',($tempdata[1],'+');
                $pairunit=$StrUnit{$tempdata[1]};
            }
        }
        if ((!defined($Ful1{$pairref}))and(!defined($STRtogo1{$pairref}))and(!defined($Repeattogo1{$pairunit}))) { ### skip if impossible to be paired with $Reads1
            next;
        }
        @Raw=&getRaw($tempdata[0],$tempdata[1],$tempdata[2],$tempdata[3],$tempdata[4],$Reads2,$RLen2);
        if ($#Raw>=0) {
            if (!defined($Index2{"$tempdata[1] $tempdata[0]"})) {
                $IndexNum2++;
                $Index2{"$tempdata[1] $tempdata[0]"}=$IndexNum2;
                $index=$IndexNum2;
                $Result2[$index][0]=$tempdata[1];
                $Result2[$index][1]=$tempdata[0];
                $Tail2[$index]=1;
            } else {
                $index=$Index2{"$tempdata[1] $tempdata[0]"};
            }
            for ($j=0;$j<=$#Raw;$j=$j+3) {
                $Result2[$index][$Tail2[$index]+1]=$Raw[$j];
                $Result2[$index][$Tail2[$index]+2]=$Raw[$j+1];
                $Result2[$index][$Tail2[$index]+3]=$Raw[$j+2];
                $Result2[$index][$Tail2[$index]+4]=0;
                $Tail2[$index]+=4;
            }
        }
    }

    for ($i=$IndexNum2;$i>=1;$i--) { ### 2nd round 1st half, prepare %Repeat2
        $strand=$Result2[$i][1];
        if ($strand eq '-') {
            $unit=$StrUnitcomple{$Result2[$i][0]};
        } else {
            $unit=$StrUnit{$Result2[$i][0]};
        }
        $dymatch_score=$dymscore[length($unit)];
        for ($j=$Tail2[$i]-3;$j>=2;$j=$j-4) {
            if (-1==$Result2[$i][$j+1]) { ### an X2 to be finished
                $tempX1=$Result2[$i][$j];
                $tempX2=$Result2[$i][$j+1];
                $tempX3=$Result2[$i][$j+2];
                if ($strand eq '-') {
                    if ($tempX1>0) {
                        $pos=$RLen2-$tempX1;
                        if (!(defined($Repeat2{$unit}))) { ### add a new [start end] to $Repeat2{$unit}
                            $upext=$pos+1-&up_extension($Reads2,$pos+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            $downext=$pos+&down_extension($Reads2,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            if ($upext<=$downext) {
                                $Repeat2{$unit}=join ' ',($upext,$downext,$RLen2,1);
                            }
                            #&getSkip(2,$unit,$dymatch_score,$mis_score,$gap_score,$RLen2);
                        } else {
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2{$unit});
                            if (0==$upext) { ### add a new [start end] to $Repeat2{$unit}, then check overlap
                                $upext=$pos+1-&up_extension($Reads2,$pos+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                $downext=$pos+&down_extension($Reads2,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                if ($upext<=$downext) {
                                    $Repeat2{$unit}=&AddRepeat($upext,$downext,$RLen2,1,$Repeat2{$unit});
                                }
                            }
                        }
                    } elsif ($tempX3>0) {
                        $pos=$tempX3+1;
                        if (!(defined($Repeat2{$unit}))) { ### add a new [start end] to $Repeat2{$unit}
                            $upext=$pos-&up_extension($Reads2,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            $downext=$pos-1+&down_extension($Reads2,$pos-1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            if ($upext<=$downext) {
                                $Repeat2{$unit}=join ' ',($upext,$downext,$RLen2,1);
                            }
                            #&getSkip(2,$unit,$dymatch_score,$mis_score,$gap_score,$RLen2);
                        } else {
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2{$unit});
                            if (0==$upext) { ### add a new [start end] to $Repeat2{$unit}, then check overlap
                                $upext=$pos-&up_extension($Reads2,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                $downext=$pos-1+&down_extension($Reads2,$pos-1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                if ($upext<=$downext) {
                                    $Repeat2{$unit}=&AddRepeat($upext,$downext,$RLen2,1,$Repeat2{$unit});
                                }
                            }
                        }
                    }
                } else {
                    if ($tempX1>0) {
                        $pos=$tempX1+1;
                        if (!(defined($Repeat2{$unit}))) { ### add a new [start end] to $Repeat2{$unit}
                            $upext=$pos-&up_extension($Reads2,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            $downext=$pos-1+&down_extension($Reads2,$pos-1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            if ($upext<=$downext) {
                                $Repeat2{$unit}=join ' ',($upext,$downext,$RLen2,1);
                            }
                            #&getSkip(2,$unit,$dymatch_score,$mis_score,$gap_score,$RLen2);
                        } else {
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2{$unit});
                            if (0==$upext) { ### add a new [start end] to $Repeat2{$unit}, then check overlap
                                $upext=$pos-&up_extension($Reads2,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                $downext=$pos-1+&down_extension($Reads2,$pos-1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                if ($upext<=$downext) {
                                    $Repeat2{$unit}=&AddRepeat($upext,$downext,$RLen2,1,$Repeat2{$unit});
                                }
                            }
                        }
                    } elsif ($tempX3>0) {
                        $pos=$RLen2-$tempX3;
                        if (!(defined($Repeat2{$unit}))) { ### add a new [start end] to $Repeat2{$unit}
                            $upext=$pos+1-&up_extension($Reads2,$pos+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            $downext=$pos+&down_extension($Reads2,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                            if ($upext<=$downext) {
                                $Repeat2{$unit}=join ' ',($upext,$downext,$RLen2,1);
                            }
                            #&getSkip(2,$unit,$dymatch_score,$mis_score,$gap_score,$RLen2);
                        } else {
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2{$unit});
                            if (0==$upext) { ### add a new [start end] to $Repeat2{$unit}, then check overlap
                                $upext=$pos+1-&up_extension($Reads2,$pos+1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                $downext=$pos+&down_extension($Reads2,$pos,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
                                if ($upext<=$downext) {
                                    $Repeat2{$unit}=&AddRepeat($upext,$downext,$RLen2,1,$Repeat2{$unit});
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for ($i=$IndexNum2;$i>=1;$i--) { ### 2nd round 2nd half, get @Result2 X2
        $strand=$Result2[$i][1];
        $ref=join '',($Result2[$i][0],$strand);
        if ($strand eq '-') {
            $unit=$StrUnitcomple{$Result2[$i][0]};
        } else {
            $unit=$StrUnit{$Result2[$i][0]};
        }
        for ($j=$Tail2[$i]-3;$j>=2;$j=$j-4) {
            if (-1==$Result2[$i][$j+1]) { ### an X2 to be finished
                $tempX1=$Result2[$i][$j];
                $tempX2=$Result2[$i][$j+1];
                $tempX3=$Result2[$i][$j+2];
                if ($strand eq '-') {
                    if ($tempX1>0) {
                        $pos=$RLen2-$tempX1;
                        ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2{$unit});
                        if (0==$upext) { ### should not happen because of $allowance
                            print ("Warning: cannot find repeats of $unit covering $pos for $Reads2\n");
                            next;
                        }
                        $tempX2=$RLen2-$upext-$tempX1+1;
                        if (($RLen2-$tempX1-int($tempX2))>$MaxGap{$Result2[$i][0]}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                            $tempX2=int($tempX2);
                        }
                        $tempX4=$uplimit;
                        if (defined($Repeat2STR{$ref})) {
                            $pos=$RLen2-$tempX1-int($tempX2)+1;
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2STR{$ref});
                            if (($upext>0)and($uplimit<$tempX4)) {
                                $tempX4=$uplimit;
                            }
                        }
                    } elsif ($tempX3>0) {
                        $pos=$tempX3+1;
                        ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2{$unit});
                        if (0==$downext) { ### should not happen because of $allowance
                            print ("Warning: cannot find repeats of $unit covering $pos for $Reads2\n");
                            next;
                        }
                        $tempX2=$downext-$tempX3;
                        if (($RLen2-$tempX3-int($tempX2))>$MaxGap{$Result2[$i][0]}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                            $tempX2=int($tempX2);
                        }
                        $tempX4=$RLen2-$downlimit+1;
                        if (defined($Repeat2STR{$ref})) {
                            $pos=$tempX3+int($tempX2);
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2STR{$ref});
                            if (($upext>0)and(($RLen2-$downlimit+1)<$tempX4)) {
                                $tempX4=$RLen2-$downlimit+1;
                            }
                        }
                    }
                } else {
                    if ($tempX1>0) {
                        $pos=$tempX1+1;
                        ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2{$unit});
                        if (0==$downext) { ### should not happen because of $allowance
                            print ("Warning: cannot find repeats of $unit covering $pos for $Reads2\n");
                            next;
                        }
                        $tempX2=$downext-$tempX1;
                        if (($RLen2-$tempX1-int($tempX2))>$MaxGap{$Result2[$i][0]}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                            $tempX2=int($tempX2);
                        }
                        $tempX4=$RLen2-$downlimit+1;
                        if (defined($Repeat2STR{$ref})) {
                            $pos=$tempX1+int($tempX2);
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2STR{$ref});
                            if (($upext>0)and(($RLen2-$downlimit+1)<$tempX4)) {
                                $tempX4=$RLen2-$downlimit+1;
                            }
                        }
                    } elsif ($tempX3>0) {
                        $pos=$RLen2-$tempX3;
                        ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2{$unit});
                        if (0==$upext) { ### should not happen because of $allowance
                            print ("Warning: cannot find repeats of $unit covering $pos for $Reads2\n");
                            next;
                        }
                        $tempX2=$RLen2-$upext-$tempX3+1;
                        if (($RLen2-$tempX3-int($tempX2))>$MaxGap{$Result2[$i][0]}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                            $tempX2=int($tempX2);
                        }
                        $tempX4=$uplimit;
                        if (defined($Repeat2STR{$ref})) {
                            $pos=$RLen2-$tempX3-int($tempX2)+1;
                            ($upext,$downext,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2STR{$ref});
                            if (($upext>0)and($uplimit<$tempX4)) {
                                $tempX4=$uplimit;
                            }
                        }
                    }
                }
                $Result2[$i][$j+1]=$tempX2;
                $Result2[$i][$j+3]=max($tempX4,1);
            }
        }
    }

    for ($i=$IndexNum2;$i>=1;$i--) { ### 3rd round 1st half, get %FlankAlign from bwa alignment
        $strand=$Result2[$i][1];
        $ref=join '',($Result2[$i][0],$strand);
        for ($j=$Tail2[$i]-3;$j>=2;$j=$j-4) {
            $tempX1=$Result2[$i][$j];
            $tempX2=$Result2[$i][$j+1];
            $tempX3=$Result2[$i][$j+2];
            if ($strand eq '-') {
                if ($tempX1>0) {
                    if (defined($FlankAlign{"2 $ref 3"})) {
                        $FlankAlign{"2 $ref 3"}.=" $tempX1";
                    } else {
                        $FlankAlign{"2 $ref 3"}=$tempX1;
                    }
                } elsif ($tempX3>0) {
                    if (defined($FlankAlign{"2 $ref 5"})) {
                        $FlankAlign{"2 $ref 5"}.=" $tempX3";
                    } else {
                        $FlankAlign{"2 $ref 5"}=$tempX3;
                    }
                }
            } else {
                if ($tempX1>0) {
                    if (defined($FlankAlign{"2 $ref 5"})) {
                        $FlankAlign{"2 $ref 5"}.=" $tempX1";
                    } else {
                        $FlankAlign{"2 $ref 5"}=$tempX1;
                    }
                } elsif ($tempX3>0) {
                    if (defined($FlankAlign{"2 $ref 3"})) {
                        $FlankAlign{"2 $ref 3"}.=" $tempX3";
                    } else {
                        $FlankAlign{"2 $ref 3"}=$tempX3;
                    }
                }
            }
        }
    }
    for ($i=$IndexNum2;$i>=1;$i--) { ### 3rd round 2nd half, get $Ful2{STRid+/-} by processing @Result2
        $strand=$Result2[$i][1];
        $ref=join '',($Result2[$i][0],$strand);
        for ($j=$Tail2[$i]-3;$j>=2;$j=$j-4) {
            $tempX1=$Result2[$i][$j];
            $tempX2=$Result2[$i][$j+1];
            $tempX3=$Result2[$i][$j+2];
            $findflank=0;
            if (-1==$tempX1) {
                if ($strand eq '-') {
                    $direct=3;
                } else {
                    $direct=5;
                }
                if (defined($FlankAlign{"2 $ref $direct"})) { ### search available X1 in %FlankAlign
                    @tempflank=split / /, $FlankAlign{"2 $ref $direct"};
                    for ($k=$#tempflank;$k>=0;$k--) {
                        if (($tempflank[$k]>$findflank)and($tempflank[$k]<=($RLen2-$tempX3))and($tempflank[$k]>=($RLen2-$tempX3-int($tempX2)))) {
                            $findflank=$tempflank[$k];
                        }
                    }
                }
                if ($findflank!=0) { ### found good X1 from %FlankAlign
                    $tempX1=$findflank;
                    $tempX2=$RLen2-$tempX1-$tempX3;
                } elsif (($RLen2-$tempX3-int($tempX2))<=$search_length) { ### if still possible, continue to obtain X1 by &flank_search
                    $end=$RLen2-$tempX3-int($tempX2);
                    if (0==$end) {
                        $end=1;
                    }
                    $start=min($search_length,$end+$UpFL{$Result2[$i][0]},$RLen2-$tempX3,$Result2[$i][$j+3]);
                    if (((int($tempX2)+$tempX3)==$RLen2)or($tempX2=~/\./)) { ### X2 replacement exists
                        $X2reach=1;
                    } else {
                        $X2reach=0;
                    }
                    while ($start>=$end) {
                        @tempflank=&flank_search(2,$ref,$direct,$start,$end,'no');
                        if (0==$#tempflank) { ### X1 not found, end search
                            last;
                        }
                        for ($k=0;$k<=$#tempflank-1;$k++) {
                            $findflankraw=$tempflank[$k];
                            $findflank=int($findflankraw+0.5);
                            if ($findflank>$flanklimit) { ### found good X1 by &flank_search, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX1=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX1=$findflank-0.1;
                                    }
                                } else {
                                    $tempX1=$findflank;
                                }
                                $tempX2=$RLen2-$findflank-$tempX3;
                                last;
                            } else { ### found not very good X1 by &flank_search
                                $reliab=&reliab_check(2,$ref,$direct,$findflank,$RLen2-$findflank-$tempX3,$tempX3);
                                if ($reliab>=$reliablimit) { ### X1 is trustful, end search
                                    if (1==$X2reach) { ### X2 replacement exists
                                        if ($findflank>$checklimit) {
                                            $tempX1=$findflankraw;
                                        } else { ### set a mark for further check in &PrepareResult
                                            $tempX1=$findflank-0.1;
                                        }
                                    } else { ### X2 cannot reach end, skip marking
                                        $tempX1=$findflank;
                                    }
                                    $tempX2=$RLen2-$findflank-$tempX3;
                                    last;
                                } else { ### X1 is not trustful, continue to search
                                    $findflank=0;
                                }
                            }
                        }
                        if ($findflank!=0) {
                            last;
                        } else {
                            $start=$tempflank[$#tempflank]-1;
                        }
                    }
                    if ((0==$findflank)and(1==$X2reach)) { ### X1 not found, but X2 reaches end or was incomplete due to reads length
                        $tempX1=0;
                        $tempX2=$RLen2-$tempX3;
                    }
                }
            } elsif (-1==$tempX3) {
                if ($strand eq '-') {
                    $direct=5;
                } else {
                    $direct=3;
                }
                if (defined($FlankAlign{"2 $ref $direct"})) { ### search available X3 in %FlankAlign
                    @tempflank=split / /, $FlankAlign{"2 $ref $direct"};
                    for ($k=$#tempflank;$k>=0;$k--) {
                        if (($tempflank[$k]>$findflank)and($tempflank[$k]<=($RLen2-$tempX1))and($tempflank[$k]>=($RLen2-$tempX1-int($tempX2)))) {
                            $findflank=$tempflank[$k];
                        }
                    }
                }
                if ($findflank!=0) { ### found good X3 from %FlankAlign
                    $tempX3=$findflank;
                    $tempX2=$RLen2-$tempX1-$tempX3;
                } elsif (($RLen2-$tempX1-int($tempX2))<=$search_length) { ### if still possible, continue to obtain X3 by &flank_search
                    $end=$RLen2-$tempX1-int($tempX2);
                    if (0==$end) {
                        $end=1;
                    }
                    $start=min($search_length,$end+$DownFL{$Result2[$i][0]},$RLen2-$tempX1,$Result2[$i][$j+3]);
                    if (((int($tempX2)+$tempX1)==$RLen2)or($tempX2=~/\./)) { ### X2 replacement exists
                        $X2reach=1;
                    } else {
                        $X2reach=0;
                    }
                    while ($start>=$end) {
                        @tempflank=&flank_search(2,$ref,$direct,$start,$end,'no');
                        if (0==$#tempflank) { ### X3 not found, end search
                            last;
                        }
                        for ($k=0;$k<=$#tempflank-1;$k++) {
                            $findflankraw=$tempflank[$k];
                            $findflank=int($findflankraw+0.5);
                            if ($findflank>$flanklimit) { ### found good X3 by &flank_search, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX3=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX3=$findflank-0.1;
                                    }
                                } else {
                                    $tempX3=$findflank;
                                }
                                $tempX2=$RLen2-$tempX1-$findflank;
                                last;
                            } else { ### found not very good X3 by &flank_search
                                $reliab=&reliab_check(2,$ref,$direct,$tempX1,$RLen2-$findflank-$tempX1,$findflank);
                                if ($reliab>=$reliablimit) { ### X3 is trustful, end search
                                    if (1==$X2reach) { ### X2 replacement exists
                                        if ($findflank>$checklimit) {
                                            $tempX3=$findflankraw;
                                        } else { ### set a mark for further check in &PrepareResult
                                            $tempX3=$findflank-0.1;
                                        }
                                    } else { ### X2 cannot reach end, skip marking
                                        $tempX3=$findflank;
                                    }
                                    $tempX2=$RLen2-$tempX1-$findflank;
                                    last;
                                } else { ### X3 is not trustful, continue to search
                                    $findflank=0;
                                }
                            }
                        }
                        if ($findflank!=0) {
                            last;
                        } else {
                            $start=$tempflank[$#tempflank]-1;
                        }
                    }
                    if ((0==$findflank)and(1==$X2reach)) { ### X3 not found, but X2 reaches end or was incomplete due to reads length
                        $tempX3=0;
                        $tempX2=$RLen2-$tempX1;
                    }
                }
            }
            if (($tempX1!=-1)and($tempX2!=-1)and($tempX3!=-1)) {
                if (defined($Ful2{$ref})) {
                    $Ful2{$ref}.=" $tempX1 $tempX2 $tempX3";
                } else {
                    $Ful2{$ref}="$tempX1 $tempX2 $tempX3";
                }
                if (($tempX2>0)or(($tempX1>0)and($tempX3>0))) {
                    $Ful2STR{$ref}=1;
                }
            }
        }
    }

    my $fakeX1=-1;
    my $fakeX2=-1;
    my $fakeX3=-1;
    my $realref='';
    my $start5=0;
    my $start3=0;
    my $end5=0;
    my $end3=0;
    my $sametotal=0;
    my $uptemp=0;
    my $downtemp=0;
    foreach $unit (keys(%Repeattogo1)) { ### process if $Reads1 is repeat enriched, and finish %Ful1
        $pairunit=$unit;
        if (($pairtype eq 'fr')or($pairtype eq 'rf')) {
            $pairunit=&complementary($pairunit);
        }
        $pairunit=&TransformSTR($pairunit);
        if ($Repeattogo1{$unit} ne 'y') {
            ($upext,$downext,$uplimit,$downlimit)=split / /, $Repeattogo1{$unit};
            $start5=max(min($search_length,$uplimit),1);
            $start3=max(min($search_length,$RLen1-$downlimit+1),1);
        } else { ### extend $Reads1 from $seedpos1 for possible STR $unit
            $dymatch_score=$dymscore[length($unit)];
            $upext=$seedpos1-&up_extension($Reads1,$seedpos1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
            $downext=$seedpos1+&down_extension($Reads1,$seedpos1,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
            $start5=min($search_length,$seedpos1-1);
            $start3=min($search_length,$RLen1-$seedpos1);
        }
        $end5=int($upext+0.5)-1;
        $end3=$RLen1-int($downext+0.5);
        if (0==$end5) {
            $end5=1;
        }
        if (0==$end3) {
            $end3=1;
        }
        $index=$SameUnit{$unit};
        $sametotal=$SameUnitNum[$index];
        if (($flankseed>0)and($sametotal>20)) {
            &getSeed($Reads1,$RLen1,$start5,$end5,$start3,$end3);
        }
        for ($i=$sametotal;$i>=1;$i--) { ### get all STRid+/- for this STR unit
            $ref=$SameUnitCont[$index][$i];
            $realref=substr($ref,0,length($ref)-1);
            if (($pairtype eq 'ff')or($pairtype eq 'rr')) {
                $pairref=$ref;
            } elsif (($pairtype eq 'fr')or($pairtype eq 'rf')) {
                $pairref=join '',($realref,'+');
                if ($pairref eq $ref) {
                    $pairref=join '',($realref,'-');
                }
            }
            if ((!defined($Ful2{$pairref}))and(!defined($STRtogo2{$pairref}))and(!defined($Repeattogo2{$pairunit}))) { ### skip if impossible to be paired with $Reads2
                next;
            }
            if (!(defined($Ful1{$ref}))) { ### skip if realignment for this STRid+/- already finished
                $pos=int(($upext+$downext)/2);
                $uplimit=$RLen1;
                $downlimit=1;
                if (defined($Repeat1STR{$ref})) {
                    ($uptemp,$downtemp,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat1STR{$ref});
                }
                $tempX1=-1;
                $findflank=0;
                if ($ref=~/-$/) {
                    $direct=3;
                    $fakeX3=$pos;
                    $fakeX2=$downext-$fakeX3;
                    $start=max(min($start3,$RLen1-$downlimit+1),1);
                    $end=$end3;
                } else {
                    $direct=5;
                    $fakeX3=$RLen1-$pos+1;
                    $fakeX2=$RLen1+1-$fakeX3-$upext;
                    $start=max(min($start5,$uplimit),1);
                    $end=$end5;
                }
                if (($RLen1-$fakeX3-int($fakeX2))>$MaxGap{$realref}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                    $fakeX2=int($fakeX2);
                }
                if (defined($FlankAlign{"1 $ref $direct"})) { ### search available X1 in %FlankAlign
                    @tempflank=split / /, $FlankAlign{"1 $ref $direct"};
                    for ($j=$#tempflank;$j>=0;$j--) {
                        if (($tempflank[$j]>$findflank)and($tempflank[$j]<=($RLen1-$fakeX3))and($tempflank[$j]>=($RLen1-$fakeX3-int($fakeX2)))) {
                            $findflank=$tempflank[$j];
                        }
                    }
                }
                if ($findflank!=0) { ### found good X1 from %FlankAlign
                    $tempX1=$findflank;
                } elsif (($RLen1-$fakeX3-int($fakeX2))<=$search_length) { ### if still possible, continue to obtain X1 by &flank_search
                    if ($start>($end+$UpFL{$realref})) {
                        $start=$end+$UpFL{$realref};
                    }
                    if ((($fakeX3+int($fakeX2))==$RLen1)or($fakeX2=~/\./)) { ### X2 replacement exists
                        $X2reach=1;
                    } else {
                        $X2reach=0;
                    }
                    while ($start>=$end) {
                        if (($flankseed>0)and($sametotal>20)) {
                            @tempflank=&flank_search(1,$ref,$direct,$start,$end,'yes');
                        } else {
                            @tempflank=&flank_search(1,$ref,$direct,$start,$end,'no');
                        }
                        if (0==$#tempflank) { ### X1 not found, end search
                            last;
                        }
                        for ($j=0;$j<=$#tempflank-1;$j++) {
                            $findflankraw=$tempflank[$j];
                            $findflank=int($findflankraw+0.5);
                            if ($findflank>$flanklimit) { ### found good X1 by &flank_search, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX1=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX1=$findflank-0.1;
                                    }
                                } else {
                                    $tempX1=$findflank;
                                }
                                last;
                            } else { ### found not very good X1 by &flank_search
                                $reliab=&reliab_check(1,$ref,$direct,$findflank,$RLen1-$findflank-$fakeX3,$fakeX3);
                                if ($reliab>=$reliablimit) { ### X1 is trustful, end search
                                    if (1==$X2reach) { ### X2 replacement exists
                                        if ($findflank>$checklimit) {
                                            $tempX1=$findflankraw;
                                        } else { ### set a mark for further check in &PrepareResult
                                            $tempX1=$findflank-0.1;
                                        }
                                    } else { ### X2 cannot reach end, skip marking
                                        $tempX1=$findflank;
                                    }
                                    last;
                                } else { ### X1 is not trustful, continue to search
                                    $findflank=0;
                                }
                            }
                        }
                        if ($findflank!=0) {
                            last;
                        } else {
                            $start=$tempflank[$#tempflank]-1;
                        }
                    }
                    if ((0==$findflank)and(1==$X2reach)) { ### X1 not found, but X2 reaches end or was incomplete due to reads length
                        $tempX1=0;
                    } elsif (0==$findflank) {
                        next;
                    }
                }

                $tempX3=-1;
                $findflank=0;
                if ($ref=~/-$/) {
                    $direct=5;
                    $fakeX1=$RLen1-$pos+1;
                    $fakeX2=$RLen1+1-$fakeX1-$upext;
                    $start=max(min($start5,$uplimit),1);
                    $end=$end5;
                } else {
                    $direct=3;
                    $fakeX1=$pos;
                    $fakeX2=$downext-$fakeX1;
                    $start=max(min($start3,$RLen1-$downlimit+1),1);
                    $end=$end3;
                }
                if (($RLen1-$fakeX1-int($fakeX2))>$MaxGap{$realref}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                    $fakeX2=int($fakeX2);
                }
                if (defined($FlankAlign{"1 $ref $direct"})) { ### search available X3 in %FlankAlign
                    @tempflank=split / /, $FlankAlign{"1 $ref $direct"};
                    for ($j=$#tempflank;$j>=0;$j--) {
                        if (($tempflank[$j]>$findflank)and($tempflank[$j]<=($RLen1-$fakeX1))and($tempflank[$j]>=($RLen1-$fakeX1-int($fakeX2)))) {
                            $findflank=$tempflank[$j];
                        }
                    }
                }
                if ($findflank!=0) { ### found good X3 from %FlankAlign
                    $tempX3=$findflank;
                } elsif (($RLen1-$fakeX1-int($fakeX2))<=$search_length) { ### if still possible, continue to obtain X3 by &flank_search
                    if ($start>($end+$DownFL{$realref})) {
                        $start=$end+$DownFL{$realref};
                    }
                    if ((($fakeX1+int($fakeX2))==$RLen1)or($fakeX2=~/\./)) { ### X2 replacement exists
                        $X2reach=1;
                    } else {
                        $X2reach=0;
                    }
                    while ($start>=$end) {
                        if (($flankseed>0)and($sametotal>20)) {
                            @tempflank=&flank_search(1,$ref,$direct,$start,$end,'yes');
                        } else {
                            @tempflank=&flank_search(1,$ref,$direct,$start,$end,'no');
                        }
                        if (0==$#tempflank) { ### X3 not found, end search
                            last;
                        }
                        for ($j=0;$j<=$#tempflank-1;$j++) {
                            $findflankraw=$tempflank[$j];
                            $findflank=int($findflankraw+0.5);
                            if ($findflank>$flanklimit) { ### found good X3 by &flank_search, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX3=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX3=$findflank-0.1;
                                    }
                                } else {
                                    $tempX3=$findflank;
                                }
                                last;
                            } else { ### found not very good X3 by &flank_search
                                $reliab=&reliab_check(1,$ref,$direct,$fakeX1,$RLen1-$findflank-$fakeX1,$findflank);
                                if ($reliab>=$reliablimit) { ### X3 is trustful, end search
                                    if (1==$X2reach) { ### X2 replacement exists
                                        if ($findflank>$checklimit) {
                                            $tempX3=$findflankraw;
                                        } else { ### set a mark for further check in &PrepareResult
                                            $tempX3=$findflank-0.1;
                                        }
                                    } else { ### X2 cannot reach end, skip marking
                                        $tempX3=$findflank;
                                    }
                                    last;
                                } else { ### X3 is not trustful, continue to search
                                    $findflank=0;
                                }
                            }
                        }
                        if ($findflank!=0) {
                            last;
                        } else {
                            $start=$tempflank[$#tempflank]-1;
                        }
                    }
                    if ((0==$findflank)and(1==$X2reach)) { ### X3 not found, but X2 reaches end or was incomplete due to reads length
                        $tempX3=0;
                    } elsif (0==$findflank) {
                        next;
                    }
                }

                if (($tempX1!=-1)and($tempX3!=-1)) {
                    $Ful1{$ref}=join ' ',($tempX1,$RLen1-int($tempX1+0.5)-int($tempX3+0.5),$tempX3);
                    $Ful1STR{$ref}=1;
                }
            }
        }
    }

    foreach $ref (keys(%STRtogo1)) { ### process if $Reads1 is repeat enriched for certain STRid+/-, and finish %Ful1
        if (!(defined($Ful1{$ref}))) { ### skip if realignment for this STRid+/- already finished
            $realref=substr($ref,0,length($ref)-1);
            if (($pairtype eq 'ff')or($pairtype eq 'rr')) {
                $pairref=$ref;
                if ($ref=~/-$/) {
                    $pairunit=$StrUnitcomple{$realref};
                } else {
                    $pairunit=$StrUnit{$realref};
                }
            } elsif (($pairtype eq 'fr')or($pairtype eq 'rf')) {
                $pairref=join '',($realref,'+');
                if ($pairref eq $ref) {
                    $pairref=join '',($realref,'-');
                    $pairunit=$StrUnitcomple{$realref};
                } else {
                    $pairunit=$StrUnit{$realref};
                }
            }
            if ((!defined($Ful2{$pairref}))and(!defined($STRtogo2{$pairref}))and(!defined($Repeattogo2{$pairunit}))) { ### skip if impossible to be paired with $Reads2
                next;
            }
            ($upext,$downext,$uplimit,$downlimit)=split / /, $STRtogo1{$ref};
            $start5=max(min($search_length,$uplimit),1);
            $start3=max(min($search_length,$RLen1-$downlimit+1),1);
            $end5=int($upext+0.5)-1;
            $end3=$RLen1-int($downext+0.5);
            if (0==$end5) {
                $end5=1;
            }
            if (0==$end3) {
                $end3=1;
            }

            $tempX1=-1;
            $findflank=0;
            if ($ref=~/-$/) {
                $direct=3;
                $fakeX3=int(($upext+$downext)/2);
                $fakeX2=$downext-$fakeX3;
                $start=$start3;
                $end=$end3;
            } else {
                $direct=5;
                $fakeX3=$RLen1-int(($upext+$downext)/2)+1;
                $fakeX2=$RLen1+1-$fakeX3-$upext;
                $start=$start5;
                $end=$end5;
            }
            if (($RLen1-$fakeX3-int($fakeX2))>$MaxGap{$realref}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                $fakeX2=int($fakeX2);
            }
            if (defined($FlankAlign{"1 $ref $direct"})) { ### search available X1 in %FlankAlign
                @tempflank=split / /, $FlankAlign{"1 $ref $direct"};
                for ($j=$#tempflank;$j>=0;$j--) {
                    if (($tempflank[$j]>$findflank)and($tempflank[$j]<=($RLen1-$fakeX3))and($tempflank[$j]>=($RLen1-$fakeX3-int($fakeX2)))) {
                        $findflank=$tempflank[$j];
                    }
                }
            }
            if ($findflank!=0) { ### found good X1 from %FlankAlign
                $tempX1=$findflank;
            } elsif (($RLen1-$fakeX3-int($fakeX2))<=$search_length) { ### if still possible, continue to obtain X1 by &flank_search
                if ($start>($end+$UpFL{$realref})) {
                    $start=$end+$UpFL{$realref};
                }
                if ((($fakeX3+int($fakeX2))==$RLen1)or($fakeX2=~/\./)) { ### X2 replacement exists
                    $X2reach=1;
                } else {
                    $X2reach=0;
                }
                while ($start>=$end) {
                    @tempflank=&flank_search(1,$ref,$direct,$start,$end,'no');
                    if (0==$#tempflank) { ### X1 not found, end search
                        last;
                    }
                    for ($j=0;$j<=$#tempflank-1;$j++) {
                        $findflankraw=$tempflank[$j];
                        $findflank=int($findflankraw+0.5);
                        if ($findflank>$flanklimit) { ### found good X1 by &flank_search, end search
                            if (1==$X2reach) { ### X2 replacement exists
                                if ($findflank>$checklimit) {
                                    $tempX1=$findflankraw;
                                } else { ### set a mark for further check in &PrepareResult
                                    $tempX1=$findflank-0.1;
                                }
                            } else {
                                $tempX1=$findflank;
                            }
                            last;
                        } else { ### found not very good X1 by &flank_search
                            $reliab=&reliab_check(1,$ref,$direct,$findflank,$RLen1-$findflank-$fakeX3,$fakeX3);
                            if ($reliab>=$reliablimit) { ### X1 is trustful, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX1=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX1=$findflank-0.1;
                                    }
                                } else { ### X2 cannot reach end, skip marking
                                    $tempX1=$findflank;
                                }
                                last;
                            } else { ### X1 is not trustful, continue to search
                                $findflank=0;
                            }
                        }
                    }
                    if ($findflank!=0) {
                        last;
                    } else {
                        $start=$tempflank[$#tempflank]-1;
                    }
                }
                if ((0==$findflank)and(1==$X2reach)) { ### X1 not found, but X2 reaches end or was incomplete due to reads length
                    $tempX1=0;
                } elsif (0==$findflank) {
                    next;
                }
            }

            $tempX3=-1;
            $findflank=0;
            if ($ref=~/-$/) {
                $direct=5;
                $fakeX1=$RLen1-int(($upext+$downext)/2)+1;
                $fakeX2=$RLen1+1-$fakeX1-$upext;
                $start=$start5;
                $end=$end5;
            } else {
                $direct=3;
                $fakeX1=int(($upext+$downext)/2);
                $fakeX2=$downext-$fakeX1;
                $start=$start3;
                $end=$end3;
            }
            if (($RLen1-$fakeX1-int($fakeX2))>$MaxGap{$realref}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                $fakeX2=int($fakeX2);
            }
            if (defined($FlankAlign{"1 $ref $direct"})) { ### search available X3 in %FlankAlign
                @tempflank=split / /, $FlankAlign{"1 $ref $direct"};
                for ($j=$#tempflank;$j>=0;$j--) {
                    if (($tempflank[$j]>$findflank)and($tempflank[$j]<=($RLen1-$fakeX1))and($tempflank[$j]>=($RLen1-$fakeX1-int($fakeX2)))) {
                        $findflank=$tempflank[$j];
                    }
                }
            }
            if ($findflank!=0) { ### found good X3 from %FlankAlign
                $tempX3=$findflank;
            } elsif (($RLen1-$fakeX1-int($fakeX2))<=$search_length) { ### if still possible, continue to obtain X3 by &flank_search
                if ($start>($end+$DownFL{$realref})) {
                    $start=$end+$DownFL{$realref};
                }
                if ((($fakeX1+int($fakeX2))==$RLen1)or($fakeX2=~/\./)) { ### X2 replacement exists
                    $X2reach=1;
                } else {
                    $X2reach=0;
                }
                while ($start>=$end) {
                    @tempflank=&flank_search(1,$ref,$direct,$start,$end,'no');
                    if (0==$#tempflank) { ### X3 not found, end search
                        last;
                    }
                    for ($j=0;$j<=$#tempflank-1;$j++) {
                        $findflankraw=$tempflank[$j];
                        $findflank=int($findflankraw+0.5);
                        if ($findflank>$flanklimit) { ### found good X3 by &flank_search, end search
                            if (1==$X2reach) { ### X2 replacement exists
                                if ($findflank>$checklimit) {
                                    $tempX3=$findflankraw;
                                } else { ### set a mark for further check in &PrepareResult
                                    $tempX3=$findflank-0.1;
                                }
                            } else {
                                $tempX3=$findflank;
                            }
                            last;
                        } else { ### found not very good X3 by &flank_search
                            $reliab=&reliab_check(1,$ref,$direct,$fakeX1,$RLen1-$findflank-$fakeX1,$findflank);
                            if ($reliab>=$reliablimit) { ### X3 is trustful, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX3=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX3=$findflank-0.1;
                                    }
                                } else { ### X2 cannot reach end, skip marking
                                    $tempX3=$findflank;
                                }
                                last;
                            } else { ### X3 is not trustful, continue to search
                                $findflank=0;
                            }
                        }
                    }
                    if ($findflank!=0) {
                        last;
                    } else {
                        $start=$tempflank[$#tempflank]-1;
                    }
                }
                if ((0==$findflank)and(1==$X2reach)) { ### X3 not found, but X2 reaches end or was incomplete due to reads length
                    $tempX3=0;
                } elsif (0==$findflank) {
                    next;
                }
            }

            if (($tempX1!=-1)and($tempX3!=-1)) {
                $Ful1{$ref}=join ' ',($tempX1,$RLen1-int($tempX1+0.5)-int($tempX3+0.5),$tempX3);
                $Ful1STR{$ref}=1;
            }
        }
    }

    foreach $unit (keys(%Repeattogo2)) { ### process if $Reads2 is repeat enriched, and finish %Ful2
        if ($Repeattogo2{$unit} ne 'y') {
            ($upext,$downext,$uplimit,$downlimit)=split / /, $Repeattogo2{$unit};
            $start5=max(min($search_length,$uplimit),1);
            $start3=max(min($search_length,$RLen2-$downlimit+1),1);
        } else { ### extend $Reads2 from $seedpos2 for possible STR $unit
            $dymatch_score=$dymscore[length($unit)];
            $upext=$seedpos2-&up_extension($Reads2,$seedpos2,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
            $downext=$seedpos2+&down_extension($Reads2,$seedpos2,$unit,$dymatch_score,$mis_score,$gap_score,$limit_score,$allowance_score-$dymatch_score);
            $start5=min($search_length,$seedpos2-1);
            $start3=min($search_length,$RLen2-$seedpos2);
        }
        $end5=int($upext+0.5)-1;
        $end3=$RLen2-int($downext+0.5);
        if (0==$end5) {
            $end5=1;
        }
        if (0==$end3) {
            $end3=1;
        }
        $index=$SameUnit{$unit};
        $sametotal=$SameUnitNum[$index];
        if (($flankseed>0)and($sametotal>20)) {
            &getSeed($Reads2,$RLen2,$start5,$end5,$start3,$end3);
        }
        for ($i=$sametotal;$i>=1;$i--) { ### get all STRid+/- for this STR unit
            $ref=$SameUnitCont[$index][$i];
            $realref=substr($ref,0,length($ref)-1);
            if (($pairtype eq 'ff')or($pairtype eq 'rr')) {
                $pairref=$ref;
            } elsif (($pairtype eq 'fr')or($pairtype eq 'rf')) {
                $pairref=join '',($realref,'+');
                if ($pairref eq $ref) {
                    $pairref=join '',($realref,'-');
                }
            }
            if (!(defined($Ful1{$pairref}))) { ### skip if impossible to be paired with $Reads1
                next;
            }
            if (!(defined($Ful2{$ref}))) { ### skip if realignment for this STRid+/- already finished
                $pos=int(($upext+$downext)/2);
                $uplimit=$RLen2;
                $downlimit=1;
                if (defined($Repeat2STR{$ref})) {
                    ($uptemp,$downtemp,$uplimit,$downlimit)=&SearchRepeat($pos,$Repeat2STR{$ref});
                }
                $tempX1=-1;
                $findflank=0;
                if ($ref=~/-$/) {
                    $direct=3;
                    $fakeX3=$pos;
                    $fakeX2=$downext-$fakeX3;
                    $start=max(min($start3,$RLen2-$downlimit+1),1);
                    $end=$end3;
                } else {
                    $direct=5;
                    $fakeX3=$RLen2-$pos+1;
                    $fakeX2=$RLen2+1-$fakeX3-$upext;
                    $start=max(min($start5,$uplimit),1);
                    $end=$end5;
                }
                if (($RLen2-$fakeX3-int($fakeX2))>$MaxGap{$realref}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                    $fakeX2=int($fakeX2);
                }
                if (defined($FlankAlign{"2 $ref $direct"})) { ### search available X1 in %FlankAlign
                    @tempflank=split / /, $FlankAlign{"2 $ref $direct"};
                    for ($j=$#tempflank;$j>=0;$j--) {
                        if (($tempflank[$j]>$findflank)and($tempflank[$j]<=($RLen2-$fakeX3))and($tempflank[$j]>=($RLen2-$fakeX3-int($fakeX2)))) {
                            $findflank=$tempflank[$j];
                        }
                    }
                }
                if ($findflank!=0) { ### found good X1 from %FlankAlign
                    $tempX1=$findflank;
                } elsif (($RLen2-$fakeX3-int($fakeX2))<=$search_length) { ### if still possible, continue to obtain X1 by &flank_search
                    if ($start>($end+$UpFL{$realref})) {
                        $start=$end+$UpFL{$realref};
                    }
                    if ((($fakeX3+int($fakeX2))==$RLen2)or($fakeX2=~/\./)) { ### X2 replacement exists
                        $X2reach=1;
                    } else {
                        $X2reach=0;
                    }
                    while ($start>=$end) {
                        if (($flankseed>0)and($sametotal>20)) {
                            @tempflank=&flank_search(2,$ref,$direct,$start,$end,'yes');
                        } else {
                            @tempflank=&flank_search(2,$ref,$direct,$start,$end,'no');
                        }
                        if (0==$#tempflank) { ### X1 not found, end search
                            last;
                        }
                        for ($j=0;$j<=$#tempflank-1;$j++) {
                            $findflankraw=$tempflank[$j];
                            $findflank=int($findflankraw+0.5);
                            if ($findflank>$flanklimit) { ### found good X1 by &flank_search, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX1=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX1=$findflank-0.1;
                                    }
                                } else {
                                    $tempX1=$findflank;
                                }
                                last;
                            } else { ### found not very good X1 by &flank_search
                                $reliab=&reliab_check(2,$ref,$direct,$findflank,$RLen2-$findflank-$fakeX3,$fakeX3);
                                if ($reliab>=$reliablimit) { ### X1 is trustful, end search
                                    if (1==$X2reach) { ### X2 replacement exists
                                        if ($findflank>$checklimit) {
                                            $tempX1=$findflankraw;
                                        } else { ### set a mark for further check in &PrepareResult
                                            $tempX1=$findflank-0.1;
                                        }
                                    } else { ### X2 cannot reach end, skip marking
                                        $tempX1=$findflank;
                                    }
                                    last;
                                } else { ### X1 is not trustful, continue to search
                                    $findflank=0;
                                }
                            }
                        }
                        if ($findflank!=0) {
                            last;
                        } else {
                            $start=$tempflank[$#tempflank]-1;
                        }
                    }
                    if ((0==$findflank)and(1==$X2reach)) { ### X1 not found, but X2 reaches end or was incomplete due to reads length
                        $tempX1=0;
                    } elsif (0==$findflank) {
                        next;
                    }
                }

                $tempX3=-1;
                $findflank=0;
                if ($ref=~/-$/) {
                    $direct=5;
                    $fakeX1=$RLen2-$pos+1;
                    $fakeX2=$RLen2+1-$fakeX1-$upext;
                    $start=max(min($start5,$uplimit),1);
                    $end=$end5;
                } else {
                    $direct=3;
                    $fakeX1=$pos;
                    $fakeX2=$downext-$fakeX1;
                    $start=max(min($start3,$RLen2-$downlimit+1),1);
                    $end=$end3;
                }
                if (($RLen2-$fakeX1-int($fakeX2))>$MaxGap{$realref}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                    $fakeX2=int($fakeX2);
                }
                if (defined($FlankAlign{"2 $ref $direct"})) { ### search available X3 in %FlankAlign
                    @tempflank=split / /, $FlankAlign{"2 $ref $direct"};
                    for ($j=$#tempflank;$j>=0;$j--) {
                        if (($tempflank[$j]>$findflank)and($tempflank[$j]<=($RLen2-$fakeX1))and($tempflank[$j]>=($RLen2-$fakeX1-int($fakeX2)))) {
                            $findflank=$tempflank[$j];
                        }
                    }
                }
                if ($findflank!=0) { ### found good X3 from %FlankAlign
                    $tempX3=$findflank;
                } elsif (($RLen2-$fakeX1-int($fakeX2))<=$search_length) { ### if still possible, continue to obtain X3 by &flank_search
                    if ($start>($end+$DownFL{$realref})) {
                        $start=$end+$DownFL{$realref};
                    }
                    if ((($fakeX1+int($fakeX2))==$RLen2)or($fakeX2=~/\./)) { ### X2 replacement exists
                        $X2reach=1;
                    } else {
                        $X2reach=0;
                    }
                    while ($start>=$end) {
                        if (($flankseed>0)and($sametotal>20)) {
                            @tempflank=&flank_search(2,$ref,$direct,$start,$end,'yes');
                        } else {
                            @tempflank=&flank_search(2,$ref,$direct,$start,$end,'no');
                        }
                        if (0==$#tempflank) { ### X3 not found, end search
                            last;
                        }
                        for ($j=0;$j<=$#tempflank-1;$j++) {
                            $findflankraw=$tempflank[$j];
                            $findflank=int($findflankraw+0.5);
                            if ($findflank>$flanklimit) { ### found good X3 by &flank_search, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX3=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX3=$findflank-0.1;
                                    }
                                } else {
                                    $tempX3=$findflank;
                                }
                                last;
                            } else { ### found not very good X3 by &flank_search
                                $reliab=&reliab_check(2,$ref,$direct,$fakeX1,$RLen2-$findflank-$fakeX1,$findflank);
                                if ($reliab>=$reliablimit) { ### X3 is trustful, end search
                                    if (1==$X2reach) { ### X2 replacement exists
                                        if ($findflank>$checklimit) {
                                            $tempX3=$findflankraw;
                                        } else { ### set a mark for further check in &PrepareResult
                                            $tempX3=$findflank-0.1;
                                        }
                                    } else { ### X2 cannot reach end, skip marking
                                        $tempX3=$findflank;
                                    }
                                    last;
                                } else { ### X3 is not trustful, continue to search
                                    $findflank=0;
                                }
                            }
                        }
                        if ($findflank!=0) {
                            last;
                        } else {
                            $start=$tempflank[$#tempflank]-1;
                        }
                    }
                    if ((0==$findflank)and(1==$X2reach)) { ### X3 not found, but X2 reaches end or was incomplete due to reads length
                        $tempX3=0;
                    } elsif (0==$findflank) {
                        next;
                    }
                }

                if (($tempX1!=-1)and($tempX3!=-1)) {
                    $Ful2{$ref}=join ' ',($tempX1,$RLen2-int($tempX1+0.5)-int($tempX3+0.5),$tempX3);
                    $Ful2STR{$ref}=1;
                }
            }
        }
    }

    foreach $ref (keys(%STRtogo2)) { ### process if $Reads2 is repeat enriched for certain STRid+/-, and finish %Ful2
        if (!(defined($Ful2{$ref}))) { ### skip if realignment for this STRid+/- already finished
            $realref=substr($ref,0,length($ref)-1);
            if (($pairtype eq 'ff')or($pairtype eq 'rr')) {
                $pairref=$ref;
            } elsif (($pairtype eq 'fr')or($pairtype eq 'rf')) {
                $pairref=join '',($realref,'+');
                if ($pairref eq $ref) {
                    $pairref=join '',($realref,'-');
                }
            }
            if (!(defined($Ful1{$pairref}))) { ### skip if impossible to be paired with $Reads1
                next;
            }
            ($upext,$downext,$uplimit,$downlimit)=split / /, $STRtogo2{$ref};
            $start5=max(min($search_length,$uplimit),1);
            $start3=max(min($search_length,$RLen2-$downlimit+1),1);
            $end5=int($upext+0.5)-1;
            $end3=$RLen2-int($downext+0.5);
            if (0==$end5) {
                $end5=1;
            }
            if (0==$end3) {
                $end3=1;
            }

            $tempX1=-1;
            $findflank=0;
            if ($ref=~/-$/) {
                $direct=3;
                $fakeX3=int(($upext+$downext)/2);
                $fakeX2=$downext-$fakeX3;
                $start=$start3;
                $end=$end3;
            } else {
                $direct=5;
                $fakeX3=$RLen2-int(($upext+$downext)/2)+1;
                $fakeX2=$RLen2+1-$fakeX3-$upext;
                $start=$start5;
                $end=$end5;
            }
            if (($RLen2-$fakeX3-int($fakeX2))>$MaxGap{$realref}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                $fakeX2=int($fakeX2);
            }
            if (defined($FlankAlign{"2 $ref $direct"})) { ### search available X1 in %FlankAlign
                @tempflank=split / /, $FlankAlign{"2 $ref $direct"};
                for ($j=$#tempflank;$j>=0;$j--) {
                    if (($tempflank[$j]>$findflank)and($tempflank[$j]<=($RLen2-$fakeX3))and($tempflank[$j]>=($RLen2-$fakeX3-int($fakeX2)))) {
                        $findflank=$tempflank[$j];
                    }
                }
            }
            if ($findflank!=0) { ### found good X1 from %FlankAlign
                $tempX1=$findflank;
            } elsif (($RLen2-$fakeX3-int($fakeX2))<=$search_length) { ### if still possible, continue to obtain X1 by &flank_search
                if ($start>($end+$UpFL{$realref})) {
                    $start=$end+$UpFL{$realref};
                }
                if ((($fakeX3+int($fakeX2))==$RLen2)or($fakeX2=~/\./)) { ### X2 replacement exists
                    $X2reach=1;
                } else {
                    $X2reach=0;
                }
                while ($start>=$end) {
                    @tempflank=&flank_search(2,$ref,$direct,$start,$end,'no');
                    if (0==$#tempflank) { ### X1 not found, end search
                        last;
                    }
                    for ($j=0;$j<=$#tempflank-1;$j++) {
                        $findflankraw=$tempflank[$j];
                        $findflank=int($findflankraw+0.5);
                        if ($findflank>$flanklimit) { ### found good X1 by &flank_search, end search
                            if (1==$X2reach) { ### X2 replacement exists
                                if ($findflank>$checklimit) {
                                    $tempX1=$findflankraw;
                                } else { ### set a mark for further check in &PrepareResult
                                    $tempX1=$findflank-0.1;
                                }
                            } else {
                                $tempX1=$findflank;
                            }
                            last;
                        } else { ### found not very good X1 by &flank_search
                            $reliab=&reliab_check(2,$ref,$direct,$findflank,$RLen2-$findflank-$fakeX3,$fakeX3);
                            if ($reliab>=$reliablimit) { ### X1 is trustful, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX1=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX1=$findflank-0.1;
                                    }
                                } else { ### X2 cannot reach end, skip marking
                                    $tempX1=$findflank;
                                }
                                last;
                            } else { ### X1 is not trustful, continue to search
                                $findflank=0;
                            }
                        }
                    }
                    if ($findflank!=0) {
                        last;
                    } else {
                        $start=$tempflank[$#tempflank]-1;
                    }
                }
                if ((0==$findflank)and(1==$X2reach)) { ### X1 not found, but X2 reaches end or was incomplete due to reads length
                    $tempX1=0;
                } elsif (0==$findflank) {
                    next;
                }
            }

            $tempX3=-1;
            $findflank=0;
            if ($ref=~/-$/) {
                $direct=5;
                $fakeX1=$RLen2-int(($upext+$downext)/2)+1;
                $fakeX2=$RLen2+1-$fakeX1-$upext;
                $start=$start5;
                $end=$end5;
            } else {
                $direct=3;
                $fakeX1=int(($upext+$downext)/2);
                $fakeX2=$downext-$fakeX1;
                $start=$start3;
                $end=$end3;
            }
            if (($RLen2-$fakeX1-int($fakeX2))>$MaxGap{$realref}) { ### remove the incomplete marker of X2 when remaining undetermined part is too long
                $fakeX2=int($fakeX2);
            }
            if (defined($FlankAlign{"2 $ref $direct"})) { ### search available X3 in %FlankAlign
                @tempflank=split / /, $FlankAlign{"2 $ref $direct"};
                for ($j=$#tempflank;$j>=0;$j--) {
                    if (($tempflank[$j]>$findflank)and($tempflank[$j]<=($RLen2-$fakeX1))and($tempflank[$j]>=($RLen2-$fakeX1-int($fakeX2)))) {
                        $findflank=$tempflank[$j];
                    }
                }
            }
            if ($findflank!=0) { ### found good X3 from %FlankAlign
                $tempX3=$findflank;
            } elsif (($RLen2-$fakeX1-int($fakeX2))<=$search_length) { ### if still possible, continue to obtain X3 by &flank_search
                if ($start>($end+$DownFL{$realref})) {
                    $start=$end+$DownFL{$realref};
                }
                if ((($fakeX1+int($fakeX2))==$RLen2)or($fakeX2=~/\./)) { ### X2 replacement exists
                    $X2reach=1;
                } else {
                    $X2reach=0;
                }
                while ($start>=$end) {
                    @tempflank=&flank_search(2,$ref,$direct,$start,$end,'no');
                    if (0==$#tempflank) { ### X3 not found, end search
                        last;
                    }
                    for ($j=0;$j<=$#tempflank-1;$j++) {
                        $findflankraw=$tempflank[$j];
                        $findflank=int($findflankraw+0.5);
                        if ($findflank>$flanklimit) { ### found good X3 by &flank_search, end search
                            if (1==$X2reach) { ### X2 replacement exists
                                if ($findflank>$checklimit) {
                                    $tempX3=$findflankraw;
                                } else { ### set a mark for further check in &PrepareResult
                                    $tempX3=$findflank-0.1;
                                }
                            } else {
                                $tempX3=$findflank;
                            }
                            last;
                        } else { ### found not very good X3 by &flank_search
                            $reliab=&reliab_check(2,$ref,$direct,$fakeX1,$RLen2-$findflank-$fakeX1,$findflank);
                            if ($reliab>=$reliablimit) { ### X3 is trustful, end search
                                if (1==$X2reach) { ### X2 replacement exists
                                    if ($findflank>$checklimit) {
                                        $tempX3=$findflankraw;
                                    } else { ### set a mark for further check in &PrepareResult
                                        $tempX3=$findflank-0.1;
                                    }
                                } else { ### X2 cannot reach end, skip marking
                                    $tempX3=$findflank;
                                }
                                last;
                            } else { ### X3 is not trustful, continue to search
                                $findflank=0;
                            }
                        }
                    }
                    if ($findflank!=0) {
                        last;
                    } else {
                        $start=$tempflank[$#tempflank]-1;
                    }
                }
                if ((0==$findflank)and(1==$X2reach)) { ### X3 not found, but X2 reaches end or was incomplete due to reads length
                    $tempX3=0;
                } elsif (0==$findflank) {
                    next;
                }
            }

            if (($tempX1!=-1)and($tempX3!=-1)) {
                $Ful2{$ref}=join ' ',($tempX1,$RLen2-int($tempX1+0.5)-int($tempX3+0.5),$tempX3);
                $Ful2STR{$ref}=1;
            }
        }
    }

    my $temp1='';
    my $temp2='';
    foreach $ref (keys(%Ful2)) { ### store paired $Ful1/2 into @RawResult
        $realref=substr($ref,0,length($ref)-1);
        if (($pairtype eq 'ff')or($pairtype eq 'rr')) {
            $pairref=$ref;
        } elsif (($pairtype eq 'fr')or($pairtype eq 'rf')) {
            $pairref=join '',($realref,'+');
            if ($pairref eq $ref) {
                $pairref=join '',($realref,'-');
            }
        }
        if ((defined($Ful1{$pairref}))and((defined($Ful1STR{$pairref}))or(defined($Ful2STR{$ref})))) {
            $temp1=substr($pairref,-1,1);
            $temp2=substr($ref,-1,1);
            $index=$RawResultIndex{$realref};
            $RawResultMax[$index]++;
            $RawResult[$index][$RawResultMax[$index]]="$temp1 $_[0] 1";
            $RawResultMax[$index]++;
            $RawResult[$index][$RawResultMax[$index]]=$Ful1{$pairref};
            $RawResultMax[$index]++;
            $RawResult[$index][$RawResultMax[$index]]="$temp2 $_[0] 2";
            $RawResultMax[$index]++;
            $RawResult[$index][$RawResultMax[$index]]=$Ful2{$ref};
        }
    }

}

sub OutputResult { ### first 2 lines for general info, then 1 empty line, then for each STR: one line for STR info, one line for how many pairs realigned
    ### following lines for each mapped pair: "pairID {X1,X2,X3}{X1,X2,X3} {X1,X2,X3}{X1,X2,X3} ..", last line ""
    ### each {}{} means a realignment for this pair, X1=length mapped to upflank, X2=length of repeat sequence, X3=length mapped to downflank
    ### first {X1 X2 X3} from reads1, second {X1 X2 X3} from reads2
    print ("Filtering out unreliable & overlapped results...\n");
    my $i=0;
    my $j=0;
    my $k=0;
    my $count=0;
    my %check=();
    my $temp='';
    my $filehandle='';
    open (OUTPUT_RESULT, ">$_[0]") or die "Couldn't open: $!";
    print OUTPUT_RESULT ("STR source: $Para[0] built by $reffile\n");
    print OUTPUT_RESULT (".sam/.bam source: $Para[1]\n\n");
    for ($i=0;$i<=$#StrInfo;$i++) {
        print OUTPUT_RESULT ("$StrInfo[$i]\n");
        @RawResult=();
        $RawResultMax[$i]=-1;
        for ($j=1;$j<=$totalnumber;$j++) {
            $filehandle=$FH[$j];
            chomp ($temp=<$filehandle>);
            $temp=~/ /;
            if ($`!=$i) {
                print ("Error reading sub results\n");
            } else {
                for ($k=0;$k<=$';$k++) {
                    chomp ($temp=<$filehandle>);
                    $RawResultMax[$`]++;
                    $RawResult[$`][$RawResultMax[$`]]=$temp;
                }
            }
        }
        &PrepareResult($i);
        if (-1==$#AlignResult) {
            print OUTPUT_RESULT ("0 pairs realigned\n");
        } else {
            %check=();
            $count=0;
            for ($j=$#AlignResult-6;$j>=0;$j=$j-7) {
                if (!(defined($check{$AlignResult[$j]}))) {
                    $check{$AlignResult[$j]}=1;
                    $count++;
                } else {
                    $check{$AlignResult[$j]}++;
                }
            }
            print OUTPUT_RESULT ("$count pairs realigned\n");
            $j=$#AlignResult-6;
            while ($j>=0) {
                $temp=$AlignResult[$j];
                for ($k=$check{$AlignResult[$j]};$k>=1;$k--) {
                    $temp=join '',($temp,' {',$AlignResult[$j+1],',',$AlignResult[$j+2],',',$AlignResult[$j+3],'}{',$AlignResult[$j+4],',',$AlignResult[$j+5],',',$AlignResult[$j+6],'}');
                    $j=$j-7;
                }
                print OUTPUT_RESULT ("$temp\n");
            }
        }
        print OUTPUT_RESULT ("\n");
    }
    close OUTPUT_RESULT;
    print ("All realignment results written in $_[0]\n");
}

sub PossibleRepeat { ### return X, means there might be a repeat region in reads $_[0] covering Xth nt
### use k-mer to check each position, if it contributes a new k-mer (1) or not (0)
### find the longest 0..0, if it passes the length test, return the middle position
### if no likely repeat region is found, return 0
### add in %Repeattogo1list ($_[1]=1) or %Repeattogo2list ($_[1]=2) for possible STR units
    my @start=();
    my @end=();
    my @head=();
    my $i=0;
    my $kmer=3; ### k-mer length
    my $range=3*$max_ul;
    my $limit=int($RepeatPercentLimit*length($_[0]));
    my $temp='';
    my %check=();
    my $record='no';
    for ($i=0;$i<=length($_[0])-$kmer;$i++) {
        $temp=substr($_[0],$i,$kmer);
        if (defined($check{$temp})) {
            if (($i-$check{$temp})<=$range) {
                if ($record eq 'no') {
                    push @start, ($i+1);
                    push @head, ($check{$temp}+1);
                    $record='yes';
                } elsif (($check{$temp}+1)<$head[$#head]) {
                    $head[$#head]=$check{$temp}+1;
                }
            } else {
                if ($record eq 'yes') {
                    push @end, ($i);
                    $record='no';
                }
            }
        } else {
            if ($record eq 'yes') {
                push @end, ($i);
                $record='no';
            }
        }
        $check{$temp}=$i;
    }
    if ($record eq 'yes') {
        push @end, (length($_[0])-$kmer+1);
    }
    if ($#start!=$#end) {
        print ("Error: Unknown\n");
    }

    my $filter=10; ### 0s shorter than this will not be used
    my $long0s=0;
    my $max0s=0;
    my $result=int(length($_[0])/2); ### this value will be used only when $limit=0 & $#start=-1
    for ($i=0;$i<=$#start;$i++) { ### get total length of 0s and middle position of longest 0s
        if (($end[$i]-$start[$i]+1)>=$filter) {
            $long0s=$long0s+$end[$i]-$start[$i]+1;
            if (($end[$i]-$start[$i]+1)>$max0s) {
                $max0s=$end[$i]-$start[$i]+1;
                $result=int(($end[$i]+$start[$i])/2);
            }
        }
    }
    if ($long0s<$limit) {
        return (0);
    }

    my $ul=0;
    my $tempseq='';
    my %sum=();
    my %mark=();
    my $maxsum=0;
    my $j=0;
    my $transtemp='';
    my $threshold=0;
    for ($ul=$min_ul;$ul<=$max_ul;$ul++) { ### process each $ul length
        %sum=();
        $threshold=0;
        for ($j=0;$j<=$#start;$j++) {
            if (($end[$j]-$start[$j]+1)>=$filter) {
                $tempseq=substr($_[0],$head[$j]-1,$end[$j]+$kmer-$head[$j]);
                %check=();
                %mark=();
                for ($i=0;$i<=length($tempseq)-$ul;$i++) {
                    $temp=substr($tempseq,$i,$ul);
                    if (defined($check{$temp})) {
                        if ((($i-$check{$temp})==$ul)or(($i-$check{$temp})==($ul+1))or(($i-$check{$temp})==(2*$ul))or(($i-$check{$temp})==(2*$ul+1))or(($i-$check{$temp})==(2*$ul-1))) {
                            $transtemp=&TransformSTR($temp);
                            if ($transtemp ne '') {
                                if (defined($sum{$transtemp})) {
                                    $sum{$transtemp}++;
                                } else {
                                    $sum{$transtemp}=1;
                                }
                                if (0==$mark{$temp}) {
                                    $sum{$transtemp}++;
                                }
                            }
                            $mark{$temp}=1;
                        } else {
                            $mark{$temp}=0;
                        }
                    } else {
                        $mark{$temp}=0;
                    }
                    $check{$temp}=$i;
                }
                $threshold=$threshold+$i*(0.99-0.14*$ul); ### require 85% for 1nt STR, 15% for 6nt STR
            }
        }
        $maxsum=0;
        foreach (keys(%sum)) { ### find the count of best repeat unit of $ul length
            if (($sum{$_}>$maxsum)and($sum{$_}>=$threshold)) {
                $maxsum=$sum{$_};
            }
        }
        if ($maxsum>0) {
            $maxsum=$maxsum-$search_depth*$ul;
            foreach (keys(%sum)) {
                if ($sum{$_}>=$maxsum) {
                    if (1==$_[1]) {
                        $Repeattogo1list{$_}=1;
                    } elsif (2==$_[1]) {
                        $Repeattogo2list{$_}=1;
                    }
                }
            }
        }
    }
    return ($result);
}

sub reliab_check { ### return best alignment score of last query nt
    ### [0]=1/2 to use full Reads1/2 as query string [1]=STRid+/- to get STR/STR+flank sequence as reference
    ### [2]=5 (check requested by a potential flank at 5' end) or 3 (check requested by a potential flank at 3' end)
    ### [3][4][5]=combination of X1,X2,X3 to check reliability
    my @queryseq=(); ### store nt processed from $_[0]
    my @refseq=(); ### store nt processed from $_[1] and $_[2]
    my $result=0;
    if (($reliablimit<=0)or($reliablimit>$flankcheck_length)) { ### skip reliability check, as the result will be 0~$flankcheck_length
        $result=$flankcheck_length/2;
        return($result);
    }
    if (((1!=$_[0])and(2!=$_[0]))or((5!=$_[2])and(3!=$_[2]))or($#_<5)) {
        print ("Warning: Invalid parameters for &reliab_search, skipped\n");
        return($result);
    }
    my $qstring='';
    if (1==$_[0]) {
        $qstring=$Reads1;
    } else {
        $qstring=$Reads2;
    }
    my $qlength=$_[3]+$_[4]+$_[5];
    my $rstring='';
    my $rlength=0;
    my $ref=$_[1];
    my $temp=substr($ref,0,length($ref)-1);
    my $wantX1=0;
    my $wantX2=0;
    my $wantX3=0;
    my $X1string='';
    my $X2string='';
    my $X3string='';
    my $i=0;
    if ($ref=~/-$/) {
        if (5==$_[2]) { ### this check is for a potential X3
            $qstring=substr($qstring,$_[5],$flankcheck_length);
            $qlength=length($qstring);
            if ($_[4]>=($qlength*2)) {
                $wantX2=$qlength*2;
            } else {
                $wantX2=$_[4];
            }
            $wantX1=$qlength*2-$wantX2;
            if (($wantX2>0)and(defined($StrDown_comple{$temp}))) {
                $X2string=substr($StrDown_comple{$temp},0,$wantX2);
            } else {
                $X2string='';
            }
            if (($wantX1>0)and(defined($UpStream_comple{$temp}))) {
                $X1string=substr($UpStream_comple{$temp},0,$wantX1);
            } else {
                $X1string='';
            }
            $rstring=join '',($X2string,$X1string);
            $rlength=length($rstring);
            for ($i=1;$i<=$qlength;$i++) { ### generate @queryseq in forward direction
                $queryseq[$i]=substr($qstring,$i-1,1);
            }
            for ($i=1;$i<=$rlength;$i++) { ### generate @refseq in forward direction
                $refseq[$i]=substr($rstring,$i-1,1);
            }
        } else { ### this check is for a potential X1
            if (($qlength-$_[3])>=$flankcheck_length) {
                $qstring=substr($qstring,$qlength-$_[3]-$flankcheck_length,$flankcheck_length);
            } else {
                $qstring=substr($qstring,0,$qlength-$_[3]);
            }
            $qlength=length($qstring);
            if ($_[4]>=($qlength*2)) {
                $wantX2=$qlength*2;
            } else {
                $wantX2=$_[4];
            }
            $wantX3=$qlength*2-$wantX2;
            if (($wantX2>0)and(defined($StrUp{$temp}))) {
                $X2string=substr($StrUp{$temp},0,$wantX2);
            } else {
                $X2string='';
            }
            if (($wantX3>0)and(defined($DownStream{$temp}))) {
                $X3string=substr($DownStream{$temp},0,$wantX3);
            } else {
                $X3string='';
            }
            $rstring=join '',($X2string,$X3string);
            $rstring=&complementary($rstring);
            $rlength=length($rstring);
            for ($i=1;$i<=$qlength;$i++) { ### generate @queryseq in reverse direction
                $queryseq[$i]=substr($qstring,$qlength-$i,1);
            }
            for ($i=1;$i<=$rlength;$i++) { ### generate @refseq in reverse direction
                $refseq[$i]=substr($rstring,$rlength-$i,1);
            }
        }
    } else {
        if (5==$_[2]) { ### this check is for a potential X1
            $qstring=substr($qstring,$_[3],$flankcheck_length);
            $qlength=length($qstring);
            if ($_[4]>=($qlength*2)) {
                $wantX2=$qlength*2;
            } else {
                $wantX2=$_[4];
            }
            $wantX3=$qlength*2-$wantX2;
            if (($wantX2>0)and(defined($StrUp{$temp}))) {
                $X2string=substr($StrUp{$temp},0,$wantX2);
            } else {
                $X2string='';
            }
            if (($wantX3>0)and(defined($DownStream{$temp}))) {
                $X3string=substr($DownStream{$temp},0,$wantX3);
            } else {
                $X3string='';
            }
            $rstring=join '',($X2string,$X3string);
            $rlength=length($rstring);
            for ($i=1;$i<=$qlength;$i++) { ### generate @queryseq in forward direction
                $queryseq[$i]=substr($qstring,$i-1,1);
            }
            for ($i=1;$i<=$rlength;$i++) { ### generate @refseq in forward direction
                $refseq[$i]=substr($rstring,$i-1,1);
            }
        } else { ### this check is for a potential X3
            if (($qlength-$_[5])>=$flankcheck_length) {
                $qstring=substr($qstring,$qlength-$_[5]-$flankcheck_length,$flankcheck_length);
            } else {
                $qstring=substr($qstring,0,$qlength-$_[5]);
            }
            $qlength=length($qstring);
            if ($_[4]>=($qlength*2)) {
                $wantX2=$qlength*2;
            } else {
                $wantX2=$_[4];
            }
            $wantX1=$qlength*2-$wantX2;
            if (($wantX2>0)and(defined($StrDown_comple{$temp}))) {
                $X2string=substr($StrDown_comple{$temp},0,$wantX2);
            } else {
                $X2string='';
            }
            if (($wantX1>0)and(defined($UpStream_comple{$temp}))) {
                $X1string=substr($UpStream_comple{$temp},0,$wantX1);
            } else {
                $X1string='';
            }
            $rstring=join '',($X2string,$X1string);
            $rstring=&complementary($rstring);
            $rlength=length($rstring);
            for ($i=1;$i<=$qlength;$i++) { ### generate @queryseq in reverse direction
                $queryseq[$i]=substr($qstring,$qlength-$i,1);
            }
            for ($i=1;$i<=$rlength;$i++) { ### generate @refseq in reverse direction
                $refseq[$i]=substr($rstring,$rlength-$i,1);
            }
        }
    }

    my $gapscore=-2;
    my $misscore=-1;
    my $matchscore=0; ### N from @refseq will = $matchscore, otherwise N from @queryseq will = $misscore*1/2
    my $x=0;
    my $y=0;
    my @alignscore=();
    my $maxy=0;
    $alignscore[0][0]=min($qlength,$rlength);
    for ($y=1;$y<=$rlength;$y++) {
        $alignscore[0][$y]=$alignscore[0][$y-1]+$gapscore;
        if ($alignscore[0][$y]<$reliablimit) {
            $maxy=$y-1;
            last;
        }
    }
    for ($x=1;$x<=$qlength;$x++) {
        $alignscore[$x][0]=$alignscore[$x-1][0]+$gapscore;
        for ($y=1;$y<=$maxy;$y++) {
            if (($queryseq[$x] eq $refseq[$y])or($refseq[$y] eq 'N')) {
                $alignscore[$x][$y]=$alignscore[$x-1][$y-1]+$matchscore;
            } elsif ($queryseq[$x] eq 'N') {
                $alignscore[$x][$y]=$alignscore[$x-1][$y-1]+$misscore*0.5;
            } else {
                $alignscore[$x][$y]=$alignscore[$x-1][$y-1]+$misscore;
            }
            if (($alignscore[$x][$y-1]+$gapscore)>$alignscore[$x][$y]) {
                $alignscore[$x][$y]=$alignscore[$x][$y-1]+$gapscore;
            }
            if (($alignscore[$x-1][$y]+$gapscore)>$alignscore[$x][$y]) {
                $alignscore[$x][$y]=$alignscore[$x-1][$y]+$gapscore;
            }
        }
        if ($y<=$rlength) {
            if (($queryseq[$x] eq $refseq[$y])or($refseq[$y] eq 'N')) {
                $alignscore[$x][$y]=$alignscore[$x-1][$y-1]+$matchscore;
            } elsif ($queryseq[$x] eq 'N') {
                $alignscore[$x][$y]=$alignscore[$x-1][$y-1]+$misscore*0.5;
            } else {
                $alignscore[$x][$y]=$alignscore[$x-1][$y-1]+$misscore;
            }
            if (($alignscore[$x][$y-1]+$gapscore)>$alignscore[$x][$y]) {
                $alignscore[$x][$y]=$alignscore[$x][$y-1]+$gapscore;
            }
        } else {
            $y--;
        }
        while ($alignscore[$x][$y]>=$reliablimit) {
            $alignscore[$x][$y+1]=$alignscore[$x][$y]+$gapscore;
            $y++;
        }
        $maxy=0;
        for ($i=$y;$i>=0;$i--) {
            if ($alignscore[$x][$i]>=$reliablimit) {
                $maxy=$i;
                last;
            }
        }
        if ($maxy>$rlength) {
            $maxy=$rlength;
        }
    }
    $x--;
    for ($y=0;$y<=$maxy;$y++) {
        if ($alignscore[$x][$y]>$result) {
            $result=$alignscore[$x][$y];
        }
    }
    return ($result);
}

sub PrepareResult { ### process @RawResult[$_[0]] and generate @AlignResult
    @AlignResult=(); ### (pairid1,X1,X2,X3,X1,X2,X3,pairid2,X1,X2,X3,X1,X2,X3, ..)
    ### pairid can duplicate, but same pairids will be together
    ### X1=length mapped to upflank, X2=length of repeat sequence, X3=length mapped to downflank
    ### first [X1 X2 X3] from reads1, second [X1 X2 X3] from reads2

    my @Sol=(); ### $Reads1/2 info processed from @RawResult[$_[0]], by replacing low reability X1/3s with X2s
    my $allowrange=2;
    my @temptrust=(); ### store all trustful X2 length found by all reads for this STRid
    my %toolong=(); ### when repeat length is too long (determined by $toolongminlen), trust when total unique number >=$toolongminnum
    my %toolongcheck=();
    my $toolongminlen=max(3*$checklimit,15);
    my $toolongminnum=2;
    my @tempdata=();
    my $i=0;
    my $j=0;
    my $k=0;
    my $m=0;
    my $trust='no';
    for ($j=$RawResultMax[$_[0]];$j>=1;$j=$j-2) { ### get @temptrust from @RawResult[$_[0]]
        @tempdata=split / /, $RawResult[$_[0]][$j];
        for ($i=$#tempdata-2;$i>=0;$i=$i-3) {
            if (($tempdata[$i]>0)and($tempdata[$i+2]>0)) {
                if ((!($tempdata[$i]=~/\./))and(!($tempdata[$i+2]=~/\./))) {
                    push @temptrust, ($tempdata[$i+1]);
                } else {
                    $k=int($tempdata[$i]+0.5);
                    $m=int($tempdata[$i+2]+0.5);
                    if (($k+$m)<=$toolongminlen) {
                        if (!(defined($toolongcheck{"$k $m"}))) {
                            $toolongcheck{"$k $m"}=1;
                            if (!(defined($toolong{$tempdata[$i+1]}))) {
                                $toolong{$tempdata[$i+1]}=1;
                            } else {
                                $toolong{$tempdata[$i+1]}++;
                            }
                        }
                    }
                }
            }
        }
    }
    foreach (keys(%toolong)) { ### add more elements to @temptrust for those with too long repeats to have long flank support
        if ($toolong{$_}>=$toolongminnum) {
            push @temptrust, ($_);
        }
    }
    for ($j=$RawResultMax[$_[0]];$j>=1;$j=$j-2) { ### get @Sol from @RawResult[$_[0]]
        push @Sol, ($RawResult[$_[0]][$j-1]);
        @tempdata=split / /, $RawResult[$_[0]][$j];
        for ($i=$#tempdata-2;$i>=0;$i=$i-3) {
            if (($tempdata[$i]=~/\./)or($tempdata[$i+2]=~/\./)) { ### mark found, check for support from other reads
                if ((0==$tempdata[$i])or(0==$tempdata[$i+2])) { ### open-end cannot have support, skip check
                    $trust='skip'; ### "yes" to accept original, "no" to replace with X2, "skip" to wait for further check by %Open
                } else {
                    $trust='no';
                    for ($k=$#temptrust;$k>=0;$k--) {
                        if ((abs($tempdata[$i+1]-$temptrust[$k]))<=$allowrange) {
                            $trust='yes';
                            last;
                        }
                    }
                }
                if ($trust eq 'no') { ### no support found, replace marked X1/3 with X2
                    if ($tempdata[$i]=~/\./) {
                        $tempdata[$i+1]=$tempdata[$i+1]+int($tempdata[$i]+0.5);
                        $tempdata[$i]=0;
                    }
                    if ($tempdata[$i+2]=~/\./) {
                        $tempdata[$i+1]=$tempdata[$i+1]+int($tempdata[$i+2]+0.5);
                        $tempdata[$i+2]=0;
                    }
                } elsif ($trust eq 'yes') { ### support found, accept original X1/3
                    $tempdata[$i]=int($tempdata[$i]+0.5);
                    $tempdata[$i+2]=int($tempdata[$i+2]+0.5);
                }
            }
        }
        push @Sol, ($tempdata[0]);
        for ($i=1;$i<=$#tempdata;$i++) {
            $Sol[$#Sol].=" $tempdata[$i]";
        }
    }

    my @Tru=(); ### $Reads1/2 info processed from @Sol, adding header 'X ' to realignment, by discarding X2 only results without open end support
    my $Open=0; ### value=1 means found open end aligments to support X2 only results & open-end results with short flanks
    for ($j=$#Sol;$j>=1;$j=$j-2) { ### get $Open from @Sol
        @tempdata=split / /, $Sol[$j];
        for ($i=$#tempdata-2;$i>=0;$i=$i-3) {
            if (($tempdata[$i+1]>0)and(($tempdata[$i]+$tempdata[$i+2])>0)and(0==($tempdata[$i]*$tempdata[$i+2]))) {
                if (($tempdata[$i]=~/\./)or($tempdata[$i+2]=~/\./)) { ### untrustful open reads don't count for support
                    next;
                }
                if (($tempdata[$i]>0)and($tempdata[$i]>$openlimit1)and($tempdata[$i]<=$openlimit2)) {
                    $Open=1;
                } elsif (($tempdata[$i+2]>0)and($tempdata[$i+2]>$openlimit1)and($tempdata[$i+2]<=$openlimit2)) {
                    $Open=1;
                }
            }
        }
    }
    for ($j=$#Sol;$j>=1;$j=$j-2) { ### get @Tru from @Sol
        push @Tru, ($Sol[$j-1]);
        push @Tru, ('X');
        @tempdata=split / /, $Sol[$j];
        for ($i=$#tempdata-2;$i>=0;$i=$i-3) {
            if (($tempdata[$i]=~/\./)or($tempdata[$i+2]=~/\./)) { ### if open-end results were skipped, finish here
                if (0==$Open) { ### no support found, replace marked X1/3 with X2
                    if ($tempdata[$i]=~/\./) {
                        $tempdata[$i+1]=$tempdata[$i+1]+int($tempdata[$i]+0.5);
                        $tempdata[$i]=0;
                    }
                    if ($tempdata[$i+2]=~/\./) {
                        $tempdata[$i+1]=$tempdata[$i+1]+int($tempdata[$i+2]+0.5);
                        $tempdata[$i+2]=0;
                    }
                } else { ### support found, accept original X1/3
                    $tempdata[$i]=int($tempdata[$i]+0.5);
                    $tempdata[$i+2]=int($tempdata[$i+2]+0.5);
                }
            }
            if (($tempdata[$i]>0)and(0==$tempdata[$i+2])and($tempdata[$i]<=$force_replace)) {
                $tempdata[$i+1]+=$tempdata[$i];
                $tempdata[$i]=0;
            } elsif (($tempdata[$i+2]>0)and(0==$tempdata[$i])and($tempdata[$i+2]<=$force_replace)) {
                $tempdata[$i+1]+=$tempdata[$i+2];
                $tempdata[$i+2]=0;
            }
            $Tru[$#Tru].=" $tempdata[$i] $tempdata[$i+1] $tempdata[$i+2]"; ### keep X2 only results even if not supported by other open end results
        }
    }

    my %Map1=(); ### $Reads1 info processed from @Tru, keep header 'X ' to realignment, by discarding [X1 X2 X3]s within any other [X1 X2 X3], and indexed by {pairid}{+/-}
    my %Map2=(); ### $Reads2 info processed from @Tru, keep header 'X ' to realignment, by discarding [X1 X2 X3]s within any other [X1 X2 X3], and indexed by {pairid}{+/-}
    my @tempindex=();
    my $temp='';
    for ($j=$#Tru;$j>=1;$j=$j-2) { ### get %Map1/2 from @Tru
        if ($Tru[$j] ne 'X') {
            @tempdata=split / /, $Tru[$j];
            for ($i=$#tempdata-2;$i>=1;$i=$i-3) { ### get rid of included [X1 X2 X3]s
                for ($k=$#tempdata-2;$k>=1;$k=$k-3) {
                    if (($k!=$i)and($tempdata[$i]<=$tempdata[$k])and($tempdata[$i+2]<=$tempdata[$k+2])) {
                        $tempdata[$i]=-1;
                        last;
                    }
                }
            }
            $temp='X';
            for ($i=$#tempdata-2;$i>=1;$i=$i-3) {
                if ($tempdata[$i]!=-1) {
                    $temp.=" $tempdata[$i] $tempdata[$i+1] $tempdata[$i+2]";
                }
            }
            @tempindex=split / /, $Tru[$j-1];
            if (1==$tempindex[2]) {
                $Map1{$tempindex[1]}{$tempindex[0]}=$temp;
            } else {
                $Map2{$tempindex[1]}{$tempindex[0]}=$temp;
            }
        }
    }

    my $ref1='';
    my $ref2='';
    my $pairid='';
    my @tempdata1=();
    my @tempdata2=();
    my $pureflank1=0;
    my $pureflank2=0;
    foreach $pairid (keys(%Map1)) { ### find paired realignment and save results to %AlignResult
        foreach $ref1 (keys(%{$Map1{$pairid}})) {
            if (($pairtype eq 'fr')or($pairtype eq 'rf')) {
                if ($ref1 eq '-') {
                    $ref2='+';
                } else {
                    $ref2='-';
                }
            } else {
                $ref2=$ref1;
            }
            if (defined($Map2{$pairid}{$ref2})) {
                @tempdata1=split / /, $Map1{$pairid}{$ref1};
                @tempdata2=split / /, $Map2{$pairid}{$ref2};
                for ($i=$#tempdata1-2;$i>=1;$i=$i-3) {
                    if ((0==$tempdata1[$i+1])and((0==$tempdata1[$i])or(0==$tempdata1[$i+2]))) {
                        $pureflank1=1;
                    } else {
                        $pureflank1=0;
                    }
                    for ($j=$#tempdata2-2;$j>=1;$j=$j-3) {
                        if ((0==$tempdata2[$j+1])and((0==$tempdata2[$j])or(0==$tempdata2[$j+2]))) {
                            $pureflank2=1;
                        } else {
                            $pureflank2=0;
                        }
                        if ((0==$pureflank1)or(0==$pureflank2)) { ### skip if both fully mapped to flank regions
                            $temp=&getPairtype($ref1,$tempdata1[$i],$tempdata1[$i+1],$tempdata1[$i+2],$ref2,$tempdata2[$j],$tempdata2[$j+1],$tempdata2[$j+2]);
                            if ($temp=~/\Q$pairtype\E/) {
                                push @AlignResult, ($pairid,$tempdata1[$i],$tempdata1[$i+1],$tempdata1[$i+2],$tempdata2[$j],$tempdata2[$j+1],$tempdata2[$j+2]);
                            }
                        }
                    }
                }
            }
        }
    }
}

### read parameters ###
my $i=0;
@Para=();
### Para[0] main result file from _findseq.pl
### Para[1] unsorted .sam or .bam (arranged by reads ID, not by mapping position) from bwa using flanking as reference
### Para[2] output file
### Para[3] maximum length to search for flanks, suggestion -T/-A +20 (-T and -A are bwa parameters with default or customized values) for minimum value
### However it may miss some reads due to bwa output (tend to softclip mismatches near reads ends). Set it to very large value will enhance the accuracy, but lower running speed. May update in the future
### Para[4] flank <= this will go into reliability check when generating %Ful1/2, others will be trusted, suggestion =15, 0 means all will be trusted
### Para[5] if there's also available X2 replacement when generating %Ful1/2, reliability >= this will be accepted or marked if flank shorter than Para[11], < this will be replaced, suggestion 0.8, <=0 means all will be trusted, 1 means very strict, >1 means none will be trusted. Not apply when no X2 replacement and will accept anyway
### Para[6] pairtype "fr"/"rr"/"rf"/"ff"
### Para[7] length of fake softclip at repeat side, as bwa tends to softclip if a mismatch happens near reads ends, suggestion =10, correlated with bwa settings (10 when default settings is a good option)
### Para[8] length of fake softclip at flank side, suggestion =5
### Para[9] ratio limit to search for possible repeat by k-mer, suggestion = 0.8-2*(-T/-A)/len, <=0 all reads will be possible, >1 no reads will be possible
### Para[10] allow search by k-mer for more possible units with counts lower than best by this*unitlength, suggestion =0
### Para[11] flank <= this will be further examined, may be replaced with X2 if not trustful when replacement is available (suggestion =3)
### a marked flank will be replaced by X2 in &PrepareResult if: (1)it's not open end, and (2)not supported by other unmarked results
### Para[12] if available, open-end results with flank <= this length will be forcely replaced by X2 (suggestion =0)
### Para[13] seed length for flank search, allowing 1 mis/gap but still may reduce accuracy, <=0 will skip this (suggestion =7)
### Para[14] ratio (repeat/fullreadslength) limit to initiate repeat-guided realignment from _Repeat/_Fake mapping, suggestion = 1-2*(-T/-A)/len
### Para[15] ratio modification for Para[14], suggestion =0.25
### _Repeat/_Fake mapping ratio will use Para[14]-Para[15], _Repeat/_Fake mapping + SW extension ratio will use Para[14]
### Para[16] nt allowance to search flank from _Fake mapping ends, suggestion = 5
### Para[17] "yes"/"no" to apply k-mer method to search more repeat candidate, "no" will disable Para[9]&[10]
### Para[18] upstream SW extension limit for _Repeat/_Fake mapped reads, suggestion =(-T/-A)+20
### Para[19] downstream SW extension limit for _Repeat/_Fake mapped reads, suggestion =reads length to allow 3' misreading, or (-T/-A)+20 to discard misread reads
### Para[20] total process number
### Para[21] the number of this run
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-h')or($ARGV[$i] eq '--help')) {
        print "Instruction:\n";
        print "\tThe Realigner module realigns the remapped pairs to each STR locus, and generates a plain text file indicating how each pair is realigned.\n";
        print "\tThe running speed is affected by both library size and number of STRs. For heavy duty work, the Realigner module allows multi-processing mode.\n";
        print "\tSee README.txt Step (1.2) and (2) for input requirements, and Step (3) for more details of the output.\n";
        print "Usage Sample (single process, read length = ~ 150bp):\n";
        print "\tperl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt>\n";
        print "Usage Sample (single process, read length = ~ 100bp):\n";
        print "\tperl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt> --repeatm 0.15 --repeate 0.4 -d 100\n";
        print "Usage Sample (multiple parallel processes):\n";
        print "\tperl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt> -m 3 -n 1\n";
        print "\tperl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt> -m 3 -n 2\n";
        print "\tperl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt> -m 3 -n 3\n";
        print "Options:\n";
        print "\t--help/-h\t\tPrint usage instructions\n";
        print "\t--str/-s <file>\t\tThe STR sequences obtained by the Finder module. See format details in README.txt Step (1.2)\n";
        print "\t--input/-i <file>\tUnsorted .sam/.bam obtained in the remapping step. See details in README.txt Step (2)\n";
        print "\t--output/-o <file>\tOutput file for realignment results. See format details in README.txt Step (3)\n";
        print "\t--fr/--rr/--rf/--ff\tIndicate how reads are paired if presented (forward-reversely/reverse-reversely/reverse-forwardly/forward-forwardly)\n";
        print "\t\t\t\t(Default = --fr)\n";
        print "\t--flank/-f <value>\tMaximum length to search for flankings\n";
        print "\t\t\t\tA setting of the minimum mapping length used in previous step plus 20nt or more is recommended\n";
        print "\t\t\t\tBy default bwa mem requires ~30nt as the minimum mapping length to export the hits\n";
        print "\t\t\t\t(Default = 50)\n";
        print "\t--seed <value>\t\tSeed length to search for flankings. A setting of 7 is recommended for best performance (Default = 7)\n";
        print "\t--untrust1 <value>\tLength threshold of short flankings to enter further check by looking into nearby repeat sequences (Default = 15)\n";
        print "\t--reliab/-r <value>\tReliability score threshold (0~1) to decide whether to remark short flankings as repeats\n";
        print "\t\t\t\tEffective for --untrust1 and remark only when acceptable\n";
        print "\t\t\t\t(Default = 0.8)\n";
        print "\t--untrust2 <value>\tLength threshold of short flankings to enter further check by looking into other supportive pairs\n";
        print "\t\t\t\tUnsupported short flankings may be remarked as repeats only when acceptable\n";
        print "\t\t\t\t(Default = 3)\n";
        print "\t--untrust3 <value>\tLength threshold to forcely remark single-side short flankings as repeats (Default = 0)\n";
        print "\t--clipo <value>\t\tLength threshold to ignore outwards (away from repeats) short clips on the flankings (Default = 5)\n";
        print "\t--clipi <value>\t\tLength threshold to ignore inwards (towards the repeats) short clips on the flankings\n";
        print "\t\t\t\tFor bwa mem default settings, 10 is recommended\n";
        print "\t\t\t\t(Default = 10)\n";
        print "\t--repeatm <value>\tThreshold of repeat ratio to read length (0~1) to initiate extension from remapped hits\n";
        print "\t\t\t\tA setting = 0.75-2*(minimum mapping length)/readlength is recommended\n";
        print "\t\t\t\t(Default = 0.35)\n";
        print "\t--repeate <value>\tThreshold of repeat ratio to read length (0~1) to accept the repeats\n";
        print "\t\t\t\tA setting = 1-2*(minimum mapping length)/readlength is recommended\n";
        print "\t\t\t\t(Default = 0.6)\n";
        print "\t--back/-b <value>\tLength to step back from remapped ends of reads hitting on repeat references, to search for flankings (Default = 5)\n";
        print "\t--up/-u <value>\t\tMaximum length of periodic Smith-Waterman extension of repeats towards 5' end of reads\n";
        print "\t\t\t\tA setting of the minimum mapping length used in previous step plus 20nt is recommended\n";
        print "\t\t\t\t(Default = 50)\n";
        print "\t--down/-d <value>\tMaximum length of periodic Smith-Waterman extension of repeats towards 3' end of reads\n";
        print "\t\t\t\tA setting of read length is recommended considering sequencing errors\n";
        print "\t\t\t\t(Default = 150)\n";
        print "\t--kratio <value>\tThreshold of repeat ratio to read length (0~1) to initiate search for potential repeat units by k-mer method\n";
        print "\t\t\t\tThe k-mer search will cost big running time. Apply it when the STR list is small and a deep realignment is required\n";
        print "\t\t\t\tTo apply, a setting = 0.8-2*(minimum mapping length)/readlength is recommended\n";
        print "\t\t\t\t(Default = 0.4, but will be disabled if neither --kratio nor --kcandi is present)\n";
        print "\t--kcandi <value>\tBesides the best potential repeat unit for each size, consider more repeat units as candidates\n";
        print "\t\t\t\tFor each unit size, candidates with frequencies close to the best within a range determined by this setting will be considered\n";
        print "\t\t\t\t(Default = 0, but will be disabled if neither --kratio nor --kcandi is present)\n";
        print "\t--multi/-m <value>\tTotal number of processes to run parallel (Default = 1)\n";
        print "\t--num/-n <value>\tThe order of this run in the parallel processes (Default = 1, but will be disabled if --multi/-m is not present)\n";
        exit;
    }
}
$Para[0]='';
$Para[1]='';
$Para[2]='';
$Para[3]=50;
$Para[4]=15;
$Para[5]=0.8;
$Para[6]='fr';
$Para[7]=10;
$Para[8]=5;
$Para[9]=0.4;
$Para[10]=0;
$Para[11]=3;
$Para[12]=0;
$Para[13]=7;
$Para[14]=0.6;
$Para[15]=0.35;
$Para[16]=5;
$Para[17]='no';
$Para[18]=50;
$Para[19]=150;
$Para[20]=1;
$Para[21]=1;
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-s')or($ARGV[$i] eq '--str')) {
        $i++;
        $Para[0]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-i')or($ARGV[$i] eq '--input')) {
        $i++;
        $Para[1]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-o')or($ARGV[$i] eq '--output')) {
        $i++;
        $Para[2]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-f')or($ARGV[$i] eq '--flank')) {
        $i++;
        $Para[3]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--seed') {
        $i++;
        $Para[13]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--fr') {
        $Para[6]='fr';
    } elsif ($ARGV[$i] eq '--rf') {
        $Para[6]='rf';
    } elsif ($ARGV[$i] eq '--rr') {
        $Para[6]='rr';
    } elsif ($ARGV[$i] eq '--ff') {
        $Para[6]='ff';
    } elsif ($ARGV[$i] eq '--untrust1') {
        $i++;
        $Para[4]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-r')or($ARGV[$i] eq '--reliab')) {
        $i++;
        $Para[5]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--untrust2') {
        $i++;
        $Para[11]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--untrust3') {
        $i++;
        $Para[12]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--clipo') {
        $i++;
        $Para[8]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--clipi') {
        $i++;
        $Para[7]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--repeatm') {
        $i++;
        $Para[15]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--repeate') {
        $i++;
        $Para[14]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-b')or($ARGV[$i] eq '--back')) {
        $i++;
        $Para[16]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-u')or($ARGV[$i] eq '--up')) {
        $i++;
        $Para[18]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-d')or($ARGV[$i] eq '--down')) {
        $i++;
        $Para[19]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-m')or($ARGV[$i] eq '--multi')) {
        $i++;
        $Para[20]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-n')or($ARGV[$i] eq '--num')) {
        $i++;
        $Para[21]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--kratio') {
        $i++;
        $Para[9]=$ARGV[$i];
        $Para[17]='yes';
    } elsif ($ARGV[$i] eq '--kcandi') {
        $i++;
        $Para[10]=$ARGV[$i];
        $Para[17]='yes';
    }
}
$Para[15]=$Para[14]-$Para[15];
if (1==$Para[20]) {
    $Para[21]=1;
}
if (($Para[0] eq '')or($Para[1] eq '')or($Para[2] eq '')) {
    print ("STR/Input/Output not appointed. Stop running.\n");
    print "Please use --help or -h to see instructions.\n";
    exit;
}

### main ###
use List::Util qw/max min/;
$search_length=max(int($Para[3]),0);
$flankcheck_length=15;
$flankgapscore=-5; ### used in &PrepareSTR and &flank_extension
$flankmisscore=-3; ### used in &PrepareSTR and &flank_extension
$flankmatchscore=0; ### used in &flank_extension
$flanklimit=int($Para[4]); ### used when generating %Ful1/2, =0 means no limit
$reliablimit=int($Para[5]*$flankcheck_length*10+0.5)/10; ### used when generating %Ful1/2, =0 means no limit
$pairtype=$Para[6]; ### "fr"/"ff"/"rf"/"rr"
$fakesoft=int($Para[7]); ### softclip <= this will be considered as fake
$fakesoftin=int($Para[8]); ### softclip <= this will be considered as fake
$RepeatPercentLimit=$Para[9]; ### used in &PossibleRepeat
$search_depth=max(int($Para[10]),0); ### used in &PossibleRepeat
$checklimit=max(int($Para[11]),0); ## used when generating %Ful1/2, =0 means no limit
$force_replace=int($Para[12]); ### used in &PrepareResult
$flankseed=int($Para[13]); ### used for repeat-guide realignment when same unit flank candidate number >20
$RepeatPercentLimitMap=$Para[14]; ### used in &getRawRepeat
$Ratiomodify=$Para[15];
$flanklimit_bymapping=max(int($Para[16]),0); ### used when generating %Ful1/2
$untrustmap=max(10,$flanklimit_bymapping); ### used in &getRawRepeat
$kmer_switch=$Para[17]; ### decide whether to use &PossibleRepeat
$upexlimit=int($Para[18]); ### used in &up_extension_fix
$downexlimit=int($Para[19]); ### used in &down_extension_fix
$totalnumber=$Para[20];
$thisnumber=$Para[21];
&PrepareSTR($Para[0]);
&PrepareFlag;
&PrepareDymscore;
$allowance_score=$mis_score; ### will -$dymatch_score when use, to allow sequencing errors or SNPs in repeats from reads
my $line='';
my @tempdata=();
my $IDnow='';
my $countlimit=10000000;
my $count=$countlimit;
my $totalcount=0;
my $temp='';
my $temppri='';
my $temppair='';
my $tempstrand='';
my $tempmap='';
@Lines1=(); ### incomplete lines of the same reads from .sam, processed by ReAligner, excluding lines for repeat mapping/unmapping information
@Lines2=(); ### incomplete lines of the same reads paired with @Lines1
### @Lines1/2: multiple (+/-, refID, U/D, mapposition, CIGAR), all copied/processed from .sam (U=ref_Upstream, D=ref_Downstream)
@Lines1Repeat=(); ### incomplete lines of the same reads from .sam, processed by ReAligner, only lines for _Repeat mapping
@Lines2Repeat=(); ### incomplete lines of the same reads paired with @Lines1
### @Lines1/2Repeat: multiple (+/-, refID, CIGAR), all copied/processed from .sam
@Lines1Fake=(); ### incomplete lines of the same reads from .sam, processed by ReAligner, only lines for _Fake mapping
@Lines2Fake=(); ### incomplete lines of the same reads paired with @Lines1
### @Lines1/2Fake: multiple (+/-, transformed unit of reads, CIGAR), all copied/processed from .sam
$Reads1=''; ### the original reads sequence of @Lines1, processed from .sam, will be processed by ReAligner
$Reads2=''; ### the original reads sequence of @Lines2
$Reads1rev=''; ### reverse of $Reads1, used in &getRawRepeat
$Reads2rev=''; ### reverse of $Reads2, used in &getRawRepeat
%RefLength=(); ### length of each reference, copy from @SQ lines
if ($Para[1]=~/\.sam$/) {
    open (INPUT_DATA, "<$Para[1]") or die "Couldn't open: $!";
    if (-s INPUT_DATA) {
        while (1) {
            chomp ($line=<INPUT_DATA>);
            if ($line=~/^\@SQ/) {
                @tempdata=split /\s+/, $line;
                $tempdata[1]=~/SN:/;
                $temp=$';
                if ($temp=~/Stream$/) {
                    $tempdata[2]=~/LN:/;
                    $RefLength{$temp}=$';
                }
            } elsif (!($line=~/^\@/)) {
                @tempdata=split /\s+/, $line;
                if ($tempdata[0] ne $IDnow) {
                    if ($IDnow ne '') {
                        if (($Reads1 eq '')or($Reads2 eq '')) {
                            print ("Warning: $IDnow reads not successfully recovered\n");
                        } else {
                            $totalcount++;
                            $count--;
                            if ((($totalcount%$totalnumber)==$thisnumber)or((($totalcount%$totalnumber)==0)and($thisnumber==$totalnumber))) {
                                &ReAligner($IDnow);
                            }
                            if (0==$count) {
                                print ("$totalcount pairs processed...\n");
                                $count=$countlimit;
                            }
                        }
                    }
                    @Lines1=();
                    @Lines2=();
                    @Lines1Repeat=();
                    @Lines2Repeat=();
                    @Lines1Fake=();
                    @Lines2Fake=();
                    $Reads1='';
                    $Reads2='';
                    $IDnow=$tempdata[0];
                }
                ($temppri,$temppair,$tempstrand,$tempmap)=($FlagPriority[$tempdata[1]],$FlagPair[$tempdata[1]],$FlagStrand[$tempdata[1]],$FlagMap[$tempdata[1]]);
                if (1==$temppair) {
                    if ($tempmap eq 'y') {
                        if ($tempdata[2]=~/_UpStream$/) {
                            if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                                push @Lines1, ($tempstrand,$`,'U',$tempdata[3],$tempdata[5]);
                            }
                        } elsif ($tempdata[2]=~/_DownStream$/) {
                            if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                                push @Lines1, ($tempstrand,$`,'D',$tempdata[3],$tempdata[5]);
                            }
                        } elsif ($tempdata[2]=~/_Repeat$/) {
                            if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                                push @Lines1Repeat, ($tempstrand,$`,$tempdata[5]);
                            }
                        } elsif ($tempdata[2]=~/_Fake$/) {
                            if ($tempstrand eq '-') {
                                $temp=&complementary($`);
                            } else {
                                $temp=$`;
                            }
                            $temp=&TransformSTR($temp);
                            if ($temp ne '') { ### this allows processing sub STRlist from full mapped .sam/.bam
                                push @Lines1Fake, ($tempstrand,$temp,$tempdata[5]);
                            }
                        } else {
                            print ("Warning: $tempdata[2] unexpected reference name from mapping results\n");
                        }
                    }
                    if ($temppri eq 'p') {
                        if ($Reads1 ne '') {
                            print ("Warning: $IDnow has multiple primary alignment\n");
                        } else {
                            $Reads1=&standardizeDNA($tempdata[9]);
                            if ($tempstrand eq '-') {
                                $Reads1=&complementary($Reads1);
                            }
                        }
                    }
                } elsif (2==$temppair) {
                    if ($tempmap eq 'y') {
                        if ($tempdata[2]=~/_UpStream$/) {
                            if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                                push @Lines2, ($tempstrand,$`,'U',$tempdata[3],$tempdata[5]);
                            }
                        } elsif ($tempdata[2]=~/_DownStream$/) {
                            if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                                push @Lines2, ($tempstrand,$`,'D',$tempdata[3],$tempdata[5]);
                            }
                        } elsif ($tempdata[2]=~/_Repeat$/) {
                            if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                                push @Lines2Repeat, ($tempstrand,$`,$tempdata[5]);
                            }
                        } elsif ($tempdata[2]=~/_Fake$/) {
                            if ($tempstrand eq '-') {
                                $temp=&complementary($`);
                            } else {
                                $temp=$`;
                            }
                            $temp=&TransformSTR($temp);
                            if ($temp ne '') { ### this allows processing sub STRlist from full mapped .sam/.bam
                                push @Lines2Fake, ($tempstrand,$temp,$tempdata[5]);
                            }
                        } else {
                            print ("Warning: $tempdata[2] unexpected reference name from mapping results\n");
                        }
                    }
                    if ($temppri eq 'p') {
                        if ($Reads2 ne '') {
                            print ("Warning: $IDnow has multiple primary alignment\n");
                        } else {
                            $Reads2=&standardizeDNA($tempdata[9]);
                            if ($tempstrand eq '-') {
                                $Reads2=&complementary($Reads2);
                            }
                        }
                    }
                } else {
                    print ("Warning: $IDnow unpaired, possibly because alignment was not done in pair-end mode\n");
                }
            }
            if (eof INPUT_DATA) {
                last;
            }
        }
    }
    if (($Reads1 eq '')or($Reads2 eq '')) {
        print ("Warning: $IDnow reads not successfully recovered\n");
    } else {
        $totalcount++;
        $count--;
        if ((($totalcount%$totalnumber)==$thisnumber)or((($totalcount%$totalnumber)==0)and($thisnumber==$totalnumber))) {
            &ReAligner($IDnow);
        }
        if (0==$count) {
            print ("$totalcount pairs processed...\n");
        }
    }
    close (INPUT_DATA);
    print ("Total $totalcount pairs processed\n");
} elsif ($Para[1]=~/\.bam$/) {
    open INPUT_DATA, "samtools view -H $Para[1] |";
    while (<INPUT_DATA>) {
        chomp;
        @tempdata=split /\s+/;
        if ($tempdata[0]=~/^\@SQ/) {
            $tempdata[1]=~/SN:/;
            $temp=$';
            if ($temp=~/Stream$/) {
                $tempdata[2]=~/LN:/;
                $RefLength{$temp}=$';
            }
        }
    }
    close (INPUT_DATA);
    open INPUT_DATA, "samtools view $Para[1] |";
    while (<INPUT_DATA>) {
        chomp;
        @tempdata=split /\s+/;
        if ($tempdata[0] ne $IDnow) {
            if ($IDnow ne '') {
                if (($Reads1 eq '')or($Reads2 eq '')) {
                    print ("Warning: $IDnow reads not successfully recovered\n");
                } else {
                    $totalcount++;
                    $count--;
                    if ((($totalcount%$totalnumber)==$thisnumber)or((($totalcount%$totalnumber)==0)and($thisnumber==$totalnumber))) {
                        &ReAligner($IDnow);
                    }
                    if (0==$count) {
                        print ("$totalcount pairs processed...\n");
                        $count=$countlimit;
                    }
                }
            }
            @Lines1=();
            @Lines2=();
            @Lines1Repeat=();
            @Lines2Repeat=();
            @Lines1Fake=();
            @Lines2Fake=();
            $Reads1='';
            $Reads2='';
            $IDnow=$tempdata[0];
        }
        ($temppri,$temppair,$tempstrand,$tempmap)=($FlagPriority[$tempdata[1]],$FlagPair[$tempdata[1]],$FlagStrand[$tempdata[1]],$FlagMap[$tempdata[1]]);
        if (1==$temppair) {
            if ($tempmap eq 'y') {
                if ($tempdata[2]=~/_UpStream$/) {
                    if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                        push @Lines1, ($tempstrand,$`,'U',$tempdata[3],$tempdata[5]);
                    }
                } elsif ($tempdata[2]=~/_DownStream$/) {
                    if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                        push @Lines1, ($tempstrand,$`,'D',$tempdata[3],$tempdata[5]);
                    }
                } elsif ($tempdata[2]=~/_Repeat$/) {
                    if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                        push @Lines1Repeat, ($tempstrand,$`,$tempdata[5]);
                    }
                } elsif ($tempdata[2]=~/_Fake$/) {
                    if ($tempstrand eq '-') {
                        $temp=&complementary($`);
                    } else {
                        $temp=$`;
                    }
                    $temp=&TransformSTR($temp);
                    if ($temp ne '') { ### this allows processing sub STRlist from full mapped .sam/.bam
                        push @Lines1Fake, ($tempstrand,$temp,$tempdata[5]);
                    }
                } else {
                    print ("Warning: $tempdata[2] unexpected reference name from mapping results\n");
                }
            }
            if ($temppri eq 'p') {
                if ($Reads1 ne '') {
                    print ("Warning: $IDnow has multiple primary alignment\n");
                } else {
                    $Reads1=&standardizeDNA($tempdata[9]);
                    if ($tempstrand eq '-') {
                        $Reads1=&complementary($Reads1);
                    }
                }
            }
        } elsif (2==$temppair) {
            if ($tempmap eq 'y') {
                if ($tempdata[2]=~/_UpStream$/) {
                    if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                        push @Lines2, ($tempstrand,$`,'U',$tempdata[3],$tempdata[5]);
                    }
                } elsif ($tempdata[2]=~/_DownStream$/) {
                    if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                        push @Lines2, ($tempstrand,$`,'D',$tempdata[3],$tempdata[5]);
                    }
                } elsif ($tempdata[2]=~/_Repeat$/) {
                    if (defined($StrUnit{$`})) { ### this allows processing sub STRlist from full mapped .sam/.bam
                        push @Lines2Repeat, ($tempstrand,$`,$tempdata[5]);
                    }
                } elsif ($tempdata[2]=~/_Fake$/) {
                    if ($tempstrand eq '-') {
                        $temp=&complementary($`);
                    } else {
                        $temp=$`;
                    }
                    $temp=&TransformSTR($temp);
                    if ($temp ne '') { ### this allows processing sub STRlist from full mapped .sam/.bam
                        push @Lines2Fake, ($tempstrand,$temp,$tempdata[5]);
                    }
                } else {
                    print ("Warning: $tempdata[2] unexpected reference name from mapping results\n");
                }
            }
            if ($temppri eq 'p') {
                if ($Reads2 ne '') {
                    print ("Warning: $IDnow has multiple primary alignment\n");
                } else {
                    $Reads2=&standardizeDNA($tempdata[9]);
                    if ($tempstrand eq '-') {
                        $Reads2=&complementary($Reads2);
                    }
                }
            }
        } else {
            print ("Warning: $IDnow unpaired, possibly because alignment was not done in pair-end mode\n");
        }
    }
    if (($Reads1 eq '')or($Reads2 eq '')) {
        print ("Warning: $IDnow reads not successfully recovered\n");
    } else {
        $totalcount++;
        $count--;
        if ((($totalcount%$totalnumber)==$thisnumber)or((($totalcount%$totalnumber)==0)and($thisnumber==$totalnumber))) {
            &ReAligner($IDnow);
        }
        if (0==$count) {
            print ("$totalcount pairs processed...\n");
        }
    }
    close (INPUT_DATA);
    print ("Total $totalcount pairs processed\n");
} else {
    print ("Error: $Para[1] is not .sam or .bam file\n");
}

my $outpath="/";
my $outfile='';
my $j=0;
@tempdata=split /\//,$Para[2];
for ($i=1;$i<$#tempdata;$i++) {
    $outpath.="$tempdata[$i]/";
}
if ($tempdata[$#tempdata]=~/\./) {
    $outfile=$`;
} else {
    $outfile=$tempdata[$#tempdata];
}
$temp=join '',($outpath,$outfile,'_hash',$thisnumber);
open (OUTPUT_RESULT, ">$temp") or die "Couldn't open: $!"; ### store sub result into file
for ($i=0;$i<=$#StrInfo;$i++) {
    print OUTPUT_RESULT ("$i $RawResultMax[$i]\n");
    for ($j=0;$j<=$RawResultMax[$i];$j++) {
        print OUTPUT_RESULT ("$RawResult[$i][$j]\n");
    }
}
close (OUTPUT_RESULT);
$temp=join '',($outpath,$outfile,'_done',$thisnumber);
open (OUTPUT_RESULT, ">$temp") or die "Couldn't open: $!"; ### generate done marker file
close (OUTPUT_RESULT);
my $donenumber=0;
my $index=0;
for ($i=1;$i<=$totalnumber;$i++) { ### check if other jobs are finished
    $temp=join '',($outpath,$outfile,'_done',$i);
    if (-e $temp) {
        $donenumber++;
    }
}
@FH=();
if ($donenumber==$totalnumber) { ### combine results to process and output, delete done marker & sub result files
    for ($i=1;$i<=$totalnumber;$i++) {
        $temp=join '',($outpath,$outfile,'_hash',$i);
        open ($FH[$i], "<$temp") or die "Couldn't open: $!";
    }
    &OutputResult($Para[2]); ### final output
    for ($i=1;$i<=$totalnumber;$i++) {
        close ($FH[$i]);
        $temp=join '',($outpath,$outfile,'_done',$i);
        system ("rm $temp");
        $temp=join '',($outpath,$outfile,'_hash',$i);
        system ("rm $temp");
    }
}
############

exit;

