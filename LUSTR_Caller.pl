#!/usr/local/bin/perl -w
### Process realignment pairs to get genotype of each STR # Bookmarks: 0,0 0,0 0,0 112,1484; CollapsedSubs: getReadsLength  shortflank_filter  forceflank_filter  reliab_fix  collect  estimate

sub getReadsLength { ### get reads length from first $_[0] pairs
    my $count=$_[0];
    my $line='';
    my %RL=();
    my $i=0;
    my @tempdata=();
    my @tempXs=();
    if (-s INPUT_DATA) {
        chomp ($line=<INPUT_DATA>);
        chomp ($line=<INPUT_DATA>);
        chomp ($line=<INPUT_DATA>);
        while (1) {
            chomp ($line=<INPUT_DATA>);
            chomp ($line=<INPUT_DATA>);
            $line=~/ /;
            for ($i=$`;$i>=1;$i--) {
                chomp ($line=<INPUT_DATA>);
                @tempdata=split / /, $line;
                $count--;
                $tempdata[1]=~s/}\{/,/;
                $tempdata[1]=~s/[{}]//g;
                @tempXs=split /,/, $tempdata[1];
                if (defined($RL{$tempXs[0]+$tempXs[1]+$tempXs[2]})) {
                    $RL{$tempXs[0]+$tempXs[1]+$tempXs[2]}++;
                } else {
                    $RL{$tempXs[0]+$tempXs[1]+$tempXs[2]}=1;
                }
                if (defined($RL{$tempXs[3]+$tempXs[4]+$tempXs[5]})) {
                    $RL{$tempXs[3]+$tempXs[4]+$tempXs[5]}++;
                } else {
                    $RL{$tempXs[3]+$tempXs[4]+$tempXs[5]}=1;
                }
                if (0==$count) {
                    last;
                }
            }
            chomp ($line=<INPUT_DATA>);
            if ((eof INPUT_DATA)or(0==$count)) {
                last;
            }
        }
    }
    seek (INPUT_DATA,0,0);
    my $result=0;
    my $maxcount=0;
    foreach (keys(%RL)) {
        if ($RL{$_}>$maxcount) {
            $maxcount=$RL{$_};
            $result=$_;
        }
    }
    return($result);
}

sub shortflank_filter { ### get %Mark to skip elements into @TopFix & modify @Realign & @Pattern to reduce noises from too short flanks
    my $i=0;
    my %Fixtotal=();
    my %Fixshort=();
    my %check=(); ### short flanks of same length will be counted only once
    my $index='';
    my $checkresult=0;
    my $addvalue=0;
    for ($i=$#Pattern;$i>=0;$i--) { ### collect information to mark fix length reads
        if (($Pattern[$i] eq '111')or($Pattern[$i] eq '101')) {
            if ((1==$Realign[$i*3])or(1==$Realign[$i*3+2])) { ### 1nt flanks likely to be noises
                $addvalue=3;
            } else {
                $addvalue=10;
            }
            if (defined($Fixtotal{$Realign[$i*3+1]})) {
                $Fixtotal{$Realign[$i*3+1]}+=$addvalue;
            } else {
                $Fixtotal{$Realign[$i*3+1]}=$addvalue;
                $Fixshort{$Realign[$i*3+1]}=0;
            }
            if (($Realign[$i*3]<=$fixflank_replace)or($Realign[$i*3+2]<=$fixflank_replace)) {
                $checkresult=0;
                if ($Realign[$i*3]<=$fixflank_replace) {
                    $index="$Realign[$i*3+1] 5 $Realign[$i*3]";
                    if (defined($check{$index})) {
                        $checkresult++;
                    }
                    $check{$index}=1;
                }
                if ($Realign[$i*3+2]<=$fixflank_replace) {
                    $index="$Realign[$i*3+1] 3 $Realign[$i*3+2]";
                    if (defined($check{$index})) {
                        $checkresult++;
                    }
                    $check{$index}=1;
                }
                if ($checkresult>0) {
                    $Fixtotal{$Realign[$i*3+1]}-=$addvalue;
                } else {
                    $Fixshort{$Realign[$i*3+1]}+=$addvalue;
                }
            }
        }
    }
    my $largest=0;
    my $numlimit=0;
    foreach (keys(%Fixtotal)) { ### mark fix length reads
        if ($_>=($ReadsLength-15)) {
            $numlimit=11;
        } elsif ($_<($ReadsLength-3*$fixflank_replace)) {
            $numlimit=10*$fixflank_number;
        } elsif ($_<($ReadsLength-2*$fixflank_replace)) {
            $numlimit=20;
        } else {
            $numlimit=11;
        }
        if (($Fixtotal{$_}==$Fixshort{$_})and($Fixtotal{$_}<$numlimit)) {
            $Mark{$_}=1;
        } elsif ($_>$largest) {
            $largest=$_;
        }
    }
    my $largestplus=$largest+$untrustexp; ### expansion possibility from reads with repeat length close to largest fix length is low

    my $largersupport=0;
    %check=();
    my $judge='';
    for ($i=$#Pattern;$i>=0;$i=$i-2) { ### check if need to modify open-end reads
        if (($Pattern[$i] eq '111')or($Pattern[$i] eq '101')or($Pattern[$i-1] eq '111')or($Pattern[$i-1] eq '101')) { ### skip if paired with fix length reads
            next;
        }
        if (($Realign[$i*3-3]<=$expflank_replace)and($Realign[$i*3-1]<=$expflank_replace)) {
            if ($Realign[$i*3-2]>$largestplus) {
                $judge='r';
            } else {
                $judge='u';
            }
        } elsif ($Realign[$i*3-2]>$largestplus) {
            if ($Pattern[$i-1] eq '110') {
                $index="5 $Realign[$i*3-3]";
            } else {
                $index="3 $Realign[$i*3-1]";
            }
            if (defined($check{$index})) {
                $judge='f';
            } else {
                $judge='e';
            }
            $check{$index}=1;
        } else {
            $judge='f';
        }
        if (($Realign[$i*3]<=$expflank_replace)and($Realign[$i*3+2]<=$expflank_replace)) {
            if ($Realign[$i*3+1]>$largestplus) {
                $judge.='r';
            } else {
                $judge.='u';
            }
        } elsif ($Realign[$i*3+1]>$largestplus) {
            if ($Pattern[$i] eq '110') {
                $index="5 $Realign[$i*3]";
            } else {
                $index="3 $Realign[$i*3+2]";
            }
            if (defined($check{$index})) {
                $judge.='f';
            } else {
                $judge.='e';
            }
            $check{$index}=1;
        } else {
            $judge.='f';
        }
        if (($judge=~/e/)or($judge eq 'fr')or($judge eq 'rf')) { ### support from long flank long expansion open-end reads or short flank open-end reads paired with trustful flanking sequences
            $largersupport++;
        }
    }

    my $target=0;
    my $j=0;
    my $judge1='';
    my $judge2='';
    if ($largersupport<$expflank_number) { ### modify untrustful open-end reads to fix length or repeat-only reads
        for ($i=$#Pattern;$i>=0;$i=$i-2) {
            if (($Realign[$i*3-2]<=$largest)or($Pattern[$i-1] eq '010')) { ### no change
                $judge1='u';
            } elsif ($Realign[$i*3-2]<=$largestplus) { ### turn to fix length
                $judge1='f';
            } elsif (($Realign[$i*3-3]<=$expflank_replace)and($Realign[$i*3-1]<=$expflank_replace)) { ### short flank open-end reads turn to repeat-only reads
                $judge1='r';
            } else { ### long flank open-end reads, no change
                $judge1='u';
            }
            if (($Realign[$i*3+1]<=$largest)or($Pattern[$i] eq '010')) { ### no change
                $judge2='u';
            } elsif ($Realign[$i*3+1]<=$largestplus) { ### turn to fix length
                $judge2='f';
            } elsif (($Realign[$i*3]<=$expflank_replace)and($Realign[$i*3+2]<=$expflank_replace)) { ### short flank open-end reads turn to repeat-only reads
                $judge2='r';
            } else { ### long flank open-end reads, no change
                $judge2='u';
            }
            if ($judge1 eq 'f') {
                $target=$largest;
                for ($j=$largest-1;$j>=$largest-2;$j--) {
                    if ((defined($Fixtotal{$j}))and(!(defined($Mark{$j})))) {
                        if ((($Realign[$i*3-2]-$j)<=$untrustexp)and($Fixtotal{$j}>$Fixtotal{$target})) {
                            $target=$j;
                        }
                    }
                }
                if ($Pattern[$i-1] eq '110') {
                    $Realign[$i*3-1]=$Realign[$i*3-2]-$target;
                    $Realign[$i*3-2]=$target;
                    if (0==$target) {
                        $Pattern[$i-1]='101';
                    } else {
                        $Pattern[$i-1]='111';
                    }
                } else {
                    $Realign[$i*3-3]=$Realign[$i*3-2]-$target;
                    $Realign[$i*3-2]=$target;
                    if (0==$target) {
                        $Pattern[$i-1]='101';
                    } else {
                        $Pattern[$i-1]='111';
                    }
                }
            }
            if ($judge2 eq 'f') {
                $target=$largest;
                for ($j=$largest-1;$j>=$largest-2;$j--) {
                    if ((defined($Fixtotal{$j}))and(!(defined($Mark{$j})))) {
                        if ((($Realign[$i*3+1]-$j)<=$untrustexp)and($Fixtotal{$j}>$Fixtotal{$target})) {
                            $target=$j;
                        }
                    }
                }
                if ($Pattern[$i] eq '110') {
                    $Realign[$i*3+2]=$Realign[$i*3+1]-$target;
                    $Realign[$i*3+1]=$target;
                    if (0==$target) {
                        $Pattern[$i]='101';
                    } else {
                        $Pattern[$i]='111';
                    }
                } else {
                    $Realign[$i*3]=$Realign[$i*3+1]-$target;
                    $Realign[$i*3+1]=$target;
                    if (0==$target) {
                        $Pattern[$i]='101';
                    } else {
                        $Pattern[$i]='111';
                    }
                }
            }
            if (($judge1 eq 'r')and($judge2 eq 'r')) { ### turn to repeat-only reads only when both together
                $Realign[$i*3-2]=$Realign[$i*3-2]+$Realign[$i*3-3]+$Realign[$i*3-1];
                $Realign[$i*3-3]=0;
                $Realign[$i*3-1]=0;
                $Pattern[$i-1]='010';
                $Realign[$i*3+1]=$Realign[$i*3+1]+$Realign[$i*3]+$Realign[$i*3+2];
                $Realign[$i*3]=0;
                $Realign[$i*3+2]=0;
                $Pattern[$i]='010';
            }
        }
    }
}

sub forceflank_filter { ### get rid of flanks <=$_[0] to further reduce noises
    my %Fixtotal=();
    my $temp='';
    for ($i=$#Pattern;$i>=0;$i--) { ### collect information to mark fix length reads to skip this filter
        if (($Pattern[$i] eq '111')or($Pattern[$i] eq '101')) {
            if (($Realign[$i*3]>$_[0])and($Realign[$i*3+2]>$_[0])) {
                $Fixtotal{$Realign[$i*3+1]}=1;
            }
        }
    }
    for ($i=$#Pattern;$i>=0;$i--) {
        if (!(defined($Fixtotal{$Realign[$i*3+1]}))) {
            if ($Realign[$i*3]<=$_[0]) {
                $Realign[$i*3+1]+=$Realign[$i*3];
                $Realign[$i*3]=0;
            }
            if ($Realign[$i*3+2]<=$_[0]) {
                $Realign[$i*3+1]+=$Realign[$i*3+2];
                $Realign[$i*3+2]=0;
            }
            $temp='';
            if ($Realign[$i*3]>0) {
                $temp.='1';
            } else {
                $temp.='0';
            }
            if ($Realign[$i*3+1]>0) {
                $temp.='1';
            } else {
                $temp.='0';
            }
            if ($Realign[$i*3+2]>0) {
                $temp.='1';
            } else {
                $temp.='0';
            }
            $Pattern[$i]=$temp;
        }
    }
}

sub reliab_fix { ### get %Reliab
    my $i=0;
    my %minratio=();
    my %nearmid=();
    my %maxmin=();
    my $temp=0;
    my $index=0;
    for ($i=$#Pattern;$i>=0;$i--) {
        if (($Pattern[$i] eq '111')or($Pattern[$i] eq '101')) {
            $index=$Realign[$i*3+1];
            $temp=($ReadsLength-$index)/2;
            if (((abs($Realign[$i*3]-$temp))<=$_[0])or((abs($Realign[$i*3+2]-$temp))<=$_[0])) {
                $nearmid{$index}=1;
            }
            $temp=(max($Realign[$i*3],$Realign[$i*3+2]))/(min($Realign[$i*3],$Realign[$i*3+2]));
            if (!(defined($minratio{$index}))) {
                $minratio{$index}=$temp;
                $maxmin{$index}=min($Realign[$i*3],$Realign[$i*3+2]);
            } else {
                $minratio{$index}=min($minratio{$index},$temp);
                $maxmin{$index}=max($maxmin{$index},min($Realign[$i*3],$Realign[$i*3+2]));
            }
        }
    }
    my $result=0;
    foreach (keys(%minratio)) {
        $result=3;
        if (!(defined($nearmid{$_}))) {
            $result--;
        }
        if ($minratio{$_}>$_[1]) {
            $result--;
        }
        if ($maxmin{$_}<$_[2]) {
            $result--;
        }
        if (3==$result) {
            $Reliab{$_}='High';
        } elsif ($result>0) {
            $Reliab{$_}='Medium';
        } else {
            $Reliab{$_}='Low';
        }
    }
}

sub collect { ### get $totalpair $totalrepeat @TopFix @Fixreads @Openreads @FixOpenreads $expjudge
    my $i=0;
    my $j=0;
    my $index=0;
    $totalpair=0; ### total pair number of Fix/Open reads
    $totalrepeat=0; ### reads number of repeat-only reads from not repeat-only pairs
    my @Fix=(); ### pair number contributing to each fixed STR length (=X.5 if a pair has 2 STR length, rarely happen)
    @Fixreads=(); ### reads number with fixed STR length from pairs contributing to $Fix[]
    @TopFix=(); ### all effective lengths of @Fix, sorted by small to large length
    $expjudge='Low'; ### open reads with flank length between 1 $_[1]/$_[1] $_[2] are considered proof 'High'/'Medium' of long expansion
    $expjudge2='Low'; ### reliability of super long expansion including repeat-only pairs based on $expjudge
    @Openreads=(); ### reads number with open-end from pairs contributing to $Open[]
    @FixOpenreads=(); ### reads number with open-end from pairs contributing to $Fix[]
    my $temp=0; ### for last $Openreads[]
    my $short1=0;
    my $short2=0;
    my $totalrepeattrust=0;
    my %check=();
    my $checkindex='';

    for ($i=$#Pattern;$i>=0;$i=$i-2) { ### get @Fix @Fixreads $expjudge $totalrepeat, prepare @Openreads
        if (($Pattern[$i] eq '111')or($Pattern[$i-1] eq '111')or($Pattern[$i] eq '101')or($Pattern[$i-1] eq '101')) { ### Fix pair
            $totalpair++;
            if ((($Pattern[$i] eq '111')or($Pattern[$i] eq '101'))and(($Pattern[$i-1] eq '111')or($Pattern[$i-1] eq '101'))) {
                $index=$Realign[$i*3+1];
                if (defined($Fix[$index])) {
                    $Fix[$index]+=0.5;
                    $Fixreads[$index]++;
                } else {
                    $Fix[$index]=0.5;
                    $Fixreads[$index]=1;
                }
                $index=$Realign[$i*3-2];
                if (defined($Fix[$index])) {
                    $Fix[$index]+=0.5;
                    $Fixreads[$index]++;
                } else {
                    $Fix[$index]=0.5;
                    $Fixreads[$index]=1;
                }
            } else {
                if (($Pattern[$i] eq '111')or($Pattern[$i] eq '101')) {
                    $index=$Realign[$i*3+1];
                } else {
                    $index=$Realign[$i*3-2];
                }
                if (defined($Fix[$index])) {
                    $Fix[$index]++;
                    $Fixreads[$index]++;
                } else {
                    $Fix[$index]=1;
                    $Fixreads[$index]=1;
                }
                if (($Pattern[$i] eq '110')or($Pattern[$i] eq '011')or($Pattern[$i-1] eq '110')or($Pattern[$i-1] eq '011')) {
                    if (($Realign[$i*3+1]<$index)or($Realign[$i*3-2]<$index)) { ### abnormal open reads can be due to flank mismatches
                        if (defined($FixOpenreads[$index])) {
                            $FixOpenreads[$index]++;
                        } else {
                            $FixOpenreads[$index]=1;
                        }
                    } else {
                        $Fixreads[$index]++;
                    }
                }
            }
        } elsif (((($Pattern[$i] eq '110')or($Pattern[$i] eq '011'))and($Pattern[$i-1] ne '111')and($Pattern[$i-1] ne '101'))or((($Pattern[$i-1] eq '110')or($Pattern[$i-1] eq '011'))and($Pattern[$i] ne '111')and($Pattern[$i] ne '101'))) { ### Open pair
            $totalpair++;
            if (($Pattern[$i] ne '010')and($Pattern[$i-1] ne '010')) {
                $index=max($Realign[$i*3+1],$Realign[$i*3-2]);
                if ((($Pattern[$i] eq '110')or($Pattern[$i] eq '011'))and(($Pattern[$i-1] eq '110')or($Pattern[$i-1] eq '011'))) {
                    if (defined($Openreads[$index])) {
                        $Openreads[$index]+=2;
                    } else {
                        $Openreads[$index]=2;
                    }
                } else {
                    if (defined($Openreads[$index])) {
                        $Openreads[$index]++;
                    } else {
                        $Openreads[$index]=1;
                    }
                }
            } else {
                $totalrepeat++;
                $temp++;
            }
            if (($Realign[$i*3]>$_[0])or($Realign[$i*3+2]>$_[0])or($Realign[$i*3-3]>$_[0])or($Realign[$i*3-1]>$_[0])) { ### exclude likely short flank noises from repeat-only pairs
                if ($Realign[$i*3]>0) {
                    $checkindex="5 $Realign[$i*3]";
                    if (!(defined($check{$checkindex}))) {
                        if ($Pattern[$i-1] eq '010') {
                            $totalrepeattrust++;
                        }
                        if ($Realign[$i*3]<=$_[1]) {
                            $short1++;
                        } elsif ($Realign[$i*3]<=$_[2]) {
                            $short2++;
                        }
                        $check{$checkindex}=1;
                    }
                }
                if ($Realign[$i*3+2]>0) {
                    $checkindex="3 $Realign[$i*3+2]";
                    if (!(defined($check{$checkindex}))) {
                        if ($Pattern[$i-1] eq '010') {
                            $totalrepeattrust++;
                        }
                        if ($Realign[$i*3+2]<=$_[1]) {
                            $short1++;
                        } elsif ($Realign[$i*3+2]<=$_[2]) {
                            $short2++;
                        }
                        $check{$checkindex}=1;
                    }
                }
                if ($Realign[$i*3-3]>0) {
                    $checkindex="5 $Realign[$i*3-3]";
                    if (!(defined($check{$checkindex}))) {
                        if ($Pattern[$i] eq '010') {
                            $totalrepeattrust++;
                        }
                        if ($Realign[$i*3-3]<=$_[1]) {
                            $short1++;
                        } elsif ($Realign[$i*3-3]<=$_[2]) {
                            $short2++;
                        }
                        $check{$checkindex}=1;
                    }
                }
                if ($Realign[$i*3-1]>0) {
                    $checkindex="3 $Realign[$i*3-1]";
                    if (!(defined($check{$checkindex}))) {
                        if ($Pattern[$i] eq '010') {
                            $totalrepeattrust++;
                        }
                        if ($Realign[$i*3-1]<=$_[1]) {
                            $short1++;
                        } elsif ($Realign[$i*3-1]<=$_[2]) {
                            $short2++;
                        }
                        $check{$checkindex}=1;
                    }
                }
            }
        } elsif ((($Pattern[$i] eq '010')and(($Pattern[$i-1] eq '100')or($Pattern[$i-1] eq '001')))or(($Pattern[$i-1] eq '010')and(($Pattern[$i] eq '100')or($Pattern[$i] eq '001')))) { ### repeat only + flank only pair
            $totalpair++;
            $totalrepeat++;
            $totalrepeattrust++;
        }
    }

    if ($short1>=$_[3]) {
        $expjudge='High';
    } elsif ($short2>=$_[3]) {
        $expjudge='Medium';
    } elsif (($short1>0)and($short2>0)) {
        $expjudge='Medium';
    }
    if (($expjudge eq 'Medium')and($totalrepeattrust>=$_[3])) {
        $expjudge='High';
    } elsif (($expjudge eq 'Low')and($totalrepeattrust>=$_[3])) {
        $expjudge='Medium';
    }
    $expjudge2=$expjudge;
    if (0==$totalrepeattrust) {
        $expjudge2='Low';
    } elsif ($totalrepeattrust<$_[3]) {
        if ($expjudge2 eq 'High') {
            $expjudge2='Medium';
        }
    }

    $Openreads[0]=0; ### finish @Openreads
    if ($#Openreads>=$#Fix) {
        push @Openreads, ($temp);
    } else {
        $Openreads[$#Fix+1]=$temp;
    }
    for ($i=1;$i<=$#Openreads;$i++) {
        if (defined($Openreads[$i])) {
            $Openreads[$i]+=$Openreads[$i-1];
        } else {
            $Openreads[$i]=$Openreads[$i-1];
        }
    }

    for ($i=$#Fix;$i>=0;$i--) { ### finish undefined $FixOpenreads[]=0
        if (!(defined($FixOpenreads[$i]))) {
            $FixOpenreads[$i]=0;
        }
    }

    my %Merge=(); ### merge certain fix lengths to nearby lengths (2nt range)
    my $temp1='';
    my $temp2='';
    for ($i=$#Fix;$i>=0;$i--) {
        if (defined($Fix[$i])) {
            $temp1='';
            if (defined($Fix[$i-1])) {
                if ($Fix[$i]<=($nearby*$Fix[$i-1])) {
                    $temp1='-';
                }
            }
            if (defined($Fix[$i-2])) {
                if ($Fix[$i]<=($nearby*$Fix[$i-2])) {
                    if (defined($Fix[$i-1])) {
                        if ($Fix[$i-1]<=($nearby*$Fix[$i-2])) {
                            $temp1='--';
                        }
                    } else {
                        $temp1='--';
                    }
                }
            }
            $temp2='';
            if (defined($Fix[$i+1])) {
                if ($Fix[$i]<=($nearby*$Fix[$i+1])) {
                    $temp2='+';
                }
            }
            if (defined($Fix[$i+2])) {
                if ($Fix[$i]<=($nearby*$Fix[$i+2])) {
                    if (defined($Fix[$i+1])) {
                        if ($Fix[$i+1]<=($nearby*$Fix[$i+2])) {
                            $temp2='++';
                        }
                    } else {
                        $temp2='++';
                    }
                }
            }
            $Merge{$i}=join '',($temp1,$temp2);
        }
    }
    foreach $i (keys(%Merge)) {
        if ($Merge{$i} eq '') {
            next;
        } elsif ($Merge{$i} eq '-') {
            $Fixreads[$i-1]+=$Fixreads[$i];
            $Fixreads[$i]=0;
            $FixOpenreads[$i-1]+=$FixOpenreads[$i];
            $FixOpenreads[$i]=0;
            $Openreads[$i-1]=$Openreads[$i];
        } elsif ($Merge{$i} eq '+') {
            $Fixreads[$i+1]+=$Fixreads[$i];
            $Fixreads[$i]=0;
            $FixOpenreads[$i+1]+=$FixOpenreads[$i];
            $FixOpenreads[$i]=0;
        } elsif ($Merge{$i} eq '-+') {
            $Fixreads[$i-1]+=$Fixreads[$i]*$Fix[$i-1]/($Fix[$i-1]+$Fix[$i+1]);
            $Fixreads[$i+1]+=$Fixreads[$i]*$Fix[$i+1]/($Fix[$i-1]+$Fix[$i+1]);
            $Fixreads[$i]=0;
            $FixOpenreads[$i-1]+=$FixOpenreads[$i]*$Fix[$i-1]/($Fix[$i-1]+$Fix[$i+1]);
            $FixOpenreads[$i+1]+=$FixOpenreads[$i]*$Fix[$i+1]/($Fix[$i-1]+$Fix[$i+1]);
            $FixOpenreads[$i]=0;
            $Openreads[$i-1]+=($Openreads[$i]-$Openreads[$i-1])*$Fix[$i-1]/($Fix[$i-1]+$Fix[$i+1]);
        } elsif ($Merge{$i} eq '--') {
            $Fixreads[$i-2]+=$Fixreads[$i];
            $Fixreads[$i]=0;
            $FixOpenreads[$i-2]+=$FixOpenreads[$i];
            $FixOpenreads[$i]=0;
            $Openreads[$i-2]=$Openreads[$i];
        } elsif ($Merge{$i} eq '++') {
            $Fixreads[$i+2]+=$Fixreads[$i];
            $Fixreads[$i]=0;
            $FixOpenreads[$i+2]+=$FixOpenreads[$i];
            $FixOpenreads[$i]=0;
        } elsif ($Merge{$i} eq '--++') {
            $Fixreads[$i-2]+=$Fixreads[$i]*$Fix[$i-2]/($Fix[$i-2]+$Fix[$i+2]);
            $Fixreads[$i+2]+=$Fixreads[$i]*$Fix[$i+2]/($Fix[$i-2]+$Fix[$i+2]);
            $Fixreads[$i]=0;
            $FixOpenreads[$i-2]+=$FixOpenreads[$i]*$Fix[$i-2]/($Fix[$i-2]+$Fix[$i+2]);
            $FixOpenreads[$i+2]+=$FixOpenreads[$i]*$Fix[$i+2]/($Fix[$i-2]+$Fix[$i+2]);
            $FixOpenreads[$i]=0;
            $Openreads[$i-2]+=($Openreads[$i]-$Openreads[$i-2])*$Fix[$i-2]/($Fix[$i-2]+$Fix[$i+2]);
        } elsif ($Merge{$i} eq '-++') {
            $Fixreads[$i-1]+=$Fixreads[$i]*$Fix[$i-1]/($Fix[$i-1]+$Fix[$i+2]);
            $Fixreads[$i+2]+=$Fixreads[$i]*$Fix[$i+2]/($Fix[$i-1]+$Fix[$i+2]);
            $Fixreads[$i]=0;
            $FixOpenreads[$i-1]+=$FixOpenreads[$i]*$Fix[$i-1]/($Fix[$i-1]+$Fix[$i+2]);
            $FixOpenreads[$i+2]+=$FixOpenreads[$i]*$Fix[$i+2]/($Fix[$i-1]+$Fix[$i+2]);
            $FixOpenreads[$i]=0;
            $Openreads[$i-1]+=($Openreads[$i]-$Openreads[$i-1])*$Fix[$i-1]/($Fix[$i-1]+$Fix[$i+2]);
        } elsif ($Merge{$i} eq '--+') {
            $Fixreads[$i-2]+=$Fixreads[$i]*$Fix[$i-2]/($Fix[$i-2]+$Fix[$i+1]);
            $Fixreads[$i+1]+=$Fixreads[$i]*$Fix[$i+1]/($Fix[$i-2]+$Fix[$i+1]);
            $Fixreads[$i]=0;
            $FixOpenreads[$i-2]+=$FixOpenreads[$i]*$Fix[$i-2]/($Fix[$i-2]+$Fix[$i+1]);
            $FixOpenreads[$i+1]+=$FixOpenreads[$i]*$Fix[$i+1]/($Fix[$i-2]+$Fix[$i+1]);
            $FixOpenreads[$i]=0;
            $Openreads[$i-2]+=($Openreads[$i]-$Openreads[$i-2])*$Fix[$i-2]/($Fix[$i-2]+$Fix[$i+1]);
        }
    }

    for ($i=$#Fix;$i>=0;$i--) { ### prepare @TopFix
        if (defined($Fix[$i])) {
            if ((!(defined($Mark{$i})))and($Merge{$i} eq '')) { ### count into @TopFix if length not merged & trustful without too short flanks
                push @TopFix, ($i);
            }
        }
    }
    for ($i=0;$i<$#TopFix;$i++) { ### finish @TopFix by sorting index
        for ($j=$i+1;$j<=$#TopFix;$j++) {
            if ($TopFix[$i]>$TopFix[$j]) {
                ($TopFix[$i],$TopFix[$j])=($TopFix[$j],$TopFix[$i]);
            }
        }
    }
}

sub estimate { ### $_[0] repeat unit length
    @call=();
    my $temp='';
    my $fraE=0;
    my $sumFO=0;
    my $sumY=0;
    my $sumYmin=0;
    my $sumYmax=0;
    my $i=0;
    my $j=0;
    my @frac=();
    my $len=0;

    my $overfix=0;
    for ($i=0;$i<=$#TopFix;$i++) {
        $len=$TopFix[$i];
        if ($len>=$ReadsLength) {
            $readslengthwarning='y';
            return('');
        }
        $sumYmin+=max(2*$len*($Fixreads[$len]-$buffer)/($ReadsLength-$len),0);
        $sumFO+=$FixOpenreads[$len];
        if ((int($sumYmin*(1-$bufferratio)+0.5))>($Openreads[$len]+$sumFO+$buffer)) {
            $overfix=1;
        }
        $sumYmax+=2*$len*($Fixreads[$len]+$buffer)/($ReadsLength-$len);
        $sumY+=max(2*$len*$Fixreads[$len]/($ReadsLength-$len)-$FixOpenreads[$len],0);
    }

    my $flank=$Openreads[$#Openreads]-$sumY;
    my $judgeE='n';
    my $mosaic=0;
    my $mosaic1=0;
    my $mosaic2=0;
    my $exp=0;
    my $bigchange=0;
    if (((int($sumYmax*(1+$bufferratio)+0.5))<($Openreads[$len]+$sumFO-$buffer))or($totalrepeat>0)or($Openreads[$#Openreads]>$Openreads[$len])) {
        if (($flank<=0)and(0==$totalrepeat)) {
            $judgeE='n';
        } elsif ($flank<=0) {
            if ($flank>(0-$buffer)) { ### persume one flank reads lost
                $flank=1;
                $judgeE='y';
            } else {
                $judgeE='n';
            }
        } else {
            $judgeE='y';
        }
        if (((int($sumYmax*(1+$bufferratio)+0.5))>=($Openreads[$len]+$sumFO-$buffer))and(0==$totalrepeat)) { ### too few remaining open-end reads may be noise
            if (($Openreads[$#Openreads]-$Openreads[$len])<=$buffer) {
                $judgeE='n';
            }
        }
    }
    if ($judgeE eq 'y') {
        if ($#TopFix>=0) {
            $temp=0;
            $bigchange=0;
            for ($i=0;$i<=$#TopFix;$i++) {
                $len=$TopFix[$i];
                $temp+=$Fixreads[$len]*2*$ReadsLength/$flank/($ReadsLength-$len);
                if ((($Fixreads[$len]-(max($Fixreads[$len]-$buffer,0)))/($ReadsLength-$len))>$bigchange) {
                    $j=$i;
                    $bigchange=($Fixreads[$len]-(max($Fixreads[$len]-$buffer,0)))/($ReadsLength-$len);
                }
            }
            $fraE=100/(1+$temp);
            if ($flank<=$buffer) {
                $mosaic1=0;
            } else {
                $temp=0;
                for ($i=0;$i<$#TopFix;$i++) {
                    $len=$TopFix[$i];
                    $temp+=$Fixreads[$len]*2*$ReadsLength/($flank-$buffer)/($ReadsLength-$len);
                }
                $len=$TopFix[$i];
                $temp+=($Fixreads[$len]+$buffer)*2*$ReadsLength/($flank-$buffer)/($ReadsLength-$len);
                $mosaic1=int(100/(1+$temp)+0.5);
            }
            $temp=0;
            for ($i=0;$i<=$#TopFix;$i++) {
                $len=$TopFix[$i];
                if ($i!=$j) {
                    $temp+=$Fixreads[$len]*2*$ReadsLength/($flank+$buffer)/($ReadsLength-$len);
                } else {
                    $temp+=(max($Fixreads[$len]-$buffer,0))*2*$ReadsLength/($flank+$buffer)/($ReadsLength-$len);
                }
            }
            $mosaic2=int(100/(1+$temp)+0.5);
        } else {
            $fraE=100;
            $mosaic1=100;
            $mosaic2=100;
        }
        $mosaic=int($fraE+0.5);
        push @frac, ($mosaic-0.1);
        if ($repeat_only>0) {
            push @frac, ($mosaic-0.2);
        }
        $temp=max($mosaic-$mosaic1,$mosaic2-$mosaic);
        $mosaic=join '+-',($mosaic,$temp);
        $exp=2*$ReadsLength*$totalrepeat/$flank+$ReadsLength;
        $temp=join "\t",((int(10*$exp/$_[0]+0.5))/10,'I',$expjudge,$mosaic);
        push @call, ($temp);
        if ($repeat_only>0) {
            $exp=2*$ReadsLength*($totalrepeat+2*$repeat_only)/$flank+$ReadsLength;
            $temp=join "\t",((int(10*$exp/$_[0]+0.5))/10,'I+',$expjudge2,$mosaic);
            push @call, ($temp);
        }
    }

    my $i2=0;
    my $len2=0;
    for ($i=0;$i<=$#TopFix;$i++) {
        $len=$TopFix[$i];
        if ($#TopFix>=1) {
            $temp=0;
            $bigchange=0;
            for ($i2=0;$i2<=$#TopFix;$i2++) {
                $len2=$TopFix[$i2];
                $temp+=($ReadsLength-$len)*$Fixreads[$len2]/($ReadsLength-$len2)/$Fixreads[$len];
                if (($i2!=$i)and((($Fixreads[$len2]-(max($Fixreads[$len2]-$buffer,0)))/($ReadsLength-$len2))>$bigchange)) {
                    $j=$i2;
                    $bigchange=($Fixreads[$len2]-(max($Fixreads[$len2]-$buffer,0)))/($ReadsLength-$len2);
                }
            }
            $mosaic=int((100-$fraE)/$temp+0.5);
            if ($Fixreads[$len]<=$buffer) {
                $mosaic1=0;
            } else {
                $temp=0;
                if ($i!=$#TopFix) {
                    for ($i2=0;$i2<$#TopFix;$i2++) {
                        if ($i2==$i) {
                            $temp+=1;
                        } else {
                            $len2=$TopFix[$i2];
                            $temp+=($ReadsLength-$len)*$Fixreads[$len2]/($ReadsLength-$len2)/($Fixreads[$len]-$buffer);
                        }
                    }
                    $len2=$TopFix[$i2];
                    $temp+=($ReadsLength-$len)*($Fixreads[$len2]+$buffer)/($ReadsLength-$len2)/($Fixreads[$len]-$buffer);
                } else {
                    for ($i2=0;$i2<($#TopFix-1);$i2++) {
                        $len2=$TopFix[$i2];
                        $temp+=($ReadsLength-$len)*$Fixreads[$len2]/($ReadsLength-$len2)/($Fixreads[$len]-$buffer);
                    }
                    $len2=$TopFix[$i2];
                    $temp+=($ReadsLength-$len)*($Fixreads[$len2]+$buffer)/($ReadsLength-$len2)/($Fixreads[$len]-$buffer);
                    $temp+=1;
                }
                $mosaic1=int((100-$fraE)/$temp+0.5);
            }
            $temp=0;
            for ($i2=0;$i2<=$#TopFix;$i2++) {
                $len2=$TopFix[$i2];
                if ($i2!=$j) {
                    if ($i2!=$i) {
                        $temp+=($ReadsLength-$len)*$Fixreads[$len2]/($ReadsLength-$len2)/($Fixreads[$len]+$buffer);
                    } else {
                        $temp+=1;
                    }
                } else {
                    $temp+=($ReadsLength-$len)*(max($Fixreads[$len2]-$buffer,0))/($ReadsLength-$len2)/($Fixreads[$len]+$buffer);
                }
            }
            $mosaic2=int((100-$fraE)/$temp+0.5);
        } elsif ($judgeE eq 'n') {
            $mosaic=100;
            $mosaic1=100;
            $mosaic2=100;
        } else {
            $temp=$Fixreads[$len]*2*$ReadsLength/$flank/($ReadsLength-$len);
            $mosaic=int(100*$temp/(1+$temp)+0.5);
            $temp=max(($Fixreads[$len]-$buffer)*2*$ReadsLength/($flank+$buffer)/($ReadsLength-$len),0);
            $mosaic1=int(100*$temp/(1+$temp)+0.5);
            if ($flank>$buffer) {
                $temp=($Fixreads[$len]+$buffer)*2*$ReadsLength/($flank-$buffer)/($ReadsLength-$len);
                $mosaic2=int(100*$temp/(1+$temp)+0.5);
            } else {
                $mosaic2=100;
            }
        }
        push @frac, ($mosaic);
        $temp=max($mosaic-$mosaic1,$mosaic2-$mosaic);
        $mosaic=join '+-',($mosaic,$temp);
        if ((1==$overfix)and($Reliab{$len} eq 'High')) {
            $temp=join "\t",((int(10*$len/$_[0]+0.5))/10,'D','Medium',$mosaic);
        } else {
            $temp=join "\t",((int(10*$len/$_[0]+0.5))/10,'D',$Reliab{$len},$mosaic);
        }
        push @call, ($temp);
    }

    for ($i=0;$i<$#call;$i++) { ### sort results by large to small fractions
        for ($j=$i+1;$j<=$#call;$j++) {
            if ($frac[$i]<$frac[$j]) {
                ($call[$i],$call[$j])=($call[$j],$call[$i]);
                ($frac[$i],$frac[$j])=($frac[$j],$frac[$i]);
            }
        }
    }
}

### read parameters ###
my $i=0;
@Para=();
### Para[0] result file from _realigner.pl
### Para[1] output file for all STRs no matter variant/non-variant
### Para[2] fix length supported only by reads with either flank <= this length MAY be marked/converted (suggestion = 3)
### Para[3] reads category of Para[2] with total reads number >= this will NOT be marked/converted (suggestion = 3)
### when repeat length is too long to allow long enough flanks, this limitation will be modified
### Para[4] open-end reads with flank <= this will be considered possible noise from repeat-only pairs (suggestion = 10)
### Para[5] if support pairs for expansion >= this, $Para[4] will be considered not noise (suggestion = 1, 0 will actually mask $Para[4])
### Para[6] reads number threshold to detect flank bias (suggestion = 20)
### Para[7] ratio threshold to determine flank bias (suggestion = 0.8)
### Para[8] "on"/"off" switch to apply single side analysis if flank bias detected (suggestion ="off" but for cases may worth the risk to turn on)
### Para[9] reads number shift limit for mosaic calculation & ratio check & indirect allele judge (suggestion = 3)
### Para[10] known value ratio shift limit for ratio check (suggestion = 0.3)
### Para[11] range for direct estimation reliability (suggestion = 15)
### Para[12] ratio threshold for direct estimation reliability (suggestion = 7)
### Para[13] length threshold for direct estimation reliability (suggestion = 10)
### Para[14~15] two ranges for expansion reliability (suggestion = 20 40)
### Para[16] minimum unique short flank & repeat-only reads number for expansion reliability (suggestion = 2, <=0 means no limit)
### Para[17] optional & skip if omitted, repeat number change (Xr/X) or nt change (Xn) threshold consider as variant into a separate file, suggestion = 100r
### Para[18] optional when $Para[17] is on, combination of 'Cnumber'+'Anumber'+'F/R'+'H/M/L'+'Y/N'+'E/B', suggestion = default = 'C10A5RMNB', order doesn't matter
### 'Cnumber' skip STR results if non repeat-only pairs are fewer than the number
### 'Anumber' skip STR results if allele fraction is lower than the number%
### 'F' apply $Para[17] to results by reads with flank only, ''='R' apply $Para[17] to results by reads including repeat-only pairs
### 'H/M/L' apply $Para[17] to results by reads including repeat-only pairs only with reliability no lower High/Medium/Low, ''='M'
### ''='N' skip results with alleles closed to reads length as not trustful, 'Y' apply $Para[17] to them
### ''='B' apply $Para[17] to both expansions and extractions, 'E' apply $Para[17] to only expansions
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-h')or($ARGV[$i] eq '--help')) {
        print "Instruction:\n";
        print "\tThe Caller module processes the realignment results, and genotypes all possible alleles at each STR.\n";
        print "\tSee README.txt Step (3) for input requirements, and Step (4) for more details of the output.\n";
        print "Usage Sample (without sub output of candidate STRs):\n";
        print "\tperl LUSTR_Caller.pl -i <realign.txt> -o <result.txt> --offtarget\n";
        print "Usage Sample (with sub output of candidate STRs):\n";
        print "\tperl LUSTR_Caller.pl -i <realign.txt> -o <result.txt> --offtarget -s 100n --subthres FH\n";
        print "Options:\n";
        print "\t--help/-h\t\tPrint usage instructions\n";
        print "\t--input/-i <file>\tRealignment results obtained in the previous step. See details in README.txt Step (3)\n";
        print "\t--output/-o <file>\tOutput file for STR genotyping results. See format details in README.txt Step (4)\n";
        print "\t--untrust1 <value>\tLength threshold of short flankings in direct size calling to enter further check (Default = 3)\n";
        print "\t--release1 <value>\tReads number threshold to trust short flankings in direct size calling\n";
        print "\t\t\t\tEffective for --untrust1, and will be automatically modified when repeat size is too long to apply long flankings\n";
        print "\t\t\t\t(Default = 3)\n";
        print "\t--untrust2 <value>\tLength threshold of short flankings from reads with single flanking to enter further check (Default = 10)\n";
        print "\t--release2 <value>\tPair number threshold to trust short flankings from reads with single flanking\n";
        print "\t\t\t\tIf untrusted, the short flankings are considered potential noise from repeat-only pairs\n";
        print "\t\t\t\t(Default = 1)\n";
        print "\t--num/-n <value>\tAllowance of reads number variation for fraction calculation and allele possibility check (Default = 3)\n";
        print "\t--ratio/-r <value>\tRatio variation allowance for allele possibility check (Default = 0.3)\n";
        print "\t--reliabdir1 <value>\tAllowance of position variation to determine reliability for repeats shorter than read length (Default = 15)\n";
        print "\t--reliabdir2 <value>\tSetting of flanking ratio threshold to determine reliability for repeats shorter than read length (Default = 7)\n";
        print "\t--reliabdir3 <value>\tSetting of flanking length threshold to determine reliability for repeats shorter than read length (Default = 10)\n";
        print "\t--reliabest1 <value>\tSetting of the shorter flanking length threshold to determine reliability for repeats longer than read length (Default = 20)\n";
        print "\t--reliabest2 <value>\tSetting of the longer flanking length threshold to determine reliability for repeats longer than read length (Default = 40)\n";
        print "\t--reliabest3 <value>\tSetting of the reads number threshold to determine reliability for alleles longer than read length (Default = 2)\n";
        print "\t--biasn <value>\t\tReads number threshold for flank bias detection (Default = 20)\n";
        print "\t--biasr <value>\t\tCalculation ratio threshold (0~1, no bias to complete bias) for flank bias detection (Default = 0.8)\n";
        print "\t--offtarget\t\tSwitch on single side analysis if flank bias detected (Default = off, but recommended for STRs with high risk of offtargets)\n";
        print "\t--sub/-s <string>\tSwitch on a sub output for STR candidates with variations beyond the provided threshold\n";
        print "\t\t\t\tUse a number or number+\"r\" to set threshold for repeat number variation, or number+\"n\" for nucleotide length variation\n";
        print "\t\t\t\t(Default = off, recommended setting = 100n or 100r if a sub output is wanted)\n";
        print "\t--subthres <string>\tModification of the sub output effective for --sub/-s\n";
        print "\t\t\t\tThe string is a combination by \"Cnumber\"+\"Apercent\"+\"F/R\"+\"H/M/L\"+\"Y/N\"+\"E/B\" in any order\n";
        print "\t\t\t\t\t\"Cnumber\" - skip candidate STRs not meeting this minimum coverage from non-repeat-only pairs (Default = C10)\n";
        print "\t\t\t\t\t\"Apercent\" - skip candidate STRs when the varied allele has a fraction lower than this minimum percent % (Default = A5)\n";
        print "\t\t\t\t\t\"F\" - only consider varied alleles estimated by non-repeat-only pairs, when \"R\" consider all estimations (Default = R)\n";
        print "\t\t\t\t\t\"H\"/\"M\"/\"L\" - minimum reliability requirement (\"High\"/\"Medium\"/\"Low\") for the alleles with variations (Default = M)\n";
        print "\t\t\t\t\t\"N\" - skip alleles with repeat length close to read length, when \"Y\" will keep them (Default = N)\n";
        print "\t\t\t\t\t\"E\" - only consider expanded alleles compared to reference, when \"B\" consider both expansion/contraction (Default = B)\n";
        print "\t\t\t\t(Default = C10A5RMNB, but will be disabled if --sub/-s is not present)\n";
        exit;
    }
}
$Para[0]='';
$Para[1]='';
$Para[2]=3;
$Para[3]=3;
$Para[4]=10;
$Para[5]=1;
$Para[6]=20;
$Para[7]=0.8;
$Para[8]='off';
$Para[9]=3;
$Para[10]=0.3;
$Para[11]=15;
$Para[12]=7;
$Para[13]=10;
$Para[14]=20;
$Para[15]=40;
$Para[16]=2;
for ($i=0;$i<=$#ARGV;$i++) {
    if (($ARGV[$i] eq '-i')or($ARGV[$i] eq '--input')) {
        $i++;
        $Para[0]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-o')or($ARGV[$i] eq '--output')) {
        $i++;
        $Para[1]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--untrust1') {
        $i++;
        $Para[2]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--release1') {
        $i++;
        $Para[3]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--untrust2') {
        $i++;
        $Para[4]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--release2') {
        $i++;
        $Para[5]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--biasn') {
        $i++;
        $Para[6]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--biasr') {
        $i++;
        $Para[7]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--offtarget') {
        $Para[8]='on';
    } elsif (($ARGV[$i] eq '-n')or($ARGV[$i] eq '--num')) {
        $i++;
        $Para[9]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-r')or($ARGV[$i] eq '--ratio')) {
        $i++;
        $Para[10]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--reliabdir1') {
        $i++;
        $Para[11]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--reliabdir2') {
        $i++;
        $Para[12]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--reliabdir3') {
        $i++;
        $Para[13]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--reliabest1') {
        $i++;
        $Para[14]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--reliabest2') {
        $i++;
        $Para[15]=$ARGV[$i];
    } elsif ($ARGV[$i] eq '--reliabest3') {
        $i++;
        $Para[16]=$ARGV[$i];
    } elsif (($ARGV[$i] eq '-s')or($ARGV[$i] eq '--sub')) {
        $i++;
        $Para[17]=$ARGV[$i];
    }
}
if (($Para[0] eq '')or($Para[1] eq '')) {
    print ("Input/Output not appointed. Stop running.\n");
    print "Please use --help or -h to see instructions.\n";
    exit;
}
if (defined($Para[17])) {
    for ($i=0;$i<=$#ARGV;$i++) {
        if ($ARGV[$i] eq '--subthres') {
            $i++;
            $Para[18]=$ARGV[$i];
        }
    }
}

### main ###
use List::Util qw/max min/;
open (INPUT_DATA, "<$Para[0]") or die "Couldn't open: $!";
$ReadsLength=&getReadsLength(100);
$fixflank_replace=max(int($Para[2]),0); ### used in &shortflank_filter
$fixflank_number=max(int($Para[3]),0); ### used in &shortflank_filter
$untrustexp=max(10,$fixflank_replace); ### used in &shortflank_filter
$expflank_replace=max(int($Para[4]),0); ### used in &shortflank_filter, also in &collect for $expjudge & $expjudge2
$expflank_number=max(int($Para[5]),0); ### used in &shortflank_filter
$biasmin=max(int($Para[6]),0);
$biaslimit=max($Para[7],0);
$SSswitch=$Para[8];
$buffer=max($Para[9],0); ### used for mosaic calculation & ratio check in &estimate
$bufferratio=max($Para[10],0); ### used for ratio check in &estimate
@fixreli=(max($Para[11],0),max($Para[12],0),max($Para[13],0)); ### used in &reliab_fix for %Reliab
@expreli=($expflank_replace,max($Para[14],0),max($Para[15],0),max($Para[16],0)); ### used in &collect for $expjudge & $expjudge2
@expreliSS=($expflank_replace,max($Para[14],0),max($Para[15],0),max($Para[16]*1.5,0)); ### used when offtarget modification applies
$nearby=0.4; ### combine reads within 2nt length range if number smaller than this ratio to the main one
my $warnfile='';
if ($Para[1]=~/\..*$/) {
    $warnfile=join '',($`,"_warning",$&);
} else {
    $warnfile=join '',($Para[1],"_warning");
}
open (OUTPUT_WARN, ">$warnfile") or die "Couldn't open: $!";
print OUTPUT_WARN ("H5 = strong evidence for offtarget upstream flank sequences, or too short downstream flanks\n");
print OUTPUT_WARN ("H3 = strong evidence for offtarget downstream flank sequences, or too short upstream flanks\n");
print OUTPUT_WARN ("V5 = strong evidence for complex variations at upstream flank near STR\n");
print OUTPUT_WARN ("V3 = strong evidence for complex variations at downstream flank near STR\n");
print OUTPUT_WARN ("NA = no/weak evidence for offtarget/variations\n");
print OUTPUT_WARN ("Following: distant flanking hits upstream/downstream, nearby flanking hits upstream/downstream\n");
print OUTPUT_WARN ("See last parts for overall duplicate filter warning and running log\n\n");
open (OUTPUT_RESULT, ">$Para[1]") or die "Couldn't open: $!";
my $line='';
my $id='';
my $j=0;
my $totalcount=0;
my $zerocount=0;
my $warncount=0;
my $callcount=0;
@call=(); ### each element = 'repeat    D/I/I+    High/Medium/Low/NA    fraction+-uncertainty' separated by '\t'
my $order=0;
my $pairusage='';
my @tempdata=();
my @tempresult=();
my @tempXs=();
my $chi=0;
my $warn='';
my $origin=0;
my $temp='';
my $ul=0;
my $homo=0;
my $flank_up=0; ### reads number, for homologous warning
my $flank_down=0; ### reads number, for homologous warning
my $flankSTR_up=0; ### reads number, for mutation warning
my $flankSTR_down=0; ### reads number, for mutation warning
$repeat_only=0; ### pair number of repeat only ones, used for duplicate filter warning
my $repeatonly_equal=0; ### number of STRs with repeat_only pairs equal to n*n (n=STR unit length)
my $repeatonly_large=0; ### number of STRs with repeat_only pairs larger than n*n (n=STR unit length)
@Realign=(); ### groups of 6 elements for X1 X2 X3 X1 X2 X3 (reads in pair)
@Pattern=(); ### groups of 2 elements for pattern of X1 X2 X3 X1 X2 X3 (reads in pair)
%Mark=(); ### length from too short flanks may be skipped for @TopFix/@Fix in &collect
%Reliab=(); ### reliability of fix lengths
$readslengthwarning='n'; ### gives warning if unexpected reads length detected in &estimate_
my $reffile='Unknown';
if (-s INPUT_DATA) {
    chomp ($line=<INPUT_DATA>);
    @tempdata=split /\s+/, $line;
    if (($tempdata[$#tempdata-2] eq 'built')and($tempdata[$#tempdata-1] eq 'by')) {
        $reffile=$tempdata[$#tempdata];
    }
    print OUTPUT_RESULT ("\# Header Line for Each STR:\n");
    print OUTPUT_RESULT ("\# Column 1 = STR ID\n");
    print OUTPUT_RESULT ("\# Column 2 = Repeat Unit\n");
    print OUTPUT_RESULT ("\# Column 3 = Reference Position (REF = $reffile)\n");
    print OUTPUT_RESULT ("\# Column 4 = Reference Repeat Number\n");
    print OUTPUT_RESULT ("\# Column 5 = Flanking+Repeat-only Realigned Pairs\n");
    print OUTPUT_RESULT ("\# Column 6 = Counts of ALL 3 Types of Estimated Alleles, 0 = Unable to Estimate\n");
    print OUTPUT_RESULT ("\# Allele Lines Sorted by Fractions for Each STR:\n");
    print OUTPUT_RESULT ("\# Column 1 = Allele Repeat Number\n");
    print OUTPUT_RESULT ("\# Column 2 = Allele Type (D/I/I+ = Direct/Indirect/Indirect+)\n");
    print OUTPUT_RESULT ("\# Column 3 = Reliability\n");
    print OUTPUT_RESULT ("\# Column 4 = Allele Fraction +- Uncertain Range, Rounded to Nearest 1%\n");
    print OUTPUT_RESULT ("\# Note: I+ Type is Alternative Estimation of I Type Including Repeat-only Pairs\n\n");
    chomp ($line=<INPUT_DATA>);
    chomp ($line=<INPUT_DATA>);
    while (1) {
        chomp ($id=<INPUT_DATA>);
        @tempdata=split /\s+/, $id;
        $temp=join "\t",($tempdata[0],$tempdata[1],$tempdata[2]);
        $tempdata[2]=~/:/;
        $tempdata[2]=$';
        $tempdata[2]=~/-/;
        $ul=length($tempdata[1]);
        $origin=int(($'-$`+1)/$ul*10+0.5)/10;
        $id=join "\t",($temp,$origin);
        chomp ($line=<INPUT_DATA>);
        if (!($line=~/^0/)) {
            $line=~/\s+/;
            $temp=$`;
            @Realign=(); ### get @Realign
            for ($i=$temp;$i>=1;$i--) {
                chomp ($line=<INPUT_DATA>);
                @tempresult=split /\s+/, $line;
                for ($j=$#tempresult;$j>=1;$j--) {
                    $tempresult[$j]=~s/}\{/,/;
                    $tempresult[$j]=~s/^{//;
                    $tempresult[$j]=~s/}$//;
                    @tempXs=split /,/, $tempresult[$j];
                    push @Realign, (@tempXs);
                }
            }
            @Pattern=(); ### get @Pattern
            for ($i=0;$i<=$#Realign;$i=$i+3) {
                $temp='';
                if ($Realign[$i]>0) {
                    $temp.='1';
                } else {
                    $temp.='0';
                }
                if ($Realign[$i+1]>0) {
                    $temp.='1';
                } else {
                    $temp.='0';
                }
                if ($Realign[$i+2]>0) {
                    $temp.='1';
                } else {
                    $temp.='0';
                }
                push @Pattern, ($temp);
            }
            %Mark=();
            &shortflank_filter; ### get %Mark, modify @Realign & @Pattern
            &forceflank_filter(1); ### get rid of super short flanks, modify @Realign & @Pattern
            &reliab_fix(@fixreli); ### get %Reliab
            $flank_up=0; ### get homologous/mutation/duplicate filter warning
            $flank_down=0;
            $flankSTR_up=0;
            $flankSTR_down=0;
            $repeat_only=0;
            for ($i=$#Pattern;$i>=0;$i=$i-2) {
                if (($Pattern[$i] eq '100')and($Realign[3*$i-1]<=10)) {
                    $flank_up++;
                } elsif (($Pattern[$i] eq '110')and($Realign[3*$i-1]<=10)) {
                    $flankSTR_up++;
                } elsif (($Pattern[$i] eq '001')and($Realign[3*$i-3]<=10)) {
                    $flank_down++;
                } elsif (($Pattern[$i] eq '011')and($Realign[3*$i-3]<=10)) {
                    $flankSTR_down++;
                } elsif (($Pattern[$i-1] eq '100')and($Realign[3*$i+2]<=10)) {
                    $flank_up++;
                } elsif (($Pattern[$i-1] eq '110')and($Realign[3*$i+2]<=10)) {
                    $flankSTR_up++;
                } elsif (($Pattern[$i-1] eq '001')and($Realign[3*$i]<=10)) {
                    $flank_down++;
                } elsif (($Pattern[$i-1] eq '011')and($Realign[3*$i]<=10)) {
                    $flankSTR_down++;
                }
                if (($Pattern[$i] eq '010')and($Pattern[$i-1] eq '010')) {
                    $repeat_only++;
                }
            }
            $warn='';
            if (($flank_up<$biasmin)and($flank_down<$biasmin)) { ### test of 2 values for equality
                $chi=0;
            } else {
                $chi=(abs($flank_up-$flank_down))/($flank_up+$flank_down);
            }
            if ($chi>=$biaslimit) { ### evidence for homologous flank
                if ($flank_up>$flank_down) {
                    $warn='H5';
                    $homo=5;
                } else {
                    $warn='H3';
                    $homo=3;
                }
            } else {
                $homo=0;
            }
            if (($flankSTR_up<$biasmin)and($flankSTR_down<$biasmin)) { ### test of 2 values for equality
                $chi=0;
            } else {
                $chi=(abs($flankSTR_up-$flankSTR_down))/($flankSTR_up+$flankSTR_down);
            }
            if ($chi>=$biaslimit) { ### evidence for mutation near flank-STR
                if ($flankSTR_up>$flankSTR_down) {
                    if ($warn eq '') {
                        $warn='V3';
                        $homo=3;
                    } elsif ($warn ne 'H5') {
                        $warn.=' V3';
                        $homo=3;
                    }
                } else {
                    if ($warn eq '') {
                        $warn='V5';
                        $homo=5;
                    } elsif ($warn ne 'H3') {
                        $warn.=' V5';
                        $homo=5;
                    }
                }
            }
            if ($warn ne '') {
                $warncount++;
            } else {
                $warn='NA';
            }
            print OUTPUT_WARN ("$id\t$warn\t$flank_up/$flank_down\t$flankSTR_up/$flankSTR_down\n");
            if ($repeat_only>(($ul+1)**2)) {
                $repeatonly_large++;
            } elsif ($repeat_only>=($ul**2-$ul+1)) {
                $repeatonly_equal++;
            }
            if ($SSswitch eq 'on') { ### modify @Realign and @Pattern to eliminate influence from offtarget homologous flanks
                if (5==$homo) {
                    for ($i=$#Pattern;$i>=0;$i=$i-2) {
                        if (($Realign[3*$i-1]<=10)and($Realign[3*$i+2]<=10)and(($Realign[3*$i-3]>10)or($Realign[3*$i]>10))) { ### delete pairs with upstream flanks but no meaningful downstream flanks
                            $Pattern[$i]='000';
                            $Pattern[$i-1]='000';
                        } elsif (($Realign[3*$i-3]<=10)and($Realign[3*$i]<=10)and(($Realign[3*$i-1]>10)or($Realign[3*$i+2]>10))) { ### duplicate pairs to mirror those misfired
                            push @Realign, ($Realign[3*$i-3],$Realign[3*$i-2],$Realign[3*$i-1],$Realign[3*$i],$Realign[3*$i+1],$Realign[3*$i+2]);
                            push @Pattern, ($Pattern[$i-1],$Pattern[$i]);
                        }
                    }
                } elsif (3==$homo) {
                    for ($i=$#Pattern;$i>=0;$i=$i-2) {
                        if (($Realign[3*$i-3]<=10)and($Realign[3*$i]<=10)and(($Realign[3*$i-1]>10)or($Realign[3*$i+2]>10))) { ### delete pairs with downstream flanks but no meaningful upstream flanks
                            $Pattern[$i]='000';
                            $Pattern[$i-1]='000';
                        } elsif (($Realign[3*$i-1]<=10)and($Realign[3*$i+2]<=10)and(($Realign[3*$i-3]>10)or($Realign[3*$i]>10))) { ### duplicate pairs to mirror those misfired
                            push @Realign, ($Realign[3*$i-3],$Realign[3*$i-2],$Realign[3*$i-1],$Realign[3*$i],$Realign[3*$i+1],$Realign[3*$i+2]);
                            push @Pattern, ($Pattern[$i-1],$Pattern[$i]);
                        }
                    }
                }
            }
            if (($SSswitch eq 'on')and($homo>0)) {
                &collect(@expreliSS);
            } else {
                &collect(@expreli); ### collect data from @Realign & @Pattern, get $totalpair $totalrepeat @TopFix @Fixreads @Openreads @FixOpenreads $expjudge
            }
            $pairusage=join '+',($totalpair,$repeat_only);
            &estimate($ul); ### estimation of all possible alleles
            if (-1==$#call) {
                print OUTPUT_RESULT ("$id\t$pairusage\t0\n\n");
            } else {
                $callcount++;
                $temp=$#call+1;
                print OUTPUT_RESULT ("$id\t$pairusage\t$temp\n");
                for ($i=0;$i<=$#call;$i++) {
                    print OUTPUT_RESULT ("$call[$i]\n");
                }
                print OUTPUT_RESULT ("\n");
            }
        } else {
            $zerocount++;
            print OUTPUT_RESULT ("$id\t0+0\t0\n\n");
        }
        chomp ($line=<INPUT_DATA>);
        $totalcount++;
        if (0==($totalcount%10000)) {
            print ("$totalcount STR processed...\n");
        }
        if (eof INPUT_DATA) {
            last;
        }
    }
}
if (0!=($totalcount%10000)) {
    print ("$totalcount STR processed...\n");
}
print ("Done\n");
print OUTPUT_WARN ("\n");
if ($readslengthwarning eq 'y') {
    print ("WARNING: Unable to determine reads length using header of $Para[0]\n");
}
if (($repeatonly_equal>=10)and($repeatonly_large<=10)and($repeatonly_large<(0.1*$repeatonly_equal))) { ### overall duplication filter warning
    print ("WARNING: positive evidence for duplication filter detected (Evidence unlike/likely: $repeatonly_large/$repeatonly_equal)\n");
    print OUTPUT_WARN ("Positive evidence for duplication filter (Evidence unlike/likely: $repeatonly_large/$repeatonly_equal)\n");
} else {
    print OUTPUT_WARN ("Negative evidence for duplication filter (Evidence unlike/likely: $repeatonly_large/$repeatonly_equal)\n");
}
close (INPUT_DATA);
close (OUTPUT_RESULT);

my $varcount=0; ### make sub output if required
my $varfile='';
my $threshold=0;
my $thresholdtype='';
my $checkcolumn='';
my $minpair=10;
my $minfrac=5;
my $usage='R';
my $reliab='M';
my $direction='B';
my $RLallele='N';
my $pass='';
my $discard1=0;
my $discard2=0;
my @ALine=();
if ($#Para>=17) {
    print ("Generating sub output file...\n");
    $thresholdtype=$Para[17];
    if ($thresholdtype=~/[Rr]$/) {
        $threshold=$`;
        $thresholdtype='r';
    } elsif ($thresholdtype=~/[Nn]$/) {
        $threshold=$`;
        $thresholdtype='n';
    } else {
        $threshold=$thresholdtype;
        $thresholdtype='r';
    }
    if ($#Para>=18) {
        $checkcolumn=$Para[18];
    }
    if (($checkcolumn=~/C\d+/)or($checkcolumn=~/c\d+/)) {
        $temp=$&;
        $temp=~/[Cc]/;
        $minpair=$';
    }
    if (($checkcolumn=~/A\d+/)or($checkcolumn=~/a\d+/)) {
        $temp=$&;
        $temp=~/[Aa]/;
        $minfrac=$';
    }
    if (($checkcolumn=~/F/)or($checkcolumn=~/f/)) {
        $usage='F';
    }
    if (($checkcolumn=~/H/)or($checkcolumn=~/h/)) {
        $reliab='H';
    } elsif (($checkcolumn=~/L/)or($checkcolumn=~/l/)) {
        $reliab='L';
    }
    if (($checkcolumn=~/E/)or($checkcolumn=~/e/)) {
        $direction='E';
    }
    if (($checkcolumn=~/Y/)or($checkcolumn=~/y/)) {
        $RLallele='Y';
    }
    if ($Para[1]=~/\..*$/) {
        $varfile=join '',($`,'_change',$threshold,$thresholdtype,'_filterC',$minpair,'A',$minfrac,$usage,$direction,$reliab,$RLallele,$&);
    } else {
        $varfile=join '',($Para[1],'_change',$threshold,$thresholdtype,'_filterC',$minpair,'A',$minfrac,$usage,$direction,$reliab,$RLallele);
    }
    open (INPUT_DATA, "<$Para[1]") or die "Couldn't open: $!";
    open (OUTPUT_RESULT, ">$varfile") or die "Couldn't open: $!";
    for ($i=1;$i<=14;$i++) { ### copy header
        chomp ($line=<INPUT_DATA>);
        print OUTPUT_RESULT ("$line\n");
    }
    if (!(eof INPUT_DATA)) {
        while (1) {
            @ALine=();
            $pass='n';
            chomp ($line=<INPUT_DATA>);
            push @ALine, ($line);
            @tempdata=split /\t/, $line;
            for ($i=$tempdata[$#tempdata];$i>=1;$i--) {
                chomp ($line=<INPUT_DATA>);
                push @ALine, ($line);
            }
            chomp ($line=<INPUT_DATA>);
            $tempdata[4]=~/\+/;
            $origin=$tempdata[3];
            $ul=length($tempdata[1]);
            if ($`>=$minpair) { ### check pair usage filter
                if ($RLallele eq 'N') { ### set filter for alleles closed to reads length
                    $discard1=($ReadsLength-1)/(length($tempdata[1]));
                    $discard2=($ReadsLength+1)/(length($tempdata[1]));
                } else {
                    $discard1=-1;
                    $discard2=-1;
                }
                for ($i=$#ALine;$i>=1;$i--) {
                    @tempdata=split /\t/, $ALine[$i];
                    if (($usage eq 'F')and($tempdata[1] eq 'I+')) { ### skip if results including repeat-only pairs not accepted
                        next;
                    }
                    if ((($tempdata[2] eq 'Medium')and($reliab eq 'H'))or(($tempdata[2] eq 'Low')and($reliab ne 'L'))) { ### check reliability filter
                        next;
                    }
                    if (($tempdata[0]>=$discard1)and($tempdata[0]<=$discard2)) { ### skip if alleles closed to reads length not accepted
                        next;
                    }
                    $tempdata[3]=~/\+-/;
                    if ($`<$minfrac) { ### check allele fraction filter
                        next;
                    }
                    if ($thresholdtype eq 'r') {
                        if ((abs($tempdata[0]-$origin))>=$threshold) { ### pass change filter
                            if (($direction eq 'B')or($tempdata[0]>=$origin)) {
                                $pass='y';
                                last;
                            }
                        }
                    } else {
                        if (((abs($tempdata[0]-$origin))*$ul)>=$threshold) { ### pass change filter
                            if (($direction eq 'B')or($tempdata[0]>=$origin)) {
                                $pass='y';
                                last;
                            }
                        }
                    }
                }
            }
            if ($pass eq 'y') {
                $varcount++;
                for ($i=0;$i<=$#ALine;$i++) {
                    print OUTPUT_RESULT ("$ALine[$i]\n");
                }
                print OUTPUT_RESULT ("\n");
            }
            if (eof INPUT_DATA) {
                last;
            }
        }
    }
    close (INPUT_DATA);
    close (OUTPUT_RESULT);
    print ("Done\n");
}

my $uncallcount=$totalcount-$zerocount-$callcount;
print OUTPUT_WARN ("\n");
print ("Total $totalcount STRs processed\n");
print ("$zerocount STRs have no realignments, skipped\n");
print ("$uncallcount STRs fail to estimate genotypes\n");
print ("$callcount STRs succeed to estimate genotypes\n");
print ("--- Check results in $Para[1]\n");
print OUTPUT_WARN ("Total $totalcount STRs processed\n");
print OUTPUT_WARN ("$zerocount STRs have no realignments, skipped\n");
print OUTPUT_WARN ("$uncallcount STRs fail to estimate genotypes\n");
print OUTPUT_WARN ("$callcount STRs succeed to estimate genotypes\n");
if ($varfile ne '') {
    print ("$varcount STRs exhibited potential expansions/truncations\n");
    print OUTPUT_WARN ("$varcount STRs exhibited potential expansions/truncations\n");
    print ("--- Check results in $varfile\n");
}
print ("$warncount STRs show evidence for offtarget flanks or complex variations\n");
print OUTPUT_WARN ("$warncount STRs show evidence for offtarget flanks or complex variations\n");
print ("--- Check results in $warnfile\n");
close (OUTPUT_WARN);

exit;








