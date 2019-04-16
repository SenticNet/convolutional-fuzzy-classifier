# readme clues

$cnt = 1;
open(FILE, $ARGV[0]);
while($line = <FILE>)
{
 if(index($line,"positive")>0 || index($line,"negative")>0){
 @list = ();
 @list = split " ",$line;
 ($word,$pos)=$line=~/word1=(.*?) pos1=(.*?) /;
# print "$word $pos \n"; 
 if($pos eq "adj"){$pos = "jj";}
 if($pos eq "noun"){$pos = "nn";}
 if($pos eq "verb"){$pos = "vb";}
 if($pos eq "adverb"){$pos = "rb";}
 if($pos eq "anypos"){$pos = "xx";}

 # add to hash table
 $name = lc($word)." ".lc($pos);
 $clues{$name}=$cnt;
 $cnts3[$cnt]=$name;
 $cnt++;
 }
}
close(FILE);


open(FILE, $ARGV[1]);
while($line = <FILE>)
{
 chomp($line);
 @list = ();
 @list = split " ",$line;
 $name = lc($list[0])." ".lc($list[1]);
 $cluind = $clues{$name};

 if($cluind > 1){}
 else{ $name = lc($list[0])." xx";
 $cluind = $clues{$name};
 }
 if($cluind > 1){}else{$cluid = 0;}

# print $name." ".$cluind." ".$clues{$name}."\n";

# if($cluind > 1 && $list[0] eq "O"){ $cnts[$cluind]++; }
# if($cluind > 1 && $list[0] eq "S"){ $cnts2[$cluind]++; }
  if($cluind > 1){ $cnts[$cluind]++;} 
}
close(FILE);


for($j=0; $j<$cnt; $j++){$ans{$cnts3[$j]}=$cnts[$j];}

 foreach $key (sort { $ans{$b} <=> $ans{$a} } keys %ans) {
      printf "%4d\t%s\n", $ans{$key}, $key;
    }

