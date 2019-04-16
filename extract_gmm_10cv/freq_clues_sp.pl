# top words , obj sent, 


$cnt = 1;
open(FILE, $ARGV[0]);
while($line = <FILE>)
{
 chomp($line);
 @list = ();
 @list = split " ",$line;
 $pos{lc($list[1])}=$cnt;
 $cnt++;
 if($cnt>20){last;}
}
close(FILE);

open(FILE, $ARGV[1]);
while($line = <FILE>)
{
 chomp($line);
 @list = ();
 @list = split " ",$line;

 for($i=0;$i<scalar(@list);$i++){
   if($pos{lc($list[$i])}>0)
   {
    $cntp[$pos{lc($list[$i])}][$d]=1;      
   }
 }

# if($list[0] eq "O"){$cntp[$cnt][$d]=1;}
# else{$cntp[$cnt][$d]=2;}

 $d++;

}
close(FILE);

for($i=1;$i<$cnt;$i++)
{
 for($j=0;$j<$d;$j++){
  if($cntp[$i][$j]==1){}else{$cntp[$i][$j]=0;}
  print $cntp[$i][$j]." ";
 }
 print "\n";
}

#for($j=0;$j<$d;$j++){print $cntp[$cnt][$j]." ";}
#print "\n";
