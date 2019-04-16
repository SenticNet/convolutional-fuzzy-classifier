

# topwords highml

$cnt = 1;
open (FILE, $ARGV[0]);
while($line = <FILE>)
{
 chomp($line);
 @list = ();
 @list = split " ",$line;
 $topword{$cnt}=$list[1];
 $cnt++;
}
close(FILE);


open(FILE, $ARGV[1]);
while($line = <FILE>)
{
 $check = 0;
 @list = ();
 @list = split ",",$line;
 $node = $list[0];
 $nodename = $topword{$node};
# print "$node $nodename"."\n";

 for($i=1;$i<21;$i++){
  if($list[$i]>0 && $node!=$i){
   $node2 = $i;
   $nodename.=" ".$topword{$node2};
   if($node2 < 50){$check=1;} 
   #print $nodename." ".$nodename2."\n";
  }
 }

 if($check==1){print $nodename."\n";}

}
close(FILE);
