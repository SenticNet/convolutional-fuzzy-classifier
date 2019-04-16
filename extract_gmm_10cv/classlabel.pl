
open(FILE, $ARGV[0]);
while($line = <FILE>)
{
  @list = ();
  @list = split "\t",$line;
  #print $list[0]."iti\n";
  if($list[0] > 0){$class = 1;}
  if($list[0] < 0){$class = 2;}
  if($list[0] == 0){$class = 3;}  
  
  $line2 = substr($line,index($line," ")+3);
  print "$class\t$line2";

}
close(FILE);