

open(FILE, $ARGV[0]);
while($line = <FILE>)
{
 @list = ();
 @list = split "\t",$line;
 if($list[0] eq "1"){$class = 0;}
 if($list[0] eq "2"){$class = 1;}
 if($list[0] eq "3"){$class = 2;}
 if($list[0] eq "4"){$class = 3;} 
 print "$class\n";
}
close(FILE);


