
for($k=0;$k<1;$k++){
$name = "fold2/trainc".$k.".txt";
$name2 = "fold1/trainc".$k.".txt";
open(FILE,$name);
open(OUT,">$name2");
while($line=<FILE>)
{
  if(length($line)>300){
    $line2=substr($line,0,299);
	$line=$line2."\n";
   }
  print OUT $line;
}
close(FILE);
close(OUT);
system("cp fold2/trainc".$k."_y fold1/");
}
