# clues stops input

$cnt = 1;
open(FILE, $ARGV[0]);
while($line = <FILE>)
{
  ($word)=$line=~/word1\=(.*?) /;
  $clues{$word}=$cnt;
  $cnt++;
  # print $word."\n";
}
close(FILE);

$cnt2 = 1;
open(FILE, $ARGV[1]);
while($line = <FILE>)
{
  chomp($line);
  if($clues{$line}>0){}
  else{
   $stop{$line}=$cnt2;
   $cnt2++;
#   print $line."\n";
  }
  # print $word."\n";
}
close(FILE);


open(FILE, $ARGV[2]);
$cnt = 0;

while($line = <FILE>)
{
 chomp($line);
 $check=0;
 @list2 = ();
 @list2 = split "\t",$line;
 

 @list = ();
 @list =  split " ",$line;
 $line2 = "";
 for($i=1;$i<scalar(@list);$i++)
 {
  if($stop{lc($list[$i])}>0){}
  else{
  $line2.=" ".$list[$i];
  $len+=length($list[$i]);
  }

  if($len>15000 || $i==scalar(@list)-1){
  # $line2 =~ s/[?;:!,.'"&#]//g;
  $line2 =~ s/[^a-zA-Z#0-9 _-]//g;
  print $list2[0]."\t".$line2."\n"; 
  $line2 = ""; $len = 0; $check=1;}
 }
 if($check==0){print "1\tempty\n";}
 $cnt++;
}
close(FILE);
