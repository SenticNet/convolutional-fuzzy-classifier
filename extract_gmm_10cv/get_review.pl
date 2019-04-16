use utf8;

use Unicode::Normalize;

$cnt=0;
open(FILE, $ARGV[0]);
while($line = <FILE>)
{
 chomp($line);
 if(index($line,"\<rating\>")>=0){
  $cnt++;
  ($rating)=$line=~/rating>(.*?)</;
  $review="";
 }

 if(index($line,"\<\/text\>")>=0 && index($line,"\<text")>=0){
   ($review)=$line=~/text>(.*?)</;
   if(length($review)>2){
   $review = NFD($review);
   $review=~ s/\pM//g;
   print $rating."\t".$review."\n"; $in=0;}
 }

 if($in==1 && index($line,"\<\/text\>")<0)
 {
  $review.=$line;
 }

 if(index($line,"\<text")>=0 && index($line,"\<\/text\>")<0){
    $review=substr($line,index($line,"\<text")+6);
    $in = 1;
 }

 if($in==1 && index($line,"\<\/text\>")>=0)
 {
  $review.=substr($line,0,index($line,"\<\/text"));
  if(length($review)>2){
  $review = NFD($review);
  $review=~ s/\pM//g;
  print $rating."\t".$review."\n"; $in=0;} 
 }
 
}
close(FILE);
