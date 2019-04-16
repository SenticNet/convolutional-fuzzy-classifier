

$delay = $ARGV[2];

open(FILE, $ARGV[0]);
while($line = <FILE>)
{
  if($cncpt < 6000){

  chomp($line);
  @list = ();
  @list = split " ",$line;
  
  # sentences with these words
  open(FILE2, $ARGV[1]);
  $line2o = <FILE2>; $cntk=0;
  while($line2 = <FILE2>)
  {
   $cnt = 0; $cntlist = ();
   for($i=0; $i<scalar(@list); $i++){
     if(index(lc($line2),lc($list[$i]))>0){
        $cntlist.=" ".$list[$i];
        $cnt++;
     }
   }

 if($cnt > 1){print $line2;}

 # time delayed features 
 if($delay == 1){   
      $cnt = 0; $cntlist = ();
      for($i=0; $i<scalar(@list); $i++){
      if(index(lc($line2),lc($list[$i]))>0){
        $cntlist.=" ".$list[$i];
        $cnt++;
      } 
      }
      #if($cnt > 1){print "1 ".$line2o; print "2 ".$line2;}
      if($cnt > 1){print $line2o; print $line2;}
  }

  $line2o = $line2;
  $cntk++;
  
  }
  close(FILE2);

   $cncpt++;

  }
}
close(FILE);
