
open(FILE,$ARGV[0]);
@lines = <FILE>;
close(FILE);

open(FILE,"../count");
$train_cnt=<FILE>;
chomp($train_cnt);
close(FILE);
$test_cnt = 4000;


$k = scalar(@lines)/5;

for($i=0;$i<1;$i++)
{
  open(TRAIN, ">fold/train$i.txt");
  open(TEST, ">fold/test$i.txt");

  $starti = $k*$i;
  $endi = $k*($i+1);

  $endi = $train_cnt;

  for($j=0;$j<scalar(@lines);$j++){

  if($j>=$starti && $j<$endi){print TRAIN $lines[$j];}
  else{print TEST $lines[$j];}

  }

  close(TRAIN);
  close(TEST);

}
