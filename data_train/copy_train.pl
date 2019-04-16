$dom1="books";
$dom2="dvd";
$dom3="electronics";
$dom4="kitchen";
$src=$dom1;$trg=$dom3;$srcpath="/home/senticteam/online/multidomain/exp1/data/$src";
$trgpath="/home/senticteam/online/multidomain/exp1/data/$trg";

system("cp $srcpath/positive2.review source/");
system("cp $srcpath/positive5.review source/");

system("cp $srcpath/negative2.review source/");
system("cp $srcpath/negative5.review source/");

system("cp $trgpath/positive2.review target/");
system("cp $trgpath/positive5.review target/");

system("cp $trgpath/negative2.review target/");
system("cp $trgpath/negative5.review target/");
