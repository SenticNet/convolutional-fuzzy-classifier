system("perl reduce_len.pl");

for($k=0;$k<1;$k++){

system("cp fold1/trainc$k.txt train_data1");
system("cp fold1/trainc".$k."_y label_y");
system("cp fold2/testb".$k."_y tests_y");

#system("cnn_run.sh");

system("perl cnn_all.pl");
system("matlab -r change_y");

sleep(20);

$chk=0;
while($chk==0)
{
 if(-e "test1_y"){$chk=1;}
 sleep(20);
}

system("cp test1_y fold/test".$k."_y");
system("cp test1_x fold/test".$k."_x");
system("cp val1_y fold/val".$k."_y");
system("cp val1_x fold/val".$k."_x");
system("cp train1_y fold/train".$k."_y");
system("cp train1_x fold/train".$k."_x");

sleep(10);

}

system("python2 pack_Data_fold.py");
