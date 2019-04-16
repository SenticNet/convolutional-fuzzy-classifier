system("python2 process_data.py GoogleNews-vectors-negative300.bin"); 
#note maximum length 80
system("python2 conv_net_sentence.py -static -word2vec");

system("rm test1_x");
system("rm train1_x");
system("rm val1_x");
system("rm test1_y");
system("rm train1_y");
system("rm val1_y");


