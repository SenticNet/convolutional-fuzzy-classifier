system("matlab -r format_lib");
system("pytho csv2lib.py train train.lib 0 False");
system("pytho csv2lib.py test test.lib 0 False");
system("pytho csv2lib.py train2 train2.lib 0 False");
system("pytho csv2lib.py test2 test2.lib 0 False");
system("pytho csv2lib.py train3 train3.lib 0 False");
system("pytho csv2lib.py test3 test3.lib 0 False");
system("pytho csv2lib.py train4 train4.lib 0 False");
system("pytho csv2lib.py test4 test4.lib 0 False");

