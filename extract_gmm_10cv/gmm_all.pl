

system("rm spanish_pos");
system("rm spanish_pos_cnt");
system("rm targetes.data");

system("perl split_length_word.pl subj_clue_es.txt stop_word_es.txt train2.review > train2.txt");

system("perl extract_pos.pl train2_tagged.review > spanish_pos");

system("perl cnt_clues.pl subj_clue_es.txt spanish_pos > spanish_pos_cnt");

system("perl freq_clues_sp.pl spanish_pos_cnt train2.txt >targetes.data");

