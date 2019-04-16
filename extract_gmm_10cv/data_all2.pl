

system("perl split_train.pl train2.review");

system("perl get_gmm_sent1.pl spanish_pos_cnt target_highml_es.para > gmm_concepts_es");

# divide train2.txt into 10 fold

for($k=0;$k<1;$k++){
 system("perl get_gmm_sent2.pl gmm_concepts_es fold/train$k\.txt > orig_es_gmm_1_id");
 system("cat orig_es_gmm_1_id fold/train$k\.txt  > fold/spanish_gmm$k");
 system("perl get_y.pl fold/spanish_gmm$k > fold/trainb".$k."_y");
 system("perl split_length_word2.pl subj_clue_es\.txt stop_word_es\.txt fold/spanish_gmm$k > fold/trainb$k\.txt");
 system("perl get_y.pl fold/test$k\.txt > fold/testb".$k."_y");
 system("perl split_length_word2.pl subj_clue_es\.txt stop_word_es\.txt fold/test$k\.txt > fold/testb$k\.txt");
 system("cat fold/trainb$k\.txt fold/testb$k\.txt > fold/trainc$k\.txt");
 system("cat fold/trainb".$k."_y fold/testb".$k."_y > fold/trainc".$k."_y");
 system("cp fold/* ../cnn_data/fold2/"); 
}

