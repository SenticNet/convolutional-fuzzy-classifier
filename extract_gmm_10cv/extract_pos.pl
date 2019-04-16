

#words_pos

open(FILE, $ARGV[0]);
while($line = <FILE>)
{
 chomp($line);
 @list = split " ",$line;
# $type = $list[0];
 for($i=0;$i<scalar(@list);$i++){
 if(index($type,"<")<0){

  if(index($list[$i],"_NOUN")>0 || index($list[$i],"_ADJ")>0 || index($list[$i],"_VERB")>0 || index($list[$i],"_ADV")>0) 
  {
#    if(index($type,"O")==0){$type="O";}
#    else{$type="S";}
    
    ($word)=$list[$i]=~/(.*?)\_/;
    if(index($list[$i],"_NOUN")>0){$type2 = "NN";}
    if(index($list[$i],"_ADJ")>0){$type2 = "JJ";}
    if(index($list[$i],"_VERB")>0){$type2 = "VB";}
    if(index($list[$i],"_ADV")>0){$type2 = "RB";}
 
   print lc($word)." $type2\n";}
 }

 }
}
close(FILE);
