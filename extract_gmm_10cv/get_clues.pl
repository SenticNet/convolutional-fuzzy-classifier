

open(FILE, $ARGV[0]);
while($line = <FILE>)
{
 @list = ();
 @list = split " ",$line;
 ($word,$pos)=$line=~/word1=(.*?) pos1=(.*?) /;
 print "$word\n"; 
}
close(FILE);