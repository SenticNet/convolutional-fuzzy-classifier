warning off all
tic

filename =sprintf('target1.para');

if exist(filename) == 2
           
data = importdata(filename);

nonodes = size(data,2)-2;
notime = size(data,1);

data1=sortrows(data,-1*size(data,2));
data = data1;

num = int32(size(data,1))*0.2;

filename2 =sprintf('target_highml_es.para');

dlmwrite(filename2,data(1:num,1:size(data,2)));

end
toc

exit;
