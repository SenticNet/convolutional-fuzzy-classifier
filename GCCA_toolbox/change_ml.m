warning off all
tic
for i=1:96

filename =sprintf('data_jul25/target%d.para',i);

if exist(filename) == 2
           
data = importdata(filename);

data2 = [data(:,1:end-2) data(:,end)];

dlmwrite(filename,data2);

end

end
toc