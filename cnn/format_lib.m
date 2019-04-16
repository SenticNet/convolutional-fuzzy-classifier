data = importdata('data_svm2_0.txt');

load ../count;
train_cnt = 12000-count;
train_cnt = 12000;

data = data(:,1:end-train_cnt);
labely = [];

for i=1:size(data,2)

if data(end,i)==1
labely(i)=-1;
end

if data(end,i)~=1
labely(i)=1;
end

end

data2 = data(1:end-1,:);
data3 = labely;

data4 = [data3' data2'];

data4 = [data4(1,:);data4];

dlmwrite('train',data4);

data = importdata('data_svm2_0.txt');

data = data(:,end-train_cnt+1:end);

labely = [];
for i=1:size(data,2)

if data(end,i)==1
labely(i)=-1;
end

if data(end,i)~=1
labely(i)=1;
end

end

data2 = data(1:end-1,:);
data3 = labely;

data4 = [data3' data2'];

data4 = [data4(1,:);data4];

dlmwrite('test',data4);

data = importdata('data_svm2_0.txt');

load ../count;
train_cnt = 12000-count;
train_cnt = 12000;

data = data(:,1:end-train_cnt);
labely = [];
for i=1:size(data,2)

if data(end,i)==2
labely(i)=-1;
end

if data(end,i)~=2
labely(i)=1;
end

end

data2 = data(1:end-1,:);
data3 = labely;

data4 = [data3' data2'];

data4 = [data4(1,:);data4];

dlmwrite('train2',data4);

data = importdata('data_svm2_0.txt');

data = data(:,end-train_cnt+1:end);
labely = [];
for i=1:size(data,2)

if data(end,i)==2
labely(i)=-1;
end

if data(end,i)~=2
labely(i)=1;
end

end

data2 = data(1:end-1,:);
data3 = labely;

data4 = [data3' data2'];

data4 = [data4(1,:);data4];

dlmwrite('test2',data4);

data = importdata('data_svm2_0.txt');

load ../count;
train_cnt = 12000-count;
train_cnt = 12000;

data = data(:,1:end-train_cnt);
labely = [];
for i=1:size(data,2)

if data(end,i)==3
labely(i)=-1;
end

if data(end,i)~=3
labely(i)=1;
end

end

data2 = data(1:end-1,:);
data3 = labely;

data4 = [data3' data2'];

data4 = [data4(1,:);data4];

dlmwrite('train3',data4);

data = importdata('data_svm2_0.txt');

data = data(:,end-train_cnt+1:end);
labely = [];
for i=1:size(data,2)

if data(end,i)==3
labely(i)=-1;
end

if data(end,i)~=3
labely(i)=1;
end

end

data2 = data(1:end-1,:);
data3 = labely;

data4 = [data3' data2'];

data4 = [data4(1,:);data4];

dlmwrite('test3',data4);

data = importdata('data_svm2_0.txt');

load ../count;
train_cnt = 12000-count;
train_cnt = 12000;

data = data(:,1:end-train_cnt);
labely = [];
for i=1:size(data,2)

if data(end,i)==4
labely(i)=-1;
end

if data(end,i)~=4
labely(i)=1;
end

end

data2 = data(1:end-1,:);
data3 = labely;

data4 = [data3' data2'];

data4 = [data4(1,:);data4];

dlmwrite('train4',data4);

data = importdata('data_svm2_0.txt');

data = data(:,end-train_cnt+1:end);
labely = [];
for i=1:size(data,2)

if data(end,i)==4
labely(i)=-1;
end

if data(end,i)~=4
labely(i)=1;
end

end

data2 = data(1:end-1,:);
data3 = labely;

data4 = [data3' data2'];

data4 = [data4(1,:);data4];

dlmwrite('test4',data4);

exit
