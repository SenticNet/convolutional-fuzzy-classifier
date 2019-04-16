for k=0:0

filename =  sprintf('../cnn/data_svm2_%d.txt',k);
data = importdata(filename);
load senticnet4d.txt
sent = senticnet4d(:,2:5);
filename = sprintf('../cnn/set_test_rnn%d',k);
test = importdata(filename);

datax  = data(1:end-1,:);
label = data(end,:)-1;
[n1 n2]=size(datax);
[n3 n4]=size(test);

ind = n2-n4;
datax1 = datax(:,1:ind);
new_train = sent(1:n1,:)'*datax1;
datax2 = datax(:,ind+1:end);
new_test = sent(1:n1,:)'*datax2;
new =[new_train new_test];

new = [new;label];

filename =  sprintf('test_vectors%d',k);
dlmwrite(filename,new');

end
exit

