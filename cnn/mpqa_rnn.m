tic
clear

%rev_ind = importdata('trainb_ind');

nclass = 4;
allfmea=zeros(10,1);
classfmea=zeros(10,nclass);
neu = 10;
loopmax = 1;

for k=0:0

for loopk=1:loopmax

if loopk==1
 
filename = sprintf('set_train_rnn%d',k);    
train1 = importdata(filename);
filename = sprintf('set_val_rnn%d',k);    
val = importdata(filename);
filename = sprintf('set_test_rnn%d',k);    
test = importdata(filename);

n10 = floor(0.2*size(test,2));

%val = test(:,1:n10);;
% val = [test test];
% 
 [val1 val2]=size(val);
 [tes1 tes2]=size(test);
% 
% val2 = floor(val2/5);
% tes2 = floor(4*(tes2/5));
% 
% val = [val(:,1:val2) test(:,1:tes2)];
% 

val = test;

[n1 n2]=size(train1);
[n3 n4a]=size(val);
[n5 n6]=size(test);

train4 = [train1 val test];

else
    train4 = newdata;
end

n4 = n2+n4a;

[n7 n8] = size(train4);
X = train4(1:end-1,:);
T = train4(end,:);
T2=transformtarget(T,nclass);
T2 = T2';
T = T2;
net = layrecnet(1:2,neu);
net.trainFcn = 'trainbr';
%net.trainFcn = 'trainbfg';
net.trainParam.epochs = 50;
%net.trainParam.lr = 0.0001;
net.trainParam.showWindow = false;
net.trainParam.showWindow=0;
net.performParam.regularization = 0.5;
net.divideFcn = 'divideind'; 
net.performFcn = 'mse'; 
net.divideParam.trainInd = 1:n2;
net.divideParam.valInd   = n2+1:n4;
net.divideParam.testInd  = n4+1:size(train4,2);

[net,tr] = train(net,X,T);
W1 = net.LW{1,1}(:,1:neu);
W2 = net.LW{1,1}(:,neu+1:2*neu);
b = net.b{1,1};
C = net.LW{2,1};

ncase=size(X,2);
Y = net(X);

%without sdp
%%{
newy = train4(end,:);
newx = [T(:,1:n4)'*net.LW{2,1};Y(:,n4+1:end)'*net.LW{2,1}];
%newx = Y'*net.LW{2,1};
newx = newx';
newdata = [newx; newy];
%%}

if loopk==loopmax
   filename = sprintf('data_svm%d_%d.txt',loopk+1,k);
   dlmwrite(filename,newdata);
   filename = sprintf('run2/train_rnn_%d',k);
   save (filename, 'net') ;
   filename = sprintf('result%d_%d.txt',loopk+1,k);
   dlmwrite(filename,Y');
end

end

[num idx] = max(T(:,n4+1:end));
[num2 idx2] = max(Y(:,n4+1:end));

dlmwrite('pred_label.txt',idx2);

cm = confusionmat(idx,idx2);

for x=1:nclass

tp = cm(x,x);
tn = cm(1,1);
for y=2:nclass
tn = tn+cm(y,y);
end
tn = tn-cm(x,x);

fp = sum(cm(:, x))-cm(x, x);
fn = sum(cm(x, :), 2)-cm(x, x);
pre(x)=tp/(tp+fp+0.01);
rec(x)=tp/(tp+fn+0.01);
%fmea(x) = 2*pre(x)*rec(x)/(pre(x)+rec(x)+0.01);
fmea(x) = (tp+tn)/(tp+fp+tn+fn);

end


classfmea(k+1,:)=fmea;
fmea = sum(fmea)/nclass;

fmea

allfmea(k+1)=fmea;

end

allfmea

dlmwrite('allfmea.txt',allfmea);
dlmwrite('classfmea.txt',classfmea);

toc
exit;
