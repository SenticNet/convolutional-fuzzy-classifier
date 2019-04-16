% ANFIS classification on raw cancer dataset
classfmea = zeros(10,4);
allfmea = zeros(10,1);

mfarr = ['gaussmf ';'gbellmf ';'gauss2mf';'trimf   ';'trapmf  '];
mfarr = cellstr(mfarr);

for mf = 1:1

for k=0:0

filename =  sprintf('test_vectors%d',k);    
CancerTrain  = importdata(filename);
CancerTrain = CancerTrain;

filename = sprintf('../cnn/set_test_rnn%d',k);
test = importdata(filename);
n3 = 2*size(test,2);

% n3=floor(0.7*n1);
% CancerTrain = [CancerTrainx(1:n3,:),CancerTrainy(1:n3)];
% CancerTest = [CancerTrainx(n3+1:end,:),CancerTrainy(n3+1:end)];

% get class label column from training dataset
CancerTest = CancerTrain(end-n3+1:end,:);
CancerTrain = CancerTrain(1:end-n3,:);

y = CancerTrain(:, end);

% prepare training data
data_train = CancerTrain(:,1:4);   % dataset with 5 attributes
data_train = CancerTrain;        % dataset combined with class label

% prepare testing data
data_test = CancerTest(:,1:4);     % dataset with 5 attributes
test_y = CancerTest(:, end);        % prepare testing dataset's class label

% ANFIS training
numMFs = 2;                     % number of membership functions
inmftype = mfarr{mf};% 'gaussmf';           % membership function type for input
outmftype = 'linear';           % membership function type for output
n_epochs = 10;                % number of epochs

% generate Sugeno-type FIS structure from data using grid partition
a = genfis1(data_train, numMFs, inmftype, outmftype)

% start training
[fis, error] = anfis(data_train, a, n_epochs);

% to check with training dataset
data_train_chk = data_train(:,1:end-1);


% ANFIS classifying training
output_train = evalfis(data_train_chk, fis);
result_train = output_train;

for i = 1:size(output_train,1)
     if(output_train(i) >= 4.5)
         result_train(i) = 3;
     elseif output_train(i) >= 1.5
         result_train(i) = 2;
     elseif output_train(i) >= 0.5
         result_train(i) = 1;
     else
         result_train(i) = 0;
     end;
end;

% ANFIS evaluation training
 correct_train = 0;
 
 for i = 1:size(result_train,1)
     if(result_train(i) == y(i))
         correct_train = correct_train + 1;
     end;
 end;

correct_train

% confusion matrix
[trainC,trainOrder] = confusionmat(y,result_train);

trainC

% ANFIS classifying testing
[output, IRR, ORR, ARR] = evalfis(data_test, fis);
result = output;


for i = 1:size(output,1)
     if(output(i) >= 2.5)
         result(i) = 3;
     elseif output(i) >= 1.5
         result(i) = 2;
     elseif output(i) >= 0.5
         result(i) = 1;
     else
         result(i) = 0;
     end;
end;


% ANFIS evaluation
 correct_test = 0;
 
 for i = 1:size(result,1)
     if(result(i) == test_y(i))
         correct_test = correct_test + 1;
     end;
 end;
 
 correct_test

% confusion matrix
[testC,testOrder] = confusionmat(test_y,result);

testC

cm = testC;
nclass = 4;

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
fmea(x) = 2*pre(x)*rec(x)/(pre(x)+rec(x)+0.01);
%fmea(x) = (tp+tn)/(tp+fp+tn+fn);

end

classfmea(k+1,:)=fmea(:);
fmea = sum(fmea)/nclass;

fmea

allfmea(k+1)=fmea;

allfmea

filename = sprintf('allfmea%d_2.txt',mf);
dlmwrite(filename,allfmea);
filename = sprintf('classfmea%d_2.txt',mf);
dlmwrite(filename,classfmea);
%dlmwrite('output.txt', result);
%exit
end

end

exit
