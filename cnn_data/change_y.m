load wordvectors;

load label_y;
load test0;

load ../count;

count = 1000;

tests = importdata('tests_y');

test0(:,end)=label_y;
test1 = test0;

test_y = test1(:,end);
test_x = test1(:,1:end-1);

[n1 n2] = size(test_x);
%n3 = floor(0.6*(n1-100));

n4 = n1-count;
n3 = floor(0.6*n4);

k = 1;

for i=1:n1

img=wordvectors(test_x(i,:)'+1,:);
data = img(:)';

  if i <= n3
     filename = sprintf('train%d_x',k);
     filename2 = sprintf('train%d_y',k);
  end

  if i>n3 && i<=n4
     filename = sprintf('val%d_x',k);
     filename2 = sprintf('val%d_y',k);
  end

   if i > n4 
      filename = sprintf('test%d_x',k);
      filename2 = sprintf('test%d_y',k);
  end
  
  dlmwrite(filename,data,'-append');
  dlmwrite(filename2,test_y(i),'-append');
  
end

%filename = sprintf('test%d_x',k);
%filename2 = sprintf('test%d_y',k);
%filenameb = sprintf('val%d_x',k);
%filename2b = sprintf('val%d_y',k);

%copyfile(filename,filenameb)
%copyfile(filename2,filename2b)

exit
