for k=0:0
    filename = sprintf('fold/test%d_y',k);
    testy = importdata(filename);
    
    num = floor(size(testy,1)/50);
    testy2 = testy(1:num*50,:);
    filename = sprintf('cnn/run2/test%d_y',k);
    dlmwrite(filename,testy2);
    
    filename = sprintf('fold/val%d_y',k);
    testy = importdata(filename);
    
    num = floor(size(testy,1)/50);
    testy2 = testy(1:num*50,:);
    filename = sprintf('cnn/run2/val%d_y',k);
    dlmwrite(filename,testy2);
    
    
    filename = sprintf('fold/train%d_y',k);
    testy = importdata(filename);
    
    num = floor(size(testy,1)/50);
    testy2 = testy(1:num*50,:);
    filename = sprintf('cnn/run2/train%d_y',k);
    dlmwrite(filename,testy2);
end