warning off all
tic
ccaStartup;
ord = 1;
for i=1:1
filename =sprintf('../extract_gmm_10cv/targetes.data');
if exist(filename) == 2
       
 data = importdata(filename);    
 data = data(1:end-1,:);
 
 %data = data1(:,1:1000);
 
 filename2 =sprintf('../extract_gmm_10cv/target%d.para',i);

 if exist(filename2)==2
    delete(filename2);
 end
 
 learn_para_mcmcbnrc_sp1(data, filename2,ord);
    
%  filename2 =sprintf('data/target%d.parab',i);
% 
%  if exist(filename2)==2
%     delete(filename2);
%  end
%  
%  learn_para_mcmcbnrc_sp1b(data, filename2,1);
    
%  filename2 =sprintf('data/en_fr_dvd/target%d.parac',i);
% 
%  if exist(filename2)==2
%     delete(filename2);
%  end
%  
%  learn_para_mcmcbnrc_sp1c(data, filename2,ord);
    
 end
  
end     % end of for

exit;
toc
