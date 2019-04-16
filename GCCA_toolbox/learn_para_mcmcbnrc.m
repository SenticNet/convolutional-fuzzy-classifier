function learn_para_mcmcgcca(data, filename, ord)

[n ncases] = size(data);

 %check corr
xcoor = data;
th = 0.5;   
    
%%% for all combinations
arr2 = -100;
for j=1:n
 for loopi = 1:n
   for loopj = loopi+1:n  
    % for loopk = loopj+1:n 
    %   for loopm = loopk+1:n
    
 %   if abs(xcoor(j,loopi))>th && abs(xcoor(j,loopj))>th
    
   arr = zeros(1,n+3);
   new_dag = zeros(n,n);
   new_dag(j,loopi)=1;
   if loopj>0 new_dag(j,loopj)=1; end
  % if loopk>0 new_dag(j,loopk)=1; end
  % if loopm>0 new_dag(j,loopm)=1; end
   
   % check if already printed
    dag = new_dag;
    likA = 0 ; A = 0; 
    arr(1)=j;
     for g=1:size(new_dag,1)
        arr(g+1)=dag(j,g);
     end
     
    arr(size(new_dag,2)+2)=A;
    arr(size(new_dag,2)+3)=likA; 
    if arr2 == -100
        arr2 = arr;
    else
        arr2=[arr2;arr];
    end

 %   end % end of corr if
    
   %   end % end of loopm
   % end % end of for loopk
  end % end of for loopj
  %dlmwrite(filename,arr2(loopi:n,:),'-append');
end % end of for loopi
end % end of for j

% dlmwrite(filename,arr2);
dagb = arr2(:,2:end-2);
arra = zeros(1,size(arr2,1));
arrb = zeros(1,size(arr2,1));
curv=-1; cntt = 0;
for l=1:size(arr2,1)
    
    dag = zeros(n,n);
    node = arr2(l,1);
    dag(node,:)=dagb(l,:); 
    
    tran = 0;
    % check if transfer is possible to previous index
 
%%{
    if curv~=-1 
  
    n1 = size(dag,1);
    i1 = node;
    yd=data(i1,ord+1:end);
    c = 0;
    xd = 0;
    for j1=1:n1
        if(dag(i1,j1)>0)
            c = c+1;
            if c>1
            xd=[xd;data(j1,1:end-ord)];
            else
            xd = data(j1,1:end-ord);
            end
        end
    end        
    data2 = [yd;xd];
    
   % if size(p1,1)==size(data2,1)
    
    ind = randi(size(data2,2)-5);
    q1 = data2(:,ind);
%    ind = randi(size(data2,2));
    r1 = data2(:,ind+2);
    
    dist1 = (p1-q1).^2;
    sum1 = sqrt(sum(dist1));

    dist2 = (r1-q1).^2;
    sum2 = sqrt(sum(dist2));

    dist3 = (p1-q1).^2;
    sum3 = sqrt(sum(dist3));
    sumall = sum1 + sum2 + sum3;
     
    if( sumall < 2*3.14/sqrt(curv) && sum2 < 3.14/sqrt(curv) )
         tran = 1;  
    end
    
    tran
    sumall
    curv
    
   % end
      
    end
%%} 

    cntt = cntt+1;
    if cntt > 10
        cntt = 0;
        tran = 0;
    end

    trafm = 1; 
    if tran==0 || curv==-1 
    [A, likA, curv, p1, data1] =  bnrc(dag,data,node,ord); 
    else
        % transform
        A1=cov(data1');
        B=cov(data2');
        one = eig(A1); two = eig(B);
        trafm = sqrt(diag(one))/(sqrt(diag(two))); 
    end

     arra(l) = A*trafm(1,1);
     arrb(l) = likA*trafm(1,1);
    
end

for l=1:size(arr2,1)
   arr2(l,end)=arrb(l);
   arr2(l,end-1)=arra(l);
end    


dlmwrite(filename,arr2);

% for j=1:n  
%  for loopi = 1:n   
%    arr2=zeros(n,n+3);
%    for loopj = loopi:n   
%    parfor loopk = loopj:n  
%    arr = zeros(1,n+3);
%    new_dag = zeros(n,n);
%    new_dag(j,loopi)=1;
%    new_dag(j,loopj)=1;
%    new_dag(j,loopk)=1; 
% 
%    % check if already printed
%     dag = new_dag;
%     
%     [A, likA] =  bnrc(dag,data,j,ord);
%     arr(1)=j;
%      for g=1:size(new_dag,1)
%         arr(g+1)=dag(j,g);
%      end
%     arr(size(new_dag,2)+2)=A;
%     arr(size(new_dag,2)+3)=likA; 
%     
%     arr2(loopk,:)=arr;
%    end % end of for loopk
%    dlmwrite(filename,arr2(loopj:n,:),'-append');
%    %arr2(loopj,:)=arr;
%   end % end of for loopj
%   %dlmwrite(filename,arr2(loopi:n,:),'-append');
% end % end of for loopi
% end % end of for j

end

function [A, likA, curv, p1, data1] = bnrc(dag,datab,nodei,ord) 
likA = -10000;
 A = log(0.01);
 
%filename = sprintf('E:/FullBNT-1.0.4/historybnrc');
%dlmwrite('/home/iti/FullBNT-1.0.7_bnrc/FullBNT-1.0.4/dag_bnrc.txt',dag);

e = zeros(1,size(dag,1));
n = size(dag,1);
p = size(datab,2);

i = nodei;
    yd=datab(i,ord+1:end);
    c = 0;
    xd = 0;

    for j=1:n
        if(dag(i,j)>0)
            c = c+1;
            if c>1
            xd=[xd;datab(j,1:end-ord)];
            else
            xd = datab(j,1:end-ord);
            end
        end
    end

       % find new curv
    data1 = [yd;xd];  
    data1 = data1/norm(data1);
    if size(data1,1)>2
        
    
    if size(data1,1)==3
    x = data1(1:2,:);
    y = data1(2:3,:);
    z = [data1(1,:);data1(3,:)];    
    else    
    x = data1(1:end-1,:);
    y = data1(2:end,:);
    z = [data1(1,:);data1(3:end,:)];
    end
    
    f = gcurvature(x,y,z);
    
    curv = -1; 
    for f1=1:size(f,2)  
      for d1=1:size(f,1) 
         if f(d1,f1)>0 && curv == -1
             if min(f(:,f1))>=0
             curv = max(f(:,f1));
         %    get smallest > 0
             for d2=1:size(f,1)
                if curv>f(d2,f1) && f(d2,f1)>0 
                     curv = f(d2,f1);
                end
             end           
             curvi = f1;  
             p1 = data1(:,curvi);
             end
         end
      end
    end
    
    if curv==-1
        p1=0;
    end
    
    else
        curv = -1;
        p1 = 0;
    end
    
  %  datamvu = [ yd ; xd ];

  %  cd /home/iti/prib/mcmc4/sep10/nips/GCCA_toolbox/drtoolbox
   % addpath(genpathKPM(pwd));

   % [mappedX, mapping] = compute_mapping(datamvu', 'MVU', 2);
  %  datamvu = mapping.vec';

   %  yd2 = datamvu(1,:);
   %  xd2 = datamvu(2:size(datamvu,1),:);
   
   % cd /home/iti/prib/mcmc4/sep10/nips/GCCA_toolbox

[test_set_predictions, chain_stats, final_pred] = bayes_mls_gaussgc([yd' xd'],[yd' xd']);
        likA = chain_stats.LL_store(5);
   
ccaStartup;
DEMO = 1;
X = [yd; xd];

sfile = ['ccademo.',num2str(DEMO),'.net'];
PVAL = 0.01;
NLAGS  = 1;

if ord > 1 
    NLAGS = -1;
end

% find best model order
if NLAGS == -1,
%    disp('finding best model order ...');
    [bic,aic] = cca_find_model_order(X,2,5);
%    disp(['best model order by Bayesian Information Criterion = ',num2str(bic)]);
%    disp(['best model order by Aikaike Information Criterion = ',num2str(aic)]);
    NLAGS = aic;
end

[ret] = cca_granger_regressbnrc(X,NLAGS,1);

%[PR,q] = cca_findsignificance(ret,PVAL,1);
% GC = ret.gc
% GC2 = GC.*PR;
% cca_pajek(PR,GC,sfile);
   
     A = 0 ;
    for c1=2:c+1
        A=A+log(1-ret.prb(1,c1)+0.01);
    end
    

end
