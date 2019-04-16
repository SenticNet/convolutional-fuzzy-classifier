function learn_para_mcmcgcca(data, filename, ord)

[n ncases] = size(data);

%%% for all combinations
for j=1:n
 for loopi = j+1:n
   arr = zeros(1,n+2);
   for loopj = loopi+1:n
%   for loopk = 1:n  

   new_dag = zeros(n,n);
   new_dag(j,loopi)=1;
   new_dag(j,loopj)=1;
 %  new_dag(j,loopk)=1; 

   % check if already printed
    dag = new_dag;
    
    [likA] =  bnrc(dag,data,j,ord);
    arr(1)=j;
     for g=1:size(new_dag,1)
        arr(g+1)=dag(j,g);
     end
   
    arr(size(new_dag,2)+2)=likA;
    
      dlmwrite(filename,arr,'-append');
    
%  end % end of for loopk
  end % end of for loopj
end % end of for loopi
end % end of for j

end

function [likA] = bnrc(dag,datab,nodei,ord) 
likA = -10000;

e = zeros(1,size(dag,1));
n = size(dag,1);
p = size(datab,2);

i = nodei;
    yd=datab(i,1:end-ord);
    c = 0;
    xd = 0;

    for j=1:n
        if(dag(i,j)>0)
            c = c+1;
            if c>1
            xd=[xd;datab(j,1+ord:end)];
            else
            xd = datab(j,1+ord:end);
            end
        end
    end

[test_set_predictions, chain_stats, final_pred] = bayes_mls_gaussgc([yd' xd'],[yd' xd']);
        likA = chain_stats.LL_store(50);



end
