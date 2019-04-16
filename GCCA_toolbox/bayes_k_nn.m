function [pred, mc, beta_store, k_store, lik_store] = ...
   bayes_k_nn(data,test,q,costs,its,burn_in,weighted);
%
% Performs Weighted/Non-weighted Bayesian K - Nearest Neighbour model
% using q state random model
% P(Y) \propto exp( \sum_{i=1}^n \sum_{j=i+1}^n beta dirac(y_i, y_j) ... 
%							I(x_j a k nearest neighbour of x_i) )
%
% where dirac(a,b) = 1 if a==b, 0 otherwise, and I(.) is identify function
%
% based on paper by Holmes and Adams (2002) JRSSB
%
% thanks to Fangxin Hong and Jonathan Fieldsend for pointing out bugs (now corrected) in the code.
%
% Input variables:
%
% data = data set with first column integer reseponse variable Y and other 
%			columns are covariate: Y=data(:,1), X=data(:,2:size(data,2))....
%			That is, data = [Y X] where Y is a column vector with integer elements 
%			Y(i) ~ {1,2,...,q} for q class classification 
%			and X is the set of features/covariates/predictors/inputs/explanatory variables .....
% test = test set, same config as data: test = [Yt Xt] where Yt is labels of test set...
%			if you don't have these then just fill them up with dummy values, Yt = ones(nt,1) 
%			where nt = number of test points. Xt is the set of covariates for the test 
%			points. Note that the number of columns in test must equal the number of 
%			columns in data
% q = number of classes
% costs = cost vector of misclassification, e.g. costs = (1,1,2) equates....
% 			to a three class problem where missclassifying the last class is 
% 			twice as costly as misclassifying the first two. 
% its = number of samples required from MCMC
% burn_in = number of discarded burn in samples
%
%				The total run length of the algorithm is burn_in + (its * sample_rate) 
%				where sample_rate is defined below.
% weighted = {0,1} indicator. If weighted = 1 then use tricube kernel else uniform.
%
%
% Outputs:
%
% pred : probability predictions on test set. pred is an nt by q matrix where nt is the 
%			number of test points and q is the number of classes where pred(i,j) gives 
%			the probability that the ith test point is of class j.
% mc : misclassification error rate on test set
% beta_store : posterior samples of beta parameter
% k_store : posterior samples of k
% lik_store : posterior samples of likelihood values
%
% THIS PROGRAM IS FOR ACADEMIC PURPOSES ONLY. PLEASE CONTACT CHRIS HOLMES 
% IF YOU WISH TO USE IT FOR COMMERCIAL WORK.
% 

% The next lines store default parameters that can be changed by the 
% user depending on performance
sample_rate = 7; % collect every 7th sample post burn in
prior_prec = 0; % we recommend an improper prior on p(beta) ~ N(0,1/prior_prec);
k_max = min(size(data,1)-1,250); % maximum k allowed 
k=3; % starting value for k
beta = 1;  % starting value for beta
STANDARDISE = 0; % if you want to standardise X to mean 0 equal variance then set 
						% STANDARDISE = 1.
sig_beta = 0.5;  % variance of proposal distribution for new values of ....
% beta in the mcmc, beta^(t+1) = beta^(t) + sig_beta * N(0,1); 
% NOTE THAT THE ACCEPTANCE RATES ON THE MCMC NEEDS TO BE MONITORED SO YOU SHOULD 
% FINE TUNE sig_beta. A good rule of thumb due to Gareth Roberts and 
% others is to set sig_beta so acceptance rate is around 20%

% first check for sensible mcmc parameters
if its < 1000 
   warning('The number of samples specified is small');
   fprintf('\nI normaly recommend at least 1000 as a crude rule of thumb');
   dum = input('\nType any key to proceed');
end
if burn_in < 1000 
   warning('The number of burn_in samples specified is small');
   fprintf('\nI normaly recommend at least 1000 as a crude rule of thumb');
   fprintf('\nthough you need to check for convergence in likelihood of samples and k etc');
   dum = input('\nType any key to proceed');
end


% now here is the main program.....

% response should be stored in the first column of the data and test
Y = data(:,1);
Yt = test(:,1);

% convert classes to 1,...,q  if necessay
if min(Y)==0
   Y=Y+1;
   Yt=Yt+1;
end

Ytrans=Y'; % we will use Y' so sensible to store it

% see if there are unequal misclassification costs
if all(costs==costs(1)) % the equal costs
   COSTS=0;
else
   COSTS=1;
   % set up weighting matrices for predictions
   costy = costs(Yt);
   cost_mat = diag(costs);
end

if q~=max(Y)
   error('Q does not match max of Y!');
end

% store covariate/input/feature/predictor/explanatory variables
X = data(:,2:size(data,2));
Xt = test(:,2:size(test,2));

% get dimensions of data 
n = size(data,1);
nt = size(test,1);
p = size(X,2);

clear data;
clear test;

if p~=size(Xt,2)
   error('Different number of predictors in training and test sets');
end

if STANDARDISE
   % then standardise data
   for i=1:p
      mx = mean(X(:,i));
      sx = std(X(:,i));
      X(:,i) = (X(:,i)-mx) / sx;
      Xt(:,i) = (Xt(:,i)-mx) / sx;
   end
end

% set pointers for nearest neighbours 
Pointer = zeros(k_max,n); 
Test_pointer = zeros(k_max,nt); 

% if weighted you'll also need to store distances
if weighted
   Dist = zeros(k_max,n);
   Test_dist = zeros(k_max,nt);
end


% and fill pointers 
fprintf('Calculating distances on training set \n');
for i=1:n
   dis = (X - (ones(n,1)*X(i,:))).^2;
   dis =  sqrt(sum(dis'));
   [so indx]=sort(dis);
   % don't include first point (which is to itself), so go from 2:k_max+1
   Pointer(:,i) = indx(2:k_max+1)';
   if weighted
      Dist(:,i) = so(2:k_max+1)';
   end
end
fprintf('Calculating distances on test set \n');
for i=1:nt
   dis = (X - (ones(n,1)*Xt(i,:))).^2;
   dis = sqrt(sum(dis'));
   [so indx]=sort(dis);
   Test_pointer(:,i)=indx(1:k_max)';
   if weighted
      Test_dist(:,i) = so(1:k_max)';
   end
end

% now set up and calculate stores for the proportion of classes 
% neighbouring each point.
% fprintf('Calculating neighbourhood class proportions \n');

sum_y = zeros(q,n,k_max);
sum_yt = zeros(q,nt,k_max);
% sum_y(i,j,k) records the proportion of the k nearest neighbours of.... 
% the jth data point that are of class i
% first for k equals 1
for j=1:q
   sum_yt(j,:,1)=Ytrans(Test_pointer(1,:))==j;
   sum_y(j,:,1)=Ytrans(Pointer(1,:))==j;
end
% sum_y_store checks to see if you've already calculated the 
% proportions
sum_y_store = zeros(1,k_max);
sum_y_store(1) = 1; % as you've just calculated it for k=1 above 
%....now for current value of k (if k ~= 1)....
if k > 1      
   points = Pointer(1:k,:);
   Test_points = Test_pointer(1:k,:);
   for j=1:q
      if weighted
         % extend bandwidth just past k-nn
         Bandwidth = ones(k,1)*(Dist(k,:)+realmin);
         % get tricube kernel weights
         Dist_k = (70/81)*(ones(k,n) - abs( Dist(1:k,:) ./ Bandwidth ).^3).^3;
         % calculate weighted distance
         sum_y(j,:,k)=sum((Y(points)==j).*Dist_k);
         % now do same for test points
         Test_bandwidth = ones(k,1)*(Test_dist(k,:)+realmin);
         test_k = (70/81)*(ones(k,nt) - abs(Test_dist(1:k,:)./Test_bandwidth).^3).^3;
         sum_yt(j,:,k)=sum((Y(Test_points)==j).*test_k);
      else % use uniform kernel
         sum_y(j,:,k)=sum(Y(points)==j)/k;
         sum_yt(j,:,k)=sum(Y(Test_points)==j)/k;
      end
      
   end
   sum_y_store(k) = 1; % as you've just calculated it for k 
end

% prd will store probability of categories for the training data
prb=zeros(q,n);

% convert Y to pointers in prb for use in updating
Y_pointer = Y + [0:q:q*n-1]';

% for each class find the prob.
prb = exp(beta * sum_y(:,:,k));

% now normalise
csum_prb=cumsum(prb);
prb = prb ./ (ones(q,1)*csum_prb(q,:));

% get the log likelihood
log_lik = sum(log(prb(Y_pointer)));


% set up store vectors 
beta_store=zeros(1,its);
k_store=zeros(1,its);
lik_store=zeros(1,its);

% set up t_prb, prob of classes in test set
t_prb=zeros(q,nt);
store_prb=zeros(q,nt);


fprintf('Starting MCMC \n');

% kst stores values of k visited by Markov chain.
kst = zeros(k_max,1);

count=0; % the number of samples collected 
i=0; % the number of iterations 
attempt=zeros(1,2);accept=zeros(1,2); % the accept/reject rate of the mcmc

while count < its
   i=i+1;
   
   % first copy current model
   
   k_new=k;
   
   beta_new = beta;
   
   prb_new=prb;
   
   % choose to  update k or beta with prob 0.5 
   if rand < 0.5  % then update both k and beta 
      
      flag=1;
      
      attempt(flag)=attempt(flag)+1; % make note that we are attempting change to k
      
      if rand < 0.5
         k_new = k + ceil(rand*4);
      else
         k_new = k - ceil(rand*4);
      end
      
      % reflect around boundaries
      if k_new < 1
         k_new = 2 + abs(k_new);
      elseif k_new > k_max
         k_new = k_max - (k_new-k_max);
      end
      
      % now change beta 
      beta_new = beta + sig_beta * randn;
      
      % prior
      prior_pen = -prior_prec*(beta_new^2 - beta^2);
      
      
   else % just update beta
      
      flag=2;
      attempt(2)=attempt(2)+1;
      
      beta_new = beta + sig_beta * randn;
      
      % prior
      prior_pen = -prior_prec*(beta_new^2 - beta^2);
      
   end
   
   % these next lines could be included as a constraint
   % if beta_new < 0 
   %    beta_new = -beta_new;
   % end
   
   
   % see if we need to calcluate the proportions 
   if sum_y_store(k_new)==0 % then we have not seen this value of k before, hence:
      points = Pointer(1:k_new,:);
      Test_points = Test_pointer(1:k_new,:);
      for j=1:q
         if weighted
            % extend bandwidth just past k-nn
            Bandwidth = ones(k_new,1)*(Dist(k_new,:)+realmin);
            % now calculate tricube kernel weights
            Dist_k = (70/81)*(ones(k_new,n) - abs( Dist(1:k_new,:) ./ Bandwidth ).^3).^3;
            % calculate weighted distance
            sum_y(j,:,k_new)=sum((Y(points)==j).*Dist_k);
            % now for test data 
            Test_bandwidth = ones(k_new,1)*(Test_dist(k_new,:)+realmin);
            test_k = (70/81)*(ones(k_new,nt) - abs(Test_dist(1:k_new,:)./Test_bandwidth).^3).^3;
            sum_yt(j,:,k_new)=sum((Y(Test_points)==j).*test_k);
         else % use uniform kernel
            sum_y(j,:,k_new)=sum(Y(points)==j)/k_new;
            sum_yt(j,:,k_new)=sum(Y(Test_points)==j)/k_new;
         end
      end
      % make a note that we've now been here and calculated sum_y() and sum_yt() 
      % for this value of k
      sum_y_store(k_new) = 1; 
   end
   
   % get new prob of classes
   prb_new = exp(beta_new * sum_y(:,:,k_new));
   
   % normalise 
   csum_prb=cumsum(prb_new);
   prb_new = prb_new ./ (ones(q,1)*csum_prb(q,:));
   
   % get new log likelihood
   cur_lik = sum(log(prb_new(Y_pointer)));
   
   
   % see if we accept new state  
   if rand < exp(cur_lik-log_lik+prior_pen)
      k=k_new;
      beta=beta_new;
      log_lik=cur_lik;
      prb=prb_new;
      accept(flag)=accept(flag)+1;
   end
   
   
   % Display statistics of mcmc every 100 iterations  
   if rem(i,100)==0
      fprintf('its %d: k %d  Beta %.3f Log lik %.3f k accept %d/%d beta accept %d/%d posterior samples collected %d/%d\n', i, k, ...
         beta, log_lik, accept(1),attempt(1),accept(2),attempt(2),count,its);
   end
   
   % after burn in we need to start collecting samples
   if i > burn_in & rem(i,sample_rate)==0
      count = count+1;
      beta_store(count)=beta;
      k_store(count)=k;
      lik_store(count)=log_lik;
      kst(k)=kst(k)+1;
      
      % now get test predictions
      
      t_prb = exp(beta * sum_yt(:,:,k));
      
      csum_prb=cumsum(t_prb);
      t_prb = t_prb ./ (ones(q,1)*csum_prb(q,:));
      
      store_prb=store_prb+t_prb;
      
   end
   
end % end of mcmc run 

% calculate final predictions
pred = store_prb/count;
if COSTS
   norm_prb = cost_mat * pred;
   
   [dum max_pred]=max(norm_prb);
   
   mis_class_ave = 100*(sum((max_pred' ~= Yt).*costy)/nt);
   
else 
   [dum max_pred]=max(pred);
   
   mis_class_ave = 100*(sum(max_pred' ~= Yt)/nt);
end
mc = mis_class_ave;

figure
subplot(3,1,1)
plot(beta_store);
ylabel('Interaction parameter, \beta');
title('Check these plots for non-convergence, use a larger burn in if needed');
subplot(3,1,2);
plot(k_store);
ylabel('k');
subplot(3,1,3);
plot(lik_store);
ylabel('Log probability');
xlabel('Post burn-in sample');
drawnow;


fprintf('\nProgram complete.');
if COSTS
   fprintf('\nCost weighted error rate on test set is %.3f', mc);
else
   fprintf('\nError rate on test set, using equal costs, is %.3f', mc);
end

fprintf('\nCheck acceptance rates to see if sig_beta needs to be changed.');
fprintf('\nCheck sample statistics for convergence.\n');



