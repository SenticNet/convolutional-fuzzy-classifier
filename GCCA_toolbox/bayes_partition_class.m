function [test_set_predictions, chain_stats] = bayes_partition_mult(data,test,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a Bayesian Partition model for Classification problems, i.e. multinomial response data:
%		see chapters 7 in "Bayesian methods for nonlinear classification and regression".
%									(2002). Denison, Holmes, Mallick and Smith: published by Wiley. 
%
%  Version 1.0 	Date: 10th October 2002
%
%	Writen by Chris Holmes: c.holmes@ic.ac.uk, for academic purposes only.
%									please contact me if you intend to use this for commercial work
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	USAGE:
%			the program simulates a Bayesian Partition model assuming a multinomial reponse variable 
%			using Markov chain Monte Carlo.
%
%        If requested, samples from the post burn Markov chain are stored in subdirectory ./PARTITION_MULT_MCMC_samples
%			which is created by the program
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%		data - training data: first column is Y (the response) which is a class label. This should take values (elements) 
%									Y(i) in set {1,...,q} for a q-class classification problem.
% 									the remaining columns are the covariates, X
%									i.e. data = [Y X];
%		test - test data: same format as 'data'
%     options - this is an optional input object containing user defined settings for the 
%					program. You can run the program without this input if you wish to use 
%					our default settings. See below for further details.
%
% OUTPUTS:
%
%		all of the post burn in mcmc model samples can be stored (if wished) into directory 
%		./PARTITION_MULT_MCMC_samples which is greated by the program. See user set paramter "SAVE_SAMPLES" below.
%		Warning...this can take up a lot of space and slow down the program. If you simply wish to 
%		generate predictions on a set of values simply use the automatically generated outputs 
%		which are as follows:
%			
%     test_set_predictions is an object containing:
%				pred_store is mean prediction on the test set
%				credible_upper is 95% credible intervals above pred_store
%				credible_lower is 95% credible intervals below pred_store
%
%		chain_stats is an object containing statistics from the McMC simulation
%			..use these as a quick and dirty method to check for convergence of the chain:
%
%					chain_stats.log_lik - stores marginal log likelihood of the post burn-in samples
%					chain_stats.k - stores the dimension of the models (number of basis functions).
%				
%
%		 Feel free to add your own outputs (mail me if you have any nice changes to make to the code).
%
%
% The next lines contain our remmended default program settings:
% The input object "options" should contain the following information:
if nargin==3 % if then using user defined options
   STANDARDISE=options.standardise; % see below for interpretation
   k_max = options.k_max; % see below   
   Near_points = options.near_points;
   mcmc_samples = options.mcmc_samples; % see below
   burn_in = options.burn_in; % see below
   alpha=options.alpha;
   SAVE_SAMPLES = options.save;
   Class_regions = options.regions;
   misclassification_costs = options.costs; 
else % use default settings....these are our defaut settings
   STANDARDISE=1; options.standardise=1;% see below for detail
   k_max = 500; options.k_max = k_max;% see below   
   Near_points=1; options.near_points = Near_points;
   mcmc_samples = 20000; options.mcmc_samples = mcmc_samples; % see below
   burn_in = 10000; options.burn_in =burn_in; % see below
   alpha=1; options.alpha=alpha;  %see below
   SAVE_SAMPLES = 0; options.save=SAVE_SAMPLES; 
   Class_regions = 1; options.regions = Class_regions;
   misclassification_costs = [1 1]; options.costs = misclassification_costs;
end
% 
%		Here are the interpretation of the user set parameters:
%
% 		STANDARDISE - binary indicator {0,1}: indicates whether to standardise data set before calculating basis functions 
%						I RECOMMEND YOU SET STANDARDISE=1
%		Near_points - binary indicator {0,1} on whether centres are located as a Gaussian mixture model around the data points
%  	k_max - maximum number of basis functions, I've set it to 500 but feel free to change it
%		mcmc_samples - the number of mcmc post burn in samples; needs to be at least 1000 as a rule of thumb, but the higher the better.
%		burn_in  - the number of burn in iterations; should be at least 1000 as a rule of thumb, but ALLWAYS CHECK THE MODEL 
%						STATISTICS in object "chain_stats" TO ENSURE CONVERGENCE OF THE MCMC CHAIN!!!!!!!!!!!!!!
% 		alpha - prior parameters of Multinomial-Dirichlet ~ Dir(alpha(1),...,alpha(q)) for q class classification.
%     SAVE_SAMPLES - binary indicator: if 1 then program saves each model sample in a sub-directory
%							called PARTITION_MULT_MCMC_samples - defualt is 0 as saving slows the program
%		Class_regions - binary indicator {0,1}: is one if we wish to specify a set of q regions in the 
%					quantization where each class dominates, see section 7.4.1 in the book. That is, 
%					the first q partitions have highly informative Dirichlet priors, 
%					e.g. Dir(100, 1, 1,...,1) for first partition; Dir(1, 100, 1,...) for second, and so on
%					this forces the j'th partition, j=1,2,...,q to predict heavily in favor of the j'th class.
%		missclassification_costs - stores cost of misclassifying class: set to [1 1] if equal costs
%					else set to [c(1) c(2) ... c(q)] where c(i) is cost of misclassifying class i
%		
if SAVE_SAMPLES
   fprintf('Warning: program will save %d model files to directory ./PARTITION_MULT_MCMC_samples \n', mcmc_samples);
   dum = input('Enter 1 to proceed or any other key to exit ');
   if dum~=1
      return;
   end
end


% now extract information from program inputs.....

% response of category indicator should be stored in the first column of the data and test
Y = data(:,1);
Yt = test(:,1);

% convert classes into Y in {1,...,q}
if min(Y==0)
   Y=Y+1;
   Yt=Yt+1;
end

q = max(Y); % for a q class classification problem

% extract predictor variables
X = data(:,2:size(data,2));
Xt = test(:,2:size(test,2));

% get dimensions of data 
n = size(data,1); % number of rows
nt = size(test,1);
p = size(X,2); % number of columns

% make matrix indicator
Y_indic = zeros(n,q);
Yt_indic = zeros(nt,q);
for j=1:q
   indx = find(Y==j);
   Y_indic(indx,j) = 1;   % so Y(i,j) = 1 if the i'th data point is of class j; Y(i,j)=0 otherwise
   indx = find(Yt==j);
   Yt_indic(indx,j) = 1;
end

% see if there are unequal misclassification costs
if all(misclassification_costs==misclassification_costs(1)) % then equal costs
   COSTS=0;
else
   COSTS=1;
   % set up weighting matrices for predictions
   costy = misclassification_costs(Yt); % cost of misclassifying each test point
   cost_mat = diag(misclassification_costs); % diagonal cost matrix for weighting predictions
end

% this program uses many calls to the gammaln() function, Eq 6.2 in Denison, Holmes, Mallick and Smith
% hence it make sense to save the numbers in a look up table as this saves loads of time
% see function get_ml() at end of this file.
Gamma.Gamma_1 = gammaln([0:n]+alpha);
Gamma.Gamma_2 = gammaln([0:n]+(q*alpha));
Gamma.Gamma_prior_factor = gammaln(alpha*q)-(q*gammaln(alpha));
if Class_regions   
   Gamma.Gamma_3 = gammaln([0:n]+(100*alpha));
   Gamma.Gamma_4 = gammaln([0:n]+((100+(q-1))*alpha));
end



% make some space by clearing input data
clear data;
clear test;

if STANDARDISE
   % then standardise data
   for i=1:p
      mx(i) = mean(X(:,i));
      sx(i) = std(X(:,i));
      X(:,i) = (X(:,i)-mx(i)) / sx(i);
      Xt(:,i) = (Xt(:,i)-mx(i)) / sx(i);
   end
else
   mx=ones(1,p); sx=ones(1,p);
end


if SAVE_SAMPLES
   % create subdirectory if needed to store model samples from MCMC
   Files = dir(pwd); % get current directory
   % check to see if results directory exists
   if ~any(strcmp({Files.name},'PARTITION_MULT_MCMC_samples'))
      mkdir PARTITION_MULT_MCMC_samples;
   end
   
   cd PARTITION_MULT_MCMC_samples;
   fid = fopen('README','w'); % create a file called read me
   st = ['Multinomial Partition MCMC run on ' datestr(now)]; % and make a note of the date and time that the program was run
   fprintf(fid,'%s',st);
   fprintf(fid,'\n');
   fprintf(fid,'Program settings stored in options.mat');
   fclose(fid);
   save options options; % save the options used to run the program
   save number_of_categories.dat q -ascii; % save number of classes
   save mx mx; save sx sx; % and save the standardising factors from the training set
   cd ..
end

% get range of data
x_min = min(X);
x_range = max(X)-min(X);

% Theta will store location of the partition centres
Theta = zeros(k_max,p);
% Dist and Dist_t will store distance of centres to test points
Dist = zeros(k_max,n);
Dist_t = zeros(k_max,nt);
% p_j_given_x will store class conditional probs. for each region
p_j_given_x = zeros(k_max,q);

if Class_regions
   % we will set up model with q regions, each having high prior prob for q'th class
   k=q;
   k_min=q; % we will not allow them to be removed
   if Near_points
      indx = ceil(rand(q,1)*n); % choose a data point
      Theta(1:q,:) = X(indx,:) + 0.1*randn(q,p);
   else
      Theta(1:q,:) = (ones(q,1)*x_min) + (rand(q,p).*(ones(q,1)*x_range)); % get location of first partition
   end
else
   % we will start the MCMC chain using a model with just one parition
   k=1;
   if Near_points
      indx = ceil(rand*n); % choose a data point
      Theta(1,:) = X(indx,:) + 0.1*randn(1,p);
   else
      Theta(1,:) = x_min + (rand(1,p).*x_range); % get location of first partition
   end
end

% for first partition...
Dist(1:k,:) = get_dist(Theta(1:k,:),X); % get distance to data points
Dist_test(1:k,:) = get_dist(Theta(1:k,:),Xt); % get distance to test points


% get log marginal likelihood of current model and posterior class probs. in each region
[marg_lik,p_j_given_x]=get_ml(Y_indic,Dist(1:k,:),alpha,Class_regions,Gamma);

fprintf('Starting mcmcing...\n');

% we wish to store these for output 
chain_stats.LL_store=zeros(1,mcmc_samples); % the log marginal likelihood
chain_stats.k_store=zeros(1,mcmc_samples); % the number of basis functions
pred_store = zeros(nt,q); % the predictions on test
test_set_predictions.pred_store = zeros(nt,q); % the final predictions
test_set_predictions.credible_upper = zeros(nt,q); % the 95% credible interval above predictions
test_set_predictions.credible_lower = zeros(nt,q); % the 95% credible interval below predictions

% how many samples must we store for the mcmc 95% credible intervals?
cred_n = ceil(0.025*mcmc_samples);
cred_class_upper = zeros(nt,cred_n,q);
cred_class_lower = zeros(nt,cred_n,q);


% count and sample are loop counters within the MCMC
count=0; sample=0;

% prop and acc will store proposal and acceptance rates of the various MCMC moves
acc=zeros(1,4); prop=zeros(1,4);


% now for the main program body....

% while we don't have enough samples......keep looping
while sample < mcmc_samples
   
   % increment a counter
   count = count+1;
   
   % display statistics every 100 iterations
   if rem(count,100)==0
      fprintf('Its %d  Collected %d/%d Acc %.3f L %.3f k %d \n', count,...
         sample,mcmc_samples,sum(acc)/sum(prop),marg_lik,k);
   end
   
   % at each iteration: first make a copy of the current model
   p_j_given_x_prop = p_j_given_x;
   Dist_prop = Dist(1:k,:);
   Dist_test_prop = Dist_test(1:k,:);
   Theta_prop = Theta(1:k,:);
   k_prop = k;
   %.....anything with a _prop extension is used to denote it as a proposal
   
   % now choose a move
   birth=0;death=0;move=0;   % no move chosen yet
   u = rand; % uniform random variable on U(0,1)
   if u < 0.33
      % add a basis function
      birth=1; flag=1;
      % check for boundary, not allowed more than k_max
      if k==k_max;
         birth=0; move=1; flag=3; % make a "move" move instead
      end
   elseif u < 0.66
      % delete a basis function
      death=1;flag=2;
      % check for boundary, not allowed to delete the intercept
      if k==1;
         death=0;move=1;flag=3; % make a "move" move instead
      end
   else
      % move a basis function
      move=1; flag=3;
   end
   
   % store which move we are attempting
   prop(flag)=prop(flag)+1;
   
   % now depending on move type update the model
   if birth
      % we're adding a basis function
      k_prop=k+1;
      % generate a new location
      if Near_points % locate centre around a data point
         indx = ceil(rand*n); % choose a data point
         Theta_prop(k_prop,:) = X(indx,:) + 0.1*randn(1,p);
      else
         Theta_prop(k_prop,:) = x_min + (rand(1,p).*x_range);
      end
      Dist_prop(k_prop,:) = get_dist(Theta_prop(k_prop,:),X);
      Dist_test_prop(k_prop,:)=-1*ones(1,nt); % this marks Dist_test as needing updating if accepted
      
   elseif death
      % we've lost a basis function
      k_prop=k-1;
      % choose a basis from the model to delete
      indx = ceil(rand*k);
      % update distance matrix
      Dist_prop(indx,:)=[];
      Theta_prop(indx,:)=[];
      Dist_test_prop(indx,:)=[];
      
   elseif move 
      % choose a partition centre in the model to change
      indx = ceil(rand*k);
      if Near_points % locate centre around a data point
         indx2 = ceil(rand*n); % choose a data point
         Theta_prop(indx,:) = X(indx2,:) + 0.1*randn(1,p);
      else
         Theta_prop(indx,:) = x_min + (rand(1,p).*x_range);
      end
      Dist_prop(indx,:) = get_dist(Theta_prop(indx,:),X);
      Dist_test_prop(indx,:)=-1*ones(1,nt); % this marks Dist_test as needing updating if accepted
      
   end
   
   
   % get marginal log likelihood of proposed model and a draw of coefficients
   [marg_lik_prop p_j_given_x_prop]=get_ml(Y_indic,Dist_prop,alpha,Class_regions,Gamma);
   
   
   % now see if we accept the proposed change to the model using ratio of probabilities.
   % note that as we draw a new basis function from the prior we only need marginal likelihoods
   if rand < exp(marg_lik_prop-marg_lik)
      % we accept the proposed changes: hence update the state of the Markov chain
      k=k_prop;
      Theta(1:k,:)=Theta_prop;
      Dist(1:k,:)=Dist_prop;
      Dist_test(1:k,:)=Dist_test_prop;
      acc(flag)=acc(flag)+1; 
      marg_lik=marg_lik_prop;
      p_j_given_x(1:k,:)=p_j_given_x_prop;
   end
   
   
   if count>burn_in
      % start collecting samples
      sample=sample+1;
      
      if SAVE_SAMPLES
         % store current model in Model_store.mat
         model.theta = Theta(1:k,:);
         model.k = k;
         model.p_j_given_x = p_j_given_x(1:k,:);
         cd PARTITION_MULT_MCMC_samples;
         save(['mcmc_model_' int2str(sample)],'model');
         cd ..
      end
      
      % check for changes to model which need updating
      indx  = find(Dist_test(1:k,1)==-1);
      for jj=1:length(indx)
         Dist_test(indx(jj),:) = get_dist(Theta(indx(jj),:),Xt);
      end
      
      % find closest centre to points
      if k==1
         W_test = ones(1,nt); % W_test is a partition label; W(i)=j if the i'th data point is assigned to j'th partition
      else
         [dum W_test] =min(Dist_test(1:k,:)); % W is a partition label; W(i)=j if the i'th data point is assigned to j'th partition
      end
      
      % initialise predictions
      pred_sample = zeros(nt,q);
      for jj=1:k
         indx = find(W_test==jj);
         pred_sample(indx,:) = ones(length(indx),1)*p_j_given_x(jj,:);
      end
      
      pred_store = pred_store + pred_sample;
      chain_stats.k_store(sample)=k;
      chain_stats.LL_store(sample)=marg_lik;
      
      % now check and store credibles
      for j=1:q
         a = pred_sample(:,j);
         if sample <=cred_n % if we still have'nt filled the credible store
            cred_class_upper(:,sample,j)=a;
            cred_class_lower(:,sample,j)=a;
            if sample==cred_n % get min and max if this is the last sample to fill the store
               [min_class_cred_upper(:,j) min_class_indx(:,j)] = min(cred_class_upper(:,:,j),[],2);
               [max_class_cred_lower(:,j) max_class_indx(:,j)] = max(cred_class_lower(:,:,j),[],2);
            end
         else % we have filled up the credible store
            % check to see if any current predictions are in upper band
            find_upper = find(a > min_class_cred_upper(:,j));
            % if there are any we must insert them into the store
            for jj=1:length(find_upper)
               row_indx = find_upper(jj);
               cred_class_upper(row_indx,min_class_indx(row_indx,j),j) = a(row_indx);
               % recalculate the minimal upper and it's index
               [min_class_cred_upper(row_indx,j) min_class_indx(row_indx,j)] = min(cred_class_upper(row_indx,:,j));
            end
            % now do the same for the lower credibles.....
            % check to see if any in lower band
            find_lower = find(a < max_class_cred_lower(:,j));
            for jj=1:length(find_lower)
               row_indx = find_lower(jj);
               cred_class_lower(row_indx,max_class_indx(row_indx,j),j) = a(row_indx);
               [max_class_cred_lower(row_indx,j) max_class_indx(row_indx,j)] = max(cred_class_lower(row_indx,:,j));
            end   
         end
      end
      % check test error and display every 100
      if rem(sample,100)==0
         if COSTS
            norm_prb = pred_store*cost_mat; % weight the predictions
            [dum max_pred]=max(norm_prb,[],2); % find most probable
            mis_class_ave = 100*(sum((max_pred ~= Yt).*costy)/nt); % get cost weighted error
            fprintf('Test er (unequal costs) %.3f \n', mis_class_ave);
         else 
            [dum max_pred]=max(pred_store,[],2); % find most probable
            mis_class_ave = 100*(sum(max_pred ~= Yt)/nt);         
            fprintf('Test er (equal costs) %.3f \n', mis_class_ave);
         end
      end
            
   end % end the store of post burn in samples
   
   
end % end the mcmc loop

% get MCMC mean
test_set_predictions.pred_store = pred_store/sample;
% get final costed misclassification error rate
if COSTS
   norm_prb = pred_store*cost_mat; % weight the predictions
   [dum max_pred]=max(norm_prb,[],2); % find most probable
   mis_class_ave = 100*(sum((max_pred ~= Yt).*costy)/nt); % get cost weighted error
   fprintf('Final Test er (unequal costs) %.3f \n', mis_class_ave);
else 
   [dum max_pred]=max(pred_store,[],2);
   mis_class_ave = 100*(sum(max_pred ~= Yt)/nt);
   fprintf('Final Test er (equal costs) %.3f \n', mis_class_ave);
end

% calculate credibles
test_set_predictions.credible_upper = min_class_cred_upper;
test_set_predictions.credible_lower = max_class_cred_lower;

if SAVE_SAMPLES
   % save results and return to parent directory
   cd PARTITION_MCMC_samples;
   save test_set_predictions test_set_predictions;
   save chain_stats chain_stats;
   cd ..
end

% this ends main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
%
%			NEXT SECTION ARE FUNCTIONS CALLED BY THE MAIN FUNCTION.
%
%					get_ml()   -  gets marginal likelihood
%					get_dist() -  returns Euclidean distance between sets of points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [log_ML,p_q_given_x] = get_ml(Y,Dist,alpha,Class_regions,Gamma)
% 
% function to calculate marginal likelihood of Bayes multinomial-dirichlet model
%
% INPUTS	
%		Y - matrix indicator of class 
%	   Dist - the distacne from centres to points
%		alpha - prior parameter for Dirichlet
%		Class_regions - {0,1} : ==1 if first q regions are assigned to the q classes (see header in main program)
%     Gamma - object containing log gamma function look up tables
%
% OUTPUTS
%		log_ML - log marginal likelihood (up to a constant)
%		p_q_given_x -  this store posterior class probabilitis within each region

% get dimensions of model
[k n]=size(Dist);
[n q] = size(Y);

if k==1;
   % only one partition
   W=ones(1,n); % W is a partition label; W(i)=j if the i'th data point is assigned to j'th partitio
else
   % find closest centres
   [dum W]=min(Dist); % W is a partition label; W(i)=j if the i'th data point is assigned to j'th partition
end

log_ML = 0; % initialise log marginal likelihood
p_q_given_x = zeros(k,q); % initialise posterior class probabilities for the k partitions

% from Eqn. 6.2 in Denison, Holmes, Mallick and Smith
log_ML = k*Gamma.Gamma_prior_factor;


for j=1:k
   
   indx = find(W==j); % get points in the jth partition
   
   Y_j = Y(indx,:);
   n_j = length(indx); % number of points assigned to this partition
   
   if n_j==0
      % empty partition....NOT ALLOWED
      log_ML = -Inf;
      return;
   elseif n_j==1
      class_sum = Y_j;
   else    
      class_sum = sum(Y_j);
   end
   
   if Class_regions & j <=q % then assign this partition to the j'th class by upweighting prior
      a_prior = alpha*ones(1,q);
      a_prior(j) = 100*a_prior(j); % multiply j'th prior parameter by 100
      % get log marg lik
      log_ML = log_ML + sum(Gamma.Gamma_1(class_sum(1:j-1)+1)) + sum(Gamma.Gamma_1(class_sum(j+1:q)+1)) + Gamma.Gamma_3(class_sum(j)+1);
      log_ML = log_ML - Gamma.Gamma_4(n_j+1);
      p_q_given_x(j,:)=class_sum + a_prior;
   else
      p_q_given_x(j,:)=class_sum + alpha;
      log_ML = log_ML + (sum(Gamma.Gamma_1(class_sum+1)) - Gamma.Gamma_2(n_j+1));
   end
   
   
end

% this ends the function get_ml()
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist] = get_dist(centre,data)
% This calculates the distance between objects
% INPUTS
%		centre - points to calcluate distance to...
%	   data - this is original data matrix
%
%  OUTPUTS
%
%		dist - Euclidean distance between centre and data

[n p]=size(data);

[nc pc] = size(centre);

if (p ~= pc) 
   error('CENTRE DIMENSION AND INPUT DIMENSION MUST MATCH ');
end

% calculate distance matrix from knot points to data points
dist = zeros(nc,n);
x_knot = sum(centre.^2,2)*ones(1,n);
x_data = sum(data.^2,2)*ones(1,nc);
xx = centre*data';
dist = x_knot + x_data' - 2*xx;
dist_sq = dist;
dist = sqrt(dist_sq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
