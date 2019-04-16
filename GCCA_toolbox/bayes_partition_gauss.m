function [test_set_predictions, chain_stats] = bayes_partition_gauss(data,test,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a Bayesian Partition model for Gaussian response data:
%		see chapters 7 in "Bayesian methods for nonlinear classification and regression".
%									(2002). Denison, Holmes, Mallick and Smith: published by Wiley. 
%
%  Version 1.0 	Date: 1st October 2002
%
%	Writen by Chris Holmes: c.holmes@ic.ac.uk, for academic purposes only.
%									please contact me if you intend to use this for commercial work
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	USAGE:
%			the program simulates a Bayesian Partition model assuming a Gaussian reponse variable 
%			using Markov chain Monte Carlo.
%
%        If requested, samples from the post burn Markov chain are stored in subdirectory ./PARTITION_MCMC_samples
%			which is created by the program
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%		data - training data: first column is Y (the response) remaining columns are 
%									covariates, X
%		test - test data: same format as 'data'
%     options - this is an optional input object containing user defined settings for the 
%					program. You can run the program without this input if you wish to use 
%					our default settings. See below for further details.
%
% OUTPUTS:
%
%		all of the post burn in mcmc model samples can be stored (if wished) into directory 
%		./PARTITION_MCMC_samples which is greated by the program. See user set paramter "SAVE_SAMPLES" below.
%		Warning...this can take up a lot of space and slow down the program. If you simply wish to 
%		generate predictions on a set of values simply use the automatically generated outputs 
%		which are......
%			
%     test_set_predictions is an object containing:
%				pred_store is mean prediction on the test set
%				credibles is 95% credible intervals around pred_store
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
   LIN = options.linear;
   Near_points = options.near_points;
   k_max = options.k_max; % see below
   mcmc_samples = options.mcmc_samples; % see below
   burn_in = options.burn_in; % see below
   alpha_1=options.alpha_1;alpha_2=options.alpha_2;
   SAVE_SAMPLES = options.save;
else % use default settings....these are our defaut settings
   STANDARDISE=1; options.standardise=1;% see below for detail
   LIN=0; options.linear=LIN;
   Near_points=1; options.near_points = Near_points;
   k_max = 200; options.k_max = k_max;% see below
   mcmc_samples = 1000; options.mcmc_samples = mcmc_samples; % see below
   burn_in = 1000; options.burn_in =burn_in; % see below
   alpha_1=0.1;alpha_2=0.1; options.alpha_1=alpha_1; options.alpha_2=alpha_2; %see below
   SAVE_SAMPLES = 1; options.save=SAVE_SAMPLES; 
end
% 
%		Here are the interpretation of the user set parameters:
%
% 		STANDARDISE - binary indicator {0,1}: indicates whether to standardise data set before calculating basis functions 
%						I RECOMMEND YOU SET STANDARDISE=1
%		LIN - binary indicator {0,1): is one if using linear model within each partition or simply an 
%				intercept.
%		Near_points - binary indicator {0,1} on whether centres are located as a Gaussian mixture model around the data points
%  	k_max - maximum number of basis functions, I've set it to 500 but feel free to change it
%		mcmc_samples - the number of mcmc post burn in samples; needs to be at least 1000 as a rule of thumb, but the higher the better.
%		burn_in  - the number of burn in iterations; should be at least 1000 as a rule of thumb, but ALLWAYS CHECK THE MODEL 
%						STATISTICS in object "chain_stats" TO ENSURE CONVERGENCE OF THE MCMC CHAIN!!!!!!!!!!!!!!
% 		alpha_1, alpha_2 - prior parameters of noise variance, 1/sig^2 ~ Gamma(alpha_1,alpha_2)
%     SAVE_SAMPLES - binary indicator: if 1 then program saves each model sample in a sub-directory
%							called PARTITION_MCMC_samples - defualt is 0 as saving slows the program
if SAVE_SAMPLES
   fprintf('Warning: program will save %d model files to directory ./PARTITION_MCMC_samples \n', mcmc_samples);
   dum = input('Enter 1 to proceed or any other key to exit ');
   if dum~=1
      return;
   end
end


% now extract information from program inputs.....

% response should be stored in the first column of the data and test
Y = data(:,1);
Yt = test(:,1);

% extract predictor variables
X = data(:,2:size(data,2));
Xt = test(:,2:size(test,2));

% get dimensions of data 
n = size(data,1); % number of rows
nt = size(test,1);
p = size(X,2); % number of columns


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
   if ~any(strcmp({Files.name},'PARTITION_MCMC_samples'))
      mkdir PARTITION_MCMC_samples;
   end
   
   cd PARTITION_MCMC_samples;
   fid = fopen('README','w'); % create a file called read me
   st = ['MCMC run on ' datestr(now)]; % and make a note of the date and time that the program was run
   fprintf(fid,'%s',st);
   fprintf(fid,'\n');
   fprintf(fid,'Program settings stored in options.mat');
   fclose(fid);
   save options options; % save the options used to run the program
   save mx mx; save sx sx; % and save the standardising factors from the training set
   cd ..
end

% get range of data
x_min = min(X);
x_range = max(X)-min(X);


% we will start the MCMC chain using a model with just one parition
k=1;
% Theta will store location of the partition centres
Theta = zeros(k_max,p);
if Near_points
   indx = ceil(rand*n); % choose a data point
   Theta(1,:) = X(indx,:) + 0.1*randn(1,p);
else
   Theta(1,:) = x_min + (rand(1,p).*x_range); % get location of first partition
end


% Dist and Dist_t will store distance of centres to test points
Dist = zeros(k_max,n);
Dist_test = zeros(k_max,nt);
% for first partition...
Dist(1,:) = get_dist(Theta(1,:),X); % get distance to data points
Dist_test(1,:) = get_dist(Theta(1,:),Xt); % get distance to test points

% begin with a fairly flat Gaussian prior on levels in partitions ~ N(0, sig2 * prec^(-1) I)
% and noise variance
prec = 0.1; % this is the precision of the coefficient prior: will be updated as part of the program
sig2 = 1; % this is the noise variance: updated during the mcmc run.


% get log marginal likelihood of current global model and a draw of the coefficients, beta, 
% the mean of beta, and the posterior sum of squares alpha_star
[marg_lik, beta, beta_mean, alpha_star]=get_ml(X,Y,Dist(1:k,:),sig2,alpha_1,alpha_2,prec,LIN);


fprintf('Starting mcmcing...\n');

% we wish to store these for output 
chain_stats.LL_store=zeros(1,mcmc_samples); % the log marginal likelihood
chain_stats.k_store=zeros(1,mcmc_samples); % the number of basis functions
pred_store = zeros(nt,1); % the predictions on test
test_set_predictions.pred_store = zeros(nt,1); % the final predictions
test_set_predictions.credibles = zeros(nt,2); % the 95% credible interval around predictions

% how many samples must we store for the mcmc 95% credible intervals?
cred_n = ceil(0.025*mcmc_samples);
cred_upper = zeros(nt,cred_n);
cred_lower = zeros(nt,cred_n);

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
      fprintf('Its %d  Collected %d/%d Acc %.3f L %.3f k %d Prec %f Noise var %f \n', count,...
         sample,mcmc_samples,sum(acc)/sum(prop),marg_lik,k,prec,sig2);
   end
   
   % at each iteration: first make a copy of the current model
   beta_prop = beta;
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
         death=0;move=0;flag=3; % make a "move" move instead
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
   [marg_lik_prop,beta_prop,beta_mean_prop,alpha_star_prop]=get_ml(X,Y,Dist_prop,sig2,alpha_1,alpha_2,prec,LIN);
   
   
   % now see if we accept the proposed change to the model using ratio of probabilities.
   % note that as we draw a new basis function from the prior we only need marginal likelihoods
   if rand < exp(marg_lik_prop-marg_lik)
      % we accept the proposed changes: hence update the state of the Markov chain
      beta=beta_prop;
      beta_mean=beta_mean_prop;
      alpha_star = alpha_star_prop;
      k=k_prop;
      Theta(1:k,:)=Theta_prop;
      Dist(1:k,:)=Dist_prop;
      Dist_test(1:k,:)=Dist_test_prop;
      acc(flag)=acc(flag)+1; 
      marg_lik=marg_lik_prop;
   end
   
   % update prior precision on beta every 10 iterations after first 200 mcmc its
   if  rem(count,10)==0 & count > 500 & k>1
      % get sum squared value of coefficients
      sumsq = sum(sum(beta(:,1:k).^2));
      % get number of coeffs
      n_coeff = size(beta(:,1:k),1)*size(beta(:,1:k),2);
      prec = (1./(0.05+0.5*(1/sig2)*sumsq)).*randgamma_mat(0.05+0.5*n_coeff,1,1);
      % prior precision has changed and hence marginal likelihood of current model has changed, so recalculate
      [marg_lik, beta, beta_mean, alpha_star]=get_ml(X,Y,Dist(1:k,:),sig2,alpha_1,alpha_2,prec,LIN);
   end
   
   
   % draw a value for the noise variance - this is needed to draw beta in function get_ml()
   % inverse variance is Gamma
   sig2_inv = (1./(0.5*(alpha_star + alpha_1)))*randgamma_mat(0.5*(n+alpha_2),1,1);
   sig2 = 1/sig2_inv;
   
   
   if count>burn_in
      % start collecting samples
      sample=sample+1;
      
      if SAVE_SAMPLES
         % store current model in Model_store.mat
         model.theta = Theta(1:k,:);
         model.beta = beta;
         model.beta_mean = beta_mean;
         model.k = k;
         cd PARTITION_MCMC_samples;
         save(['mcmc_model_' int2str(sample)],'model');
         cd ..
      end
      
      % check for changes to model which need updating
      indx  = find(Dist_test(1:k,1)==-1);
      for jj=1:length(indx)
         Dist_test(indx(jj),:) = get_dist(Theta(indx(jj),:),Xt);
      end
      
      if k==1
         W_test = ones(1,nt); % W_test stores indictors of which points are assigned to which partitions
      else
         [dum W_test] =min(Dist_test(1:k,:)); % W_test stores indictors of which points are assigned to which partitions
      end
      
      if LIN % if linear in each partition
         xt = [ones(nt,1) Xt];
      else	% using constant
         xt = ones(nt,1);
      end
      
      aa = xt*beta_mean(:,1:k); % mean prediction of test data
      aaa = xt*beta(:,1:k); % sample prediction of test data
      pred_mean = zeros(nt,1);
      pred_sample = zeros(nt,1);
      for jj=1:k
         indx = find(W_test==jj);
         pred_mean(indx) = aa(indx,jj);
         pred_sample(indx)=aaa(indx,jj);
      end
      
      
      % store statistics
      pred_store = pred_store + pred_mean;
      chain_stats.k_store(sample)=k;
      chain_stats.LL_store(sample)=marg_lik;
      
      % store credibles
      a = pred_sample; % using draw of beta not mean of beta
      if sample <=cred_n % if we still have'nt filled the credible store
         cred_upper(:,sample)=a;
         cred_lower(:,sample)=a;
         if sample==cred_n % get min and max if this is the last sample to fill the store
            [min_cred_upper min_indx] = min(cred_upper,[],2);
            [max_cred_lower max_indx] = max(cred_lower,[],2);
         end
      else % we have filled up the credible store
         % check to see if any current predictions are in upper band
         find_upper = find(a > min_cred_upper);
         % if there are any we must insert them into the store
         for j=1:length(find_upper)
            row_indx = find_upper(j);
            cred_upper(row_indx,min_indx(row_indx)) = a(row_indx);
            % recalculate the minimal upper and it's index
            [min_cred_upper(row_indx) min_indx(row_indx)] = min(cred_upper(row_indx,:));
         end
         % now do the same for the lower credibles.....
         % check to see if any in lower band
         find_lower = find(a < max_cred_lower);
         for j=1:length(find_lower)
            row_indx = find_lower(j);
            cred_lower(row_indx,max_indx(row_indx)) = a(row_indx);
            [max_cred_lower(row_indx) max_indx(row_indx)] = max(cred_lower(row_indx,:));
         end   
      end
      
      
      % check the test error and display
      pred_t = pred_store/sample;
      test_er = sum((Yt-pred_t).^2);
      if rem(sample,100)==0
         fprintf('Test er %.3f \n', test_er);
      end
      
   end % end the store of post burn in samples
   
end % end the mcmc loop

% get MCMC mean
test_set_predictions.pred_store = pred_store/sample;

% calculate credibles
test_set_predictions.credibles = [min_cred_upper max_cred_lower];

% check the final test error and display
pred_t = pred_store/sample;
test_er = sum((Yt-pred_t).^2);
fprintf('Final Test er %.3f \n', test_er);


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
%					get_ml()   -  gets marginal likelihood and draws beta
%					randgamma_mat() - draws variates from a gamma distribution
%					get_dist() - returns Euclidean distance between objects
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [log_ML,beta,beta_mean,a_star] = get_ml(X,Y,Dist,sig2,a,b,prec,Lin)
% 
% function to calculate marginal likelihood of Bayes linear model, Y ~ N(X beta, sig2 I)
% with normal-inverse-gamma prior on beta, sig2 ~ NIG(0,prec I, a, b)
%
% INPUTS	
%		X - the design matrix
%		Y - the response
%	   Dist - the distacne from centres to points
%		sig2 - a draw from the noise variance
%		a, b - prior parameters for noise variance
%		prec - precision of normal prior on beta: beta | sig2 ~ N(0, sig2 * (1/prec) * I)
%		Lin - indicator for linear or intercept
%
% OUTPUTS
%		log_ML - log marginal likelihood (up to a constant)
%		beta - a draw from the posterior distribution of beta
%		beta_mean - the posterior mean vector for beta
%		a_star - the posterior sum_squares

[n p]=size(X);

if Lin
   X = [ones(n,1) X];
   [n p]=size(X);
else
   X=ones(n,1);
   [n p]=size(X);
end

[k n]=size(Dist);

if k==1;
   W=ones(1,n); % W indicates which points are assigened to which centres
else
   [dum W]=min(Dist); % W indicates which points are assigened to which centres
end

% make prior precision (inverse-variance) matrix......
prior_prec = prec*eye(p);

beta = zeros(p,k); % intialise partition coefficients
beta_mean = zeros(p,k);

log_ML = 0; % initialise log marginal likelihood

a_star = 0; % initialise sum sq.

sum_ss =0; % sum_ss stores sum of squares

for j=1:k % for each partition
   
   indx = find(W==j); % get points in the jth partition
   
   X_j = X(indx,:); % get subset of data
   Y_j = Y(indx);
   
   
   if isempty(indx) % then partition contains no points, NOT ALLOWED
      
      log_ML = -Inf;
      return;
      
   else % partition contains data 
      
      % calculate posterior variance covariance matrix and precision
      post_P = X_j'*X_j + prior_prec;
      post_V = inv(post_P);
      
      % get posterior mean of beta
      beta_mean(:,j) = post_V*X_j'*Y_j;
      
      % calculate log of the square root of determinant of post_V by using Cholesky decomposition
      [R]=chol(post_V);
      % this is nice as the log of square root of determinant of post_V is just the 
      % sum of the log of the diagonal elements of R, where post_V = R'*R, R is upper triangular
      half_log_det_post = sum(log(diag(R)));
      
      % now calculate log of square root of determinant of prior (this is easy as prior on beta is diagonal) 
      half_log_det_prior = -0.5*p*log(prec);
      
      % now calculate posterior sum_squares
      sum_ss = sum_ss +  beta_mean(:,j)'*post_P*beta_mean(:,j);
      
      % finally log marginal likelihood is
      log_ML = log_ML + ( half_log_det_post - half_log_det_prior );
      
   end
   
   % Now draw a value of beta from conditional posterior distribution....
   % making use of previous cholesky decomposition
   Rsig2 = sqrt(sig2)*R;
   Rsig2 = Rsig2';
   beta(:,j) =  beta_mean(:,j) + Rsig2*randn(p,1);
   
end

a_star = (Y'*Y) - sum_ss;
log_ML = log_ML - (0.5*(n+b))*log(0.5*(a+a_star));

% this ends the function get_ml()
%
%
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function gamma_store=randgamma_mat(a,N,M);
% This function is used to draw gamma random variables.
% I can't remember who I nicked this off????
% if you let me know then I'll credit them, please email me.
%
% It was addapted from........
%
%
%RANDGAMM Generates N gamma random deviates.
%	RANDGAMM(A,N) is a random deviate from the standard gamma
%	distribution with shape parameter A.
%
%	B*RANDGAMM(A,N) is a random deviate from the gamma distribution
%	with shape parameter A and scale parameter B.  The distribution
%	then has mean A*B and variance A*B^2.
%
%	See RAND.
% GKS 31 July 93
% Algorithm for A >= 1 is Best's rejection algorithm XG
% Adapted from L. Devroye, "Non-uniform random variate
% generation", Springer-Verlag, New York, 1986, p. 410.
% Algorithm for A < 1 is rejection algorithm GS from
% Ahrens, J.H. and Dieter, U. Computer methods for sampling
% from gamma, beta, Poisson and binomial distributions.
% Computing, 12 (1974), 223 - 246.  Adapted from Netlib
% Fortran routine.

gamma_store = zeros(N,M);

for count_col = 1:M
   
   for count_ro = 1:N
      
      a = a(1);
      if a < 0,
         gam = NaN;
      elseif a == 0,
         gam = 0;
      elseif a >= 1,
         b = a-1;
         c = 3*a-0.75;
         accept = 0;
         while accept == 0,
            u = rand(2,1);
            w = u(1)*(1-u(1));
            y = sqrt(c/w)*(u(1)-0.5);
            gam = b+y;
            if gam >= 0,
               z = 64*w^3*u(2)^2;
               accept = ( z<=1-2*y^2/gam );
               if accept == 0,
                  if b == 0,
                     accept = ( log(z)<=-2*y );
                  else
                     accept = ( log(z)<=2*(b*log(gam/b)-y) );
                  end;
               end;
            end;
         end;
      else
         aa = 0;
         b = 1 + .3678794*a;
         accept = 0;
         while accept == 0,
            p = b*rand(1);                        
            if p < 1, 
               gam = exp(log(p)/a);
               accept = (-log(rand(1)) >= gam);
            else
               gam = -log((b-p)/a);
               accept = (-log(rand(1)) >= (1-a)*log(gam));
            end;
         end;
      end;
      
      gamma_store(count_ro,count_col) = gam;
      
   end;
   
end;




