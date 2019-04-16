function [test_set_predictions, chain_stats] = bayes_mls_gauss(data,test,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a Bayesian Piecewise multivariate linear model for Gaussian response data:
%		see chapters 3, 4 in "Bayesian methods for nonlinear classification and regression".
%									(2002). Denison, Holmes, Mallick and Smith: published by Wiley. 
%     and Holmes and Mallick (2001). "Bayesian multivariate linear splines", JRSB.
%
%  Version 1.0 	Date: 7th October 2002
%
%	Writen by Chris Holmes: c.holmes@ic.ac.uk, for academic purposes only.
%									please contact me if you intend to use this for commercial work
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	USAGE:
%			the program simulates a Bayesian Multivariate Linear Spline (MLS) model assuming a Gaussian reponse variable 
%			using Markov chain Monte Carlo.
%
%        If requested, samples from the post burn Markov chain are stored in subdirectory ./MLS_MCMC_samples
%			which is created by the program
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
%		./MLS_MCMC_samples which is greated by the program. See user set paramter "SAVE_SAMPLES" below.
%		Warning...this can take up a lot of space and slow down the program. If you simply wish to 
%		generate predictions on a set of values simply use the automatically generated outputs 
%		which are......
%			
%     test_set_predictions is an object containing:
%				pred_store is mean prediction on the test set
%				credibles is 95% credible intervals around pred_store
%				pred_beta_store is mean of local linear effects for each covariate
%				credibles is 95% credible interval of the local linear effects
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
   interaction=options.interaction; % see below
   k_max = options.k_max; % see below
   mcmc_samples = options.mcmc_samples; % see below
   burn_in = options.burn_in; % see below
   alpha_1=options.alpha_1;alpha_2=options.alpha_2;
   SAVE_SAMPLES = options.save;
   Calculate_local_linear = options.local_beta;
else % use default settings....these are our defaut settings
   STANDARDISE=1; options.standardise=1;% see below for detail
   interaction=2; options.interaction = interaction; % see below
   k_max = 200; options.k_max = k_max;% see below
   mcmc_samples = 100; options.mcmc_samples = mcmc_samples; % see below
   burn_in = 100; options.burn_in =burn_in; % see below
   alpha_1=0.1;alpha_2=0.1; options.alpha_1=alpha_1; options.alpha_2=alpha_2; %see below
   SAVE_SAMPLES = 0; options.save=SAVE_SAMPLES; 
   Calculate_local_linear = 1; options.local_beta = 1;
end
% 
%		Here are the interpretation of the user set parameters:
%
% 		STANDARDISE - binary indicator {0,1}: indicates whether to standardise data set before calculating basis functions 
%						I RECOMMEND YOU SET STANDARDISE=1
%     interaction - set of allowed number of interactions in the spline...
%						...E.g., interaction=1 gives additive model; interaction=2 allows for bivariate interaction as well as additive struture
%  	k_max - maximum number of basis functions, I've set it to 100 but feel free to change it
%		mcmc_samples - the number of mcmc post burn in samples; needs to be at least 1000 as a rule of thumb, but the higher the better.
%		burn_in  - the number of burn in iterations; should be at least 1000 as a rule of thumb, but ALLWAYS CHECK THE MODEL 
%						STATISTICS in object "chain_stats" TO ENSURE CONVERGENCE OF THE MCMC CHAIN!!!!!!!!!!!!!!
% 		alpha_1, alpha_2 - prior parameters of noise variance, 1/sig^2 ~ Gamma(alpha_1,alpha_2)
%     SAVE_SAMPLES - binary indicator: if 1 then program saves each model sample in a sub-directory
%							called MLS_MCMC_samples - defualt is 0 as saving slows the program
%     Calculate_local_linear - binary indicator {0,1}, is set to one if you wish to calculate 
%						the local linear coefficients
if SAVE_SAMPLES
   % fprintf('Warning: program will save %d model files to directory ./MLS_MCMC_samples \n', mcmc_samples);
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

% set interaction level to be at most the number of covariates
interaction = min(p,interaction);

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
   if ~any(strcmp({Files.name},'MLS_MCMC_samples'))
      mkdir MLS_MCMC_samples;
   end
   
   cd MLS_MCMC_samples;
   fid = fopen('README','w'); % create a file called read me
   st = ['MCMC run on ' datestr(now)]; % and make a note of the date and time that the program was run
   %fprintf(fid,'%s',st);
   %fprintf(fid,'\n');
   %fprintf(fid,'Program settings stored in options.mat');
   fclose(fid);
   save options options; % save the options used to run the program
   save mx mx; save sx sx; % and save the standardising factors from the training set
   cd ..
end



% we will start the MCMC chain using a model with just an intercept (constant) term
k=1; % only one basis function in the model (the intercept)
% allocate space for the design matrix
X_mls = zeros(n,k_max);
% set first column to intercept
X_mls(:,1)=1;
% and same for the design matrix of the test data
Xt_mls = zeros(nt,k_max);
Xt_mls(:,1)=1;
beta_local = zeros(nt,p); % local linear coeffs are zero

%
% Recall.....
% the model is linear in a nonlinear Mls basis set (see chapters 3,4 in book) :
%
%						Y = X_mls * beta  + epsilon
%						epsilon ~ N(0, sig2 I)
%
%						Priors:
%
%						beta ~ N(0, sig2 * prec^{-1} I)
%						sig2 ~ Gamma(alpha,beta)
%
%
% begin with a fairly flat Gaussian prior on beta ~ N(0, sig2 * prec^(-1) I)
% and noise variance
prec = 0.1; % this is the precision of the coefficient prior: will be updated as part of the program
sig2 = 1; % this is the noise variance: updated during the mcmc run.

% we need Y'*Y repeatedly so lets store it
YtY = Y'*Y;

% get log marginal likelihood of current intercept model and a draw of the basis coefficients, beta, 
% the mean of beta, and the posterior sum of squares alpha_star
[marg_lik, beta, beta_mean, alpha_star]=get_ml(X_mls(:,1:k),Y,YtY,sig2,alpha_1,alpha_2,prec);

% set basis function parameters for intercept:
% these are just dummy values for the intercept but will take 
% realistic values for each MLS basis function
basis.inter=0;
basis.knot = zeros(interaction+1,1);
basis.var = zeros(1,interaction);
basis.lr = 0;
basis.mx = 0;
basis.sx=0;
basis.beta_local =zeros(nt,p); 
%
basis_parameters(1)=basis;


% fprintf('Starting mcmcing...\n');

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

if Calculate_local_linear
   pred_beta_store = zeros(nt,p);
   test_set_predictions.pred_beta_store = zeros(nt,p); % the mean local linear effect
   test_set_predictions.beta_credibles = zeros(nt,2*p); % the 95% credibles of local linear effect
   
   % and for the mcmc 95% credible intervals for the local linear effects
   cred_beta_upper = zeros(nt,cred_n,p);
   cred_beta_lower = zeros(nt,cred_n,p);
end

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
    %  fprintf('Its %d  Collected %d/%d Acc %.3f L %.3f k %d Prec %f \n', count,...
         sample,mcmc_samples,sum(acc)/sum(prop),marg_lik,k,prec);
   end
   
   % at each iteration: first make a copy of the current model
   beta_prop = beta;
   beta_local_prop = beta_local;
   X_prop = X_mls(:,1:k);
   Xt_prop = Xt_mls(:,1:k);
   k_prop = k;
   basis_params_prop = basis_parameters;
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
         death=0;move=0;flag=3; % note move is set to zero! as we will just re-draw beta if k==1
      end
   else
      % move a basis function
      move=1; flag=3;
      if k==1
         move=0; % just re-draw coefficient
      end
   end
   
   % store which move we are attempting
   prop(flag)=prop(flag)+1;
   
   % now depending on move type update the model
   if birth
      % we're adding a basis function
      k_prop=k+1;
      % choose a random depth for the new basis function

      indx_d = ceil(rand*interaction);
      % update design matrix with a draw from a Mls basis function
      [X_prop(:,k_prop) Xt_prop(:,k_prop) basis_prop]=gen_mls_basis(X,Xt,indx_d);
      % update basis_parameters
      basis_params_prop(k_prop)=basis_prop;
      
   elseif death
      % we've lost a basis function
      k_prop=k-1;
      % choose a basis from the model to delete, NOT THE INTERCEPT THOUGH
      indx = ceil(rand*(k-1))+1;
      % update design matrix
      X_prop(:,indx)=[];
      Xt_prop(:,indx)=[];
      basis_params_prop(indx)=[];
      
   elseif move 
      % choose a basis from the model to swap with another in dictionary, not the intercept
      indx = ceil(rand*(k-1))+1;
      % choose a depth for the new basis function
      indx_d = ceil(rand*interaction);
      % update design matrix
      [X_prop(:,indx) Xt_prop(:,indx) basis_prop]=gen_mls_basis(X,Xt,indx_d);
      % update basis function parameters
      basis_params_prop(indx)=basis_prop;
   end
   
   
   % get marginal log likelihood of proposed model and a draw of coefficients
   [marg_lik_prop,beta_prop,beta_mean_prop,alpha_star_prop]=get_ml(X_prop(:,1:k_prop),Y,YtY,sig2,alpha_1,alpha_2,prec);
   
   
   % now see if we accept the proposed change to the model using ratio of probabilities.
   % note that as we draw a new basis function from the prior we only need marginal likelihoods
   if rand < exp(marg_lik_prop-marg_lik)
      % we accept the proposed changes: hence update the state of the Markov chain
      beta=beta_prop;
      beta_mean=beta_mean_prop;
      alpha_star = alpha_star_prop;
      basis_parameters = basis_params_prop;
      k=k_prop;
      X_mls(:,1:k)=X_prop;
      Xt_mls(:,1:k) = Xt_prop;
      acc(flag)=acc(flag)+1; 
      marg_lik=marg_lik_prop;
   end
   
   % update prior precision on beta every 10 iterations after first 200 mcmc its
   if  rem(count,10)==0 & count > 200 & k>1
      % get sum squared value of coefficients
      sumsq = sum(sum(beta(2:k,:).^2));
      prec = (1./(0.05+0.5*(1/sig2)*sumsq)).*randgamma_mat(0.05+0.5*(k-1),1,1);
      % prior precision has changed and hence marginal likelihood of current model has changed, so recalculate
      [marg_lik, beta, beta_mean, alpha_star]=get_ml(X_mls(:,1:k),Y,YtY,sig2,alpha_1,alpha_2,prec);
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
         model.basis = basis_parameters;
         model.beta = beta;
         model.beta_mean = beta_mean;
         model.k = k;
         cd MLS_MCMC_samples;
         save(['mcmc_model_' int2str(sample)],'model');
         cd ..
      end
      
      % get mean predictions
      a = Xt_mls(:,1:k)*beta_mean(1:k,:); % using the posterior mean of beta
      
      % store statistics
      pred_store = pred_store + a;
      chain_stats.k_store(sample)=k;
      chain_stats.LL_store(sample)=marg_lik;
      
      % store credibles
      a = Xt_mls(:,1:k)*beta(1:k,:); % using draw of beta not mean of beta
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

      
      if Calculate_local_linear % calculating local linear coefficients
         
         % first get local linear parameters
         local_beta = zeros(nt,p);
         for j=2:k
            local_beta = local_beta + beta(j)*basis_parameters(j).beta_local;
         end
         pred_beta_store = pred_beta_store + local_beta;
         % now check and store credibles
         for j=1:p
            a = local_beta(:,j);
            if sample <=cred_n % if we still have'nt filled the credible store
               cred_beta_upper(:,sample,j)=a;
               cred_beta_lower(:,sample,j)=a;
               if sample==cred_n % get min and max if this is the last sample to fill the store
                  [min_beta_cred_upper(:,j) min_beta_indx(:,j)] = min(cred_beta_upper(:,:,j),[],2);
                  [max_beta_cred_lower(:,j) max_beta_indx(:,j)] = max(cred_beta_lower(:,:,j),[],2);
               end
            else % we have filled up the credible store
               % check to see if any current predictions are in upper band
               find_upper = find(a > min_beta_cred_upper(:,j));
               % if there are any we must insert them into the store
               for jj=1:length(find_upper)
                  row_indx = find_upper(jj);
                  cred_beta_upper(row_indx,min_beta_indx(row_indx,j),j) = a(row_indx);
                  % recalculate the minimal upper and it's index
                  [min_beta_cred_upper(row_indx,j) min_beta_indx(row_indx,j)] = min(cred_beta_upper(row_indx,:,j));
               end
               % now do the same for the lower credibles.....
               % check to see if any in lower band
               find_lower = find(a < max_beta_cred_lower(:,j));
               for jj=1:length(find_lower)
                  row_indx = find_lower(jj);
                  cred_beta_lower(row_indx,max_beta_indx(row_indx,j),j) = a(row_indx);
                  [max_beta_cred_lower(row_indx,j) max_beta_indx(row_indx,j)] = max(cred_beta_lower(row_indx,:,j));
               end   
            end
         end
      end
      
      
      % check the test error and display
      pred_t = pred_store/sample;
      test_er = sum((Yt-pred_t).^2);
      if rem(sample,100)==0
         % fprintf('Test er %.3f \n', test_er);
      end

   end % end the store of post burn in samples
   
end % end the mcmc loop

% get MCMC mean
test_set_predictions.pred_store = pred_store/sample;

% check the final test error and display
pred_t = pred_store/sample;
test_er = sum((Yt-pred_t).^2);
% fprintf('Final Test er %.3f \n', test_er);

% calculate credibles
test_set_predictions.credibles = [min_cred_upper max_cred_lower];

if Calculate_local_linear
   test_set_predictions.pred_beta_store = pred_beta_store/sample; % the mean local linear effect
   test_set_predictions.beta_credibles = [min_beta_cred_upper max_beta_cred_lower]; % the 95% credibles of local linear effect
end
if SAVE_SAMPLES
   % save results and return to parent directory
   cd MLS_MCMC_samples;
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
%					gen_mls_basis()     - generates a MLS basis function
%					randgamma_mat() - draws variates from a gamma distribution
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [log_ML,beta,beta_mean,a_star] = get_ml(X,Y,YtY,sig2,a,b,prec)
% 
% function to calculate marginal likelihood of Bayes linear model, Y ~ N(X beta, sig2 I)
% with normal-inverse-gamma prior on beta, sig2 ~ NIG(0,prec I, a, b)
%
% INPUTS	
%		X - the design matrix
%		Y - the response
%		YtY - the sum squared of response values, Y'*Y
%		sig2 - a draw from the noise variance
%		a, b - prior parameters for noise variance
%		prec - precision of normal prior on beta: beta | sig2 ~ N(0, sig2 * (1/prec) * I)
%
% OUTPUTS
%		log_ML - log marginal likelihood (up to a constant)
%		beta - a draw from the posterior distribution of beta
%		beta_mean - the posterior mean vector for beta
%		a_star - the posterior sum_squares

[n p]=size(X);

% make prior precision (inverse-variance) matrix......
prior_prec = prec*eye(p);
prior_prec(1,1)=0; % improper prior on intercept (first col of X)

% calculate posterior variance covariance matrix and precision
post_P = X'*X + prior_prec;
post_V = inv(post_P);

% get posterior mean of beta
beta_mean = post_V*X'*Y;

% calculate log of the square root of determinant of post_V by using Cholesky decomposition
[R]=chol(post_V);
% this is nice as the log of square root of determinant of post_V is just the 
% sum of the log of the diagonal elements of R, where post_V = R'*R, R is upper triangular
half_log_det_post = sum(log(diag(R)));

% now calculate log of square root of determinant of prior (this is easy as prior on beta is diagonal) 
half_log_det_prior = -0.5*(p-1)*log(prec);
%.......note that we use (p-1) as we use improper prior on intecept, beta(1) ~ N(0, infinity)
%.....this does not cause (Lyndley-Bartlett paradox) problems as we allways include an intercept in the model

% now calculate posterior sum_squares
a_star = YtY - beta_mean'*post_P*beta_mean;

% finally log marginal likelihood is
log_ML = half_log_det_post - half_log_det_prior - (0.5*(n+b))*log(0.5*(a+a_star));

% Now draw a value of beta from conditional posterior distribution....
% making use of previous cholesky decomposition
Rsig2 = sqrt(sig2)*R;
Rsig2 = Rsig2';
beta =  beta_mean + Rsig2*randn(p,1);

% this ends the function get_ml()
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
function [x, xt, basis]=gen_mls_basis(X,Xt,depth);
% generates random mls basis 
% INPUTS:
%		X - data
%		interaction - number of terms-1 in the basis function
% 		order - order of basis (typically order \in {0 1}, i.e. step or linear)
%		Xt - OPTIONAL data for prediction of basis on
%
% OUTPUTS:
%		x - basis response on X
%		xt - basis response on Xt
%		basis - stores basis parameters

% get data sizes
[n p] = size(X);
nt = size(Xt,1);


% check for integrety
if depth > p
   error('Depth of interaction greater than number of variables');
end

% knot_pos stores knot position, we need an intercept as well for each basis
knot=zeros(depth+1,1);

% lr stores whether each knot is `left' or "right" facing
lr = rand>0.5;

% reponse of basis function stored in 
x = zeros(n,1);
xt = zeros(nt,1);

% now make basis function

% repeat until we get a non-zero basis function
not_finished =1;

while not_finished
   
   % choose a data point to locate knot on
   data_indx = ceil(rand*n);
   
   
   % var stores indicator of which covariates are used
   var = zeros(depth,1);
   
   for j=1:depth
      
      % choose a variable not already used
      not_ok=1;
      while not_ok
         ind = ceil(rand*p);
         if ~ismember(ind,var(1:j))
            var(j)=ind;
            not_ok=0;
         end
      end
      
   end
   
   x_v = X(data_indx,var); % data point we are centering on
   % generate rotation on basis function
   knot(1:depth) = randn(1,depth);  
   % normalise
   knot(1:depth) = knot(1:depth) / sqrt(sum(knot.^2));
   % ensure plane passes through the chosen data point
   knot(depth+1) = -x_v*knot(1:depth);
   
   % now get response
   temp  = [X(:,var) ones(n,1)]*knot;
   temp_t = [Xt(:,var) ones(nt,1)]*knot;
   
   % truncated linear spline
   if lr==0
      temp = max(0,temp);
      temp_t = max(0,temp_t);
   else
      temp = min(0,temp);
      temp_t=min(0,temp_t);
   end
   
   x = temp;
   xt = temp_t;
   
   % null basis functions
   not_finished = ~any(x);
   
end

% make a note of which test points are covered by the basis function
indx_t = xt ~=0;

% standardise function
mx = mean(x);
stx = std(x);
x = (x-mx)/stx;
xt = (xt-mx)/stx;

% store basis function paramters
basis.inter=depth;
basis.knot = knot;
basis.var = var;
basis.lr = lr;
basis.mx = mx;
basis.sx=stx;
% get local beta's
basis.beta_local = zeros(nt,p);
for j=1:depth
   basis.beta_local(indx_t,var(j))=knot(j)/stx;
end


% end of function gen_mls_basis()
%
%
%
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




