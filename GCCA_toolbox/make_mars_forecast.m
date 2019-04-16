function [test_predictions] = make_mars_forecast(X);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is uses the stored output from "bayes_mars_gauss.m" to make predictions on a test set
%
%
%  Version 1.0 	Date: 1st October 2002
%
%	Writen by Chris Holmes: c.holmes@ic.ac.uk, for academic purposes only.
%									please contact me if you intend to use this for commercial work
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	USAGE:
%			
%        Samples from the post burn Markov chain are stored in subdirectory ./Mars_mcmc_samples
%
%			THIS PROGRAM MUST BE RUN IN THE PARENT DIRECTORY OF "./MARS_MCMC_samples"
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
%		covariates X, n by p matrix, n data points, p predictors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% OUTPUTS:
%			
%     test_predictions is an object containing:
%				pred_store is mean prediction on the test set
%				credibles is 95% credible intervals around pred_store
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get data dimensions
[n p]=size(X);

pred_store = zeros(n,1); % mean predictions
test_predictions.pred_store = zeros(n,1); % the predictions stored
test_predictions.credibles = zeros(n,2); % the 95% credible interval around predictions

% change directory to where MCMC samples are stored
cd MARS_MCMC_samples

% get options
load options;

its = options.mcmc_samples;

% get standardising factors
load mx; load sx;
for j=1:p
   X(:,j) = (X(:,j)-mx(j))/sx(j);
end


% how many samples must we store for the mcmc 95% credible intervals?
cred_n = ceil(0.025*its);
cred_upper = zeros(n,cred_n);
cred_lower = zeros(n,cred_n);

sample = 0;

for i=1:its
   
   if rem(i,100)==0
      fprintf('Calculating %d / %d samples \n', i, its);
   end
   
   sample = sample+1;
   
   str = ['mcmc_model_' int2str(i)];
   load(str);
   
   X_basis = zeros(n,model.k);
   
   % now form basis model
   X_basis(:,1) = 1; % intercept term
   
   for j=2:model.k
      X_basis(:,j) = gen_mars_response(X,model.basis(j));
   end
   
   % get mean predictions
   a = X_basis*model.beta_mean; % using the posterior mean of beta
   
   pred_store = pred_store + a;
   
   % store credibles
   a = X_basis*model.beta; % using draw of beta not mean of beta
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
   
end

% get MCMC mean
test_predictions.pred_store = pred_store/its;

% calculate credibles
test_predictions.credibles = [min_cred_upper max_cred_lower];

cd ..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_base]=gen_mars_response(X,basis);
% generates random mars basis
% INPUTS:
%		X - data
%     basis - stores basis paramters

% get data sizes
n = size(X,1);

% make temp response
temp = zeros(n,basis.inter);

for j=1:basis.inter
   
   if basis.lr(j)==0
      temp(:,j) = min(0,X(:,basis.var(j))-basis.knot(j));
   else
      temp(:,j) = max(0,X(:,basis.var(j))-basis.knot(j));
   end
   
end


if basis.order==0
   temp = temp~=0;
else
   temp=temp.^basis.order;
end


% tensor product
X_base = prod(temp,2);

X_base = (X_base-basis.mx)/basis.sx;




