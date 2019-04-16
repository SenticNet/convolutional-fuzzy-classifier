function [test_predictions] = make_rbf_forecast(X);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is uses the stored output from "bayes_rbf_gauss.m" to make predictions on a test set
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
%        Samples from the post burn Markov chain are stored in subdirectory ./RBF_MCMC_samples
%
%			THIS PROGRAM MUST BE RUN IN THE PARENT DIRECTORY OF "./RBF_MCMC_samples"
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
cd RBF_MCMC_samples

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
   if options.lin_terms==1
      X_basis(:,1:p+1)=[ones(n,1) X];
      k_min = p+1;
   else
      k_min=1;
      X_basis(:,1) = 1; % intercept term
   end
   
   for j=k_min+1:model.k
      this_basis = model.basis(j);
      Xv = X(:,this_basis.var);
      X_basis(:,j) = gen_rbf_response(Xv,this_basis.knot,options.basis_type,0,options.basis_var);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X] = gen_rbf_response(data,knots,basis,lin_terms,var)
% This is a function to calculate radial basis function design matrix from data
%
% INPUTS
%	   data - this is original data matrix
%		knots - this is the knot locations
%     basis - basis type out of 
%							1. "Gaussian" 
% 							2. "Cubic"
%							3. "Thin-plate"
% 							4. "Linear"
% 							5. "Hardy_multi-quadric"
% 							6. "Inverse_Hardy_multi-quadric"
% 							7. "Wavelets"
%
% 		lin_terms - lin_terms=1 if you want to include linear terms in the design
% 		var  - variance parameter, needed for certain basis functions
%
%  OUTPUTS
%
%		X - the basis design matrix

[dat_no in_dim]=size(data);

[knot_no knot_dim] = size(knots);

if (knot_dim ~= in_dim) 
   error('KNOT DIMENSION AND INPUT DIMENSION MUST MATCH ');
end

% calculate distance matrix from knot points to data points
dist = zeros(knot_no,dat_no);
x_knot = sum(knots.^2,2)*ones(1,dat_no);
x_data = sum(data.^2,2)*ones(1,knot_no);
xx = knots*data';
dist = x_knot + x_data' - 2*xx;
dist_sq = dist';
dist = sqrt(dist_sq);


switch lower(basis)
	case 'thin-plate'
   	fi0 = find(dist~=0);
   	dist(fi0) = log(dist(fi0))*dist_sq(fi0);
   	X = dist;
	case 'cubic'
   	X = dist.^3;
	case 'gaussian'
   	X = exp(-var*(dist_sq));
	case 'linear'
   	X = dist;
	case 'hardy_multi-quadric'
   	X = sqrt(var + dist_sq);   
	case 'inverse_hardy_multi-quadric'
   	X = 1 ./ sqrt(var + dist_sq);
	case 'wavelet' % continuous wavelet function see Holmes and Mallick (1999), IEEE Trans. Neural Nets for details
   	X = (in_dim - (var*dist_sq)) .* exp(-0.5*(var*dist_sq)); 
	otherwise
   	error(['Unrecognised basis function ' basis]);
end

% end of function rbf_response()
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

