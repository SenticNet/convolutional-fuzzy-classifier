ccaStartup
N = 2000
settleTime = 500; 
N = N + settleTime;

  N = N + settleTime;
    X = cca_normrnd(0,1,5,N);   % 5 variables
    for i=5:N,
        X(1,i) = X(1,i) + 0.6.*X(1,i-1) + 0.65.*X(2,i-2);
        X(2,i) = X(2,i) + 0.5.*X(2,i-1) - 0.3.*X(2,i-2) - 0.3.*X(3,i-4) + 0.6.*X(4,i-1);
        X(3,i) = X(3,i) + 0.8.*X(3,i-1) - 0.7.*X(3,i-2) - 0.1.*X(5,i-3);
        X(4,i) = X(4,i) + 0.5.*X(4,i-1) + 0.9.*X(3,i-2) + 0.4.*X(5,i-2);
        X(5,i) = X(5,i) + 0.7.*X(5,i-1) - 0.5.*X(5,i-2) - 0.2.*X(3,i-1);
    end
    X = X(:,settleTime+1:end);

dlmwrite('tnn5.data',X);
