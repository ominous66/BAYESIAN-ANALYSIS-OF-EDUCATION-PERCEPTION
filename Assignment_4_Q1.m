%% Data preparation
df = xlsread("EPAA_PewTeachMode.xlsx");

Response = df(:,1);	
Age100 = df(:,2);
Income100k = df(:,3);
SqIncome = df(:,4);	
OnlineCourse = df(:,5);	
Agelth65Enroll = df(:,6);	
Female = df(:,7);	
PostGrad = df(:,8);
CollGrad = df(:,9);	
BelowGrad = df(:,10);	
HSnbelow = df(:,11);	
Fulltime = df(:,12);	
Parttime = df(:,13);	
Unemployed = df(:,14);	
Black = df(:,15);	
White = df(:,16);	
OtherRace = df(:,17);	
Urban = df(:,18);	
Suburban = df(:,19);	
Rural = df(:,20);	
Northeast = df(:,21);	
Midwest = df(:,22);	
West = df(:,23);	
South = df(:,24);
intercept = zeros(size(South));

x = [intercept Age100 Income100k SqIncome OnlineCourse Agelth65Enroll Female PostGrad CollGrad BelowGrad Fulltime Parttime White Black Urban Suburban Northeast West South];
y = Response;
n = size(x,1);
k = size(x,2);
%% Priors
mcmc = 6;
burn = 2500;
beta0 = ones(k,1);
B0 = 1000*eye(k,k);
z = zeros(n,mcmc);
z(:,1) = y;

beta = ones(k,mcmc);
beta(:,1) = ones(19,1);
w = ones(n,mcmc);

p = 0.1;
%% MCMC Algorithm
h = waitbar(0, 'Starting simulations...');
tic

for j = 2:mcmc
    beta(:,j) = draw_beta(x,z,p,B0,beta0);
    tau = sqrt(2/(p*(1-p)));
    theta = (1-2*p)/(p*(1-p));
     
    for i = 1:n
        w(i) = draw_wi(x,z,0.1, beta(:,j), i);
        if y(i) == 0
            z(i) = truncated_normal(x(i,:)*beta(:,j) + theta * w(i), tau * tau * w(i), -Inf, 0);
        elseif y(i) == 1
            z(i) = truncated_normal(x(i,:)*beta(:,j) + theta * w(i), tau * tau * w(i), 0, Inf);
        end
        
        
    end
     waitbar(j/mcmc, h, sprintf('Iteration %d of %d | Elapsed: %.2f sec', ...
            j, mcmc, toc));
end
close(h);


%% Helper Functions
function beta = draw_beta(x,z,p, B0, beta0)
    theta = (1-2*p)/(p*(1-p));
    tau = sqrt(2/(p*(1-p)));
    n = size(x,1);
    k = size(x,2);
    sum1 = zeros(k);
    sum2 = zeros(k,1);

    for i = 1:n
        sum1 = sum1 + (x(i,:)' * x(i,:))/(exprnd(1) * tau * tau) ;
    end

    B = (sum1 + B0\eye(k))\eye(k);
    for i = 1:n
        sum2 = sum2 + (x(i,:)' * (z(i) - theta*exprnd(1)))/(exprnd(1) * tau * tau) ;
    end
    betaTilde = B * (sum2 + (B0\eye(k))*beta0);
    beta = mvnrnd(betaTilde, B); 
end

function w = draw_wi(x,z,p, beta, i)
    theta = (1-2*p)/(p*(1-p));
    tau = sqrt(2/(p*(1-p)));
    n = size(x,1);
    k = size(x,2);

    lambda = ((z(i) - x(i,:)*beta)/(tau))*((z(i) - x(i,:)*beta)/(tau));
    eeta = ((theta*theta)/(tau*tau) + 2);
    w = gigrnd(0.5, lambda, eeta,1);
end

function samples = draw_truncated_normal(mu, sigma, a, b)
    % Function to generate n random samples from a truncated normal distribution.
    %
    % Inputs:
    % - mu: Mean of the normal distribution.
    % - sigma: Standard deviation of the normal distribution.
    % - a: Lower bound of truncation.
    % - b: Upper bound of truncation.
    % - n: Number of samples to generate.
    %
    % Output:
    % - samples: Vector of n samples drawn from the truncated normal distribution.

    % Initialize the samples array
    samples = 0;

    % Calculate the CDF values for the bounds
    cdf_a = normcdf(a, mu, sigma);
    cdf_b = normcdf(b, mu, sigma);
    
    % Acceptance-Rejection method
    for i = 1:1
        accepted = false;
        while ~accepted
            % Step 1: Sample from the normal distribution
            x = mu + sigma * randn;  % Sample from normal distribution N(mu, sigma^2)
            
            % Step 2: Check if the sample is within the bounds [a, b]
            if x >= a && x <= b
                % Step 3: Calculate the acceptance criterion
                u = rand;  % Uniform random number
                acceptance_prob = (normcdf(x, mu, sigma) - cdf_a) / (cdf_b - cdf_a);  % Probability of acceptance
                
                % Step 4: Accept or reject the sample based on the acceptance probability
                if u < acceptance_prob
                    samples = x;  % Accept the sample
                    accepted = true;  % Exit the loop
                end
            end
        end
    end
end


















