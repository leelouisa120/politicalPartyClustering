function y = NegativeLogLikelihood(x)
clear probfunc;

%% 2012- 2 Issues
global subset covmin num_of_cluster;
dim = num_of_cluster;

mu = x(1:dim,:);
sigma = x(dim+1:end,:);
% 
% if any(any(mu > 7))
%     y = 10^10;
%     return;
% end

if any(any(sigma < covmin)) % check if there is any negative number in input variable
    y = 10^10;    % give a big value to the result
    return;                % return to fminsearch - do not execute the rest of the code
end

% Calculate Probability Density Function for Each Individual
for cc = 1:num_of_cluster
    mean = mu(cc,:);
    cov = sigma(cc,:);
    probfunc(:,cc) = mvnpdf(subset,mean,cov);
end

log_of_pdf = log((sum(probfunc,2))./num_of_cluster);
LogLikelihood = sum(log_of_pdf);
y = -LogLikelihood;
end