clear all;
close all;
% Reading Data from Excel Sheet (Party and 3 Issues)
num_of_cluster = input('Number of Clusters: ');
filename = 'ANES 2012-1990 Data.xlsx';
sheet = 4;

% Storing Range of Data Values for Each Year (Including Party)
xlRange2012 = 'E5:H5918     ';  % 2012
% xlRange2008 = 'E5920:H8241  ';  % 2008
% xlRange2004 = 'E8243:H9454  ';  % 2004
% xlRange2000 = 'E9456:H10466 ';  % 2000
% xlRange1998 = 'E10468:H11748';  % 1998
% xlRange1996 = 'E11750:H13463';  % 1996
% xlRange1994 = 'E13465:H15259';  % 1994
% xlRange1992 = 'E15261:H17745';  % 1992
% xlRange1990 = 'E17747:H19726';  % 1990

% Store each Range(char) into a cell vector (convert from char to
% str)
% xlRange_year = [xlRange2012; xlRange2008; xlRange2004; xlRange2000; xlRange1998; xlRange1996; xlRange1994; xlRange1992; xlRange1990];
% xlRange_year = [xlRange2012; xlRange2008; xlRange2004];
xlRange_year = [xlRange2012];
%-23139.5
xlRange_year = cellstr(xlRange_year);

Year = [2012,2008,2004,2000,1998,1996,1994,1992,1990];

for ii = 1:length(xlRange_year)
    xlRange = char(xlRange_year(ii,:));   % Convert from cell vector back to char
    data = xlsread(filename,sheet,xlRange);
%     m = 1;
%     n = 1;
%     
%     for jj = 1:size(data,1)
%         if data(jj,1) < 3 && data(jj,1) ~= 0  % Democrat and Strong Democrat
%             democrat(m,:) = data(jj,:);
%             m = m+1;
%         elseif data(jj,1) > 5 % Republican and Strong Republican
%             republican(n,:) = data(jj,:);
%             n = n+1;
%         end
%     end
%     mu_democrat = mean(democrat(:,2:end));
%     s_democrat(ii) = std2(democrat(:,2:end));
%     mu_republican = mean(republican(:,2:end));
%     s_republican(ii) = std2(republican(:,2:end));
    
    subset = data(:,2:end);    % All columns but first are political issues
    partyscale = data(:,1);    % first column is party affiliation
    
%     Find the Gaussian Fit Model
    options = statset('Display','final','MaxIter',1000,'TolFun',1e-9);
        obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
    
    % Show Parameter Values
    model_means = obj.mu;
    model_cov = obj.Sigma;

%     mu = repmat(model_means(1,:),size(subset,1),1);
%     sigma = repmat([model_cov(1,1,1), model_cov(2,2,1), model_cov(3,3,1)],size(subset,1),1);
%     probfunc = normpdf(subset,mu,sigma);
for cc = 1:num_of_cluster
    mu = model_means(cc,:);
    sigma = model_cov(:,:,cc);
    probfunc(:,cc) = mvnpdf(subset,mu,sigma);
%     for kk = 1:size(subset,1);
%         probfunc(kk,:) = mymvnpdf(subset(kk,:),mu,sigma);
%     end

end

log_of_pdf = log((sum(probfunc,2))./num_of_cluster);
loglikelihood = sum(log_of_pdf)
end