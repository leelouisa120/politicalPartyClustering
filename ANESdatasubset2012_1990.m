clear all;
close all;

mystring = sprintf(['Choice of Program - \n'...
    '1. Compare years\n',...
    '2. Run all\n',...
    '3. Distance btw clusters\n'...
    '4. Political party means\n'...
    '5. Variance\n'...
    '6. Vote Predictions\n'...
    '7. Distance from center): \n']);

program = input(mystring);

switch(program)
    case 1 % Compare the GMM clusters of each year between 2012-1990
        % Reading Data from Excel Sheet (Party and 3 Issues)
        num_of_cluster = input('Number of Clusters: ');
        filename = 'ANES 2012-1990 Data.xlsx';
        sheet = 4;
        
        % Storing Range of Data Vales for Each Year (Including Party)
        xlRange2012 = 'E5:H5918     ';  % 2012
        xlRange2008 = 'E5920:H8241  ';  % 2008
        xlRange2004 = 'E8243:H9454  ';  % 2004
        xlRange2000 = 'E9456:H10466 ';  % 2000
        xlRange1998 = 'E10468:H11748';  % 1998
        xlRange1996 = 'E11750:H13463';  % 1996
        xlRange1994 = 'E13465:H15259';  % 1994
        xlRange1992 = 'E15261:H17745';  % 1992
        xlRange1990 = 'E17747:H19726';  % 1990
        
        % Store each Range(char) into a cell vector (convert from char to
        % str)
        xlRange_year = [xlRange2012; xlRange2008; xlRange2004; xlRange2000; xlRange1998; xlRange1996; xlRange1994; xlRange1992; xlRange1990];
%       xlRange_year = [xlRange2012; xlRange2008; xlRange2004];
        xlRange_year = cellstr(xlRange_year);
        
        Year = [2012,2008,2004,2000,1998,1996,1994,1992,1990];
        figure;
        for ii = 1:length(xlRange_year)
            xlRange = char(xlRange_year(ii,:));   % Convert from cell vector back to char
            data = xlsread(filename,sheet,xlRange);
            
            subset = data(:,2:end);    % All columns but first are political issues
            partyscale = data(:,1);    % first column is party affiliation
            
            % Finding Correlation Matrix for Data
            % [RHO,PVAL] = corr(subset);
            
            % Find the Gaussian Fit Model
            options = statset('Display','final','MaxIter',500,'TolFun',1e-9);
        obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
            
            % Show Parameter Values
            model_means = obj.mu;
            clustermeans(:,:,ii) = obj.mu;
            model_cov = obj.Sigma;
            cluster_cov(:,ii) = {obj.Sigma};
            
            AIC(ii) = obj.AIC;
            BIC(ii) = obj.BIC;
            
            % Find the distance between each Cluster Mean
            if num_of_cluster == 2
                dist_btw_cluster(ii) = norm(model_means(1,:)-model_means(2,:));
            end
            
            % Plotting cluster means onto PCA axes
            mu = mean(subset); % find mean of data
            
            % Adding Noise to the repeated rows
            [~,unique_rows,~] = unique(subset,'rows');
            duplicate_rows = setdiff(1:size(subset,1),unique_rows);
            subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
            
            % Conducting PCA Analysis on the data
            [coefs,score] = pca(subset);
            
            h(ii) = subplot(3,3,ii);
            
            % Plotting Party Affiliation through Color
            for jj = 1:length(partyscale)
                color = mapcolor(partyscale(jj));
                plot(score(jj,1),score(jj,2),'o','Color',color);
                hold on;
            end
            title(sprintf('%d PCA',Year(ii)));
            xlabel('Component 1');
            ylabel('Component 2');
            hold on;
            
            % find model mean in the centered space
            for kk = 1:size(model_means,1)
                mean_centered(kk,:) = model_means(kk,:) - mu;
            end
            
            % now mapp centered model_mean to the pca components
            mean_mapped = mean_centered/coefs';
            plot(mean_mapped(:,1), mean_mapped(:,2), 'kx','LineWidth',4,'MarkerSize',20);
            hold on;
        end
        linkaxes(h,'xy');
        
    case 2 % GMM clustering for all data points from 2012-1990
        % Reading Data from Excel Sheet
        num_of_cluster = input('Number of Clusters: ');
        filename = 'ANES 2012-1990 Data.xlsx';
        sheet = 5;
        xlRange = 'E5:H19718';
        data = xlsread(filename,sheet,xlRange);
        
        subset = data(:,2:end); % All columns but first are political issues
        partyscale = data(:,1);    % first column is party affiliation
        
        % Finding Correlation Matrix for Data
        [RHO,PVAL] = corr(subset);
        
        % Find the Gaussian Fit Model
        options = statset('Display','final','MaxIter',500,'TolFun',1e-9);
        obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
        
        % Show Parameter Values
        model_means = obj.mu;
        model_cov = obj.Sigma;
        
        % Plotting cluster means onto PCA axes
        mu = mean(subset); % find mean of data
        
        % Adding Noise to the repeated rows
        [~,unique_rows,~] = unique(subset,'rows');
        duplicate_rows = setdiff(1:size(subset,1),unique_rows);
        subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
        
        % Conducting PCA Analysis on the data
        [coefs,score] = pca(subset);
        
        % Plotting Party Affiliation through Color
        for ii = 1:length(partyscale)
            color = mapcolor(partyscale(ii));
            plot(score(ii,1),score(ii,2),'o','Color',color);
            hold on;
        end
        
        title('2012-1990 Principle Component Analysis');
        xlabel('Component 1');
        ylabel('Component 2');
        hold on;
        
        % find model mean in the centered space
        for kk = 1:size(model_means,1)
            mean_centered(kk,:) = model_means(kk,:) - mu;
        end
        
        % now mapp centered model_mean to the pca components
        mean_mapped = mean_centered/coefs';
        plot(mean_mapped(:,1), mean_mapped(:,2), 'kx','LineWidth',4,'MarkerSize',20);
        
    case 3 % Find Distances between the GMM found clusters and Compare over time
        % Reading Data from Excel Sheet (Party and 3 Issues)
        num_of_cluster = input('Number of Clusters: ');
        filename = 'ANES 2012-1990 Data.xlsx';
        sheet = 4;
        
        % Storing Range of Data Vales for Each Year (Including Party)
        xlRange2012 = 'E5:H5918     ';  % 2012
        xlRange2008 = 'E5920:H8241  ';  % 2008
        xlRange2004 = 'E8243:H9454  ';  % 2004
        xlRange2000 = 'E9456:H10466 ';  % 2000
        xlRange1998 = 'E10468:H11748';  % 1998
        xlRange1996 = 'E11750:H13463';  % 1996
        xlRange1994 = 'E13465:H15259';  % 1994
        xlRange1992 = 'E15261:H17745';  % 1992
        xlRange1990 = 'E17747:H19726';  % 1990
        
        % Store each Range(char) into a cell vector (convert from char to
        % str)
        xlRange_year = [xlRange2012; xlRange2008; xlRange2004; xlRange2000; xlRange1998; xlRange1996; xlRange1994; xlRange1992; xlRange1990];
        %        xlRange_year = [xlRange2012; xlRange2008; xlRange2004];
        xlRange_year = cellstr(xlRange_year);
        
        Year = [2012,2008,2004,2000,1998,1996,1994,1992,1990];
        
        for ii = 1:length(xlRange_year)
            xlRange = char(xlRange_year(ii,:));   % Convert from cell vector back to char
            data = xlsread(filename,sheet,xlRange);
            
            subset = data(:,2:end);    % All columns but first are political issues
            partyscale = data(:,1);    % first column is party affiliation
            
            % Finding Correlation Matrix for Data
            % [RHO,PVAL] = corr(subset);
            
            % Find the Gaussian Fit Model
            options = statset('Display','final','MaxIter',500,'TolFun',1e-9);
        obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
            
            % Show Parameter Values
            model_means = obj.mu;
            model_cov = obj.Sigma;
            
            AIC(ii) = obj.AIC;
            BIC(ii) = obj.BIC;
            
            % Find the distance between each Cluster Mean
            if num_of_cluster == 2
                dist_btw_cluster(ii) = norm(model_means(1,:)-model_means(2,:));
            end
        end
        if num_of_cluster == 2
            plot(Year,dist_btw_cluster);
            title('Distance between Cluster Means vs. Time');
            xlabel('Year');
            ylabel('Distance');
        end
        
    case 4 % Calculating Political Party Means and Distances vs. Time for Cluster and Party Means
        % Reading Data from Excel Sheet (Party and 3 Issues)
        num_of_cluster = input('Number of Clusters: ');
        filename = 'ANES 2012-1990 Data.xlsx';
        sheet = 4;
        
        % Storing Range of Data Values for Each Year (Including Party)
        xlRange2012 = 'E5:H5918     ';  % 2012
        xlRange2008 = 'E5920:H8241  ';  % 2008
        xlRange2004 = 'E8243:H9454  ';  % 2004
        xlRange2000 = 'E9456:H10466 ';  % 2000
        xlRange1998 = 'E10468:H11748';  % 1998
        xlRange1996 = 'E11750:H13463';  % 1996
        xlRange1994 = 'E13465:H15259';  % 1994
        xlRange1992 = 'E15261:H17745';  % 1992
        xlRange1990 = 'E17747:H19726';  % 1990
        
        % Store each Range(char) into a cell vector (convert from char to
        % str)
        xlRange_year = [xlRange2012; xlRange2008; xlRange2004; xlRange2000; xlRange1998; xlRange1996; xlRange1994; xlRange1992; xlRange1990];
        %        xlRange_year = [xlRange2012; xlRange2008; xlRange2004];
        %        xlRange_year = [xlRange2012];
        
        xlRange_year = cellstr(xlRange_year);
        
        Year = [2012,2008,2004,2000,1998,1996,1994,1992,1990];
        figure;
        
        for ii = 1:length(xlRange_year)
            xlRange = char(xlRange_year(ii,:));   % Convert from cell vector back to char
            data = xlsread(filename,sheet,xlRange);
            
            % Sorting Individuals by Party Identification
            m = 1;
            n = 1;
            
            for jj = 1:size(data,1)
                if data(jj,1) < 3 && data(jj,1) ~= 0  % Democrat and Strong Democrat
                    democrat(m,:) = data(jj,:);
                    m = m+1;
                elseif data(jj,1) > 5 % Republican and Strong Republican
                    republican(n,:) = data(jj,:);
                    n = n+1;
                end
            end
            
            % Calculating party means
            mu_democrat = mean(democrat(:,2:end));
            s_democrat(ii) = std2(democrat(:,2:end));
            mu_republican = mean(republican(:,2:end));
            s_republican(ii) = std2(republican(:,2:end));
            
            subset = data(:,2:end);    % All columns but first are political issues
            partyscale = data(:,1);    % first column is party affiliation
            
            % Finding Correlation Matrix for Data
            % [RHO,PVAL] = corr(subset);
            
            % Find the Gaussian Fit Model
            options = statset('Display','final','MaxIter',500,'TolFun',1e-9);
        obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
            
            % Show Parameter Values
            model_means = obj.mu;
            clustermeans(:,:,ii) = obj.mu;
            model_cov = obj.Sigma;
            cluster_cov(:,ii) = {obj.Sigma};
            
            AIC(ii) = obj.AIC;
            BIC(ii) = obj.BIC;
            
            % cluster and posterior probablity of each instance
                % note that1: [~,clustIdx] = max(p,[],2)
                [clustInd,~,p] = cluster(obj, subset);
                tabulate(clustInd)
                
            % Find the distance between each Cluster Mean and centers of
            % opinion of democrats and republicans
            if num_of_cluster == 2
                dist_btw_cluster(ii) = norm(model_means(1,:)-model_means(2,:));
                dist_btw_party(ii) = norm(mu_democrat - mu_republican);
            end
            
            % Plotting cluster means onto PCA axes
            mu = mean(subset); % find mean of data
            
            % Adding Noise to the repeated rows
            [~,unique_rows,~] = unique(subset,'rows');
            duplicate_rows = setdiff(1:size(subset,1),unique_rows);
            subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
            
            %             for jj = 1:length(partyscale)
            %                 color = mapcolor(partyscale(jj));
            %                 partycolor(jj,:) = color;
            %
            %             end
            %
            %             % Plotting Points in 3D Scatter Plot
            %             scatter3(subset(:,1),subset(:,2),subset(:,3),25,partycolor)
            %             hold on;
            %             plot3(model_means(:,1),model_means(:,2),model_means(:,3), 'gx','LineWidth',4,'MarkerSize',20);
            %             hold on;
            %             plot3(mu_democrat(:,1), mu_democrat(:,2),mu_democrat(:,3), 'b^', 'LineWidth', 3, 'MarkerSize', 10);
            %             hold on;
            %             plot3(mu_republican(:,1), mu_republican(:,2),mu_republican(:,3), 'r^', 'LineWidth', 3, 'MarkerSize', 10);
            
            % Conducting PCA Analysis on the data
            [coefs,score] = pca(subset);
            
            h(ii) = subplot(3,3,ii);
            
            % Plotting Party Affiliation through Color
            for jj = 1:length(partyscale)
                color = mapcolor(partyscale(jj));
                plot(score(jj,1),score(jj,2),'o','Color',color);
                hold on;
            end
            title(sprintf('%d PCA',Year(ii)));
            xlabel('Component 1');
            ylabel('Component 2');
            hold on;
            
            % Find model mean in the centered space
            for kk = 1:size(model_means,1)
                mean_centered(kk,:) = model_means(kk,:) - mu;
            end
            mu_democrat_centered = mu_democrat - mu;
            mu_republican_centered = mu_republican - mu;
            
            % Mapp centered model_mean to the pca components
            mean_mapped = mean_centered/coefs';
            mu_democrat_mapped = mu_democrat_centered/coefs';
            mu_republican_mapped = mu_republican_centered/coefs';
            
            plot(mean_mapped(:,1), mean_mapped(:,2), 'kx','LineWidth',4,'MarkerSize',20);
            plot(mu_democrat_mapped(:,1), mu_democrat_mapped(:,2), 'b^', 'LineWidth', 3, 'MarkerSize', 12,'MarkerEdgeColor','k','MarkerFaceColor','b');
            plot(mu_republican_mapped(:,1), mu_republican_mapped(:,2), 'r^', 'LineWidth', 3, 'MarkerSize', 12,'MarkerEdgeColor','k','MarkerFaceColor','r');
            hold on;
        end
        linkaxes(h,'xy');
        
        if num_of_cluster == 2
            figure;
            plot(Year,dist_btw_cluster);
            title('Distance vs. Time');
            xlabel('Year');
            ylabel('Distance');
            hold on;
            plot(Year, dist_btw_party);
            legend('distance btw cluster means','distance btw party means');
        end
        
        % Calculate ratio of distance between party and stdev of each party
        % to see if distance is significantly large (a ratio > 2 is
        % considered pretty large)
        distanceratio(1,:) = dist_btw_party./s_democrat;
        distanceratio(2,:) = dist_btw_party./s_republican;
        
    case 5 % Calculating StDev of each Political Issue to see which Issue separates people more
        % Reading Data from Excel Sheet (Party and 3 Issues)
        filename = 'ANES 2012-1990 Data.xlsx';
        sheet = 4;
        
        % Storing Range of Data Vales for Each Year (Including Party)
        xlRange2012 = 'E5:H5918     ';  % 2012
        xlRange2008 = 'E5920:H8241  ';  % 2008
        xlRange2004 = 'E8243:H9454  ';  % 2004
        xlRange2000 = 'E9456:H10466 ';  % 2000
        xlRange1998 = 'E10468:H11748';  % 1998
        xlRange1996 = 'E11750:H13463';  % 1996
        xlRange1994 = 'E13465:H15259';  % 1994
        xlRange1992 = 'E15261:H17745';  % 1992
        xlRange1990 = 'E17747:H19726';  % 1990
        
        % Store each Range(char) into a cell vector (convert from char to
        % str)
        xlRange_year = [xlRange2012; xlRange2008; xlRange2004; xlRange2000; xlRange1998; xlRange1996; xlRange1994; xlRange1992; xlRange1990];
        %        xlRange_year = [xlRange2012; xlRange2008; xlRange2004];
        xlRange_year = cellstr(xlRange_year);
        
        Year = [2012,2008,2004,2000,1998,1996,1994,1992,1990];
        
        for ii = 1:length(xlRange_year)
            xlRange = char(xlRange_year(ii,:));   % Convert from cell vector back to char
            data = xlsread(filename,sheet,xlRange);
            disp(xlRange);
            
            subset = data(:,2:end);    % All columns but first are political issues
            partyscale = data(:,1);    % first column is party affiliation
            
            % Find Standard Deviation of all the individuals
            S(ii,:) = std(subset);
            
        end
        plot(Year, S(:,1));
        hold on;
        plot(Year, S(:,2));
        hold on;
        plot(Year, S(:,3));
        title('Standard Deviation vs. Year of 3 Issues');
        xlabel('Year');
        ylabel('Standard Deviation');
        legend('Issue 1: Guaranteed Jobs and Income','Issue 2: Aid to Minority', 'Issue 3: Gov. Service-Spending');
        
    case 6 % Predict which party the individual will vote for using distance from party mean
        % Reading Data from Excel Sheet (Party and 3 Issues)
        num_of_cluster = input('Number of Clusters: ');
        filename = 'ANES 2012-1990 Data.xlsx';
        sheet = 4;
        
        % Storing Range of Data Values for Each Year (Including Party)
        xlRange2012 = 'E5:H5918     ';  % 2012
        xlRange2008 = 'E5920:H8241  ';  % 2008
        xlRange2004 = 'E8243:H9454  ';  % 2004
        xlRange2000 = 'E9456:H10466 ';  % 2000
        xlRange1998 = 'E10468:H11748';  % 1998
        xlRange1996 = 'E11750:H13463';  % 1996
        xlRange1994 = 'E13465:H15259';  % 1994
        xlRange1992 = 'E15261:H17745';  % 1992
        xlRange1990 = 'E17747:H19726';  % 1990
        
        % Store each Range(char) into a cell vector (convert from char to
        % str)
%         xlRange_year = [xlRange2012; xlRange2008; xlRange2004; xlRange2000; xlRange1998; xlRange1996; xlRange1994; xlRange1992; xlRange1990];
        %        xlRange_year = [xlRange2012; xlRange2008; xlRange2004];
               xlRange_year = [xlRange2012];
        
        xlRange_year = cellstr(xlRange_year);
        
        Year = [2012,2008,2004,2000,1998,1996,1994,1992,1990];
        figure;
        
        for ii = 1:length(xlRange_year)
            disp(ii);
            xlRange = char(xlRange_year(ii,:));   % Convert from cell vector back to char
            data = xlsread(filename,sheet,xlRange);
            
            % Sorting individuals by Party Identification
            m = 1;
            n = 1;
            
            for jj = 1:size(data,1)
                if data(jj,1) < 3 && data(jj,1) ~= 0  % Democrat and Strong Democrat
                    democrat(m,:) = data(jj,:);
                    m = m+1;
                elseif data(jj,1) > 5 % Republican and Strong Republican
                    republican(n,:) = data(jj,:);
                    n = n+1;
                end
            end
            
            % Calculating the mean of each party
            mu_democrat = mean(democrat(:,2:end));
            s_democrat(ii) = std2(democrat(:,2:end));
            mu_republican = mean(republican(:,2:end));
            s_republican(ii) = std2(republican(:,2:end));
            
            subset = data(:,2:end);    % All columns but first are political issues
            partyscale = data(:,1);    % first column is party affiliation
            
            % Plotting party means onto PCA axes
            mu = mean(subset); % find mean of data
            
            % Adding Noise to the repeated rows
            [~,unique_rows,~] = unique(subset,'rows');
            duplicate_rows = setdiff(1:size(subset,1),unique_rows);
            subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
            
            % Conducting PCA Analysis on the data
            [coefs,score] = pca(subset);
            
            h(ii) = subplot(3,3,ii);
            
            % Plotting Party Affiliation through Color and Party Vote
            % through squares(vote for democrat) or diamonds(vote for
            % republican)
            o = 1;
            p = 1;
            for jj = 1:length(partyscale)
                color = mapcolor(partyscale(jj));  % returns color of party affiliation
                dist_democrat = norm(subset(jj,:)- mu_democrat);  % Calculates distance from democratic mean
                dist_republican = norm(subset(jj,:) - mu_republican); % Calculates distance from republican mean
                if dist_democrat < dist_republican
                    plot(score(jj,1),score(jj,2),'s','Color',color);  % Plot a square
                    vote_democrat(o,:,ii) = [partyscale(jj) subset(jj,:)];  % Store the individual in vote_democrat
                    o = o+1;
                else
                    plot(score(jj,1),score(jj,2),'d','Color',color); % Plot a diamond
                    vote_republican(p,:,ii) = [partyscale(jj) subset(jj,:)]; % Store the individual in vote_republican
                    p = p+1;
                end
                hold on;
            end
            title(sprintf('%d Voting Predictions',Year(ii)));
            xlabel('Component 1');
            ylabel('Component 2');
            hold on;
            
            % Find model mean in the centered space
            mu_democrat_centered = mu_democrat - mu;
            mu_republican_centered = mu_republican - mu;
            
            % Mapp centered model_mean to the pca components
            mu_democrat_mapped = mu_democrat_centered/coefs';
            mu_republican_mapped = mu_republican_centered/coefs';
            
            % Plot Party Means on PCA axes
            plot(mu_democrat_mapped(:,1), mu_democrat_mapped(:,2), 'b^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','b');
            plot(mu_republican_mapped(:,1), mu_republican_mapped(:,2), 'r^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','r');
            hold on;
        end
        linkaxes(h,'xy');
    case 7 % Calculating Distance from Center of Cluster and Party Means
        % Reading Data from Excel Sheet (Party and 3 Issues)
        num_of_cluster = input('Number of Clusters: ');
        filename = 'ANES 2012-1990 Data.xlsx';
        sheet = 4;
        
        % Storing Range of Data Values for Each Year (Including Party)
        xlRange2012 = 'E5:H5918     ';  % 2012
        xlRange2008 = 'E5920:H8241  ';  % 2008
        xlRange2004 = 'E8243:H9454  ';  % 2004
        xlRange2000 = 'E9456:H10466 ';  % 2000
        xlRange1998 = 'E10468:H11748';  % 1998
        xlRange1996 = 'E11750:H13463';  % 1996
        xlRange1994 = 'E13465:H15259';  % 1994
        xlRange1992 = 'E15261:H17745';  % 1992
        xlRange1990 = 'E17747:H19726';  % 1990
        
        % Store each Range(char) into a cell vector (convert from char to
        % str)
        xlRange_year = [xlRange2012; xlRange2008; xlRange2004; xlRange2000; xlRange1998; xlRange1996; xlRange1994; xlRange1992; xlRange1990];
        %        xlRange_year = [xlRange2012; xlRange2008; xlRange2004];
        %        xlRange_year = [xlRange2012];
        
        xlRange_year = cellstr(xlRange_year);
        
        Year = [2012,2008,2004,2000,1998,1996,1994,1992,1990];
        
        for ii = 1:length(xlRange_year)
            xlRange = char(xlRange_year(ii,:));   % Convert from cell vector back to char
            data = xlsread(filename,sheet,xlRange);
            
            % Sorting Individuals by Party Identification
            m = 1;
            n = 1;
            
            for jj = 1:size(data,1)
                if data(jj,1) < 3 && data(jj,1) ~= 0  % Democrat and Strong Democrat
                    democrat(m,:) = data(jj,:);
                    m = m+1;
                elseif data(jj,1) > 5 % Republican and Strong Republican
                    republican(n,:) = data(jj,:);
                    n = n+1;
                end
            end
            
            % Calculating party means
            mu_democrat = mean(democrat(:,2:end));
            s_democrat(ii) = std2(democrat(:,2:end));
            mu_republican = mean(republican(:,2:end));
            s_republican(ii) = std2(republican(:,2:end));
            
            subset = data(:,2:end);    % All columns but first are political issues
            partyscale = data(:,1);    % first column is party affiliation
            
            % Find the Gaussian Fit Model
            options = statset('Display','final','MaxIter',500,'TolFun',1e-9);
            obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
            
            % Show Parameter Values
            model_means = obj.mu;
            model_means = sortrows(model_means,1);
            model_cov = obj.Sigma;
            
            AIC(ii) = obj.AIC;
            BIC(ii) = obj.BIC;
            
            % cluster and posterior probablity of each instance
            % note that1: [~,clustIdx] = max(p,[],2)
            [clustInd,~,p] = cluster(obj, subset);
            tabulate(clustInd)
            
            % Find mean of all data
            mu = mean(subset); % find mean of data
            
            % Find the distance from center to each Cluster Mean and center
            % to each party mean
            if num_of_cluster == 2
                dist_from_cluster1(ii) = norm(mu - model_means(1,:));
                dist_from_cluster2(ii) = norm(mu - model_means(2,:));
                dist_from_dem_party(ii) = norm(mu - mu_democrat);
                dist_from_rep_party(ii) = norm(mu - mu_republican);
            end
            
        end

        if num_of_cluster == 2
            figure;
            plot(Year,dist_from_cluster1,'x-','LineWidth',3,'MarkerSize',10);
            hold on;
            plot(Year,dist_from_cluster2,'x-','LineWidth',3,'MarkerSize',10);
            hold on;
            plot(Year,dist_from_dem_party,'LineWidth',3);
            hold on;
            plot(Year,dist_from_rep_party,'LineWidth',3);
            hold on;
            title('Distances from Center');
            xlabel('Year');
            ylabel('Distance from Center');
            hold on;
            legend('Cluster 1','Cluster 2','Party 1','Party 2');
        end
        
        
end


