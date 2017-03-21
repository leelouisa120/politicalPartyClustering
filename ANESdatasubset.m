clear all;
close all;

program = input('1. GMM Model, 2. Setting floor on Covariance: ');
switch(program)
    case 1
        dimension = input('Number of Variables: ');
        
        switch(dimension)
            case 2
                inputclusters = input('1 cluster or 2. multiple clusters: ');
                if inputclusters == 1
                    % Reading Data from Excel Sheet
                    num_of_cluster = input('Number of Clusters: ');
                    filename = 'ANES 2012 Data.xlsx';
                    sheet = 4;
                    xlRange = 'E5:G5918';
                    data = xlsread(filename,sheet,xlRange);
                    
                    % Sort Individuals by which Party they identified with
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
                    
                    % Find mean of each party
                    mu_democrat = mean(democrat(:,2:end));
                    s_democrat = std2(democrat(:,2:end));
                    mu_republican = mean(republican(:,2:end));
                    s_republican = std2(republican(:,2:end));
                    
                    % Separate the party identifications and opinions of other
                    % issues
                    subset = data(:,2:end);    % All columns but first are political issues
                    partyscale = data(:,1);    % first column is party affiliation
                    
                    % Finding Correlation Matrix for Data
                    [RHO,PVAL] = corr(subset);
                    
                    
                    
                    % Find the Gaussian Fit Model
                    options = statset('Display','final','MaxIter',1000,'TolFun',1e-9);
                    obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
                    
                    % Save 2 Parameter Values
                    model_means = obj.mu;
                    model_cov = obj.Sigma;
                    AIC = obj.AIC;
                    BIC = obj.BIC;
                    
                    % cluster and posterior probablity of each instance
                    % note that1: [~,clustIdx] = max(p,[],2)
                    [clustInd,~,p] = cluster(obj, subset);
                    tabulate(clustInd)
                    
                    % Number of Democrats/Republicans in Each Cluster
                    for kk = 1:num_of_cluster
                        clear cluster;
                        cluster = partyscale(clustInd==kk,:);
                        for ll = 1:7
                            count(ll,kk) = sum(cluster == ll);
                        end
                    end
                    
                    % Adding Noise to the repeated rows
                    [~,unique_rows,~] = unique(subset,'rows');
                    duplicate_rows = setdiff(1:size(subset,1),unique_rows);
                    subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
                    
                    % Indicating Party Affiliation through Color
                    for ii = 1:length(partyscale)
                        color = mapcolor(partyscale(ii));
                        partycolor(ii,:) = color;
                    end
                    
                    % Plotting Data Set
                    scatter(subset(:,1),subset(:,2),25,partycolor,'filled');
                    hold on;
                    
                    % Plot the contour of the Gaussian Model
                    h = ezcontour(@(x,y)pdf(obj,[x y]),[0 8],[0 8]);
                    h.LineWidth = 2.0;
                    hold on;
                    
                    % label axis
                    xlabel('Political Issue 1');
                    ylabel('Political Issue 2');
                    title('2012 - 2 Issues, 2 Clusters');
                    
                    % Plot the Means got by GMM
                    plot(model_means(:,1),model_means(:,2),'kx','LineWidth',4,'MarkerSize',20);
                    
                    % Plot the Party Means
                    plot(mu_democrat(:,1), mu_democrat(:,2), 'b^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','b');
                    plot(mu_republican(:,1), mu_republican(:,2), 'r^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','r');
                else
                    % Reading Data from Excel Sheet
                    num_of_cluster = input('Number of Clusters to Compare: ');
                    filename = 'ANES 2012 Data.xlsx';
                    sheet = 4;
                    xlRange = 'E5:G5918';
                    data = xlsread(filename,sheet,xlRange);
                    
                    % Sort Individuals by which Party they identified with
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
                    
                    % Find mean of each party
                    mu_democrat = mean(democrat(:,2:end));
                    s_democrat = std2(democrat(:,2:end));
                    mu_republican = mean(republican(:,2:end));
                    s_republican = std2(republican(:,2:end));
                    
                    figure;
                    
                    for cc = 1:num_of_cluster
                        % Separate the party identifications and opinions of other
                        % issues
                        subset = data(:,2:end);    % All columns but first are political issues
                        partyscale = data(:,1);    % first column is party affiliation
                        
                        % Finding Correlation Matrix for Data
                        [RHO,PVAL] = corr(subset);
                        
                        % Find the Gaussian Fit Model
                        options = statset('Display','final','MaxIter',1000,'TolFun',1e-9);
                        obj = gmdistribution.fit(subset,cc,'CovType','diagonal','Options',options,'Regularize',1e-5);
                        
                        % Save 2 Parameter Values
                        model_means = obj.mu;
                        cluster_means(:,:,cc) = {obj.mu};
                        model_cov = obj.Sigma;
                        cluster_cov(:,:,cc) = {obj.Sigma};
                        AIC(cc,:) = obj.AIC;
                        BIC(cc,:) = obj.BIC;
                        
                        % cluster and posterior probablity of each instance
                        % note that1: [~,clustIdx] = max(p,[],2)
                        [clustInd,~,p] = cluster(obj, subset);
                        %             tabulate(clustInd)
                        
                        
                        h(cc) = subplot(4,4,cc);
                        
                        % Adding Noise to the repeated rows
                        [~,unique_rows,~] = unique(subset,'rows');
                        duplicate_rows = setdiff(1:size(subset,1),unique_rows);
                        subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), 0.5 , 1);
                        
                        % Indicating Party Affiliation through Color
                        for ii = 1:length(partyscale)
                            color = mapcolor(partyscale(ii));
                            partycolor(ii,:) = color;
                        end
                        
                        % Plotting Data Set
                        scatter(subset(:,1),subset(:,2),25,partycolor,'filled');
                        hold on;
                        
                        
                        % Plot the contour of the Gaussian Model
                        h = ezcontour(@(x,y)pdf(obj,[x y]),[0 10],[0 10]);
                        h.LineWidth = 2.0;
                        hold on;
                        
                        title(sprintf('%d Cluster(s)',cc));
                        xlabel('Component 1');
                        ylabel('Component 2');
                        hold on;
                        
                        % Plot the Means got by GMM
                        plot(model_means(:,1),model_means(:,2),'kx','LineWidth',4,'MarkerSize',20);
                        hold on;
                        
                        % Plot the Party Means
                        plot(mu_democrat(:,1), mu_democrat(:,2), 'b^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','b');
                        plot(mu_republican(:,1), mu_republican(:,2), 'r^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','r');
                        hold on;
                    end
                    %             linkaxes(h,'xy');
                end
                
            case 3
                % Reading Data from Excel Sheet
                num_of_cluster = input('Number of Clusters: ');
                filename = 'ANES 2012 Data.xlsx';
                sheet = 4;
                xlRange = 'G5:I5918';
                subset = xlsread(filename,sheet,xlRange);
                
                % Finding Correlation Matrix for Data
                [RHO,PVAL] = corr(subset);
                
                % Find the Gaussian Fit Model
                options = statset('Display','final');
                obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
                
                % Show Parameter Values
                model_means = obj.mu;
                model_cov = obj.Sigma;
                
                % Adding Noise to the repeated rows
                [~,unique_rows,~] = unique(subset,'rows');
                duplicate_rows = setdiff(1:size(subset,1),unique_rows);
                subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
                
                % Indicating Party Affiliation through Color
                partysheet = 3;
                Range = 'E5:E5918';
                partyscale = xlsread(filename,partysheet,Range);
                
                for ii = 1:length(partyscale)
                    color = mapcolor(partyscale(ii));
                    partycolor(ii,:) = color;
                end
                
                % Plotting Points in 3D Scatter Plot
                scatter3(subset(:,1),subset(:,2),subset(:,3),25,partycolor,'filled')
                hold on;
                plot3(model_means(:,1),model_means(:,2),model_means(:,3), 'gx','LineWidth',4,'MarkerSize',20);
                
                % Plotting Points on PCA Axes
                [coefs,score] = pca(subset);
                
                figure;
                for jj = 1:length(partycolor)
                    plot(score(jj,1),score(jj,2),'o','Color',partycolor(jj,:));
                    hold on;
                end
                title('Principle Component Analysis');
                xlabel('Component 1');
                ylabel('Component 2');
                hold on;
                
                % Plotting cluster means onto PCA axes
                mu = mean(subset); % find mean of data
                
                % find model mean in the centered space
                for kk = 1:size(model_means,1)
                    mean_centered(kk,:) = model_means(kk,:) - mu;
                end
                
                % now mapp centered model_mean to the pca components
                mean_mapped = mean_centered/coefs';
                plot(mean_mapped(:,1), mean_mapped(:,2), 'kx','LineWidth',4,'MarkerSize',20);
                
            case 14
                % Choose whether you want GMM model or Prediction of which party
                % individual will vote for
                program = input('1.GMM, 2.Predicting Voting: ');
                switch(program)
                    case 1 % GMM model for 14 Political Issues
                        % Reading Data from Excel Sheet
                        num_of_cluster = input('Number of Clusters: ');
                        filename = 'ANES 2012 Data.xlsx';
                        sheet = 4;
                        xlRange = 'E5:S5918';
                        data = xlsread(filename,sheet,xlRange);
                        
                        % Sort Individuals by which Party they identified with
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
                        
                        % Find mean of each party
                        mu_democrat = mean(democrat(:,2:end));
                        s_democrat = std2(democrat(:,2:end));
                        mu_republican = mean(republican(:,2:end));
                        s_republican = std2(republican(:,2:end));
                        
                        % Separate the party identifications and opinions of other
                        % issues
                        subset = data(:,2:end);    % All columns but first are political issues
                        partyscale = data(:,1);    % first column is party affiliation
                        
                        % Finding Correlation Matrix between Issues(columns)
                        [RHO,PVAL] = corr(subset);
                        
                        % Find the Gaussian Fit Model
                        options = statset('Display','final','MaxIter',1000,'TolFun',1e-9);
                        obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
                        
                        % Save Parameter Values
                        model_means = obj.mu;
                        model_cov = obj.Sigma;
                        
                        AIC = obj.AIC
                        BIC = obj.BIC
                        
                        % cluster and posterior probablity of each instance
                        % note that1: [~,clustIdx] = max(p,[],2)
                        [clustInd,~,p] = cluster(obj, subset);
                        tabulate(clustInd)
                        
                        % Storing mean of subset for plotting cluster means onto PCA axes
                        mu = mean(subset); % find mean of data
                        
                        % Adding Noise to repeated rows
                        [~,unique_rows,~] = unique(subset,'rows');
                        duplicate_rows = setdiff(1:size(subset,1),unique_rows);
                        subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
                        
                        % Indicating Party Affiliation through Color
                        partysheet = 3;
                        Range = 'E5:E5918';
                        partyscale = xlsread(filename,partysheet,Range);
                        
                        for ii = 1:length(partyscale)
                            color = mapcolor(partyscale(ii));
                            partycolor(ii,:) = color;
                        end
                        
                        % Plotting points onto PCA axes
                        [coefs,score] = pca(subset);
                        
                        for jj = 1:length(partycolor)
                            plot(score(jj,1),score(jj,2),'o','Color',partycolor(jj,:));
                            hold on;
                        end
                        
                        title('Principle Component Analysis');
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
                        plot(mu_democrat_mapped(:,1), mu_democrat_mapped(:,2), 'b^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','b');
                        plot(mu_republican_mapped(:,1), mu_republican_mapped(:,2), 'r^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','r');
                        hold on;
                    case 2 % Predicting Voting through Distance to party means
                        % Reading Data from Excel Sheet
                        num_of_cluster = input('Number of Clusters: ');
                        filename = 'ANES 2012 Data.xlsx';
                        sheet = 4;
                        xlRange = 'E5:S5918';
                        data = xlsread(filename,sheet,xlRange);
                        
                        % Sorting Individuals by Party identification
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
                        
                        % Find mean of each Party
                        mu_democrat = mean(democrat(:,2:end));
                        s_democrat = std2(democrat(:,2:end));
                        mu_republican = mean(republican(:,2:end));
                        s_republican = std2(republican(:,2:end));
                        
                        % Separate Party Identification with other opinions
                        subset = data(:,2:end);    % All columns but first are political issues
                        partyscale = data(:,1);    % first column is party affiliation
                        
                        % Plotting cluster means onto PCA axes
                        mu = mean(subset); % find mean of data
                        
                        % Conducting PCA Analysis on the data
                        [coefs,score] = pca(subset);
                        
                        % Adding Noise to the repeated rows
                        [~,unique_rows,~] = unique(subset,'rows');
                        duplicate_rows = setdiff(1:size(subset,1),unique_rows);
                        subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
                        
                        % Plotting Party Affiliation through Color and Party Vote
                        % through squares(vote for democrat) or diamonds(vote for
                        % republican)
                        figure;
                        o = 1;
                        p = 1;
                        for jj = 1:length(partyscale)
                            color = mapcolor(partyscale(jj));  % returns color of party affiliation
                            dist_democrat = norm(subset(jj,:)- mu_democrat);  % Calculates distance from democratic mean
                            dist_republican = norm(subset(jj,:) - mu_republican); % Calculates distance from republican mean
                            if dist_democrat < dist_republican
                                plot(score(jj,1),score(jj,2),'s','Color',color);  % Plot a square
                                vote_democrat(o,:) = [partyscale(jj) subset(jj,:)];  % Store the individual in vote_democrat
                                o = o+1;
                            else
                                plot(score(jj,1),score(jj,2),'d','Color',color); % Plot a diamond
                                vote_republican(p,:) = [partyscale(jj) subset(jj,:)]; % Store the individual in vote_republican
                                p = p+1;
                            end
                            hold on;
                        end
                        title('2012 Voting Predictions');
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
        end
    case 2
        dimension = input('Number of Variables: ');
        
        switch(dimension)
            case 2
                %                 global subset covmin num_of_cluster;
                covmin = 1;
                % Reading Data from Excel Sheet
                num_of_cluster = input('Number of Clusters: ');
                filename = 'ANES 2012 Data.xlsx';
                sheet = 4;
                xlRange = 'E5:G5918';
                data = xlsread(filename,sheet,xlRange);
                
                % Sort Individuals by which Party they identified with
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
                
                % Find mean of each party
                mu_democrat = mean(democrat(:,2:end));
                s_democrat = std2(democrat(:,2:end));
                mu_republican = mean(republican(:,2:end));
                s_republican = std2(republican(:,2:end));
                
                % Separate the party identifications and opinions of other
                % issues
                
                subset = data(:,2:end);    % All columns but first are political issues
                partyscale = data(:,1);    % first column is party affiliation
                
                %         for num_of_cluster = 1: 2
                rng('default');
                
                model_means = randi(7,num_of_cluster,size(subset,2));
                cov = ones(num_of_cluster,size(subset,2));
                options = optimset('MaxIter', 3000,'MaxFunEvals',3000,...
                    'TolFun',1, 'TolX', 0.1);
                [x,fval,exitflag] = fminsearch(@NegativeLogLikelihood,...
                    [model_means; cov],options);
                LogLikelihood = -fval;
                dim = num_of_cluster;
                mu = x(1:dim,:);
                sigma = x(dim+1:end,:);
                
                
                numParam = num_of_cluster + (size(subset,2)*num_of_cluster);
                aic = aicbic(LogLikelihood,numParam)
                %         end
                %         plot(1:10,aic,'o-','LineWidth',3)
                % title('AIC Cluster Analysis');
                % xlabel('Number of Clusters');
                % ylabel('AIC');
                
                % Adding Noise to the repeated rows
                [~,unique_rows,~] = unique(subset,'rows');
                duplicate_rows = setdiff(1:size(subset,1),unique_rows);
                subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
                
                % Indicating Party Affiliation through Color
                for ii = 1:length(partyscale)
                    color = mapcolor(partyscale(ii));
                    partycolor(ii,:) = color;
                end
                
                % Plotting Data Set
                scatter(subset(:,1),subset(:,2),25,partycolor,'filled');
                hold on;
                
                % Plot the Means got by fminsearch
                plot(mu(:,1),mu(:,2),'gx','LineWidth',4,'MarkerSize',20);
                
                % Plot the Party Means
                plot(mu_democrat(:,1), mu_democrat(:,2), 'b^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','b');
                plot(mu_republican(:,1), mu_republican(:,2), 'r^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','r');
                
                title('2012 - 2 Issues');
                xlabel('Issue 1');
                ylabel('Issue 2');
                
            case 14
                global subset covmin num_of_cluster;
                covmin = 1;
                
                % Reading Data from Excel Sheet
                num_of_cluster = input('Number of Clusters: ');
                filename = 'ANES 2012 Data.xlsx';
                sheet = 4;
                xlRange = 'E5:S5918';
                data = xlsread(filename,sheet,xlRange);
                
                % Sort Individuals by which Party they identified with
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
                
                % Find mean of each party
                mu_democrat = mean(democrat(:,2:end));
                s_democrat = std2(democrat(:,2:end));
                mu_republican = mean(republican(:,2:end));
                s_republican = std2(republican(:,2:end));
                
                % Separate the party identifications and opinions of other
                % issues
                
                subset = data(:,2:end);    % All columns but first are political issues
                partyscale = data(:,1);    % first column is party affiliation
                
                % Plotting cluster means onto PCA axes
                mu = mean(subset); % find mean of data
                
                %         for num_of_cluster = 1: 2
                rng('default');
                
                initial_means = randi(7,num_of_cluster,size(subset,2));
                cov = ones(num_of_cluster,size(subset,2));
                options = optimset('MaxIter', 3000,'MaxFunEvals',3000,...
                    'TolFun',1, 'TolX', 0.1);
                [x,fval,exitflag] = fminsearch(@NegativeLogLikelihood,...
                    [initial_means; cov],options);
                LogLikelihood = -fval;
                dim = num_of_cluster;
                
                model_means = x(1:dim,:);
                sigma = x(dim+1:end,:);
                
                
                numParam = num_of_cluster + (size(subset,2)*num_of_cluster);
                aic = aicbic(LogLikelihood,numParam)
                %         end
                %         plot(1:10,aic,'o-','LineWidth',3)
                % title('AIC Cluster Analysis');
                % xlabel('Number of Clusters');
                % ylabel('AIC');
                
                % Adding Noise to the repeated rows
                [~,unique_rows,~] = unique(subset,'rows');
                duplicate_rows = setdiff(1:size(subset,1),unique_rows);
                subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
                
                % Indicating Party Affiliation through Color
                
                
                for ii = 1:length(partyscale)
                    color = mapcolor(partyscale(ii));
                    partycolor(ii,:) = color;
                end
                
                % Plotting points onto PCA axes
                [coefs,score] = pca(subset);
                
                for jj = 1:length(partycolor)
                    plot(score(jj,1),score(jj,2),'o','Color',partycolor(jj,:));
                    hold on;
                end
                
                title('2012-14 Issues');
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
                plot(mu_democrat_mapped(:,1), mu_democrat_mapped(:,2), 'b^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','b');
                plot(mu_republican_mapped(:,1), mu_republican_mapped(:,2), 'r^', 'LineWidth', 3, 'MarkerSize', 13,'MarkerEdgeColor','k','MarkerFaceColor','r');
                hold on;
                
        end
end