clear all
close all
%% Colors string
cl = ['r','g','b','k','y','m','c'];
%% Create Clusters
% Cluster 1
N = 100; % Number of generated points
sig1 = 5; % 1st dimension variation
sig2 = 5; % 2st dimension variation
sig = [sig1, sig2];
mu = [-9, 8];% mean
Sigma =diag(sig);
R = chol(Sigma);
x1 = repmat(mu,N,1) + randn(N,2)*R; % create points
% Cluster 2
N = 100; % Number of generated points
sig1 = 5; % 1st dimension variation
sig2 = 5; % 2st dimension variation
sig = [sig1, sig2];
mu = [6, 9];% mean
Sigma =diag(sig);
R = chol(Sigma);
x2 = repmat(mu,N,1) + randn(N,2)*R; % create points
% Cluster 3
N = 100; % Number of generated points
sig1 = 5; % 1st dimension variation
sig2 = 5; % 2st dimension variation
sig = [sig1, sig2];
mu = [3, -2.5];% mean
Sigma =diag(sig);
R = chol(Sigma);
x3 = repmat(mu,N,1) + randn(N,2)*R; % create points
%% Plot the data
% figure;
% plot(x1(:,1),x1(:,2),'r.'); hold on;
% plot(x2(:,1),x2(:,2),'b.');
% plot(x3(:,1),x3(:,2),'k.');
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% axis equal;
%% Group the data
x = [x1;x2;x3];
%% Initiate cluster centers (CC)
%set number of clusters
M = 3;
for m = 1:M
    CC(m,:) = x(end-m+1,:);
end
%% Randomise data order
msize = length(x);
permidx = randperm(msize);
x = x(permidx,:);
%% First plot
% Plot the data
figure;
dp = subplot(3,2,1:4);
plot(x(:,1),x(:,2),'b.','markersize',9); hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal;
% Plot cluster centers
plot(CC(:,1),CC(:,2),'*r','markersize',9);
title('Data Space')
xp = subplot(3,2,5);
title('Learning Parameter')
xlabel ('Iteration Number')

qp = subplot(3,2,6);
title('Quantisation Error')
xlabel ('Iteration Number')

%% Start iterations
it = 1;
pt = 0.1; % pause time
%% Initialise learning parameter
tau = 1;
ksi=ones(6,1)*0;
ksi(1) = exp(-(it/tau));
Qerr=zeros(6,1);
cla(xp)
subplot(3,2,5);
bar(ksi,'r');hold on
plot([0 7.5],[0.001 0.001],'k--');
title('Learning Parameter')
xlabel ('Iteration Number')
set(gca,'YScale','log');
set(gca,'YLim',[0.0001 1]); hold off
while ksi(it) > 0.001 % iterations loop
    %% initialise data tag
    dataCflag = zeros(msize,1);
    %% Find distance of all data points to cluster centers
    dists = pdist2(x,CC);
    %% loop through data points
    for ii = 1:msize
        dp = subplot(3,2,1:4);
        plot(x(ii,1),x(ii,2),'dr','markersize',6);
        % Find nearest CC in this point
        [~,p] = min(dists(ii,:));
        plot(CC(p,1),CC(p,2),'or','markersize',9);
        %% Move CC towards this data point
        ve =  x(ii,:)-CC(p,:);%vector between them
        plot ([CC(p,1),x(ii,1)],[CC(p,2),x(ii,2)],'-g','linewidth',2)% plot vector between them
        ve = ve*ksi(it); %vector after ksi aplied
        quiver(CC(p,1),CC(p,2),ve(1),ve(2),0,'-r','MaxHeadSize',9,'linewidth',2) ;hold off;
        CC(p,:) = CC(p,:)+ve; %new CC position
        % update distances for this CC
        cd = pdist2(x,CC(p,:));
        dists(:,p) = cd;
        % update data cluster flag
        dataCflag(ii) = p;
        %% Data Iteration Plot update
        % Plot Data
        pause(pt)
        cla(dp)
        dp = subplot(3,2,1:4);
        plot(x(:,1),x(:,2),'b.','markersize',9); hold on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        axis equal;
        % Plot cluster centers
        plot(CC(:,1),CC(:,2),'*r','markersize',9);
        title('Data Space')
        
        % plot(CC(p,1),CC(p,2),'or','markersize',8);
    end
    %% Calculate quantasation error for this iteration
    for cn = 1:M
        Qerr(it) = Qerr(it) + sum(pdist2(x(dataCflag==cn,:),CC(cn,:)));
    end
    Qerr(it) = Qerr(it)/msize;
    %% Plot clusters for this iteration
    cla(dp)
    dp = subplot(3,2,1:4);
    for cn = 1:M
        plot(x(dataCflag==cn,1),x(dataCflag==cn,2),[cl(cn),'.']); hold on;
        plot(CC(cn,1),CC(cn,2),[cl(cn),'o'],'markersize',8);
        plot(CC(cn,1),CC(cn,2),[cl(cn),'x'],'markersize',8);
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        axis equal;
    end
    hold off;
    title('Data Space')
    %% Plot ksi and quantasation error
    
    
    cla(qp)
    qp = subplot(3,2,6);
    bar(Qerr);
    title(['Quantisation Error'])
    xlabel ('Iteration Number')
    pause(pt*20)
    
    % Permutate data
    permidx = randperm(msize);
    x = x(permidx,:);
    pt = pt*0.1; % increase plot speed for each iteration
    it = it + 1; % increase iteration
    %% Update Learning parameter
    ksi(it) = exp(-(it/tau));
    cla(xp)
    subplot(3,2,5);
    bar(ksi,'r');hold on
    plot([0 7.5],[0.001 0.001],'k--');
    title('Learning Parameter')
    set(gca,'YScale','log');
    xlabel ('Iteration Number')
    set(gca,'YLim',[0.0001 1]); hold off
end













