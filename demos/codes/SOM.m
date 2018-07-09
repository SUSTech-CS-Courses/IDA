%% illustration of the code
% Noting that SOM is not used for data clustering but data
% dimensionality reduction that similar with PCA does
%%
clear all
close all
%% Create data
% data 1
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [1, 2];% mean
Sigma =diag(sig);
R = chol(Sigma);
x1 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 2
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [3, 6];% mean
Sigma =diag(sig);
R = chol(Sigma);
x2 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 3
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [5, 11];% mean
Sigma =diag(sig);
R = chol(Sigma);
x3 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 4
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [9, 9];% mean
Sigma =diag(sig);
R = chol(Sigma);
x4 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 5
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [13, 6];% mean
Sigma =diag(sig);
R = chol(Sigma);
x5 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 6
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [16, 9];% mean
Sigma =diag(sig);
R = chol(Sigma);
x6 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 7
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [18, 13];% mean
Sigma =diag(sig);
R = chol(Sigma);
x7 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 8
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [20, 18];% mean
Sigma =diag(sig);
R = chol(Sigma);
x8 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 9
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [22, 23];% mean
Sigma =diag(sig);
R = chol(Sigma);
x9 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 10
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [25, 27];% mean
Sigma =diag(sig);
R = chol(Sigma);
x10 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 11
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [28, 22];% mean
Sigma =diag(sig);
R = chol(Sigma);
x11 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 12
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [31, 17];% mean
Sigma =diag(sig);
R = chol(Sigma);
x12 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 13
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [35, 21];% mean
Sigma =diag(sig);
R = chol(Sigma);
x13 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 14
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [40, 22];% mean
Sigma =diag(sig);
R = chol(Sigma);
x14 = repmat(mu,N,1) + randn(N,2)*R; % create points
% data 15
N = 50; % Number of generated points
sig1 = 2; % 1st dimension variation
sig2 = 2; % 2st dimension variation
sig = [sig1, sig2];
mu = [45, 20];% mean
Sigma =diag(sig);
R = chol(Sigma);
x15 = repmat(mu,N,1) + randn(N,2)*R; % create points
% whole dataset
x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12;x13;x14;x15];
%% Initiate nodes (Nd)
%set number of nodes
M = 20;
% randomly initiate the position of every node
% xc = 40 * rand(1,30);
% yc = 20 * rand(1,30);
% xc = sort(xc);
% yc = sort(yc);
% Nd(:,1) = xc;
% Nd(:,2) = yc;
Nd = zeros(M, 2);
for m = 1:M
    Nd(m,:) = [m*2.25 1];
end
%% Randomise data order
msize = length(x);
permidx = randperm(msize);
x = x(permidx,:);
%% First plot
% Plot the data
figure;
dp = subplot(3,2,1:4);
plot(x(:,1),x(:,2),'b.','markersize',4); hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal;
% Plot nodes
plot(Nd(:,1),Nd(:,2),'-*r','markersize',9);
title('Data Space')
xp = subplot(3,2,5);
title('Learning Parameter')
xlabel ('Iteration Number')

qp = subplot(3,2,6);
title('the value of sigma')
xlabel ('Iteration Number')
%% Start iterations
it = 1;
pt = 0.1; % pause time
%% Initialise learning parameter
tau = 1;
v = 2;
ksi=ones(6,1)*0;
sig = ones(6,1)*0;
ksi(1) = exp(-(it/tau));
sig(1) = 8*exp(-(it/v));
cla(xp)
subplot(3,2,5);
bar(ksi,'r');
hold on
plot([0 7.5],[0.001 0.001],'k--');
title('Learning Parameter')
xlabel ('Iteration Number')
set(gca,'YScale','log');
set(gca,'YLim',[0.0001 1]); hold off
cla(qp)
subplot(3,2,6);
bar(sig,'g');
hold on
title('value of sigma')
xlabel ('Iteration Number')
hold off
%% making video
vd = VideoWriter('SOM_video.avi');
open(vd);
while ksi(it) > 0.001 % iterations loop
    %% initialise data tag
    dataCflag = zeros(msize,1);
    %% Find distance of all data points to the nodes
    dists = pdist2(x,Nd);
    %% loop through data points
    for ii = 1:msize
        dp = subplot(3,2,1:4);
        plot(x(ii,1),x(ii,2),'dr','markersize',6);
        % Find nearest CC in this point
        [~,p] = min(dists(ii,:));
        plot(Nd(p,1),Nd(p,2),'or','markersize',9);
        %% Move nodes towards this data point
        ve =  x(ii,:)-Nd(p,:);%vector between them
        plot ([Nd(p,1),x(ii,1)],[Nd(p,2),x(ii,2)],'-g','linewidth',2)% plot vector between them
        for jj = 1:length(Nd)
            vect(jj,:) = ve*ksi(it)*exp(-(p-jj)^2/sig(it)^2); %vector after ksi and hfunct aplied
        end 
        quiver(Nd(p,1),Nd(p,2),vect(p,1),vect(p,2),0,'-r','MaxHeadSize',9,'linewidth',2) ;hold off;
        Nd = Nd+vect; %new nodes position
        % update distances for this node
        cd = pdist2(x,Nd);
        dists = cd;
        %% Data Iteration Plot update
        % Plot Data
        pause(pt)
        cla(dp)
        dp = subplot(3,2,1:4);
        plot(x(:,1),x(:,2),'b.','markersize',4); hold on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        axis equal;
        % Plot the nodes
        plot(Nd(:,1),Nd(:,2),'-*r','markersize',9);
        title('Data Space')
        
        % plot(Nd(p,1),Nd(p,2),'or','markersize',8);
        frame = getframe(gcf);
        writeVideo(vd, frame);
    end
    %% Plot ksi and quantasation error
    % Permutate data
    permidx = randperm(msize);
    x = x(permidx,:);
    pt = pt*0.1; % increase plot speed for each iteration
    it = it + 1; % increase iteration
     %% Update Learning parameter
    ksi(it) = exp(-(it/tau));
    sig(it) = 8*exp(-it/v);
    cla(qp)
    qp = subplot(3,2,6);
    bar(sig, 'g');
    title(['value of sigma'])
    xlabel ('Iteration Number')
    pause(pt*20)
    cla(xp)
    subplot(3,2,5);
    bar(ksi,'r');hold on
    plot([0 7.5],[0.001 0.001],'k--');
    title('Learning Parameter')
    set(gca,'YScale','log');
    xlabel ('Iteration Number')
    set(gca,'YLim',[0.0001 1]); hold off
end
close(vd);