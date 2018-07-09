clear
close all
%% generate the data
x = 0:3:100;
y = 0:3:100;
X = x .* ones(length(x), length(x));
Y = y' .* ones(length(y), length(y));
% init Z
Z = zeros(26, 26);
% init the parameters
sigma = 1000;
Ex = 50;
Ey = 50;
% calculate Z
for row = 1:1:length(X)
    for col = 1:1:length(Y)
        Z(row, col) = (x(row) - Ex).*(x(row) - Ex) + (y(col) - Ey).*(y(col) - Ey);
    end
end
Z = -Z/(2 * sigma);
Z = 100 * exp(Z) + 50 * ones(length(X), length(Y));
x0 = 10:4:90;
y0 = 10:4:90;
X0 = x0 .* ones(length(x0), length(x0));
Y0 = y0' .* ones(length(y0), length(y0));
Z0 = zeros(length(x0), length(y0));
[row, col] = size(Z);
% generate all the data pionts
data = zeros(row*col, 3);
k = 1;
for i = 1:row
    for j = 1:col
        data(k, :) = [x(i) y(j) Z(i, j)];
        k = k +1;
    end
end
%% init the nodes
[row1,col1] = size(X0);
Nd = zeros(row1*col1, 3);
k = 1;
for i = 1:row1
    for j = 1:col1
        Nd(k, :) = [x0(i) y0(j) Z0(i, j)];
        k = k +1;
    end
end
%% Randomise data order
msize = length(data);
xpermidx = randperm(msize);
data = data(xpermidx,:);
%% First plot
% Plot the data
figure;
dp = subplot(3,2,1:4);
h1 = mesh(X,Y,Z);
hold on;
set(h1,'facealpha',0.5)
plot3(X, Y, Z,'*b')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal;
% Plot nodes
mesh(X0, Y0, Z0);
plot3(X0, Y0, Z0, 'or');
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
sig(1) = 100*exp(-(it/v));
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
while ksi(it) > 0.001 % iterations loop
    %% Find distance of all data points to the nodes
    dists = pdist2(data,Nd);
    %% loop through data points
    for ii = 1:msize
        dp = subplot(3,2,1:4);
        plot3(data(ii,1),data(ii,2),data(ii,3),'dr','markersize',6);
        % Find nearest Node in this point
        [~,p] = min(dists(ii,:));
        plot3(Nd(p,1),Nd(p,2),Nd(p,3),'or','markersize',9);
        %% Move nodes towards this data point
        ve =  data(ii,:)-Nd(p,:);%vector between them
        plot3 ([Nd(p,1),data(ii,1)],[Nd(p,2),data(ii,2)],[Nd(p,3),data(ii,3)],'-g','linewidth',2)% plot vector between them
        % Move the nodes
        vect = zeros(length(Nd), 3);
        for jj = 1:length(Nd)
            vect(jj,:) = ve*ksi(it)*exp(-((Nd(p,1)-Nd(jj,1))^2 + (Nd(p,2)-Nd(jj,2))^2)/sig(it)^2); %vector after ksi and hfunct aplied
        end 
        quiver3(Nd(p,1),Nd(p,2),Nd(p,3),vect(p,1),vect(p,2),vect(p,3),0,'-r','MaxHeadSize',9,'linewidth',2) ;hold off;
        Nd = Nd + vect; %new nodes position
        % update distances for all nodes
        cd = pdist2(data,Nd);
        dists = cd;
        %% Data Iteration Plot update
        % Plot Data
        pause(pt)
        cla(dp)
        dp = subplot(3,2,1:4);
        h2 = mesh(X, Y, Z);
        hold on;
        set(h2,'facealpha',0.5)
        plot3(data(:,1),data(:,2),data(:,3),'b.','markersize',4); hold on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        axis equal;
        % Plot the nodes
        kk = 1;
        X00 = zeros(row1, col1);
        Y00 = zeros(row1, col1);
        Z00 = zeros(row1, col1);
        for i = 1:row1
            for j = 1:col1
                X00(i,j) = Nd(kk,1);
                Y00(i,j) = Nd(kk,2);
                Z00(i,j) = Nd(kk,3);
                kk = kk +1;
            end
        end
        mesh(X00,Y00,Z00);
        plot3(X00,Y00,Z00,'.r','markersize',9);
        axis([-20 130 -20 130 -20 150]);
        title('Data Space')
    end
    %% Plot ksi and quantasation error
    % Permutate data
    permidx = randperm(msize);
    data = data(permidx,:);
    pt = pt*0.1; % increase plot speed for each iteration
    it = it + 1; % increase iteration
     %% Update Learning parameter
    ksi(it) = exp(-(it/tau));
    sig(it) = 100*exp(-it/v);
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
%% plot the final result
% plot the original data points
figure(2)
h3 = mesh(X, Y, Z);
hold on;
set(h3,'facealpha',0.5)
plot3(data(:,1),data(:,2),data(:,3),'b.','markersize',4); hold on;
axis([-20 130 -20 130 -20 150]);
% plot the final SOM result
figure(3)
h4 = mesh(X00,Y00,Z00);
hold on;
set(h4,'facealpha',0.5);
plot3(X00,Y00,Z00,'.r','markersize',9);
axis([-20 130 -20 130 -20 150]);