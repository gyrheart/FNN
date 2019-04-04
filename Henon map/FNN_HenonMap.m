


%% Henon map
max_dimension = 10;
% set parameters
tau = 2; % time delay
rtol = 15;
atol = 2;


num_point = 1e3; % time points
x = zeros(num_point,1);
y = zeros(num_point,1);
x(1) = 0;
y(1) = 0;
a = 1.4;
b = 0.3;
delta =0e-3;

% generate Henon map
for k = 2:num_point
    x(k) = 1-a*x(k-1)^2+y(k-1)+delta*randn(1);
    y(k) = b*x(k-1)+delta*randn(1);
end
x1 = x(1:end-2*tau);
x2 = x(tau+1:end - tau);
x3 = x(2*tau+1:end);

% select three representative points in the trajectory
p1 = 16;
p2 = 50;
p3 = 12;

markersize = 5;
hFig3 = figure(3);
set(hFig3,'units','centimeters','position',[0,0,18,10])
subplot(2,3,1)
plot(x,y,'.','color',[0.5 0.5 0.5],'markersize',1);
hold on
plot(x(p1),y(p1),'go','markerfacecolor','g','markersize',markersize)
hold on
plot(x(p2),y(p2),'bo','markerfacecolor','b','markersize',markersize)
hold on
plot(x(p3),y(p3),'ro','markerfacecolor','r','markersize',markersize)
hold off
xlabel('x(t)')
ylabel('y(t)')

% temporal activities
subplot(2,3,2)
plot(x(1:100),'-','color',[0.5 0.5 0.5]);
xlabel('t')
ylabel('x(t)')

% 1D embedding
subplot(2,3,3)
plot(x,ones(length(x),1),'.','color',[0.5 0.5 0.5],'markersize',1);
hold on
plot(x(p1),ones(length(x),1),'go','markerfacecolor','g','markersize',markersize)
hold on
plot(x(p2),ones(length(x),1),'bo','markerfacecolor','b','markersize',markersize)
hold on
plot(x(p3),ones(length(x),1),'ro','markerfacecolor','r','markersize',markersize)
hold off
axis([-1.5 1.5 0.8 1.2])
xlabel('x(t)')
title('1-D embedding')



% 2D embedding
subplot(2,3,4)
plot(x1,x2,'.','color',[0.5 0.5 0.5],'markersize',1);
hold on
plot(x1(p1),x2(p1),'go','markerfacecolor','g','markersize',markersize)
hold on
plot(x1(p2),x2(p2),'bo','markerfacecolor','b','markersize',markersize)
hold on
plot(x1(p3),x2(p3),'ro','markerfacecolor','r','markersize',markersize)
hold off
xlabel('x(t)')
ylabel('x(t+2)')
title('2-D embedding')

% 3D embedding
subplot(2,3,5)
plot3(x1,x2,x3,'.','color',[0.5 0.5 0.5],'markersize',1);
hold on
plot3(x1(p1),x2(p1),x3(p1),'go','markerfacecolor','g','markersize',markersize)
hold on
plot3(x1(p2),x2(p2),x3(p2),'bo','markerfacecolor','b','markersize',markersize)
hold on
plot3(x1(p3),x2(p3),x3(p3),'ro','markerfacecolor','r','markersize',markersize)
hold off
xlabel('x(t)')
ylabel('x(t+2)')
zlabel('x(t+4)')
title('3-D embedding')
view(60,-40)


%% FNN embedding
fnn = f_fnn(x,tau,max_dimension,rtol,atol);
%
figure(3);
subplot(2,3,6)
plot(fnn,'k-')
set(gca,'xtick',0:2:10,'fontsize',10)

xlabel('Embedding dimension')
ylabel('False NN percentage')



