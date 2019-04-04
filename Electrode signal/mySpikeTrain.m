
% clc
% clearvars
% figure(1)
% max_dimension = 10;
% tau=2;
% rtol = 15;
% atol = 2;
% for k = 1:3
% % [signals target r1] = generatenoisysamples('Duration', 0.2, 'SampleRate', 100000, 'N_Targets',k) ;
% % subplot(2,3,k)
% % plot(signals)
% % save(['target',num2str(k),'_signal.mat'],'signals')
% load (['target',num2str(k),'_signal.mat'],'signals')
% % fnn = f_fnn(signals,tau,max_dimension,rtol,atol);
% load(['target',num2str(k),'_fnn.mat'],'fnn')
% subplot(2,3,k+3)
% plot(fnn)
% % save(['target',num2str(k),'_fnn.mat'],'fnn')
% end

%%


clc
clearvars
hFig = figure(2);
set(hFig,'units','centimeters','position',[0,0,16,10])
max_dimension = 8;
tau=2;
rtol = 15;
atol = 2;
[signals target r1] = generatenoisysamples('Duration', 0.1, 'SampleRate', 10000) ;
% save(['Duration01_signal.mat'],'signals')
load(['Duration01_signal.mat'],'signals')
subplot(2,2,1)
plot(signals,'k-')
axis([0 1000 -1 1])
set(gca,'fontsize',10,'ticklength',[0.03 0.03],'Xcolor','k')
xlabel('Sampling points')
ylabel('Electrode signal')
title('Sampling rate: 10kHz' )
% fnn = f_fnn(signals,tau,max_dimension,rtol,atol);
% save(['Duration01_fnn.mat'],'fnn')
load(['Duration01_fnn.mat'],'fnn')
subplot(2,2,2)
plot(fnn,'k-')
axis([0 8 0 100])
set(gca,'fontsize',10,'ticklength',[0.03 0.03],'Xcolor','k')
xlabel('Embedding dimension')
ylabel('FNNP')
title('Sampling rate: 10kHz' )
% [signals target1 r1] = generatenoisysamples('Duration', 1, 'SampleRate', 2000);
% save(['Duration10_signal.mat'],'signals')
load(['Duration10_signal.mat'],'signals')
subplot(2,2,3)
plot(signals,'k-')
axis([0 200 -1 1])
set(gca,'fontsize',10,'ticklength',[0.03 0.03],'Xcolor','k')
xlabel('Sampling points')
ylabel('Electrode signal')
title('Sampling rate: 2kHz' )
fnn = f_fnn(signals,tau,max_dimension,rtol,atol);
% save(['Duration10_fnn.mat'],'fnn')
load(['Duration10_fnn.mat'],'fnn')
subplot(2,2,4)
plot(fnn,'k-')
axis([0 8 0 100])
set(gca,'fontsize',10,'ticklength',[0.03 0.03],'Xcolor','k')
xlabel('Embedding dimension')
ylabel('FNNP')
title('Sampling rate: 2kHz' )
