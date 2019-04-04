
% Generate electroce signals with different sampling rates and use FNN to
% test the differences of the mininal embedding dimensions of the two
% signals.

clc
clearvars
hFig = figure(2);
set(hFig,'units','centimeters','position',[0,0,16,10])
max_dimension = 8;
tau=2;
rtol = 15;
atol = 2;
% generate signals with sampling rate of 10kHz
[signals target r1] = generatenoisysamples('Duration', 0.1, 'SampleRate', 10000) ;
subplot(2,2,1)
plot(signals,'k-')
axis([0 1000 -1 1])
set(gca,'fontsize',10,'ticklength',[0.03 0.03],'Xcolor','k')
xlabel('Sampling points')
ylabel('Electrode signal')
title('Sampling rate: 10kHz' )
% FNN embeddig
fnn = f_fnn(signals,tau,max_dimension,rtol,atol);
subplot(2,2,2)
plot(fnn,'k-')
axis([0 8 0 100])
set(gca,'fontsize',10,'ticklength',[0.03 0.03],'Xcolor','k')
xlabel('Embedding dimension')
ylabel('FNNP')
title('Sampling rate: 10kHz' )


% generate signals with sampling rate of 2kHz
[signals target1 r1] = generatenoisysamples('Duration', 1, 'SampleRate', 2000);
subplot(2,2,3)
plot(signals,'k-')
axis([0 200 -1 1])
set(gca,'fontsize',10,'ticklength',[0.03 0.03],'Xcolor','k')
xlabel('Sampling points')
ylabel('Electrode signal')
title('Sampling rate: 2kHz' )
% FNN embeddig
fnn = f_fnn(signals,tau,max_dimension,rtol,atol);
subplot(2,2,4)
plot(fnn,'k-')
axis([0 8 0 100])
set(gca,'fontsize',10,'ticklength',[0.03 0.03],'Xcolor','k')
xlabel('Embedding dimension')
ylabel('FNNP')
title('Sampling rate: 2kHz' )
