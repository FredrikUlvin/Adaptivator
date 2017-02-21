%% Init
close all;

normMat = @(k) ( k - ones(size(k))*min(k(:)) ) ./ (max(k(:))-min(k(:)));
wa = load('wada.mat');
woa = load('woada.mat');
n = 0;

SSTD  = std(mcS);
SMEAN = mean(mcS);
wincovarSTD = std(mcwincovar);
wincovarMEAN = mean(mcwincovar);
innocovarREAL = mcinnocov(1,:);

SwinERR = mean(abs(mcS-mcwincovar));
StruERR = mean(abs(mcS-mcinnocov));

innoSTD = std(abs(mcinno));
innoMEAN = mean(abs(mcinno));

alphaSTD = std(mcalpha);
alphaMEAN = mean(mcalpha);



%%
n = n + 1;
figure(n)
hold on;
errorbar(len(wincovarSTD), wincovarMEAN, wincovarSTD,'g','LineWidth',.2)
errorbar(len(SSTD)+.2*ones(1, length(SSTD)),SMEAN,SSTD, 'b', 'LineWidth', .2)
plot(len(wincovarSTD),wincovarMEAN, 'g', 'LineWidth',2)
plot(len(SSTD)+.2*ones(1, length(SSTD)),SMEAN,'b', 'LineWidth',2)
plot(len(innocovarREAL), innocovarREAL, 'r', 'LineWidth',2)
set(gca,'YScale','log')
legend('Win. std.', 'Est S std.','Win. CoVar. mean','Est S std. mean','Best ','Location','se')
title('Tracking Innovation Covariance');
xlabel('Steps')
grid on;
axis([0 400 0 30]);
hold off;

%%
n = n+1;
figure(n)
hold on;

subplot(2,1,1)
plot(len(SSTD), SSTD, len(wincovarSTD), wincovarSTD)
grid on;
legend('Estimated','Sample', 'Location', 'se')
title('S Tracking std.')

subplot(2,1,2)
plot(len(SwinERR), SwinERR, len(StruERR), StruERR)
legend('estimated','sample')
axis([0 400 0 16]);
title('mean S tracking error')
grid on;


%% INNO [TODO: ADD NORMAL KF]
n = n + 1;
figure(n)

hold on;
%errorbar(len(innoSTD), innoMEAN, innoSTD,'LineWidth',.1)
plot(len(innoMEAN), woa.innoMEAN, len(wa.innoMEAN), wa.innoMEAN)
%set(gca,'YScale','log')
grid on;
title('Innovation variance');
legend('KF with est.', 'KF wo est.', 'Location', 'se')
xlabel('Steps');
ylabel('Std.');
%axis([0 400 0 3]);
hold off;

%% ALHPA
n = n + 1;
figure(n)

hold on;
errorbar(len(alphaSTD), log2(alphaMEAN), log2(alphaSTD),'LineWidth',.1)
plot(len(alphaMEAN), log2(alphaMEAN))
plot(len(alphaMEAN), [ones(1,length(alphaMEAN)/2)*log2(Rpre) ones(1,length(alphaMEAN)/2)*log2(Rpost)], 'r', 'LineWidth',2)
%set(gca,'YScale','log')
xlabel('Steps');
ylabel('log2(\alpha^2)');
legend( '\alpha^2 std.', '\alpha^2 mean', 'Theoretical \alpha^2','Location','SouthEast')
grid on;
axis([0 400 -2 7]);
title('Scaling variable \alpha^2');
%axis([0 400 0 30]);
hold off;
