function plotESINDy(sig_nn,nWrongTerms,nWrongTermsE,modelError,modelErrorE)

figure
subplot(2,1,1)
plot(sig_nn,nWrongTerms); hold on
plot(sig_nn,nWrongTermsE);
legend({'SINDy','E-SINDy'},'Location','northwest')
xlabel('noise level')
ylabel('number of wrong terms')
subplot(2,1,2)
plot(sig_nn,modelError); hold on
plot(sig_nn,modelErrorE);
xlabel('noise level')
ylabel('model coefficient error')