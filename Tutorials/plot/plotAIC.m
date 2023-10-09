function plotAIC(aic,logL,k,lambda)

[minAIC,IminAIC] = min(aic);

figure
plot(-log(lambda),aic-minAIC,'LineWidth',1.4); hold on
plot(-log(lambda),logL-minAIC); hold on
plot(-log(lambda),2*k)
plot(-log(lambda(IminAIC)),aic(IminAIC)-minAIC,'r*')
ylabel(sprintf('AIC'),'Interpreter','latex')
xlabel(sprintf('$-\\mathrm{log}(\\lambda)$'),'Interpreter','latex')
legend({'AIC','Model prediction error: logL','Model complexity: 2k','Selected model'},'Interpreter','latex','Location','northwest')
