function plotCV(Xi,epsilon,lambda,minCV,n)

% logl = -log(lambda);

titn = ['x','y','z'];

ylims = [-12 12; -10 30; -4 2];

figure

for i = 1:n
    eps = squeeze(epsilon(i,:,:));
    ks = squeeze(Xi(:,i,:,:));
    m_eps = mean(eps,2);
    % std_eps = std(eps');
    m_ksi = mean(ks,3);

    % figure

    subplot(2,3,i)
    plot(-log(lambda),m_eps); hold on
    % errorbar(-log(lambda),m_eps,std_eps); hold on
    plot(-log(lambda(minCV)),m_eps(minCV),'r*')
    % xlabel('-log(lambda)')
    if i == 1
        ylabel(sprintf('$||\\mathbf{\\dot\\mathbf{x}-\\Theta\\xi}||_2^2$'),'Interpreter','latex')
    end
    title(titn(i),'Interpreter','latex')
    set(gca,'XTickLabel',[]);
    if i==3
        legend({'MSE','min(MSE)'},'Interpreter','latex')
    end
    % plot(lambda,m_eps)
    % set ( gca, 'xdir', 'reverse' )
    % xt = get(gca, 'xTick');
    % set(gca, 'xTickLabel',flip(xt))
    
    % figure
    subplot(2,3,i+3)
    plot(-log(lambda),m_ksi); hold on
    plot([-log(lambda(minCV)) -log(lambda(minCV))],ylims(i,:),'r--')
    xlabel(sprintf('$-\\mathrm{log}(\\lambda)$'),'Interpreter','latex')
    ylabel(sprintf('$\\mathbf{\\xi}_%d$',i),'Interpreter','latex')
    % ylim(ylims(i,:))
    % plot(lambda,m_ksi)
    % set ( gca, 'xdir', 'reverse' )
    % xt = get(gca, 'xTick');
    % set(gca, 'xTickLabel',flip(xt))
end