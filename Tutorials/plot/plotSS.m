function plotSS(XiSSE,lambda,m,n)

% logl = -log(lambda);

titn = ['x','y','z'];

% XiM = squeeze(mean(XiSS,4)); % inclusion probability   
XiEM = squeeze(mean(XiSSE,4)); % inclusion probability   

sStable = 0.8; % stable support
active = zeros(m,n);
active(2:3,1) = 1;
active([2:3 7],2) = 1;
active([4 6],3) = 1;

figure
for i = 1:n
    subplot(1,3,i)
    % plot(-log(lambda),squeeze(XiEM(:,i,:)))
    if i == 3
        plot([-1 -1],[0 0],'k--'); hold on % plot for legend
    end
    plot(-log(lambda),squeeze(XiEM(active(:,i)==1,i,:)),'b'); hold on
    plot(-log(lambda),squeeze(XiEM(active(:,i)==0,i,:)),'k--'); hold on
    plot([-log(lambda(1)) -log(lambda(end))],[sStable sStable], 'r')
    xlim([-log(lambda(1)) -log(lambda(end))])
    if i == 1
        ylabel('stability, or importance  measure','Interpreter','latex')
    end
    title(titn(i),'Interpreter','latex')
    if i==3
        legend({'non-active terms','true active terms'},'Interpreter','latex')
    end
    xlabel(sprintf('$-\\mathrm{log}(\\lambda)$'),'Interpreter','latex')
end



% ylims = [-12 12; -10 30; -4 2];
% 
% figure
% 
% for i = 1:n
%     eps = squeeze(epsilon(i,:,:));
%     ks = squeeze(Xi(:,i,:,:));
%     m_eps = mean(eps,2);
%     m_ksi = mean(ks,3);
% 
%     % figure
% 
%     subplot(2,3,i)
%     plot(-log(lambda),m_eps); hold on
%     plot(-log(lambda(minCV)),m_eps(minCV),'r*')
%     % xlabel('-log(lambda)')
%     if i == 1
%         ylabel(sprintf('$||\\mathbf{\\dot\\mathbf{x}-\\Theta\\xi}||_2^2$'),'Interpreter','latex')
%     end
%     title(titn(i),'Interpreter','latex')
%     set(gca,'XTickLabel',[]);
%     if i==3
%         legend({'MSE','min(MSE)'},'Interpreter','latex')
%     end
%     % plot(lambda,m_eps)
%     % set ( gca, 'xdir', 'reverse' )
%     % xt = get(gca, 'xTick');
%     % set(gca, 'xTickLabel',flip(xt))
% 
%     % figure
%     subplot(2,3,i+3)
%     plot(-log(lambda),m_ksi); hold on
%     plot([-log(lambda(minCV)) -log(lambda(minCV))],ylims(i,:),'r--')
%     xlabel(sprintf('$-\\mathrm{log}(\\lambda)$'),'Interpreter','latex')
%     ylabel(sprintf('$\\mathbf{\\xi}_%d$',i),'Interpreter','latex')
%     % ylim(ylims(i,:))
%     % plot(lambda,m_ksi)
%     % set ( gca, 'xdir', 'reverse' )
%     % xt = get(gca, 'xTick');
%     % set(gca, 'xTickLabel',flip(xt))
% end