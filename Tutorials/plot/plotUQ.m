function plotUQ(XiE,XiTrue,Xi,lib)

%% colors figure1
% blue: x
blue = [44,127,184]./255;
% green: y
green = [49,163,84]./255;
% orange: z
orange = [240,59,32]./255;
colorsNew = [blue; green; orange];

grey = [1 1 1]*0.7;

Ct = [80 80 80]/255;
Clb2 = [187 187 187]/255;

lw = 1.0;
lwT = 3;

fos = 10;%16; % fontsize

sizeX = 600;
sizeY = 500;

figure('Position', [10 10 sizeX sizeY])

nbins = 50; % 10;%
nP1 = size(XiE,1);
nP2 = size(XiE,2);

xmin = [-12 -3 -5];
xmax = [12 30 3];
nn = 1;
for ip = 1:nP1
    for jp = 1:nP2
        xPtest = XiE(ip,jp,:);
        xPtest = xPtest(:);
        
        Clb = colorsNew(jp,:);
        
        subplot(nP1,nP2,nn)
        if sum(xPtest~=0) > 0
            % h = histfit(xPtest,nbins,'kernel'); hold on
            % % h = histfit(xPtest,nbins); hold on
            % h(1).FaceColor = Clb2;
            % h(1).FaceAlpha = 0;
            % h(1).EdgeColor = Clb2;
            % h(1).EdgeAlpha = 0;
            h = histogram(xPtest,nbins); hold on
            h.FaceColor = 'k';
            h.FaceAlpha = 0.9;
            h.EdgeAlpha = 0;
            % h(2).Color = Clb2;
            % h(2).LineWidth = 1;
            % xlim([xmin(jp) xmax(jp)])
            % mY = max(h(2).YData);
            % area(h(2).XData,h(2).YData,'FaceColor',Clb,'FaceAlpha',0.4); hold on
            mY = max(h.Values);
            if isinf(mY) || isnan(mY)
                mY = 1000;
            end
            ylim([0 mY])
            % if (XiTrue(ip,jp)~=0)
            %     patch([xmin(jp) xmin(jp) xmax(jp) xmax(jp)],[0 mY mY 0],Clb,'FaceAlpha',0.1); hold on
            % end
            if (XiTrue(ip,jp)~=0)
                % patch([xmin(jp) xmin(jp) xmax(jp) xmax(jp)],[0 mY mY 0],Clb,'FaceAlpha',0.1); hold on
                set(gca, 'XColor',Clb, 'YColor',Clb)       
            else
                set(gca, 'XColor',[.7 .7 .7], 'YColor',grey)
            end
            if XiTrue(ip,jp) ~= 0
                plot(XiTrue(ip,jp),0,'x','Color','m','Linewidth',lw,'MarkerSize',14); hold on
            end
            if Xi(ip,jp) ~= 0
                plot([Xi(ip,jp) Xi(ip,jp)],[0 mY],'-','Color',Clb,'Linewidth',lwT); hold on
            else
                plot([mean(xPtest) mean(xPtest)],[0 mY],'-','Color',grey,'Linewidth',lw); hold on
            end
        else
            plot([0 0],[0 1],'-','Color',grey,'Linewidth',lw)
            ylim([0 1])
            set(gca, 'XColor',[.7 .7 .7], 'YColor',[.7 .7 .7])
        end
        % xlim([xmin(jp) xmax(jp)])
      
        yticks([])
        % if ip~=nP1
        %     xticks([])
        % end
        % if ip==nP1
        %     set(gca,'ticklabelinterpreter','latex','FontSize',fos)
        % end
        
        if ip == 1
        	if jp == 1
                title('xDot','interpreter','latex','FontSize',fos)
            elseif jp == 2
                title('yDot','interpreter','latex','FontSize',fos)
            else 
                title('zDot','interpreter','latex','FontSize',fos)
            end
        end
        
        if jp == 1
            ylabel(lib{ip,1},'interpreter','latex','FontSize',fos,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
        end
        
        nn = nn + 1;
    end
end

