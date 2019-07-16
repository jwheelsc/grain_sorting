function [fft_xy] = fft_ptclshape(xx,yy,mP,NLi)

    y_bar = mP(2);
    x_bar = mP(1);
    th = linspace(0,2*pi,500);
    T = th(2)-th(1);
    Fs = 1/T;
    lx = length(xx);
    jf = 0:lx-1; 
    kArr = [1:20];

    YF = fft(yy);
    P2 = abs(YF/lx);
    P1 = P2(1:lx/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    fF = Fs*(0:(lx/2))/lx;
    subplot(2,2,3)
    hold on
%     figure
%     plot(fF*2*pi,P1)


%     NLArr = [1,5,10,50]

%     colr = {'r','m','g','k'}
%     NLi = 5;
%     for NL = 1:length(NLArr)
        filg = [];
        g = [];
%         NLi = NLArr(NL)

%         Y = zeros(1,length(th))
%         for N = 1:length(kArr)
%             Y = (A_k(N)*cos(N*th) + B_k(N)*sin(N*th))+Y
%         end

        NFFT = lx;
        g = fft(yy,NFFT);
        g = g(1:length(g)/2);
        realg = abs(g);
        angleg = unwrap(angle(g));
        modes = [1:length(g)];
        NLi = 50;
        g( modes > NLi ) = 0;
        filg = ifft(g,NFFT,'symmetric');
%     f2=figure(2)
%     hold on
%     subplot(2,2,4)
%         hold on
%         p4(NL) = plot(xx,filg,'-','color',colr{NL},'linewidth',2)
%         p4 = plot(xx,filg,'-','color','r','linewidth',2)

%     subplot(2,2,3)
%     f1=figure(1)
%     hold on
    xE = filg.*cos(xx)+x_bar;
    yE = filg.*sin(xx)+y_bar;

    fft_xy = [xE',yE'];
% keyboard
% if NLi == 5
% hold on
%     p3 = plot([xE,xE(1)],[yE,yE(1)],'-','color','m','linewidth',2)
% elseif NLi == 10
%         p3 = plot([xE,xE(1)],[yE,yE(1)],'-','color','r','linewidth',2)
% end

%     axis equal
%         keyboard

%     end
% 
%     subplot(2,2,4)
%     legend([p4(1) p4(2) p4(3) p4(4)],{'k = 1','k = 5','k = 10','k = 50'},'location','southeast')

end





