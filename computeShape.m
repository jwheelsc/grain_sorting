
close all
clear all

run loadSample_specs
run mineral_colors.m

scf = 25/38;

M_Arr   = [2,   3,   4,   8,   20,  32,  38,  39,  42,  44] 
% M_Arr = 2
pgt_Arr = [0.8, 0.7, 0.8, 0.8, 0.8, 0.6, 0.8, 0.6, 0.8, 0.8]

sedA =     [1 2   3  4   5   6   7   8  9  10    11  12  13    14  15    16    17    18   19  20 ];
rockA =    [1 2   3  4   5   6   7   8  9  10    11  12  13    14  15    16    17    18   19   20  21  22  23  24  25  26 ];

msS = msS(sedA);
msR = msR(rockA);

MinDiam = 2

%%
for M = 1:length(M_Arr)

    M1 = M_Arr(M);
    pgt = pgt_Arr(M);

    %% Here are the plutonic rock samples

    folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\';
    [nmsR] = dir([folderR]);
    matNmR = {nmsR.name};
    matNmR = matNmR(3:end-1);
    matNmR = matNmR(rockA);
    matNmR = matNmR(msR==0);
    gType = 'P rock'

    for i = 1:length(matNmR)

        nmsRM = dir([folderR matNmR{i} '\*.mat']);  
        matNmRM = {nmsRM.name}

        D_elM = [];

        for k = 1:length(matNmRM)


            varSt = load([folderR matNmR{i} '\' matNmRM{k}]);   %% open the .mat file with mnrl MTx
            datIn = varSt.mnrlMtx;
            [mnrl,minsN,minNFull] = abbvMins(datIn,7);   %% concatenate some of the solution solution minerals
            mnrl(end,:) = [];
            numPix = sum(mnrl(:,2:37),2);
            mnrl = mnrl./numPix;
            D = sqrt(numPix*(scf^2)*4/pi);
            elM_L = mnrl(:,M1)>pgt  & D > MinDiam;  %% find the minerals with greater than x M
            elM = find(elM_L);       

            I_mtx = varSt.I_mtx;
            D_elM = D(elM);

            PF = [];
            AR = [];
            fft_R = [];
            pep = [];
            chp = [];
            fft_lgth = [];

            for j = 1:length((elM))    %% for each of those ptcls, grab it, and replace it in a new Mtx

                ind_j = elM(j);
                elMtx = find(I_mtx==ind_j);
                [rc] = ind2rc(size(I_mtx,1),elMtx);
                [ ptclP , cvxhP, arrN, arrNe ] = findPerims( rc );

                [cvxhP] = unique_j(cvxhP);
                cvxhP = [cvxhP;cvxhP(1,:)];

                [chp(j)]=perimLength(cvxhP(:,1),cvxhP(:,2));
                [pep(j)]=perimLength(ptclP(:,1),ptclP(:,2));
%                 PF(j) = (chp./pep).^2;

                n_p = length(arrNe);
                x_bar = sum(arrNe(:,1))/n_p; % find the center of mass of x
                y_bar = sum(arrNe(:,2))/n_p; % find the center of mass of y
                mP = [x_bar,y_bar];
                    
                close all
                f1 = figure(1)
                subplot(2,2,1)
%                 plot(arrN(:,1),arrN(:,2),'r+')
                hold on
                plot(arrNe(:,1),arrNe(:,2),'b.')
                plot(x_bar,y_bar,'o')
                plot(cvxhP(:,1),cvxhP(:,2),'k-','linewidth',2)
                plot(ptclP(:,1),ptclP(:,2),'b-+')

                [cvxhP] = interpPerim(cvxhP,360);
                figure(1)
                subplot(2,2,1)
                hold on
                plot(cvxhP(:,1),cvxhP(:,2),'k-');
                axis equal
                [ptclR, ptclA]=unrollParticle(mP,cvxhP);  %this unrols the particle, then finds the aspect ratio

                [ptclR_p, ptclA_p]=unrollParticle(mP,ptclP);  %this unrols the particle, then finds the aspect ratio


                [AR(j),R_tMaj] = aspectRatio(ptclR, ptclA, mP);

                subplot(2,2,2)
                p1 = plot(ptclR_p, ptclA_p,'b+')
                subplot(2,2,4)
                plot(ptclR_p, ptclA_p,'b+')
                subplot(2,2,2)

                xx = linspace(ptclR_p(1),ptclR_p(end),500);
                [ptclR_p,is] = unique(ptclR_p);
                yy = interp1(ptclR_p, ptclA_p(is),xx);

                hold on
                plot(xx,yy,'b-')
                p2 = plot(ptclR,ptclA,'k-','linewidth',2)

                legend([p1 p2],{'Perimeter','Convex hull'},'location','southeast')

                [fft_xy] = fft_ptclshape(xx,yy,mP);
                [fft_lgth(j)] = perimLength(fft_xy(:,1),fft_xy(:,2));
%                 fft_R(j) = (fft_lgth/pep).^2;

                subplot(2,2,3)
                hold on
                plot(ptclP(:,1),ptclP(:,2),'b-+')

                for sp = 1:2
                    subplot(2,2,(sp*2)-1)
                    grid on
                    ylabel('y (\mum)')
                    xlabel('x (\mum)')
                    set(gca,'fontsize',20)
                end
                subplot(2,2,1)
                text(5,35,'Convex hull','fontsize',20)
                
                subplot(2,2,3)
                text(5 ,35,'Fourier outline','fontsize',20)
    
                for sp = 1:2
                    subplot(2,2,(sp*2))
                    grid on
                    ylabel('Radius (\mum)')
                    xlabel('Angle from center (radians)')
                    set(gca,'fontsize',20)
                end
                
                savePDFfunction(f1,['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures_revisedCol\shape\grainOutline\grain_' num2str(ind_j)])
            return

            end
            
            ['i = ' num2str(i)]
            ['k = ' num2str(k)]
            minLB = minsN{M1}
            shapePR_Fl(i).Fd(k).(minLB) = [D_elM,AR',chp',pep',fft_lgth',elM]; 
        end
    end


    % end
return
    %%

    folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\';
    [nmsR] = dir([folderR]);
    matNmR = {nmsR.name};
    matNmR = matNmR(3:end-1);
    matNmR = matNmR(sedA);
    matNmR = matNmR(msS==0);
    gType = 'MX sed'

    for i = 1:length(matNmR)

        nmsRM = dir([folderR matNmR{i} '\*.mat']);  
        matNmRM = {nmsRM.name}

        D_elM = []

        for k = 1:length(matNmRM)


            varSt = load([folderR matNmR{i} '\' matNmRM{k}]);   %% open the .mat file with mnrl MTx
            datIn = varSt.mnrlMtx;
            [mnrl,minsN,minNFull] = abbvMins(datIn,7);   %% concatenate some of the solution solution minerals
            mnrl(end,:) = [];
            numPix = sum(mnrl(:,2:37),2);
            mnrl = mnrl./numPix;
            D = sqrt(numPix*(scf^2)*4/pi);

            elM_L = mnrl(:,M1)>pgt  & D > MinDiam;  %% find the minerals with greater than x M
            elM = find(elM_L);       

            I_mtx = varSt.I_mtx;
            D_elM = D(elM);

            PF = [];
            AR = [];
            fft_R = [];
            pep = [];
            chp = [];  
            fft_lgth = [];

            for j = 1:length((elM))    %% for each of those ptcls, grab it, and replace it in a new Mtx

                ind_j = elM(j);
                elMtx = find(I_mtx==ind_j);
                [rc] = ind2rc(size(I_mtx,1),elMtx);
                [ ptclP , cvxhP, arrN, arrNe ] = findPerims( rc );

                [cvxhP] = unique_j(cvxhP);
                cvxhP = [cvxhP;cvxhP(1,:)];

                [chp(j)]=perimLength(cvxhP(:,1),cvxhP(:,2));
                [pep(j)]=perimLength(ptclP(:,1),ptclP(:,2));
%                 PF(j) = (chp./pep).^2;

                n_p = length(arrNe);
                x_bar = sum(arrNe(:,1))/n_p; % find the center of mass of x
                y_bar = sum(arrNe(:,2))/n_p; % find the center of mass of y
                mP = [x_bar,y_bar];

                [cvxhP] = interpPerim(cvxhP,360);

                [ptclR, ptclA]=unrollParticle(mP,cvxhP);  %this unrols the particle, then finds the aspect ratio

                [ptclR_p, ptclA_p]=unrollParticle(mP,ptclP);  %this unrols the particle, then finds the aspect ratio

                [AR(j),R_tMaj] = aspectRatio(ptclR, ptclA, mP);

                xx = linspace(ptclR_p(1),ptclR_p(end),500);
                [ptclR_p,is] = unique(ptclR_p);
                yy = interp1(ptclR_p, ptclA_p(is),xx);

                [fft_xy] = fft_ptclshape(xx,yy,mP);
                [fft_lgth(j)] = perimLength(fft_xy(:,1),fft_xy(:,2));
%                 fft_R(j) = (fft_lgth/pep).^2;

            end
            ['i = ' num2str(i)]
            ['k = ' num2str(k)]
            minLB = minsN{M1}
            shapeMXS_Fl(i).Fd(k).(minLB) = [D_elM,AR',chp',pep',fft_lgth',elM]; 
        end
    end

    %%
    folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\';
    [nmsR] = dir([folderR])
    matNmR = {nmsR.name}
    matNmR = matNmR(3:end-1)
    matNmR = matNmR(rockA)
    matNmR = matNmR(msR==1)
    gType = 'MS rock'


    for i = 1:length(matNmR)

        nmsRM = dir([folderR matNmR{i} '\*.mat']);  
        matNmRM = {nmsRM.name};

        D_elM = []

        for k = 1:length(matNmRM)


            varSt = load([folderR matNmR{i} '\' matNmRM{k}]);   %% open the .mat file with mnrl MTx
            datIn = varSt.mnrlMtx;
            [mnrl,minsN,minNFull] = abbvMins(datIn,7);   %% concatenate some of the solution solution minerals
            mnrl(end,:) = [];
            numPix = sum(mnrl(:,2:37),2);
            mnrl = mnrl./numPix;
            D = sqrt(numPix*(scf^2)*4/pi);

            elM_L = mnrl(:,M1)>pgt  & D > MinDiam;  %% find the minerals with greater than x M
            elM = find(elM_L);       

            I_mtx = varSt.I_mtx;
            D_elM = D(elM);

            PF = [];
            AR = [];
            fft_R = [];
            pep = [];
            chp = [];
            fft_lgth = [];

            for j = 1:length((elM))    %% for each of those ptcls, grab it, and replace it in a new Mtx

                ind_j = elM(j);
                elMtx = find(I_mtx==ind_j);
                [rc] = ind2rc(size(I_mtx,1),elMtx);
                [ ptclP , cvxhP, arrN, arrNe ] = findPerims( rc );

                [cvxhP] = unique_j(cvxhP);
                cvxhP = [cvxhP;cvxhP(1,:)];

                [chp(j)]=perimLength(cvxhP(:,1),cvxhP(:,2));
                [pep(j)]=perimLength(ptclP(:,1),ptclP(:,2));
%                 PF(j) = (chp./pep).^2;

                n_p = length(arrNe);
                x_bar = sum(arrNe(:,1))/n_p; % find the center of mass of x
                y_bar = sum(arrNe(:,2))/n_p; % find the center of mass of y
                mP = [x_bar,y_bar];

                [cvxhP] = interpPerim(cvxhP,360);

                [ptclR, ptclA]=unrollParticle(mP,cvxhP);  %this unrols the particle, then finds the aspect ratio

                [ptclR_p, ptclA_p]=unrollParticle(mP,ptclP);  %this unrols the particle, then finds the aspect ratio

                [AR(j),R_tMaj] = aspectRatio(ptclR, ptclA, mP);

                xx = linspace(ptclR_p(1),ptclR_p(end),500);
                [ptclR_p,is] = unique(ptclR_p);
                yy = interp1(ptclR_p, ptclA_p(is),xx);

                [fft_xy] = fft_ptclshape(xx,yy,mP);
                [fft_lgth(j)] = perimLength(fft_xy(:,1),fft_xy(:,2));
%                 fft_R(j) = (fft_lgth/pep).^2;

            end
            ['i = ' num2str(i)]
            ['k = ' num2str(k)]
            minLB = minsN{M1}
            shapeMSR_Fl(i).Fd(k).(minLB) = [D_elM,AR',chp',pep',fft_lgth',elM]; 
        end
    end

    %%
    folderR = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\';
    [nmsR] = dir([folderR]);
    matNmR = {nmsR.name};
    matNmR = matNmR(3:end-1);
    matNmR = matNmR(sedA);
    matNmR = matNmR(msS==1);
    gType = 'MS sed'

    for i = 1:length(matNmR)

        nmsRM = dir([folderR matNmR{i} '\*.mat']);  
        matNmRM = {nmsRM.name}

        D_elM = [];

        for k = 1:length(matNmRM)


            varSt = load([folderR matNmR{i} '\' matNmRM{k}]);   %% open the .mat file with mnrl MTx
            datIn = varSt.mnrlMtx;
            [mnrl,minsN,minNFull] = abbvMins(datIn,7);   %% concatenate some of the solution solution minerals
            mnrl(end,:) = [];
            numPix = sum(mnrl(:,2:37),2);
            mnrl = mnrl./numPix;
            D = sqrt(numPix*(scf^2)*4/pi);

            elM_L = mnrl(:,M1)>pgt  & D > MinDiam;  %% find the minerals with greater than x M
            elM = find(elM_L);       

            I_mtx = varSt.I_mtx;
            D_elM = D(elM);

            PF = [];
            AR = [];
            fft_R = [];
            pep = [];
            chp = [];
            fft_lgth = [];

            for j = 1:length((elM))    %% for each of those ptcls, grab it, and replace it in a new Mtx

                ind_j = elM(j);
                elMtx = find(I_mtx==ind_j);
                [rc] = ind2rc(size(I_mtx,1),elMtx);
                [ ptclP , cvxhP, arrN, arrNe ] = findPerims( rc );

                [cvxhP] = unique_j(cvxhP);
                cvxhP = [cvxhP;cvxhP(1,:)];

                [chp(j)]=perimLength(cvxhP(:,1),cvxhP(:,2));
                [pep(j)]=perimLength(ptclP(:,1),ptclP(:,2));
%                 PF(j) = (chp./pep).^2;

                n_p = length(arrNe);
                x_bar = sum(arrNe(:,1))/n_p; % find the center of mass of x
                y_bar = sum(arrNe(:,2))/n_p; % find the center of mass of y
                mP = [x_bar,y_bar];

                [cvxhP] = interpPerim(cvxhP,360);

                [ptclR, ptclA]=unrollParticle(mP,cvxhP);  %this unrols the particle, then finds the aspect ratio

                [ptclR_p, ptclA_p]=unrollParticle(mP,ptclP);  %this unrols the particle, then finds the aspect ratio

                [AR(j),R_tMaj] = aspectRatio(ptclR, ptclA, mP);

                xx = linspace(ptclR_p(1),ptclR_p(end),500);
                [ptclR_p,is] = unique(ptclR_p);
                yy = interp1(ptclR_p, ptclA_p(is),xx);

                [fft_xy] = fft_ptclshape(xx,yy,mP);
                [fft_lgth(j)] = perimLength(fft_xy(:,1),fft_xy(:,2));
%                 fft_R(j) = (fft_lgth/pep).^2;

            end
            ['i = ' num2str(i)]
            ['k = ' num2str(k)]
            minLB = minsN{M1}
            shapeMSS_Fl(i).Fd(k).(minLB) = [D_elM,AR',chp',pep',fft_lgth',elM]; 
        end
    end

end

save('D:\Code\Summer_2013_data\mineral_data\qemscan_tif\dist_out_stats\shapePtcl_3.mat','shapeMSR_Fl','shapeMSS_Fl','shapePR_Fl','shapeMXS_Fl')
