clc
close all
clear variables

run mineral_colors

folder = 'D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\images\revisedColors\rock_frags\MI8152-AUG16\1 - 25 um\';
[nms] = dir([folder '\*.TIF'])
matNm = {nms.name}

for H = 1:length(matNm)
%     for H = 2
    fname = [folder matNm{H}];
    
    nm = matNm{H};
    for h = 1:length(nm)
    
        if strcmp(nm(h),'.')
            img = nm(1:h-1);
        end
        
    end
 
    [I,cmap] = imread(fname);

    Irgb = ind2rgb(I,cmap);
    figure
    imshow(Irgb)

    %% grid search for all grains
    %%% here you cruise through each pixel of interest, and see what minerlas
    %%% are touching the other minerals. Here you don't really care what
    %%% the minerals are, you're just finding isolated grains

    lenMin = length(mins);
    [num_r,num_c] = size(I);
    ar_i = [-1 0 1 1 1 0 -1 -1];
    ar_j = [-1 -1 -1 0 1 1 1 0];

    I_mtx = zeros(num_r,num_c);
    I_ack = zeros(num_r,num_c);
    count = 1;

    % in this first loop you go through and row by column and assign colored
    % pixels a tag. First, you find the colored pixel, then cirlce to see if
    % there is a neighbouring colored pixel with a previous tag, and if not, give
    % a new one
    for j = 2:num_c-1

        for i = 2:num_r-1

            rgbA = [Irgb(i,j,1),Irgb(i,j,2),Irgb(i,j,3)]; % get the rgb of the ith and jth pixel

            if sum(rgbA)<3

                for k = 1:8

                    i_in = i+ar_i(k);
                    j_in = j+ar_j(k);
                    pix(k,:) = [i_in,j_in];
                    I_mtx_k(k) = I_mtx(i_in,j_in);

                end

                elG = find(I_mtx_k>0);

                if elG
    %                 I_mtx_k(elG(1))
                    I_mtx(i,j) = I_mtx_k(elG(1));
    %                 keyboard
                elseif isempty(elG)
                    I_mtx(i,j) = count;
                end

                count = count+1;

            end
        end
    end

    num_u = []
    num_u(1) = length(unique(I_mtx));

    % because the first loop can't get it right, you go through each tagged
    % pixel to see if a neighbouring pixel has a lower value, and grab it. Do
    % this enough times such that there is no longer a change in the number of
    % tags
    d_numu(1) = 10
    u = 1
    while d_numu > 1

        for j = 2:num_c-1

            for i = 2:num_r-1

                I_ind = I_mtx(i,j);

                if I_ind~=0

                    for k = 1:8

                        i_in = i+ar_i(k);
                        j_in = j+ar_j(k);
                        pix(k,:) = [i_in,j_in];
                        a(k) = I_mtx(i_in,j_in);

                    end

                    new_ind = min(a(find(a>0)));

                    if new_ind
                        I_mtx(i,j) = new_ind;
                    end

                end
            end
        end

    num_u(u+1) = length(unique(I_mtx));
    d_numu = abs(num_u(u+1)-num_u(u));
    u = u+1;
    end

    figure
    plot(1:u,num_u,'-o')


    % here you can give each tagged pixel a sequentially lower tag with no
    % missing values
    xU = unique(I_mtx);
    xArr = 0:length(xU)-1;
    for i = 1:length(xU)
        I_mtx(I_mtx==xU(i))=xArr(i);
    end
    
    % size and mineralogy

    for i = 1:length(xU)
        sizeMtx(i) = sum(sum(I_mtx==i));
    end

    

    %% Here you repeat the same process but with just the mineral that you care about
    for m = 2:length(mins)
        M_col = min_col(m,:)
        % turn all non M minerals to white
        I_Mr = (Irgb(:,:,1)==M_col(1));
        I_Mg = (Irgb(:,:,2)==M_col(2));
        I_Mb = (Irgb(:,:,3)==M_col(3));
        I_M =  1-(I_Mr.*I_Mg.*I_Mb);
% 
%         close all
%         imshow(I_M)    
%         keyboard
        %% grid search for all grains of mineral M
        %%% here you cruise through each pixel and find the isolated clusters
        %%% of mineral M

        [num_r,num_c] = size(I);
        ar_i = [-1 0 1 1 1 0 -1 -1];
        ar_j = [-1 -1 -1 0 1 1 1 0];

        I_mtxM = zeros(num_r,num_c);
        I_ackM = zeros(num_r,num_c);
        count = 1;

        % in this first loop you go through and row by column and assign colored
        % pixels a tag. First, you find the colored pixel, then cirlce to see if
        % there is a neighbouring colored pixel with a previous tag, and if not, give
        % a new one
        for j = 2:num_c-1

            for i = 2:num_r-1

                if I_M(i,j) == 0

                    for k = 1:8

                        i_in = i+ar_i(k);
                        j_in = j+ar_j(k);
                        pix(k,:) = [i_in,j_in];
                        I_mtxM_k(k) = I_mtxM(i_in,j_in);

                    end

                    elG = find(I_mtxM_k>0);

                    if elG
        %                 I_mtx_k(elG(1))
                        I_mtxM(i,j) = I_mtxM_k(elG(1));
                    elseif isempty(elG)
        %                 keyboard
                        I_mtxM(i,j) = count;
                    end

                    count = count+1;

                end
            end
        end

        num_u = []
        num_u(1) = length(unique(I_mtxM));

        % because the first loop can't get it right, you go through each tagged
        % pixel to see if a neighbouring pixel has a lower value, and grab it. Do
        % this enough times such that there is no longer a change in the number of
        % tags
        d_numu(1) = 10
        u = 1
        while d_numu > 1

            for j = 2:num_c-1

                for i = 2:num_r-1

                    I_ind = I_mtxM(i,j);

                    if I_ind~=0

                        for k = 1:8

                            i_in = i+ar_i(k);
                            j_in = j+ar_j(k);
                            pix(k,:) = [i_in,j_in];
                            a(k) = I_mtxM(i_in,j_in);

                        end

                        new_ind = min(a(find(a>0)));

                        if new_ind
                            I_mtxM(i,j) = new_ind;
                        end

                    end
                end
            end

        num_u(u+1) = length(unique(I_mtxM));
        d_numu = abs(num_u(u+1)-num_u(u));
        u = u+1;
        end

%         figure
%         plot(1:u,num_u,'-o')


        % here you can give each tagged pixel a sequentially lower tag with no
        % missing values
        xUM = unique(I_mtxM);
        xArr = 0:length(xUM)-1;
        for i = 1:length(xUM)
            I_mtxM(I_mtxM==xUM(i))=xArr(i);
        end



        %% size and mineralogy
        sizeMtxM = []
        for i = 2:length(xUM)-1
            elS = find(I_mtxM==i);
            indS = I_mtx(elS(1));
            sizeMtxM(i,2) = sum(sum(I_mtxM==i));
            sizeMtxM(i,1) = sizeMtx(indS);

        end
%         close all
%         plot(sizeMtxM(:,1),sizeMtxM(:,2),'.')
        clustM.(mins{m}) = sizeMtxM;

    end
    %%
    return
    save(['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\mineral_clusters\rock\' img '.mat'],'clustM','I_mtx','I_mtxM')

end












