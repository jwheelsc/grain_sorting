clc
close all
clear variables

run mineral_colors

folder = 'D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\images\revisedColors\rock_frags\MI8152-AUG16\1 - 25 um\';
[nms] = dir([folder '\*.TIF'])
matNm = {nms.name}


for H = 1:length(matNm)
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

    %% grid search
    %%% here you cruise through each pixel of interest, and see what minerlas
    %%% are touching the other minerals. 

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

return
    % here you can give each tagged pixel a sequentially lower tag with no
    % missing values
    xU = unique(I_mtx);
    xArr = 0:length(xU)-1;
    for i = 1:length(xU)
        I_mtx(I_mtx==xU(i))=xArr(i);
    end

   

    %% size and mineralogy
    clc
    mnrlMtx = zeros(length(xU),length(mins));
    I_mtx_C = zeros(size(I_mtx));

    for i = 1:length(xU)
    %     figure
    %     testRGB(I_mtx,i)
        sizeMtx(i) = sum(sum(I_mtx==i));

        M = find(I_mtx==i);

        M_Arr = zeros(1,length(mins));
        for k = 1:length(M)

            Mk = M(k);
            rc = ind2rc(size(Irgb,1),Mk);
            rw = int32(rc(1));
            cl = int32(rc(2));   
            rgbAr = [Irgb(rw,cl,1),Irgb(rw,cl,2),Irgb(rw,cl,3)];
            diffMtx = min_col-rgbAr;
            el_k = find(sum(abs(diffMtx),2)==0);          
            M_Arr(el_k) = M_Arr(el_k)+1;
            I_mtx_C(rw,cl) = el_k;

        end
        mnrlMtx(i,:)=M_Arr;
 
    end

    %%
    save(['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\' img '.mat'],'I_mtx','I_mtx_C','mnrlMtx')

end












