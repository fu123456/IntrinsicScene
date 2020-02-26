function [I,R]=intrinsicScene(src,alpha,beta,lambda,vareps,K,debug)
% Notes: alpha, weight for illumination beta, weight for reflectance. If you
% want to obtain better results, you can modify the parameters: alpha, beta,
% lambda
    if (~exist('alpha','var'))
        alpha=0.1;
    end
    if (~exist('beta','var'))
        beta=0.001;
    end
    if (~exist('lambda','var'))
        lambda=0.01;
    end
    if (~exist('vareps','var'))
        vareps=0.002;
    end
    if (~exist('K','var'))
        K=30;
    end
    if (~exist('debug','var'))
        debug=true;
    end
    eps=0.0001;
    if size(src,3)==1
        S=src;
        gray=src;
    else
        hsv=rgb2hsv(src);
        S=hsv(:,:,3);
        gray=rgb2gray(src);
    end

    B=ones(size(S)).*0.5;

    % intialization of I and S
    I=S;%single channel
    I=ones(size(S));
    R=S;

    % eplison
    eplisonR_plot=[];
    eplisonI_plot=[];
    if debug==true
        fprintf('-- Stop iteration until eplison < %02f or K>%d\n',vareps,K);
    end
    for iter=1:K
        preI=I;
        preR=R;

        %% update reflectance
        R=S./I;
        R=L0Smoothing(R,beta);
        eplisonR=norm(R-preR,'fro')/norm(preR,'fro');
        eplisonR_plot=[eplisonR_plot eplisonR];

        %% update illuination
        ux=ones(size(S));
        uy=ones(size(S));
        ux(:,end)=0;
        uy(:,end)=0;
        I=solveLinearSystem(S,R,ux,uy,alpha,B,lambda);
        eplisonI=norm(I-preI,'fro')/norm(preI,'fro');
        eplisonI_plot=[eplisonI_plot eplisonI];

        %% iteration until convergence
        if debug==true
            fprintf('Iter #%d : eplisonI=%f; eplisonR=%f\n',iter,eplisonI,eplisonR);
        end
        if(eplisonI<vareps||eplisonR<vareps)
            break;
        end
    end
    % final result
    if size(src,3)==3
        hsv(:,:,3)=R;
        R=hsv2rgb(hsv);
    end
end

function dst = solveLinearSystem(s, ir, uvx, uvy, alphabet, b, lambda)
    [h, w] = size(s);
    hw = h * w;
    uvx = uvx(:);
    uvy = uvy(:);
    ux = padarray(uvx, h, 'pre'); ux = ux(1:end-h);
    uy = padarray(uvy, 1, 'pre'); uy = uy(1:end-1);
    D = uvx+ux+uvy+uy;
    T = spdiags([-uvx, -uvy],[-h,-1],hw,hw);
    MN = T + T' + spdiags(D, 0, hw, hw);
    ir2 = ir.^2;
    ir2 = spdiags(ir2(:),0,hw,hw);
    DEN = ir2 + alphabet * MN + lambda * speye(hw,hw);
    NUM = ir.*s + lambda * b;
    L = ichol(DEN,struct('michol','on'));
    [dst,~] = pcg(DEN, NUM(:), 0.01, 40, L, L');
    dst = reshape(dst, h, w);
end

%% this function is from
% Xu, L., Lu, C., Xu, Y., & Jia, J., Image smoothing via $l_0$ gradient
% minimization, Transactions on Graphics, 30(6), 174 (2011).
function S = L0Smoothing(Im, lambda, kappa)
    if ~exist('kappa','var')
        kappa = 2.0;
    end
    if ~exist('lambda','var')
        lambda = 2e-2;
    end
    S = im2double(Im);
    betamax = 1e5;
    fx = [1, -1];
    fy = [1; -1];
    [N,M,D] = size(Im);
    sizeI2D = [N,M];
    otfFx = psf2otf(fx,sizeI2D);
    otfFy = psf2otf(fy,sizeI2D);
    Normin1 = fft2(S);
    Denormin2 = abs(otfFx).^2 + abs(otfFy ).^2;
    if D>1
        Denormin2 = repmat(Denormin2,[1,1,D]);
    end
    beta = 2*lambda;
    while beta < betamax
        Denormin   = 1 + beta*Denormin2;
        % h-v subproblem
        h = [diff(S,1,2), S(:,1,:) - S(:,end,:)];
        v = [diff(S,1,1); S(1,:,:) - S(end,:,:)];
        if D==1
            t = (h.^2+v.^2)<lambda/beta;
        else
            t = sum((h.^2+v.^2),3)<lambda/beta;
            t = repmat(t,[1,1,D]);
        end
        h(t)=0; v(t)=0;
        % S subproblem
        Normin2 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)];
        Normin2 = Normin2 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)];
        FS = (Normin1 + beta*fft2(Normin2))./Denormin;
        S = real(ifft2(FS));
        beta = beta*kappa;
    end
end