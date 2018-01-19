function run_chern(blockIndex,totalBlocks,ysteps)
tic

xi0 = 15;

Nc = 8*xi0;
ksteps = 7000;

kFamin = 5.2*pi;
kFamax = 5.6*pi;


ychoice = 1; % choose the y-axis (1 for eps0, 2 for lambda)
if ychoice == 1
    
    lam = 0.05;
    eps0min = -0.3;
    eps0max = 0.3;
    
    kFa = kFamin + (blockIndex-1)*(kFamax-kFamin)/(totalBlocks-1);  % the value of kFa (the point in the x-axis)
    eps0vec = linspace(eps0min,eps0max,ysteps);
    
    if blockIndex == 1
           save('out/outc1-par','totalBlocks','ysteps','eps0min', 'eps0max', 'lam', 'xi0','Nc','ksteps');
    end
    
    
    qvec = zeros(1,ysteps); % initialize the Chern number vector
    for ind = 1:ysteps
        eps0 = eps0vec(ind);
        xi = xi0;
        qvec(ind) = chern(lam,xi,kFa,eps0,Nc,ksteps);
    end
    
elseif ychoice == 2
    
    eps0 = 0.0;
    lammin = 0;
    lammax = 0.05;
    
    kFa = kFamin + (blockIndex-1)*(kFamax-kFamin)/(totalBlocks-1);  % the value of kFa (the point in the x-axis)
    lamvec = linspace(lammin,lammax,ysteps);
    
    if blockIndex == 1
        save('out/outc1-par','totalBlocks','ysteps','lammin', 'lammax', 'eps0', 'xi0','Nc','ksteps');
    end
    
    
    qvec = zeros(1,ysteps); % initialize the Chern number vector
    for ind = 1:ysteps
        lam = lamvec(ind);
        xi = xi0*sqrt(1+lam^2);
        qvec(ind) = chern(lam,xi,kFa,eps0,Nc,ksteps);
    end
    
end


t = toc

filename=strcat('out/outc1-',int2str(blockIndex));
save(filename, 'kFa', 'qvec','t');


fprintf('SUCCESS blockIndex %d\n',blockIndex);
exit;
