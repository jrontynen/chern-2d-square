% clear('totalBlocks','eps0min', 'eps0max', 'lam', 'xi0','Nc','ksteps'); 
clear('totalBlocks','lammin', 'lammax', 'eps0', 'xi0','Nc','ksteps');

load('out/outc1-par')

kFavec = [];
time = 0;
% eps0vec = linspace(eps0min,eps0max,ysteps);
lamvec = linspace(lammin,lammax,ysteps);
qmat = zeros(ysteps,totalBlocks);

for index = [1:totalBlocks]

    % read the output from the jobs
    filename = strcat( 'out/outc1-', int2str(index) );
    load(filename);
    
    kFavec = cat(2,kFavec,kFa);
    qmat(:,index) = qvec;

    if t > time
        time = t;
    end
end

nhrs = floor(time/3600);
nmins = floor( (time - nhrs*3600)/60 );
sprintf('time = %i h %i min',nhrs,nmins)





%Select only the Chern numbers within tolerance:
tol = 0.1;
qvec = qmat(:);
qtol = qvec( abs(qvec-round(qvec)) < tol );


%Count the number of different Chern numbers within tolerance:
bins = max(-500,round(min(qtol))):min(500,round(max(qtol)));
edges = [(bins-0.5) (bins(end)+0.5)];
counts = histcounts(qtol,edges);

%Select only the Chern numbers with a minimum count:
mincount = 50;
I = find(counts>mincount);
bins_mincount = bins(I) % Chern numbers with at least the minimum count

%And make a Gaussian fit:
I2 = min(I):max(I);
F = fit(bins(I2).',counts(I2).','gauss1');
coeval = coeffvalues(F);



%Create a color scheme for the phase diagram:

caxisstep = 0.4; % give as a fraction of 10, that is 0.1, 0.2, etc
Ncolor = 2*10*caxisstep;
Nblank = 10 - Ncolor;

blank = repmat([0 0 0], Nblank, 1);
skipcolor = repmat([0 0 0], Ncolor, 1);
maincolors = [0 0 0.5; 0 0 1; 0.6 0.2 0.8; 0 1 1; 1 1 0; ...
    1 1 1; ...
    0 1 0; 1 0.1 0.6; 1 0.5 0; 1 0 0; 0.5 0 0];
shade1 = repmat([0.7 0.95 0.4], Ncolor, 1);
shade2 = repmat([0.15 0.2 0.1], Ncolor, 1);
shade3 = repmat([0.5 0.75 0.9], Ncolor, 1);
shade4 = repmat([0.9999 0.5 0.5], Ncolor, 1);
% shade4 = [0 0.3 0.3; 0 0.3 0.3];

Iup = find(bins_mincount > 5);
Idown = find(bins_mincount < -5);
inc1 = 0;
inc2 = 0;
if length(Iup) > 1
    inc1 = (shade2-shade1)/(length(Iup)-1);
end
if length(Idown) > 1
    inc2 = (shade4-shade3)/(length(Idown)-1);
end

map = [];
for integer = bins_mincount(1):bins_mincount(end)
    
    if integer < -5
        
        if any(bins_mincount == integer)
            map = [map; shade3; blank];
            shade3 = shade3 + inc2;
        else
            map = [map; skipcolor; blank];
        end
        
    end
    if (integer >= -5) && (integer <= 5)
        
        if any(bins_mincount == integer)
            map = [map; repmat(maincolors(integer+6,:),Ncolor,1); blank];
        else
            map = [map; skipcolor; blank];
        end
        
    end
    if integer > 5
        
        if any(bins_mincount == integer)
            map = [map; shade1; blank];
            shade1 = shade1 + inc1;
        else
            map = [map; skipcolor; blank];
        end
        
    end
    
end


%Plot the phase diagram:

figure
% imagesc(kFavec/pi,eps0vec,qmat)
imagesc(kFavec/pi,lamvec,qmat)
set(gca,'YDir','normal')
contourcbar;
% title(sprintf('\\lambda=%0.2f, \\xi=%0.0f, N_c=%0.0f, #k=%i^2, tol=%0.2f, cstep=%0.2f, mincount=%i'...
%     ,lam,xi0,Nc,ksteps,tol,caxisstep,mincount),'fontweight','normal','fontsize',9)
title(sprintf('\\epsilon_0=%0.2f, \\xi=%0.0f, N_c=%0.0f, #k=%i^2, tol=%0.2f, cstep=%0.2f, mincount=%i'...
    ,eps0,xi0,Nc,ksteps,tol,caxisstep,mincount),'fontweight','normal','fontsize',9)
xlabel('k_Fa/\pi')
% ylabel('\epsilon_0')
ylabel('\lambda')

map = map(1:end-Nblank,:);
caxis([bins_mincount(1) - caxisstep, bins_mincount(end) + caxisstep])
colormap(map)
    





% 
% figure
% histogram(qtol,edges)
% xlim([bins_mincount(1) bins_mincount(end)])
% hold on
% h = plot(F);
% hold off
% title(sprintf('k_Fa/\\pi=%0.2f...%0.2f, \\lambda=%0.2f, \\xi=%0.0f, N_{1c}=%0.0f'...
%     ,kFavec(1)/pi,kFavec(end)/pi,lam,xi0,Nc),'fontweight','normal','fontsize',9)
% % title(sprintf('k_Fa/\\pi=%0.2f...%0.2f, \\epsilon_0=%0.2f, \\xi=%0.0f, N_{1c}=%0.0f'...
% %     ,kFavec(1)/pi,kFavec(end)/pi,eps0,xi0,Nc),'fontweight','normal','fontsize',9)
% legend(h,sprintf('%0.0f*e^{-((x-%0.2f)/%0.2f)^2}',coeval))
% xlabel('q')
% ylabel('#q')
% F


