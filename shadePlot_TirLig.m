lonGrid_sP = unique(TUVclean.LonLat(:,1));
latGrid_sP = unique(TUVclean.LonLat(:,2));

TOT_grid = NaN.*ones(length(latGrid_sP),length(lonGrid_sP),1);

for i=1:length(TUVclean.LonLat(:,1))
    TUVclean.totVel(i) = sqrt((TUVclean.U(i))^2 + (TUVclean.V(i))^2);
end

for i=1:length(TUVclean.LonLat(:,1))
    lonGrid_idx = find(lonGrid_sP==TUVclean.LonLat(i,1));
    latGrid_idx = find(latGrid_sP==TUVclean.LonLat(i,2));
    
    if (not(isnan(TUVclean.totVel(i))))
        TOT_grid(latGrid_idx,lonGrid_idx,1) = TUVclean.totVel(i);
    end
end

% Create the lat/lon grid for pcolor drawing
[Plg,Plt]=meshgrid(lonGrid_sP,latGrid_sP);

m_pcolor(Plg, Plt, TOT_grid);
shading flat;
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box', 'fancy', 'tickdir', 'in', 'xlabeldir','end','fontsize',10);

[X,Y]=m_ll2xy(10.2373,43.8579);
line(X,Y,'marker','square','markersize',4,'color','r');
text(X,Y,' VIAR','vertical','top');

[X,Y]=m_ll2xy(9.6593000,44.1435167);
line(X,Y,'marker','square','markersize',4,'color','r');
text(X,Y,' PCOR','vertical','top');

[X,Y]=m_ll2xy(9.8492167,44.0263667);
line(X,Y,'marker','square','markersize',4,'color','r');
text(X,Y,' TINO','vertical','top');

hold on;

un=TUVclean.U./sqrt(TUVclean.U.^2+TUVclean.V.^2);
wn=TUVclean.V./sqrt(TUVclean.U.^2+TUVclean.V.^2);

m_quiver(TUVclean.LonLat(:,1),TUVclean.LonLat(:,2),un,wn, 'k');

clear Plg Plt TOT_grid un wn
