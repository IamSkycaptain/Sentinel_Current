%% Correction part
clear
load ocn_list
load('scal_corr_sw1.mat')
load('scal_corr_sw2.mat')
load('scal_corr_sw3.mat')
scal_corr_sw{1} = scal_corr_sw1;
scal_corr_sw{2} = scal_corr_sw2;
scal_corr_sw{3} = scal_corr_sw3;
load('range_corr_sw1.mat')
load('range_corr_sw2.mat')
load('range_corr_sw3.mat')
range_corr_sw{1} = range_corr_sw1;
range_corr_sw{2} = range_corr_sw2;
range_corr_sw{3} = range_corr_sw3;
corr_range_v_sw = {};
load('corr_azi.mat')
load('off_land_median.mat')
load('corr_azi_sw.mat')
load('off_land_median_sw.mat')

for p = 1:size(ocn_list,1)
    clear lon_sw lat_sw v_sw h_sw inc_sw land_sw
    ncfile = ocn_list(p,:);
    lon = ncread(ncfile,'rvlLon');
    lat = ncread(ncfile,'rvlLat');
    v = ncread(ncfile,'rvlRadVel');
    h = ncread(ncfile,'rvlHeading');
    inc = ncread(ncfile,'rvlIncidenceAngle');
    land = ncread(ncfile,'rvlLandFlag');
    [swath,row,col] = size(v);
    for sw = 1:3
        for r = 1:row
            for c = 1:col
                lat_sw{sw}(r,c) = lat(sw,r,c);
                lon_sw{sw}(r,c) = lon(sw,r,c);
                v_sw{sw}(r,c) = v(sw,r,c);
                h_sw{sw}(r,c) = h(sw,r,c);
                inc_sw{sw}(r,c) = inc(sw,r,c);
                land_sw{sw}(r,c) = land(sw,r,c);
            end
        end
    end
    clear lon lat v h r c row col inc land
    lon = [lon_sw{1}(:);lon_sw{2}(:);lon_sw{3}(:)];
    lat = [lat_sw{1}(:);lat_sw{2}(:);lat_sw{3}(:)];
    land = [land_sw{1}(:);land_sw{2}(:);land_sw{3}(:)];
    v = [v_sw{1}(:);v_sw{2}(:);v_sw{3}(:)];
    
   	ocean_flag = 1; %-------- Use only Oceanic pixels
    if ocean_flag==1
        v_sw{1}(land_sw{1}==1) = nan;
        v_sw{2}(land_sw{2}==1) = nan;
        v_sw{3}(land_sw{3}==1) = nan;
    end
    
  	%----- Swap dimension to make it coincide with the geography
    for sw = 1:3
        lon_sw{sw} = flipud(lon_sw{sw}');
        lat_sw{sw} = flipud(lat_sw{sw}');
        v_sw{sw} = flipud(v_sw{sw}');
        h_sw{sw} = flipud(h_sw{sw}');
        inc_sw{sw} = flipud(inc_sw{sw}');
        land_sw{sw} = flipud(land_sw{sw}');
    end
    v_sw{1}(:,[108,109]) = nan; %-------------------------- Check that the problem occur in every acquisition or not.
    v_sw{2}(:,[129,130]) = nan;
    v_sw{3}(:,125) = nan;
    
    %----- 1. Detect outliers
    corr_range_v_sw = {};
    for sw = 1:3
        [row, column] = size(v_sw{sw});
        clear v_ix_nan_row v_ix_nan_column
        for r = 1:row
            v_loop = v_sw{sw}(r,:);
            MAD = nanmedian(abs(v_loop-nanmedian(v_loop)));
            v_thres = MAD*1.48*3;
            v_ix_nan_row(r,:) = abs(v_loop-nanmedian(v_loop))>v_thres;
        end
        for c = 1:column
            v_loop = v_sw{sw}(:,c);
            MAD = nanmedian(abs(v_loop-nanmedian(v_loop)));
            v_thres = MAD*1.48*3;
            v_ix_nan_column(:,c) = abs(v_loop-nanmedian(v_loop))>v_thres;
        end
        v_ix_nan = v_ix_nan_row+v_ix_nan_column;
        v_sw{sw}(v_ix_nan>0) = nan;
        
        if column<130
            v_sw{sw} = [v_sw{sw} nan(row,130-column)];  %--------------- Try to make it possible to be averaged. They will be averaged many times so I think it's still OK.
        end
        
        %----- Correct Scallopping
        if row == 233
            corr_scal_v_sw = v_sw{sw}-scal_corr_sw{sw}([1:end-1],:);
        else
            corr_scal_v_sw = v_sw{sw}-scal_corr_sw{sw};
        end
        %----- Correct Antenna Electronic Mispointing (Range Variation)
        corr_temp = range_corr_sw{sw};
        corr_temp = corr_temp';
        corr_temp = repmat(corr_temp,row,1);
        corr_range_v_sw{sw} = corr_scal_v_sw-corr_temp;
    end
    %----- Correct Platforms Navigations and attitude errors (Azimuth Variation and Median data)
   	v_sw1 = [nan(14,size(corr_range_v_sw{1},2));corr_range_v_sw{1}];
    v_sw2 = [nan(7,size(corr_range_v_sw{1},2));corr_range_v_sw{2};nan(7,size(corr_range_v_sw{1},2))];
    v_sw3 = [corr_range_v_sw{3};nan(14,size(corr_range_v_sw{1},2))];
    v_sw = [v_sw1 v_sw2 v_sw3];
    corr_v_sw = v_sw-repmat(corr_azi{p,1},1,size(v_sw,2));
    corr_v_sw = corr_v_sw-off_land_median(p,1);
    d_day = ncfile(15:22);
    
    fig_flag = 0;
    if fig_flag==1
        figure
        image(corr_v_sw,'CDataMapping','scaled')
        colorbar
        grid
        pause
        close all
    end
    save(['corr_v_'+string(d_day)+'.mat'],'corr_v_sw')
    
    corr_v{1} = v_sw([15:end],[1:size(corr_range_v_sw{1},2)]);
    corr_v{2} = v_sw([8:end-7],[size(corr_range_v_sw{1},2)+1:size(corr_range_v_sw{1},2)*2]);
    corr_v{3} = v_sw([1:end-14],[(size(corr_range_v_sw{1},2)*2+1):(size(corr_range_v_sw{1},2)*3)])
    
    for sw = 1:3  %------------ ***Need to flip back because the Wind-wave correction are ordered without flip 
        corr_v{sw} = flipud(corr_v{sw})';
    end
    v = [corr_v{1}(:);corr_v{2}(:);corr_v{3}(:)];
    
  	lon = lon(find(land==0)); 
    lat = lat(find(land==0));
    v = v(find(land==0));
    
    %--------- Correct Wind-wave velocities
    load(['VV_HH_'+string(d_day)+'.mat'])
    current_v = [lon lat v-VV_HH];
    
   	fig_flag = 1;
    if fig_flag==1
        figure
        scatter3(lon,lat,current_v(:,3),4,current_v(:,3),'filled')
        colorbar
        colormap(jet(64))
        view(2)
        pause
        close all
    end
    save(['./Current_v/current_v_'+string(d_day)+'.mat'],'current_v')
end










%% Plotting for the presentation slide 
figure
image(corr_range_sw,'CDataMapping','scaled')
colorbar
grid
%----- Plotting Scalloping 
figure
plot([1:234],scal_corr_sw1,'b-','LineWidth',3)
grid
hold on
plot([1:234],scal_corr_sw2,'r-','LineWidth',3)
plot([1:234],scal_corr_sw3,'g-','LineWidth',3)
axis([0 234 -0.15 0.15])
ylabel('RVL scalloping correction (m/s)')
xlabel('azimuth direction (pix)')
set(gca,'fontsize', 20)

%------- Plotting Antenna electronic mispointing correction
x = [1:size(azi_median_med{sw},2)]';
y = azi_median_med{sw}';
ix_nan = find(isnan(y));
x(ix_nan) = [];
y(ix_nan) = [];
A = [x ones(size(x,1),1)];
m(:,sw) = lscov(A,y);
figure
plot(x,y,'bx')
hold on 
grid
plot(x,A*m(:,sw),'ro-')
sw = 2;
x = [1:size(azi_median_med{sw},2)]';
y = azi_median_med{sw}';
ix_nan = find(isnan(y));
x(ix_nan) = [];
y(ix_nan) = [];
A = [x ones(size(x,1),1)];
m(:,sw) = lscov(A,y);
plot(x+130,y,'bx')
plot(x+130,A*m(:,sw),'ro-')
sw = 3;
x = [1:size(azi_median_med{sw},2)]';
y = azi_median_med{sw}';
ix_nan = find(isnan(y));
x(ix_nan) = [];
y(ix_nan) = [];
A = [x ones(size(x,1),1)];
m(:,sw) = lscov(A,y);
plot(x+130+130,y,'bx')
plot(x+130+130,A*m(:,sw),'ro-')
ylabel('RVL (m/s)')
xlabel('range direction (pix)')
set(gca,'fontsize', 20)

%------- Plotting Platform navigation and attitude errors 
figure
hold on
grid
for p = 1:size(corr_azi,1)
    plot([1:size(corr_azi{p},1)],corr_azi{p},'-','LineWidth',3)
end
box on
ylabel('Land remainging median RVL (m/s)')
xlabel('azimuth direction (pix)')
set(gca,'fontsize', 20)

%------- Plotting Platform navigation and attitude errors 
figure
hold on
grid
plot([1:size(off_land_median,1)],off_land_median(:,1),'bo', 'MarkerFaceColor', 'b','MarkerSize',8)
box on
ylabel('Land median RVL (m/s)')
xlabel('Image number')
set(gca,'fontsize', 20)










