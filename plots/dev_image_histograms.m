
dat = niftiRead('161223_M_scan7_FULLDATA.nii.gz');

b0 = dat.data(:,:,:,1);

d1 = dat.data(:,:,:,1);
d2 = dat.data(:,:,:,2);
d3 = dat.data(:,:,:,3);



figure; hold on
[ x0, y0 ] = hist(b0(:), 100); 
[ x1, y1 ] = hist(d1(:), 100); 
[ x2, y2 ] = hist(d2(:), 100); 
[ x3, y3 ] = hist(d3(:), 100);

plot(y0, x0, 'color', [0 0 0]);
plot(y1, x1, 'color', [1 0 0]);
plot(y2, x2, 'color', [0 1 0]);
plot(y3, x3, 'color', [0 0 1]);
xlim([0 2000]);


figure; hold on;
for ii = 1:size(dat.data, 4)
    
    if ii == 1
        b0 = dat.data(:,:,:,ii);
        [ x, y ] = hist(b0(:), 100);
        plot(y, x, 'color', [0 0 0]);
        continue
    end
    
    vol = dat.data(:,:,:,ii);
    [ x, y ] = hist(vol(:), 100);
    plot(y, x, 'color', [0.9, 0.1 0]);
    
end
xlim([0 2000]);




figure; 

subplot(2, 2, 1); 
hist(log10(b0(:)), 100); 

subplot(2, 2, 2); 
hist(log10(d1(:)), 100); 

subplot(2, 2, 3); 
hist(log10(d2(:)), 100); 

subplot(2, 2, 4); 
hist(log10(d3(:)), 100);





