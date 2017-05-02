function [liverCoords, nonliverCoords] = FindLiver(PETdata, CTdata, Model)
% Function that finds a point within the liver after cropping to be ran
% through region growing algorithm

%% 
%Crop image around body (using PET image because CT has bed in it)
mipxz = max(PETdata,[],1);
mipyz = max(PETdata,[],2);

mipxz = squeeze(mipxz);
mipyz = squeeze(mipyz);

[xs,zs1] = find(mipxz > 1);
[ys,zs2] = find(mipyz > 1);

CTcropped = CTdata(min(ys):max(ys),min(xs):max(xs),min(zs1):max(zs1));
PETcropped = PETdata(min(ys):max(ys),min(xs):max(xs),min(zs1):max(zs1));

%Calculate gradient in x- and y- directions
xgradient_img = zeros(size(CTcropped));
ygradient_img = zeros(size(CTcropped));
for sli = 1:size(CTcropped,3)
    [xgradient_img(:,:,sli), ygradient_img(:,:,sli)] = imgradientxy(CTcropped(:,:,sli));
end

%Run SIFT on cropped CT and extract points of interest
keys = detectSift3D(CTcropped);
[desc, coords] = extractSift3D(keys);


%specify patch size (number of neighbors in x,y direction and number of
%slices in z-direction
numneighbs = 7;
numSlice = 3;

%iterate through points of interest, calculate feature in patch, run
%through SVM model

for i = 1:size(coords,1)
    coord = coords(i,:);
    
    %generate patch around point of interest
    xregionOI = coord(1)-numneighbs:coord(1)+numneighbs;
    yregionOI = coord(2)-numneighbs:coord(2)+numneighbs;
    zregionOI = coord(3)-numSlice:coord(3)+numSlice;
    
    xregionOI(xregionOI < 1) = 1; xregionOI(xregionOI > size(CTcropped,1)) = size(CTcropped,1);
    yregionOI(yregionOI < 1) = 1; yregionOI(yregionOI > size(CTcropped,2)) = size(CTcropped,2);
    zregionOI(zregionOI < 1) = 1; zregionOI(zregionOI > size(CTcropped,3)) = size(CTcropped,3);
    
    
    %Calculate features
    needCT = CTcropped(xregionOI, yregionOI, zregionOI);
    needPET = PETcropped(xregionOI, yregionOI, zregionOI);
    xgradOI = xgradient_img(xregionOI, yregionOI, zregionOI);
    ygradOI = ygradient_img(xregionOI, yregionOI, zregionOI);
    Xnew = [median(needCT(:)),mad(needCT(:),1),max(needCT(:)),min(needCT(:)),entropy(needCT),...
          median(needPET(:)),mad(needPET(:),1),max(needPET(:)),min(needPET(:)),entropy(needPET)];
    
    %Calculate label using SVM
    [newLabel,score] = predict(Model, Xnew);

    coords(i,4) = newLabel;
    coords(i,5:6) = score;
end

% add back sections that were cropped
ind = coords(:,4) == 1;
liverCoords = coords(ind,:);
liverCoords(:,1) = liverCoords(:,1) + min(ys);
liverCoords(:,2) = liverCoords(:,2) + min(xs);
liverCoords(:,3) = liverCoords(:,3) + min(zs1);

ind2 = coords(:,4) == 0;
nonliverCoords = coords(ind2,:);
nonliverCoords(:,1) = nonliverCoords(:,1) + min(ys);
nonliverCoords(:,2) = nonliverCoords(:,2) + min(xs);
nonliverCoords(:,3) = nonliverCoords(:,3) + min(zs1);

% Plot the results
fg = figure();
sliceOI = round(mean(liverCoords(:,1)));
imshow(squeeze(CTdata(sliceOI + min(zs1), :, :)), [0,200]);
hold on
plot(liverCoords(:,3), liverCoords(:,2), 'g.', 'MarkerSize', 10);
plot(nonliverCoords(:,3), nonliverCoords(:,2), 'r.', 'MarkerSize', 10);