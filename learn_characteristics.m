%% Script that generates features and classifier for all training patients
%

%list of training patients
patients = {'DK1', 'DK2','DK3','DK4','DK8', 'DLBCL1','DLBCL2','DLBCL4','DLBCL5','MSSM_NHL1'};

% How many samples do you want per patient?
Nsamp = 1000;

characteristics = zeros([Nsamp*size(patients,2), 7]);
for pat = 1
    
    mydir = 'E:\IGTResearch\CS766 Project\Train\';
    ct = loadAmOrNrrd([mydir, patients{pat}, '_ct_resamp.nrrd']);
    pet = loadAmOrNrrd([mydir, patients{pat}, '_pet.nrrd']);
    liverMask = loadAmOrNrrd([mydir, patients{pat}, '_liver_resamp.nrrd']);
    
    PETdata = pet.data;
    PETdata = double(PETdata);
    CTdata = ct.data;
    CTdata = double(CTdata);
    liverMask = liverMask.data;
    liverMask = double(liverMask);
    
    PETdata = imgaussfilt3(PETdata,1);
    CTdata = imgaussfilt3(CTdata,1);
    
    %Crop image around body (using PET image because CT has bed in it)
    mipxz = max(PETdata,[],1);
    mipyz = max(PETdata,[],2);
    
    mipxz = squeeze(mipxz);
    mipyz = squeeze(mipyz);
    
    [xs,zs1] = find(mipxz > 1);
    [ys,zs2] = find(mipyz > 1);
    
    CTcropped = CTdata(min(ys):max(ys),min(xs):max(xs),min(zs1):max(zs1));
    PETcropped = PETdata(min(ys):max(ys),min(xs):max(xs),min(zs1):max(zs1));
    liverMask = liverMask(min(ys):max(ys),min(xs):max(xs),min(zs1):max(zs1));
    
 
    keys = detectSift3D(CTcropped);
    [desc, coords] = extractSift3D(keys);
    
    %find how many coordinates are in the liver
    liverCoords = [];
    nonLiverCoords = [];
    for pt = 1:size(coords,1)
        if liverMask(coords(pt,1), coords(pt,2), coords(pt, 3)) == 1
            liverCoords(size(liverCoords,1)+1,:) = [coords(pt,1), coords(pt,2), coords(pt, 3)];
        else
            nonLiverCoords(size(nonLiverCoords,1)+1,:) = [coords(pt,1), coords(pt,2), coords(pt, 3)];
        end
    end
    
    % if number of points in liver is above 150, sample 150, else sample
    % all points
    if size(liverCoords,1) > 200
        Nliver = size(liverCoords,1);
        Nnonliver = size(nonLiverCoords,1);
        sampleLiverCoords = liverCoords(randsample(1:Nliver, 200),:);
        sampleNonLiverCoords = nonLiverCoords(randsample(1:Nnonliver, Nsamp-200),:);
    else
        Nliver = size(liverCoords,1);
        Nnonliver = size(nonLiverCoords,1);
        sampleLiverCoords = liverCoords;
        sampleNonLiverCoords = nonLiverCoords(randsample(1:Nnonliver, Nsamp-Nliver),:);
    end
    allSamples = [sampleLiverCoords; sampleNonLiverCoords];
    numneighbs = 7;
    numSlices = 3;
    for i=1:Nsamp
        coord = allSamples(i, :);
        xregionOI = coord(1)-numneighbs:coord(1)+numneighbs;
        yregionOI = coord(2)-numneighbs:coord(2)+numneighbs;
        zregionOI = coord(3)-numSlices:coord(3)+numSlices;
        
        xregionOI(xregionOI < 1) = 1; xregionOI(xregionOI > size(CTcropped,1)) = size(CTcropped,1);
        yregionOI(yregionOI < 1) = 1; yregionOI(yregionOI > size(CTcropped,2)) = size(CTcropped,2);
        zregionOI(zregionOI < 1) = 1; zregionOI(zregionOI > size(CTcropped,3)) = size(CTcropped,3);

        needCT = CTcropped(xregionOI, yregionOI, zregionOI);


        %find parameters of guess
        inputrow = Nsamp*(pat-1)+i;
        characteristics(inputrow, 1) = median(needCT(:));
        characteristics(inputrow, 2) = mad(needCT(:),1);
        characteristics(inputrow, 3) = max(needCT(:));
        characteristics(inputrow, 4) = min(needCT(:));
        characteristics(inputrow, 5) = entropy(needCT);
        
        characteristics(inputrow, 6) = median(needPET(:));
        characteristics(inputrow, 7) = mad(needPET(:),1);
        characteristics(inputrow, 8) = max(needPET(:));
        characteristics(inputrow, 9) = min(needPET(:));
        characteristics(inputrow, 10) = entropy(needPET);


        if liverMask(coord(1), coord(2), coord(3)) > 0
            characteristics(inputrow, 11) = 1;
        else
            characteristics(inputrow, 11) = 0;
        end
    end
end


save('E:\IGTResearch\CS766 Project\Train\training_data_7neighbs_3slices_PETres', 'characteristics')
