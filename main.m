%% 
% This code is the main code for locating and segmenting the liver 
%
%
addpath('SIFT_matlab')
%% Run SVM model on given dataset
mydir = 'C:\Users\aweisman\Box Sync\Amy_Documents\2ndYear_Courses\CS766\Data\';
characteristics = load([mydir, '\training_data_7neighbs_3slices_PETres']);
X = characteristics.characteristics;
training = X(:, 1:10);
group = X(:,11);

%learning from CT using PET resolution, 7 neighbors, 3 slices
SVMModel = fitcsvm(training, group, 'Standardize', 'on', 'KernelFunction', 'polynomial', 'PolynomialOrder', 2, 'BoxConstraint', 0.36, 'Cost', [0 10; 1 0]);

CVSVMModel = crossval(SVMModel);
kfoldLoss(CVSVMModel)

%%
%list of train patients and test patients, only DK5 and DK6 used for
%example patients
trainpatients = {'DK1', 'DK2','DK3','DK4','DK8', 'DLBCL1','DLBCL2','DLBCL4','DLBCL5','MSSM_NHL1'};
testpatients = {'MSSM4', 'MSSM9', 'MSSM14', 'DK5', 'DK6', 'DK7', 'DLBCL3', 'MSSM_NHL3', 'MSSM1', 'MSSM2'};

%If you're interested in calculating sensitivity/specificity/positive
%predictive value and negative predictive value
all_sens = zeros([size(testpatients,2),1]);
all_spec = zeros([size(testpatients,2),1]);
all_ppv = zeros([size(testpatients,2),1]);
all_npv = zeros([size(testpatients,2),1]);
pt_inliver = zeros([size(testpatients,2),1]);

%run through all patients
for pat = 4:5
    
    %Load ct, pet, and liver mask (for testing accuracy)
    ct = loadAmOrNrrd([mydir, testpatients{pat}, '_ct.nrrd']);
    pet = loadAmOrNrrd([mydir, testpatients{pat}, '_pet.nrrd']);
    liver = loadAmOrNrrd([mydir, testpatients{pat}, '_liverMask.nrrd']);
    
    PETdata = pet.data;
    CTdata = ct.data;
    liverMask = liver.data;
    
    PETdata = double(PETdata);
    CTdata = double(CTdata);
    liverMask = double(liverMask);
    
    %Smooth data with gaussian filter
    PETdata = imgaussfilt3(PETdata,1);
    CTdata = imgaussfilt3(CTdata,1);
    
    % Run code to find the liver and non liver points of interest
    [liverCoords, nonliverCoords] = FindLiver(PETdata,CTdata, SVMModel);
    
    %Calculate sens, spec, ppv, and npv
    allpts = size(liverCoords,1) + size(nonliverCoords,1);
    tp = 0; %true positives
    fp = 0; % false positives
    tn = 0; %true negatives
    fn = 0; %false negatives
    for pt = 1:size(liverCoords,1)
        ptOI = [liverCoords(pt,1), liverCoords(pt,2), liverCoords(pt,3)];
        if liverMask(ptOI(1), ptOI(2), ptOI(3)) > 0
            tp = tp+1;
        else
            fp = fp+1;
        end
    end
    for pt = 1:size(nonliverCoords,1)
        ptOI = [nonliverCoords(pt,1), nonliverCoords(pt,2), nonliverCoords(pt,3)];
        if liverMask(ptOI(1), ptOI(2), ptOI(3)) == 0
            tn = tn+1;
        else
            fn = fn+1;
        end
    end
    sens = tp/(tp+fn);
    spec = tn/(tn+fp);
    ppv = tp/(tp+fp);
    npv = tn/(tn+fn);
    
    all_sens(pat) = sens; all_spec(pat) = spec; all_ppv(pat) = ppv; all_npv(pat) = npv;
    
    % Calculate median of x,y, and z coordinates from findLiver
    z_coord = round(median(liverCoords(:,3)));
    row_coord = round(median(liverCoords(:,1)));
    col_coord = round(median(liverCoords(:,2)));
    
    %Test if the point is in the liver?
    pt_inliver(pat) = liverMask(row_coord, col_coord, z_coord);
    
    if pt_inliver(pat) == 1
        %if the point is in the liver, segment liver and save mask 
        newMask = liverSegment(CTdata, PETdata, row_coord, col_coord, z_coord);
        newMask = int16(newMask);
        nrrdWriter([mydir, testpatients{pat}, '_liverSegmented.nrrd'], newMask, pet.voxel_size, ct.start, 'raw');
    
        %Generate liver sphere and calculate mean in liver
        ind = find(newMask == 1);
        [r,c,s] = ind2sub(size(newMask), ind);
        mid_s = median(s);
        mid_r = median(r); 
        side_c = quantile(c, 0.75);
        
        diam_cm = 3;
        voxelsizes = ct.voxel_size;
        xy_size = voxelsizes(1);
        z_size = voxelsizes(3);
        xy_rad = (diam_cm/xy_size)/2;
        z_rad = (diam_cm/z_size)/2;
        
        
        [rr, cc, zz] = meshgrid(1:size(CTdata,1), 1:size(CTdata,2), 1:size(CTdata,3));
        S = sqrt( ((rr-side_c).^2 + (cc-mid_r).^2)./(xy_rad^2) + ((zz-mid_s).^2./(z_rad^2)))<=1;
        ind_sphere = find(S>0);
        
       
        liver_sphere = zeros(size(CTdata));
        liver_sphere(ind_sphere)=1;
        liver_sphere = int16(liver_sphere);
        nrrdWriter([mydir, testpatients{pat}, '_liverSphere.nrrd'], liver_sphere, pet.voxel_size, ct.start, 'raw');
        
        liverMean = mean(PETdata(ind_sphere));
 
    else 
        msg = 'Point is not located in the liver';
        error(msg);
    end
    
end


