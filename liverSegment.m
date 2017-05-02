%% This function segments the liver using a region growing algorithm
% 
% Inputs are CT data, PET data, and x,y,z points of initial guess
%
% Returns a binary mask of final segmentation the size of the CT data
%%

function mask = liverSegment(CTdata, PETdata, row_coord, col_coord, z_coord)

% Set threshold parameters
CTgradThresh = 0.015; 

% preallocate array
mask = zeros(size(CTdata));

%Make an initial guess about PET and CT area using a small window 
PEToi = PETdata(row_coord-5:row_coord+5, col_coord-5:col_coord+5, z_coord-3:z_coord+3);
PETmax = max(PEToi(:));
PETmin = min(PEToi(:));

CToi = CTdata(row_coord-5:row_coord+5, col_coord-5:col_coord+5, z_coord-3:z_coord+3);
CTmax = max(CToi(:));
CTmin = min(CToi(:));

if CTmin>0
    fac = 1.5;
else
    fac = 0.5;
end

%Create x-, y-, and z-direction gradients for CT
% (not done for PET because the image quality is so bad)
[CTxgradient_img, CTygradient_img, CTzgradient_img] = imgradientxyz(CTdata);


% normalize gradients and take absolute value
CTxgradient_img = abs(CTxgradient_img);
CTxgradient_img = CTxgradient_img./max(CTxgradient_img(:));

CTygradient_img = abs(CTygradient_img);
CTygradient_img = CTygradient_img./max(CTygradient_img(:));

CTzgradient_img = abs(CTzgradient_img);
CTzgradient_img = CTzgradient_img./max(CTzgradient_img(:));


% add the initial pixel to the queue
guesses = [row_coord, col_coord, z_coord];
[nRow, nCol, nSli] = size(CTdata);

%%% START OF REGION GROWING ALGORITHM
while size(guesses, 1)
  % the first queue position determines the new values
  xnew = guesses(1,1);
  ynew = guesses(1,2);
  znew = guesses(1,3);
  
  % delete the first queue position
  guesses(1,:) = [];
    
  % check the neighbors for the current position
  for xx = -1:1
    for yy = -1:1
      for zz = -1:1
            
        if xnew+xx > 0  &&  xnew+xx <= nRow &&...          % within the x-bounds?
           ynew+yy > 0  &&  ynew+yy <= nCol &&...          % within the y-bounds?          
           znew+zz > 0  &&  znew+zz <= nSli &&...          % within the z-bounds?
           any([xx, yy, zz])       &&...      % i/j/k of (0/0/0) is redundant
           ~mask(xnew+xx, ynew+yy, znew+zz) &&...          % pixelposition already set?
           CTxgradient_img(xnew+xx, ynew+yy, znew+zz) < CTgradThresh && ... 
           CTygradient_img(xnew+xx, ynew+yy, znew+zz) < CTgradThresh && ...
           CTzgradient_img(xnew+xx, ynew+yy, znew+zz) < CTgradThresh && ...
           CTdata(xnew+xx, ynew+yy, znew+zz) < 2.5*CTmax && ... 
           CTdata(xnew+xx, ynew+yy, znew+zz) > fac*CTmin && ...
           PETdata(xnew+xx, ynew+yy, znew+zz) < 2.5*PETmax && ...
           PETdata(xnew+xx, ynew+yy, znew+zz) > 0.5*PETmin
       
           mask(xnew+xx, ynew+yy, znew+zz) = 1; 

           % add the current pixel to the computation queue (recursive)
           guesses(end+1,:) = [xnew+xx, ynew+yy, znew+zz];
    
        end        
      end
    end  
  end
end

%smooth edges using imerode and imdilate
se = strel('sphere', 3);
mask = imerode(mask, se);
mask = imdilate(mask, se);

% Get rid of small islands and fill holes
for slice = 1:size(CTdata,3)
    maskSlice = mask(:,:,slice);
    maskSlice = bwareaopen(maskSlice, 150);
    maskSlice = imfill(maskSlice, 'holes');
    mask(:,:,slice) = maskSlice;
end

% find liver component if there were broken pieces
maskLabelled = bwlabeln(mask);
liverLabel = maskLabelled(row_coord, col_coord, z_coord);
ind = maskLabelled == liverLabel;
finalMask = zeros(size(mask));
finalMask(ind) = 1;

mask = finalMask;

