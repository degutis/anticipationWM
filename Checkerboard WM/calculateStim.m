function [params, oneStim, retro] = calculateStim(params)


% 4 parts of a stimulus
% 12 levels of contrast
% same total level of luminance in each stim

contrasts_stim = [0:1/(length(params.contrasts)-1):1];

contrasts_stim_matrix = [contrasts_stim(1:3:end);contrasts_stim(2:3:end);contrasts_stim(3:3:end)];

%things that should be equal in one matrix:
numCues = 1:2;
numLocations = 1:params.numSegments;
numContrastLevels = contrasts_stim_matrix(1,:); 
numContrastBins = 3;
retro.singleBlock = [repelem(numCues,length(numLocations)*length(numContrastLevels));repmat(numLocations,1,length(numCues)*length(numContrastLevels));repelem(numContrastLevels,length(numCues)*length(numLocations))];

iA = 1:numContrastBins:length(retro.singleBlock);
iB = 2:numContrastBins:length(retro.singleBlock);
iC = 3:numContrastBins:length(retro.singleBlock);

% each block 

retro.Block{1} = retro.singleBlock;
retro.Block{2} = retro.singleBlock;
retro.Block{3} = retro.singleBlock;

mixBlocks = randperm(params.numBlocksPer);

retro.Block{mixBlocks(1)}(3,iB) =retro.Block{mixBlocks(1)}(3,iB)+contrasts_stim(2);
retro.Block{mixBlocks(1)}(3,iC) =retro.Block{mixBlocks(1)}(3,iC)+contrasts_stim(2)*2;

retro.Block{mixBlocks(2)}(3,iA) =retro.Block{mixBlocks(2)}(3,iA)+contrasts_stim(2)*2;
retro.Block{mixBlocks(2)}(3,iC) =retro.Block{mixBlocks(2)}(3,iC)+contrasts_stim(2);

retro.Block{mixBlocks(3)}(3,iA) =retro.Block{mixBlocks(3)}(3,iA)+contrasts_stim(2);
retro.Block{mixBlocks(3)}(3,iB) =retro.Block{mixBlocks(3)}(3,iB)+contrasts_stim(2)*2; 

retro.Block{mixBlocks(1)} = retro.Block{mixBlocks(1)}(:,randperm(size(retro.Block{mixBlocks(1)},2)));
retro.Block{mixBlocks(2)} = retro.Block{mixBlocks(2)}(:,randperm(size(retro.Block{mixBlocks(2)},2)));
retro.Block{mixBlocks(3)} = retro.Block{mixBlocks(3)}(:,randperm(size(retro.Block{mixBlocks(3)},2)));

% again generates blocks with 32 trials
oneStim.singleBlock =  [repmat(numLocations,1,length(numContrastLevels)*2);repelem(numContrastLevels,length(numLocations)*2)];

iA = 1:numContrastBins:length(oneStim.singleBlock);
iB = 2:numContrastBins:length(oneStim.singleBlock);
iC = 3:numContrastBins:length(oneStim.singleBlock);

oneStim.Block{1} = oneStim.singleBlock;
oneStim.Block{2} = oneStim.singleBlock;
oneStim.Block{3} = oneStim.singleBlock;

mixBlocks = randperm(params.numBlocksPer);

oneStim.Block{mixBlocks(1)}(2,iB) =oneStim.Block{mixBlocks(1)}(2,iB)+contrasts_stim(2);
oneStim.Block{mixBlocks(1)}(2,iC) =oneStim.Block{mixBlocks(1)}(2,iC)+contrasts_stim(2)*2;

oneStim.Block{mixBlocks(2)}(2,iA) =oneStim.Block{mixBlocks(2)}(2,iA)+contrasts_stim(2)*2;
oneStim.Block{mixBlocks(2)}(2,iC) =oneStim.Block{mixBlocks(2)}(2,iC)+contrasts_stim(2);

oneStim.Block{mixBlocks(3)}(2,iA) =oneStim.Block{mixBlocks(3)}(2,iA)+contrasts_stim(2);
oneStim.Block{mixBlocks(3)}(2,iB) =oneStim.Block{mixBlocks(3)}(2,iB)+contrasts_stim(2)*2; 

oneStim.Block{mixBlocks(1)} = oneStim.Block{mixBlocks(1)}(:,randperm(size(oneStim.Block{mixBlocks(1)},2)));
oneStim.Block{mixBlocks(2)} = oneStim.Block{mixBlocks(2)}(:,randperm(size(oneStim.Block{mixBlocks(2)},2)));
oneStim.Block{mixBlocks(3)} = oneStim.Block{mixBlocks(3)}(:,randperm(size(oneStim.Block{mixBlocks(3)},2)));

params.cue = cell(1,length(params.Blocks));
block=1;
for i = 1:length(params.Blocks)
    if params.Blocks(i)==1
        params.cue{i} = retro.Block{block}(1,:); 
        block=block+1;
    end
end

params.phase = {}; 

for bl0 = 1:length(params.Blocks)
    if params.Blocks(bl0) == 1
        params.phase{bl0} = repmat([1;-1],1,params.numTrials);
        params.phase{bl0} = Shuffle(params.phase{bl0});

    elseif params.Blocks(bl0) == 2
        params.phase{bl0} = repmat([1,-1],1,params.numTrials/2);
        params.phase{bl0} = params.phase{bl0}(randperm(length(params.phase{bl0}))); 
    end
end

params.probeFull = {};
block1 = 1;
block2 = 1;
for bl = 1:length(params.Blocks)
    if params.Blocks(bl)==1
        params.probeFull{bl} = retro.Block{block1}(2,:);
        block1 = block1+1;
        
    elseif params.Blocks(bl)==2     
        params.probeFull{bl} = retro.Block{block2}(2,:);
        block2 = block2+1; 
    end
end
       
params.Contrasts_matrix_stim1 = cell(1,length(params.Blocks));
params.Contrasts_matrix_stim2 = cell(1,length(params.Blocks));


block1 = 1;
block2 = 1;

for bl3 = 1:length(params.Blocks)     
    if params.Blocks(bl3)==1
        for stim = 1:params.numTrials
            stimulus = retro.Block{block1}(:,stim);
            if stimulus(1) == 1 % retrocue either 1 or 2
                contrast_select = sum(sum((round(stimulus(3),3)==round(contrasts_stim_matrix,3)),2).*contrasts_stim_matrix);
                contrast_full1 = NaN(1,params.numSegments);
                contrast_full1(stimulus(2)) = stimulus(3);
                contrast_left = contrast_select(round(contrast_select,3)~=round(stimulus(3),3));
                contrast_left = contrast_left(randperm(length(contrast_left)));
                contrast_full1(isnan(contrast_full1)) = contrast_left; 
                params.Contrasts_matrix_stim1{bl3}{stim} = contrast_full1;
                
                while 1 % make sure that contrast of stim1(A) ~= stim2(A).
                    contrast_full2 = contrast_full1(randperm(length(contrast_full1)));
                    if contrast_full2 ~= contrast_full1
                        params.Contrasts_matrix_stim2{bl3}{stim} = contrast_full2; 
                        break
                    end
                end
                
            elseif stimulus(1) ==2
                contrast_select = sum(sum((round(stimulus(3),3)==round(contrasts_stim_matrix,3)),2).*contrasts_stim_matrix);
                contrast_full2 = NaN(1,params.numSegments);
                contrast_full2(stimulus(2)) = stimulus(3);
                contrast_left = contrast_select(round(contrast_select,3)~=round(stimulus(3),3));
                contrast_left = contrast_left(randperm(length(contrast_left)));
                contrast_full2(isnan(contrast_full2)) = contrast_left; 
                params.Contrasts_matrix_stim2{bl3}{stim} = contrast_full2;
                
                while 1 % make sure that contrast of stim1(A) ~= stim2(A).
                    contrast_full1 = contrast_full2(randperm(length(contrast_full2)));
                    if contrast_full1 ~= contrast_full2
                        params.Contrasts_matrix_stim1{bl3}{stim} = contrast_full1; 
                        break
                    end
                end
            end
        end
        block1 = block1+1;  
                
    elseif params.Blocks(bl3)==2
        for stim = 1:params.numTrials
            stimulus = retro.Block{block2}(:,stim);
            contrast_select = sum(sum((round(stimulus(3),2)==round(contrasts_stim_matrix,2)),2).*contrasts_stim_matrix);
            contrast_full1 = NaN(1,params.numSegments);
            contrast_full1(stimulus(2)) = stimulus(3);
            contrast_left = contrast_select(round(contrast_select,3)~=round(stimulus(3),3));
            contrast_left = contrast_left(randperm(length(contrast_left)));
            contrast_full1(isnan(contrast_full1)) = contrast_left; 
            params.Contrasts_matrix_stim1{bl3}{stim} = contrast_full1;
        end
        block2 = block2+1; 
    end
end
end
            
