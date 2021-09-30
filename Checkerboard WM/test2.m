Contrasts_matrix_stim1 = cell(1,length(params.Blocks));
Contrasts_matrix_stim2 = cell(1,length(params.Blocks));


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
                Contrasts_matrix_stim1{bl3}{stim} = contrast_full1;
                
                while 1 % make sure that contrast of stim1(A) ~= stim2(A).
                    contrast_full2 = contrast_full1(randperm(length(contrast_full1)));
                    if contrast_full2 ~= contrast_full1
                        Contrasts_matrix_stim2{bl3}{stim} = contrast_full2; 
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
                Contrasts_matrix_stim2{bl3}{stim} = contrast_full2;
                
                while 1 % make sure that contrast of stim1(A) ~= stim2(A).
                    contrast_full1 = contrast_full2(randperm(length(contrast_full2)));
                    if contrast_full1 ~= contrast_full2
                        Contrasts_matrix_stim1{bl3}{stim} = contrast_full1; 
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
            Contrasts_matrix_stim1{bl3}{stim} = contrast_full1;
        end
        block2 = block2+1; 
    end
end
