    close=0; 
    contrast_number = 2;
    startTime = GetSecs;
    time.answer = 5;
    Probe_stim_texture = [36 37 38 39];
    response = [];
    expectedKeys = {'a','l','space'};
    ListenChar(2)
    KbQueueCreate;
    KbQueueStart;

    while strcmp(response,'m')==0 
        response = [];
        Probe_stim_texture(contrast_number)
        if GetSecs > startTime + time.answer
            response = 'm'; %m for miss
            break
        else
            response = [];
            [~,keyCode]=KbQueueCheck;            
            if keyCode == 0
                continue
            else
                response = KbName(keyCode);
                % Check if the key pressed was one we expected
                rep = find(strcmp(response,expectedKeys) == 1, 1);
                % If the response was not one of the keys we were expecting, then make
                % response empty again and wait for another response. 
                if isempty(rep)
                    response = [];
                end
            end
            
            if response == 'a'
                if contrast_number == 1
                    contrast_number = length(Probe_stim_texture);
                else
                    contrast_number = contrast_number-1;
                end
            elseif response == 'l'
                if contrast_number == length(Probe_stim_texture)
                    contrast_number = 1;
                else
                    contrast_number = contrast_number+1;
                end
            elseif strcmp(response,'space')
                close=1;
                break
            end
            
%             rt = secs - startTime;
        end
        

%         Screen('DrawTexture', window, Probe_stim_texture(contrast_number)); 
%         Screen('Flip', window);       
    end
        ListenChar(0)

    
    %     expectedTime = expectedTime+time.answer;