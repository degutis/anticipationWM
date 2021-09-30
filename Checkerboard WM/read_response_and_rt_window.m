function [response, rt] = read_response_and_rt_window(responseWindow)

startTime  = GetSecs;

% Create a cell array containing all the key presses we are expecting - any
% key press that is not stored in the array will be ignored
expectedKeys = {'a','l','space'};


% response empty to start
response = [];

% While we haven't received a response...
while isempty(response)
    
    % If we have been going on for too long then break
    if GetSecs > startTime + responseWindow
        response = 'm'; %m for miss
        break
    end

    % Get the response and timing of th response
%     [~, secs, keyCode] = KbCheck;
    [keyCode,secs]= getResponse(0.1,true);

    
    % Translate keyCode into the key pressed
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
       
end

% Calculate reaction time
rt = secs - startTime;

end