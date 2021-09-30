function [response, timestamp]= getResponse(timeout, return_after_response)
%     KbEventFlush
    response = 0;

    % start stopwatch
    tic
    while toc < timeout
        % poll the state of the keyboard
        [keyisdown, when, keyCode] = KbCheck;

        % if there wasn't a response before and there is a
        % keyboard press available
        if response == 0 && keyisdown
            timestamp = when;% - startTime;
            response = find(keyCode);
            if return_after_response
                break;
            end
        end
    end

        % if no response yet after timeout
        if (response == 0)
            timestamp = GetSecs;
        end
end