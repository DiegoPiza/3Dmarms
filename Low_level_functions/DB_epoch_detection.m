function epochs = DB_epoch_detection(data, threshold, sign)
    % Detects epochs in data based on given threshold and sign
    
    % Input:
    % data - Input data
    % threshold - Threshold value for epoch detection
    % sign - Sign flag (true for positive threshold, false for negative threshold)
    
    % Output:
    % epochs - Structure array containing detected epochs
    
    if sign
        epoch = data > threshold;
    else
        epoch = data < threshold;
    end
    
    % Initialize variables to store epoch information
    starttime = [];
    endtime = [];
    epochlength = [];
    epochcount = 0;
    epochdetect = true;
    ii = 1;

    while epochdetect
        if ii == length(data)
            epochdetect = false;
        end
        if epoch(ii)
            epochcount = epochcount + 1;
            epochs(epochcount).start_time = ii;
            isepoch = true;
            while isepoch
                if epoch(ii) || epoch(ii + 1)
                    if ii == length(data)
                        isepoch = false;
                    end
                    isepoch = true;
                    ii = ii + 1;
                else
                    isepoch = false;
                end
            end
            epochs(epochcount).end_time = ii;
            epochs(epochcount).epoch_length = epochs(epochcount).end_time - epochs(epochcount).start_time;

        else
            ii = ii + 1;
        end
    end
end
