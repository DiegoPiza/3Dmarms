function chan=unit2ch(unitID,subject)
%outputs channel number from unit ID in character form and subject ID ie.
%'C' or 'PB'

neuron=unitID;

if subject=='C'
    if (neuron(12))=='_'   
        chan=(neuron(11));
    else
        chan=(neuron(11:12));
    end
    
    elseif subject=='PB'
    if (neuron(13))=='_'   
        chan=(neuron(12));
    else
        chan=(neuron(12:13));
    end
   
end
chan=double(string(chan));
end
