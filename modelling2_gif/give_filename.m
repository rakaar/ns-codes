function fname = give_filename(c_num,off_on,training,stim,stage)
if c_num == 10
    c_num = '10_12';
elseif c_num == 9
    c_num = '9_13';
end

if stage == 'i'
    stage = 'initial';
elseif stage == 'f'
    stage = 'final';
end
    fname = strcat('c21_',c_num,'_',off_on,'_',training,'_trained_',stim,'_',stage,'.mat');
end