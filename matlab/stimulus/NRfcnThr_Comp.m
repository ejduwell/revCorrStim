function parVal = NRfcnThr_Comp(perc_thr, fcn)
    % This function takes an input Naka Rushton function (fcn) obtained from
    % QUEST procedure model fitting, and solves it for an arbitrary %
    % accuracy threshold (perc_thr).
    
    % get the matlab release runnung on this machine
    [mlb_v, mlb_d] = version;
    expression = '(R[0-9][0-9][0-9][0-9].)';
    mlb_version = regexp(mlb_v,expression,'match');
    mlb_year=str2double(mlb_version{1,1}(2:end-1));
    
    % assign function input to f
    f = fcn;
    
    % set up symbolic math variables
    if mlb_year < 2020
        syms cf real
    else
    syms cf [1 1] real
    end
    
    % Solve for the psychophysical parameter value at the specified
    % percent accuracy threshold
    parVal = double(max(solve(f.nr_rmax*((cf^f.nr_n)/(f.nr_c50^f.nr_n + cf^f.nr_n)) + f.nr_b == perc_thr,cf)));

end