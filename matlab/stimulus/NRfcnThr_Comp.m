function parVal = NRfcnThr_Comp(perc_thr, fcn)
    % This function takes an input Naka Rushton function (fcn) obtained from
    % QUEST procedure model fitting, and solves it for an arbitrary %
    % accuracy threshold (perc_thr).

    % assign function input to f
    f = fcn;
    
    % set up symbolic math variables
    syms cf [1 1] real
    
    % Solve for the psychophysical parameter value at the specified
    % percent accuracy threshold
    parVal = double(max(solve(f.nr_rmax*((cf^f.nr_n)/(f.nr_c50^f.nr_n + cf^f.nr_n)) + f.nr_b == perc_thr,cf)));

end