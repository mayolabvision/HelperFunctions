function [dat_new,these_alignTimes,goodFlag_new] = fix_specificSessions(session_name,np_mask,ripple_mask,alignTimes,dat,goodFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if isequal(session_name,'kendra_scrappy_0136a_g0') 
        if goodFlag % only used for kendra_scrappy_0136a_g0
            if sum(np_mask) >= numel(ripple_mask)
                these_alignTimes = alignTimes(np_mask);
                goodFlag_new = true;

                if numel(ripple_mask) ~= sum(ripple_mask)
                    dat_new = dat(ripple_mask);
                    fprintf('\n dat has %d rows', numel(dat_new))
                end
            elseif sum(np_mask) < numel(ripple_mask) % only used for kendra_scrappy_0136a_g0
                first_block_start = find(np_mask, 1, 'first');
                first_block_end = first_block_start + find(~np_mask(first_block_start:end), 1, 'first') - 2;
                first_zero_index = first_block_end + 1;
 
                good_alignTimes1 = alignTimes(first_block_start:first_block_end);
                remaining_alignTimes = alignTimes(first_zero_index:end);
                good_alignTimes2 = remaining_alignTimes(1:696);
    
                these_alignTimes = [good_alignTimes1; good_alignTimes2];
                dat_new(772:869) = [];
    
                goodFlag_new = false;
            end
        else
            these_alignTimes = alignTimes(end-313:end);
        end
    end
end