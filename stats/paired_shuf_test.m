function pvals_shuff = paired_shuf_test(base_fr,resp_fr,nShuff)

pvals_shuff = nan(size(base_fr,2),1);
for n = 1:size(base_fr,2)
    base = base_fr(:,n);
    resp  = resp_fr(:,n);

    % Skip neurons with no variance
    if std(resp - base) == 0
        pvals_shuff(n) = 1;
        continue
    end

    Treal = mean(resp - base);

    Tshuff = zeros(nShuff,1);
    for k = 1:nShuff
        flip = rand(size(base)) > 0.5;   % random trial-wise swaps

        resp_sh = resp;
        base_sh = base;

        resp_sh(flip)  = base(flip);
        base_sh(flip) = resp(flip);

        Tshuff(k) = mean(resp_sh - base_sh);
    end

    % Two-sided p-value
    pvals_shuff(n) = mean(abs(Tshuff) >= abs(Treal));
end

end