function r = tfm_rnd_norm_lim(lb,ub,sig,avg, nx,ny)
    if nargin < 5
        nx = 1;
        ny = 1;
    end
    r = randn(ny, nx).*sig+avg;
    b_o = r<lb | r>ub;
    while sum(b_o(:)) > 0
        r(b_o) = randn(1, sum(b_o(:))).*sig+avg;
        b_o = r<lb | r>ub;
    end
%     hist(r(:))
end