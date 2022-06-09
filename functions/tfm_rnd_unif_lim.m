function r = tfm_rnd_unif_lim(lb,ub,nx,ny)
    if nargin < 3
        nx = 1;
        ny = 1;
    end
    r = (ub-lb).*rand(ny, nx) + lb;
end