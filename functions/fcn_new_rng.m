function rnd_gen = fcn_new_rng()
    [h, m, sec] = hms(datetime);
    M = month(datetime); 
    d = day(datetime);
    seed = str2double(sprintf('%d', M,d,h,m,round(sec)));
    rng(seed);
    rnd_gen = rng;
end