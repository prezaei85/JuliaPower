# a 2-bus test case
ps = Caseps()

ps.baseMVA = 100.0

ps.bus = [
#ID type Pd Qd Gs Bs area Vmag Vang basekV Vmax Vmin lam P lam Q mu Vx mu Vn locX locY 
 1 3 0 0 0 0 1 1 0 230 1 1 1.05 0 0 0 0 1.5 1;
 2 1 0 0 0 0 1 1 0 230 1 1 1.05 0 0 0 0 2 1;
];

ps.branch = [
#from to R X B rateA rateB rateC tap shift status 
 1 2 0 0.1 0 1000 1000 1000 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
];

ps.gen = [
#bus Pg Qg Qmax Qmin Vsp mBase status Pmax Pmin mu_Px mu_Pn mu_Qx mu_Qn type cost 
 1 100 0 500 -500 1.00 100 1 1000 0 0 0 0 0 3 11.669;
];

ps.shunt = [
#bus P Q frac_S frac_Z status type value 
 2 400. 200. 1 0 1 1 100000 0;
];

ps.areas = [
 1 1;
];

ps.gencost = [
 2 0 0 3 0.00533 11.669 213.1;
];

ps = updateps(ps)
save("data/case2_ps.jld", "ps", ps) # needs: using Makeps, HDF5, JLD
