
clear

load mergers.mat

clf

subplot(2,2,1)
hist([mergers.dist],10); % luminosity distance
xlabel("luminosity distance Mpc")
title(sprintf("min dist=%3.1f Mpc ; max dist=%3.1f Mpc",min([mergers.dist]),max([mergers.dist])))

subplot(2,2,2)
plot([mergers.RA],[mergers.dec],".","MarkerSize",5);
xlabel("RA (deg)")
xlabel("dec (deg)")
title(sprintf("num of mergers: %d -- %d yrs",length(mergers),duration))

subplot(2,2,3)
hist([mergers.mass1],10);
xlabel("mass 1 M_{sun}")

subplot(2,2,4)
hist([mergers.mass2],10);
xlabel("mass 2 M_{sun}")

print("-dpng","-FHelvetica:14","mergers.png")