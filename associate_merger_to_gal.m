clear

% load ad-hoc galaxy catalog
load("galaxies.mat");
N_gal = rows(galaxies);

% generate merger per galaxies count in observable volume
duration=1; % yr
R_merger=23.5; % Myr^-3 per MW galaxy from Table 2 of Dominik et al 2012. For model A, with solar metal.
N_merger=poissrnd(R_merger*duration/1e6*N_gal);

% randomly select galaxies that are candidates for hosting merger
selected_idx = unidrnd(N_gal,N_merger,1);
selected_gal = galaxies(selected_idx,:);

% load binary population
cbc_pop = load("ANSNS.02_proc.dat");
N_cbc = rows(cbc_pop);

% random association of mergers with galaxies
selected_idx = unidrnd(N_cbc,N_merger,1);
selected_cbc = cbc_pop(selected_idx,:);

printf("num of mergers: %d for %f yr\n",N_merger,duration)

% generate associate sky location: random or from the Glade catalog
load("glade_catalog.mat");
completeness_distance = 200; % Mpc
RA = unifrnd(0,360,N_merger,1);
dec = acosd(unifrnd(-1,+1,N_merger,1));
idx = find(selected_gal(:,3) < completeness_distance);

% replace by coordinate of closest galaxy if within completeness_distance
for n=1:length(idx)
    [_,candidate_ix]=min(abs(selected_gal(idx(n),3)-[catalog.dist]));
    RA(n) = catalog(candidate_ix).RA;
    dec(n)= catalog(candidate_ix).dec;
endfor

% merge catalog
merger_pop = [selected_gal selected_cbc RA dec];

% make plots
subplot(2,2,1)
hist(merger_pop(:,3),10); % luminosity distance
xlabel("luminosity distance Mpc")
title(sprintf("min dist=%3.1f Mpc ; max dist=%3.1f Mpc",min(merger_pop(:,3)),max(merger_pop(:,3))))

subplot(2,2,2)
hist(merger_pop(:,4),10);
xlabel("mass 1 M_{sun}")
title(sprintf("num of mergers: %d for %1.1f yr",N_merger,duration))

% save results
fid=fopen("mergers.txt","w");
fprintf(fid,"## luminosity distance (Mpc) redshift mass1 (Msun) mass2 (Msun) RA dec\n");
for n=1:length(merger_pop)
    fprintf(fid,"%4.4f %1.4f %2.4f %2.4f %4.4f %4.4f\n",merger_pop(n,3),merger_pop(n,1),...
    		       merger_pop(n,4),merger_pop(n,5),merger_pop(n,7),merger_pop(n,8));
endfor
fclose(fid);

mergers = readcatalog("mergers.txt",{"dist","z","mass1","mass2","RA","dec"}," ")
save -V7 mergers.mat mergers