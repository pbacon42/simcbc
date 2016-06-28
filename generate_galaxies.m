clear

% compute horizon distance and observable volume
z_limit=.12;
dHpc=comoving_dist(z_limit); % pc
dH=dHpc/1e6; % Mpc
VH=(2*dH)^3; % Mpc^3

% generate galaxies coordinates
R_gal=0.0116; % Mpc^-3 
N_gal=poissrnd(R_gal*VH);

x=dH*rand(N_gal,1);
y=dH*rand(N_gal,1);
z=dH*rand(N_gal,1);

d_comoving=sqrt(x.^2+y.^2+z.^2); % Mpc
idx=find(d_comoving<=dH);
d_comoving=d_comoving(idx);
z=inv_comov_dist(d_comoving);
dL=lum_dist(z)/1e6; % Mpc

subplot(2,2,1)
hist(d_comoving,50);

subplot(2,2,2)
hist(z,50);

subplot(2,2,3)
hist(dL,50)

fid=fopen("galaxies.txt","w");
fprintf(fid,"## z comoving distance (Mpc) luminosity distance (Mpc)\n");
for n=1:length(z)
    fprintf(fid,"%2.4f %4.4f %4.4f\n",z(n),d_comoving(n),dL(n));
endfor
fclose(fid);