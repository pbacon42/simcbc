function Z=inv_comov_dist(D_comov,CosmoPars);
%---------------------------------------------------------------------------
% inv_comov_dist function                                          Cosmology
% Description: Given the comoving distance
%              to calculate the corresponding redshift.
% Input  : - Vector of comoving distance.
%            'D_comov'  - luminosity distance (default).
%          - Cosmological parameters : [H0, \Omega_{m}, \Omega_{\Lambda}],
%            or cosmological parameters structure, or a string containing
%            parameters source name, default is 'wmap3' (see cosmo_pars.m).
% Output : - Redshift.
% Reference : Perlmutter et al. 1997 ApJ, 483, 565
%             Oke & Sandage 1968 ApJ, 154, 21
%             Peterson, B.M., 1997, AGN, p.165
% Tested : Matlab 5.1
%     By : Eran O. Ofek                  March 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Z=inv_lum_comov(35,'dm');   
%
% Reliable: 2
%---------------------------------------------------------------------------

if (nargin==1),
   CosmoPars = "wmap3";
end

if (ischar(CosmoPars)==0 && isstruct(CosmoPars)==0),
   % do nothing
else
   Par = cosmo_pars(CosmoPars);
   CosmoPars = [Par.H0, Par.OmegaM, Par.OmegaL, Par.OmegaRad];
end
H0       = CosmoPars(1);
OmegaM   = CosmoPars(2);
OmegaL   = CosmoPars(3);
if (length(CosmoPars)==3),
   OmegaRad = 0;
else
   OmegaRad = CosmoPars(4);
end

Ndm     = length(D_comov);
Z       = zeros(Ndm,1);

% for Idm=1:1:Ndm,
%   #if mod(Idm,round(Ndm/10))==0
%   #printf("* "); fflush(stdout);
%   #endif
%   #Z(Idm) = fun_binsearch(@lum_dist,LumDist(Idm),[1e-15 100],1e-3,CosmoPars);
%   #fun=inline(sprintf("comoving_dist(x,params)-%g",1e6*D_comov(Idm)),"x","params");
%   fun_str=sprintf("comoving_dist(x,\"wmap3\")-%g",1e6*D_comov(Idm));
%   fun=inline(fun_str,"x");
%   #Z(Idm)=bisection(fun,1e-15,100,1e-3,CosmoPars);
%   Z(Idm)=fzero(fun,1);
% end

fun_str= arrayfun(@(x) sprintf("comoving_dist(x,\"wmap3\")-%g",x), 1e6*D_comov,'UniformOutput', false);
solver = @(str) fzero(inline(str,"x"),1);
Z = cellfun(solver,fun_str,'UniformOutput', true);
