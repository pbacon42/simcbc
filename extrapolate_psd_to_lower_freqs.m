% Extend PSD model to lower frequency using an ad-hoc extrapolation

clear

psd_infile = "aLIGO_80Mpc_PSD.txt";
[~,file_name,file_ext]=fileparts(psd_infile);

data = load(psd_infile);

freq = data(:,1);
psd = data(:,2);

logfreq = log10(freq);
df = logfreq(2)-logfreq(1);

flow_cut = 15; % Hz
flow_extend = .2; % Hz
psdflow_extend = 1e-42; % 1/Hz

N_extend = length(freq) + round(log10(flow_extend/freq(1))/df);

[~, idx] = min(abs(freq - flow_cut));
psdflow_cut = psd(idx);

freq_extend = logspace(log10(flow_extend), logfreq(end), N_extend);
psd_extend = ones(1,N_extend);

[~, idx] = min(abs(freq_extend - flow_cut));

psd_extend(idx:end) = interp1(freq, psd, freq_extend(idx:end));
psd_extend(end)=psd(end); % XXX dirty fix of a bug in interp1 XXX
psd_extend(1:idx) = logspace(log10(psdflow_extend),log10(psdflow_cut),idx);

psd_outfile = [file_name,"_extended",file_ext];

model_extended = [freq_extend;psd_extend]';

save("-ascii",psd_outfile,"model_extended")

clf
loglog(freq,psd,"b;original;");
hold on
loglog(freq_extend,psd_extend,"r;extended;");
xlabel("frequency (Hz)");
ylabel("power spec density (1/Hz)");




