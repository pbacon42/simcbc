
mergers = readcatalog("mergers.txt",{"dist","z","mass1","mass2","RA","dec"}," ")
save -V7 mergers.mat mergers