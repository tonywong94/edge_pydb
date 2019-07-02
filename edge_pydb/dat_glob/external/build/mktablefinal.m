fp2=fopen('dearray.csv','r');
C=textscan(fp2,'%s %f %f','Delimiter',',','CommentStyle','#','Headerlines',1);
fclose(fp2);
oname=C{1};
vsys=C{2};
snr=C{3};
ix=strcmp(oname,'NGC4211N');
oname{ix}='NGC4211NED02';
ix=strcmp(oname,'UGC05498');
oname{ix}='UGC05498NED01';

clear p doc;
p.Name=oname;
doc.Name='Standard CALIFA name for galaxy';

% LEDA table of properties
leda=parseLedaTable('hyperleda_complete.txt');

for i=1:length(p.Name),
    gname=p.Name{i};
    if (strcmp(gname,'ARP220')),
        lname='ARP220';
    elseif (strcmp(gname,'IC2247')),
        lname='UGC04299';
    else
        lname=gname;
    end
    x=strcmp(lname,deblank(leda.name));
    j=find(x);
    if (isempty(j)),
        fprintf(1,'Name %s (%s) not found in LEDA table\n',lname,gname);
    else
        p.ledaRA(i)=leda.al2000(j);
        p.ledaDE(i)=leda.de2000(j);
        p.ledaA_Bgal(i)=leda.ag(j);
        p.ledaType(i)=leda.t(j);
        p.ledaD25(i)=0.1*10.^leda.logd25(j);
        p.ledaAxRatio(i)=10.^leda.logr25(j);
        p.ledaPA(i)=leda.pa(j);
        p.ledaIncl(i)=leda.incl(j);
        p.ledaVrad(i)=leda.vrad(j);
        p.ledaVmaxg(i)=leda.vmaxg(j);
        p.ledaVrot(i)=leda.vrot(j);
        p.ledaMorph{i}=leda.type(j);
        p.ledaBar{i}=leda.bar(j);
        p.ledaRing{i}=leda.ring(j);
        p.ledaMultiple{i}=(leda.multiple(j));
        p.ledaBt(i)=leda.bt(j);
        p.ledaIt(i)=leda.it(j);
        p.ledaMfir(i)=leda.mfir(j);
        p.ledaM21(i)=leda.m21(j);
        p.ledaVvir(i)=leda.vvir(j);
        p.ledaModz(i)=leda.modz(j);
    end
end
p.ledaModz(strcmp('NGC2880',p.Name))=31.71; % There seems to be a gross error in LEDA for this galaxy, substituting with NED
p.ledaModz(strcmp('NGC4211NED02',p.Name))=34.85; % Missing from LEDA, substituting with NED
p.ledaRA(strcmp('NGC4211NED02',p.Name))=183.905511/15;
p.ledaDE(strcmp('NGC4211NED02',p.Name))=28.169634;
p.ledaA_Bgal(strcmp('NGC4211NED02',p.Name))=0.099;
p.ledaModz(strcmp('UGC05498NED01',p.Name))=34.72; % Missing from LEDA, substituting with NED
p.ledaRA(strcmp('UGC05498NED01',p.Name))=153.015250/15;
p.ledaDE(strcmp('UGC05498NED01',p.Name))=23.085447;
p.ledaA_Bgal(strcmp('UGC05498NED01',p.Name))=0.114;
p.ledaDistMpc=10.^(p.ledaModz/5-5); % Computing Mpc distance from distance modulus
doc.ledaRA='RA J2000 from LEDA (hours).';
doc.ledaDE='DEC J2000 from LEDA (degrees).';
doc.ledaA_Bgal='Galactic extinction in B band from LEDA (mag). See doc for ag in LEDA. Based on Schlegel et al. (1998).';
doc.ledaType='Morphological code from LEDA (-5 to 10). See doc for t in LEDA.';
doc.ledaD25='D25 B-band apparent diameter according to LEDA, linearized to be in arcmin. See doc for logd25 in LEDA.';
doc.ledaAxRatio='Axis ratio from LEDA, linearized. See doc for logr25 in LEDA.';
doc.ledaPA='PA from LEDA (degrees). See doc for pa in LEDA.';
doc.ledaIncl='Inclination from LEDA (degrees). See doc for incl in LEDA.';
doc.ledaVrad='Heliocentric radial velocity from radio according to LEDA (km/s). It really corresponds to cz, so it is in the optical convention. Convert to radio convention using cz=v/(1-v/c). See doc for vrad in LEDA.';
doc.ledaVmaxg='Apparent maximum rotation velocity of the HI uncorrected by inclination (km/s). See doc for vmaxg in LEDA.';
doc.ledaVrot='Maximum rotation velocity according to LEDA (km/s). See doc for vrot in LEDA.';
doc.ledaMorph='Morphology according to LEDA. See doc for type in LEDA.';
doc.ledaBar='Detection of a bar reported (B). See doc for bar in LEDA.';
doc.ledaRing='Detection of a ring reported (R). See doc for ring in LEDA.';
doc.ledaMultiple='Object belongs to a multiple system (M). See doc for multiple in LEDA.';
doc.ledaBt='Apparent B total magnitude from LEDA (mag). See doc for bt in LEDA.';
doc.ledaIt='Apparent I total magnitude from LEDA (mag). See doc for it in LEDA.';
doc.ledaMfir='Far infrared flux from LEDA (mag). Computed as mfir=-2.5 log(2.58*f60+f100)+14.75 with fluxes in Jy. See doc for mfir in LEDA.';
doc.ledaM21='HI line flux from LEDA (mag). Computed as m21=2.5 log f + 17.40. See doc for m21 in LEDA.';
doc.ledaVvir='Galaxy systemic velocity corrected by Virgo Infall from LEDA. See doc for vvir in LEDA.';
doc.ledaModz='Kinematic distance modulus from LEDA. Uses modz=5 log(D_L)+25, with D_L luminosity distance from cosmology with Ho=70, Om=0.27, and Ol=0.73 based on Virgo infall corrected velocity. See doc for modz in LEDA. NGC2880 has a very wrong value in LEDA, substituted with NED scaled to same cosmology. NGC4211 and UGC05498 are both missing from LEDA and were substituted in the same way.';
doc.ledaDistMpc='Distance in Mpc corresponding to ledaModz kinematic distance modulus.';

% FIR photometry from the measurements on the images
% by Theodoros Bitsakis, compiling Vega magnitudes. 
% For two galaxies those are missing, and the measurements come from the 
% allWISE catalog
fp=fopen('apertures_WISE.dat','r');
C=textscan(fp,'%s %s %f %f %f','CommentStyle','#');
fclose(fp);
fid=fopen('wise_photometry_extended.dat','r');
F=textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
fclose(fid);

for i=1:length(p.Name),
    j=strcmp(p.Name{i},C{2});
    ix=find(j);
    if isempty(ix),
        disp(p.Name{i});
    else
        cn=C{1}(ix);
        j=strcmp(cn,F{1});
        ix=find(j);
        p.W1(i)=F{2}(ix);
        p.eW1(i)=F{3}(ix);
        p.W2(i)=F{4}(ix);
        p.eW2(i)=F{5}(ix);
        p.W3(i)=F{6}(ix);
        p.eW3(i)=F{7}(ix);
        p.W4(i)=F{8}(ix);
        p.eW4(i)=F{9}(ix);
    end
end
doc.W1='W1 Vega magnitude from image photometry from table by Th. Bitsakis. Two galaxies (NGC0598 and NGC4676A) are not in the tables, and W1 is from w1mpro in the allWISE catalog in IPAC.';
doc.eW1='Error in W1 in Vega magnitudes from table by Th. Bitsakis. For NGC0598 and NGC4676A is computed as w1mpro+2.5*log(w1snr) from AllWISE catalog.';
doc.W2='W2 see W1 doc';
doc.eW2='Error in W2 see W1 doc';
doc.W3='W3 see W1 doc';
doc.eW3='Error in W3 see W1 doc';
doc.W4='W4 see W1 doc';
doc.eW4='Error in W4 see W1 doc';
t=parseIpacTable('califa_wisesearch.txt');
missing={'NGC0528','NGC4676A'};
for ix=1:2,
    i=find(strcmp(missing{ix},t.name_01));
    j=find(strcmp(missing{ix},p.Name));
    p.W1(j)=t.w1mpro(i);
    p.eW1(j)=t.w1mpro(i)+2.5*log10(t.w1snr(i));
    p.W2(j)=t.w2mpro(i);
    p.eW2(j)=t.w2mpro(i)+2.5*log10(t.w2snr(i));
    p.W3(j)=t.w3mpro(i);
    p.eW3(j)=t.w3mpro(i)+2.5*log10(t.w3snr(i));
    p.W4(j)=t.w4mpro(i);
    p.eW4(j)=t.w4mpro(i)+2.5*log10(t.w4snr(i));
end


% Use Tony de20_flux.csv fluxes
tt=parseCsvTable('de20_flux.csv');
for i=1:length(p.Name),
    ix=find(strcmp(p.Name{i},tt.Galaxy));
    if isempty(ix),
        if strcmp(p.Name{i},'NGC4211NED02'), 
            ix=find(strcmp('NGC4211N',tt.Galaxy));
        elseif strcmp(p.Name{i},'UGC05498NED01'), 
            ix=find(strcmp('UGC05498',tt.Galaxy)); 
        end
    end
    p.coNomask(i)=tt.Nomask(ix);
    p.coeNomask(i)=tt.Nom_err(ix);
    p.coDilated(i)=tt.Dilated(ix);
    p.coeDilated(i)=tt.Dil_err(ix);
    p.coSmooth(i)=tt.Smooth(ix);
    p.coeSmooth(i)=tt.Smo_err(ix);
    p.coMask2d(i)=tt.Mask2d(ix);
    p.coeMask2d(i)=tt.Mask2d_err(ix);
    r=rfits(['DEdata/DEarray20/' tt.Galaxy{ix} '.co.de20_str.mom0.fits'],'head');
    p.coBmaj(i)=3600*r.bmaj;
    p.coBmin(i)=3600*r.bmin;
    p.coBpa(i)=r.bpa;
    p.coDvhel(i)=r.dvhel;
end
doc.coNomask='CO flux from unmasked EDGE cube (Jy km/s).';
doc.coeNomask='Error in CO flux from unmasked EDGE cube (Jy km/s).';
doc.coDilated='CO flux from masked EDGE cube using Dilated method (Jy km/s).';
doc.coeNomask='Error in CO flux from masked EDGE cube using Dilated method (Jy km/s).';
doc.coSmooth='CO flux from masked EDGE cube using Smooth method (Jy km/s).';
doc.coeSmooth='Error in CO flux from masked EDGE cube using Smooth method (Jy km/s).';
doc.coMask2d='CO flux from masked EDGE cube using fixed 2D mask method (Jy km/s).';
doc.coeMask2d='Error in CO flux from masked EDGE cube using fixed 2D mask method (Jy km/s).';
doc.coBmaj='Beam major axis from CO map (arcsec).';
doc.coBmin='Beam minor axis from CO map (arcsec).';
doc.coBpa='Beam PA from CO map (degrees).';
doc.coDvhel='Correction from LSR to Heliocentric (km/s)';

% Use get_mag_cubes_v2.2.csv for the SDSS photometry
mp=parseSebastianTable('get_mag_cubes_v2.2.csv');
mg=parseSebastianTable('Pipe3D_NSA_CALIFA-DR3_candidates.csv'); % Parameters up to column 138 are from the NASA-Sloan Atlas
me=parseSebastianTable('get_proc_elines_CALIFA.csv');

for i=1:length(p.Name),
    j=find(strcmp(p.Name{i},mg.CALIFANAME));
    if isempty(j), fprintf(1,'Galaxy %s not found in Pipe3D_NSA_CALIFA-DR3_candidates.csv\n',p.Name{i}); continue; end
    p.nsaZdist(i)=mg.ZDIST(j);
    p.nsaAu(i)=mg.EXTINCTION_u(j);
    p.nsaAg(i)=mg.EXTINCTION_g(j);
    p.nsaAr(i)=mg.EXTINCTION_r(j);
    p.nsaAi(i)=mg.EXTINCTION_i(j);
    p.nsaAz(i)=mg.EXTINCTION_z(j);
    p.caMass(i)=mg.logMassPipe3D(j);
    p.caeMass(i)=mg.errorlogMassPipe3D(j);
    p.caSFR(i)=mg.logSFRPipe3D(j);
    p.caeSFR(i)=mg.errorlogSFRPipe3D(j);
    p.caOH(i)=mg.logOHPipe3D(j);
    p.caeOH(i)=mg.errorlogOHPipe3D(j);
    p.caAvgas(i)=mg.Av_gasPipe3D(j);
    p.caeAvgas(i)=mg.errorAv_gasPipe3D(j);
    p.caAvstars(i)=mg.Av_starsPipe3D(j);
    p.caeAvstars(i)=mg.errorAv_starsPipe3D(j);
end
doc.nsaZdist='From NSA Atlas, distance estimate using pecular velocity model of Willick et al. (1997); multiply by c/H0 for Mpc.';
doc.nsaAu='From NSA Atlas, Galactic extinction in u from Schlegel, Finkbeiner & Davis (1997).';
doc.nsaAg='From NSA Atlas, Galactic extinction in g from Schlegel, Finkbeiner & Davis (1997).';
doc.nsaAr='From NSA Atlas, Galactic extinction in r from Schlegel, Finkbeiner & Davis (1997).';
doc.nsaAi='From NSA Atlas, Galactic extinction in i from Schlegel, Finkbeiner & Davis (1997).';
doc.nsaAz='From NSA Atlas, Galactic extinction in z from Schlegel, Finkbeiner & Davis (1997).';
doc.caMass='Stellar mass as log(Mass/Msun) from column 149 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';
doc.caeMass='Error in stellar mass as log(Mass/Msun) from column 149 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';
doc.caSFR='SFR as log(SFR/[Msun/yr]) from column 151 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';
doc.caeSFR='Error in SFR as log(SFR/[Msun/yr]) from column 151 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';
doc.caOH='Oxygen abundance as 12+log(O/H) from column 153 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';
doc.caeOH='Error in oxygen abundance as 12+log(O/H) from column 153 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';
doc.caAvgas='Nebular extinction as Av from column 155 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';
doc.caeAvgas='Error in nebular extinction as Av from column 155 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';
doc.caAvstars='Stellar extinction as Av from column 157 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';
doc.caeAvstars='Error in stellar extinction as Av from column 157 in Pipe3D_NSA_CALIFA-DR3_candidates.csv';

for i=1:length(p.Name),
    j=find(strcmp(p.Name{i},mp.name_obj));
    if isempty(j), fprintf(1,'Galaxy %s not found in get_mag_cubes_v2.2.csv\n',p.Name{i}); continue; end
    p.Su(i)=mp.ubandmag(j);
    p.Sg(i)=mp.gbandmag(j);
    p.Sr(i)=mp.rbandmag(j);
    p.Si(i)=mp.ibandmag(j);
%    p.Sz(i)=mp.zbandmag(j);
    p.caB(i)=mp.Bbandmag(j);
    p.caV(i)=mp.Vbandmag(j);
    p.caR(i)=mp.Rbandmag(j);
    p.caRe(i)=mp.Rearcsec(j);
    p.caeRe(i)=mp.errorRearcsec(j);
    p.caEllipticity(i)=mp.ellipticy(j);
    p.caPA(i)=mp.Padeg(j);
    p.caR50(i)=mp.R50arcsec(j);
    p.caeR50(i)=mp.errorR50arcsec(j);
    p.caR90(i)=mp.R90arcsec(j);
    p.caeR90(i)=mp.errorR90arcsec(j);
end
doc.Su='SDSS u magnitude from CALIFA synthetic photometry, from ubandmag in get_mag_cubes_v2.2.csv. Foreground extinction corrected.';
doc.Sg='SDSS g magnitude from CALIFA synthetic photometry, from gbandmag in get_mag_cubes_v2.2.csv. Foreground extinction corrected.';
doc.Sr='SDSS r magnitude from CALIFA synthetic photometry, from rbandmag in get_mag_cubes_v2.2.csv. Foreground extinction corrected.';
doc.Si='SDSS i magnitude from CALIFA synthetic photometry, from ibandmag in get_mag_cubes_v2.2.csv. Foreground extinction corrected.';
%doc.Sz='SDSS z magnitude from zbandmag in get_mag_cubes_v2.2.csv';
doc.caB='B magnitude from Bbandmag in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction.';
doc.caV='V magnitude from Vbandmag in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction.';
doc.caR='R magnitude from Rbandmag in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction.';
doc.caRe='Equivalent radius Re (arcsec) from Re in get_mag_cubes_v2.2.csv';
doc.caeRe='Error in equivalent radius Re (arcsec) from error Re in get_mag_cubes_v2.2.csv';
doc.caEllipticity='Ellipticity sqrt(1-b^2/a^2) from ellipticy in get_mag_cubes_v2.2.csv';
doc.caPA='PA (degrees) from Pa in get_mag_cubes_v2.2.csv';
doc.caR50='R50 (arcsec) from R50 in get_mag_cubes_v2.2.csv';
doc.caeR50='Error in R50 (arcsec) from R50 in get_mag_cubes_v2.2.csv';
doc.caR90='R90 (arcsec) from R90 in get_mag_cubes_v2.2.csv';
doc.caeR90='Error in R90 (arcsec) from R90 in get_mag_cubes_v2.2.csv';

for i=1:length(p.Name),
    j=find(strcmp(p.Name{i},me.name));
    if isempty(j), fprintf(1,'Galaxy %s not found in get_proc_elines_CALIFA.csv\n',p.Name{i}); continue; end
    p.caOH_O3N2(i)=me.OH_O3N2(j);
    p.caZgas(i)=me.z_gas(j);
    p.caZstars(i)=me.z_stars(j);
    p.caAge(i)=me.log_age_mean_LW(j);
    p.caeAge(i)=me.s_log_age_mean_LW(j);
end
doc.caOH_O3N2='CALIFA metallicity from O3N2, from OH_O3N2 column in get_proc_elines_CALIFA.csv';
doc.caZgas='CALIFA redshift for gas lines, from z_gas column in get_proc_elines_CALIFA.csv';
doc.caZstars='CALIFA redshift for stars, from z_stars column in get_proc_elines_CALIFA.csv';
doc.caAge='CALIFA mean age, from log_age_mean_LW column in get_proc_elines_CALIFA.csv';
doc.caeAge='Error in CALIFA mean age, from s_log_age_mean_LW column in get_proc_elines_CALIFA.csv';

p.caDistMpc=luminosity_distance(p.caZgas);
doc.caDistMpc='Luminosity distance in Mpc computed from caZgas assuming Ho=70, Om=0.27, Ol=0.73';
p.coMmol=(1.05e4*p.coMask2d.*p.caDistMpc.^2./(1+p.caZgas));
doc.coMmol='Molecular mass from coMask2d flux and caDistMpc distance and redshift using Eq. 3 in Bolatto, Wolfire, & Leroy (2013).';
p.coeMmol=(1.05e4*p.coeMask2d.*p.caDistMpc.^2./(1+p.caZgas));
doc.coeMmol='Error in molecular mass from coeMask2d flux and caDistMpc distance and redshift.';
wisezp=[5.4188e-15,2.5172e-15,3.5878e-16,2.0876e-17]*1e7; % photometric zero points erg/s/cm2 Jarrett + 2011
wisenu=[8.856e13,6.45e13,2.675e13,1.346e13]; % iso frequencies of wise filters from Jarrett
wisednu=[1.751e13,1.465e13,1.133e13,2.496e12]; % filter widths from Jarrett
p.W1lum=wisenu(1)/wisednu(1)*wisezp(1)*10.^(-p.W1/2.5)*4*pi.*(p.caDistMpc*1e6*3.09e18).^2;
p.W2lum=wisenu(2)/wisednu(2)*wisezp(2)*10.^(-p.W2/2.5)*4*pi.*(p.caDistMpc*1e6*3.09e18).^2;
p.W3lum=wisenu(3)/wisednu(3)*wisezp(3)*10.^(-p.W3/2.5)*4*pi.*(p.caDistMpc*1e6*3.09e18).^2;
p.W4lum=wisenu(4)/wisednu(4)*wisezp(4)*10.^(-p.W4/2.5)*4*pi.*(p.caDistMpc*1e6*3.09e18).^2;
doc.W1lum='Luminosity in W1 (erg/s) from caDistMpc and W1 magnitude using the zero point, frequency, and bandwidth from Jarrett et al. (2011).';
doc.W2lum='Luminosity in W2 (erg/s) from caDistMpc and W2 magnitude using the zero point, frequency, and bandwidth from Jarrett et al. (2011).';
doc.W3lum='Luminosity in W3 (erg/s) from caDistMpc and W3 magnitude using the zero point, frequency, and bandwidth from Jarrett et al. (2011).';
doc.W4lum='Luminosity in W4 (erg/s) from caDistMpc and W4 magnitude using the zero point, frequency, and bandwidth from Jarrett et al. (2011).';
p.W4SFR=2.1e-43*p.W4lum;
doc.W4SFR='SFR (Msun/yr) from W4lum using Catalan-Torrecilla et al. (2014) calibration.';

becca=parseCsvTable('EDGE_srcparameters.csv');
for i=1:length(p.Name),
    j=find(strcmp(p.Name{i},becca.Source));
    if isempty(j), fprintf(1,'Galaxy %s not found in EDGE_srcparameters.csv\n',p.Name{i}); end
    p.coVsys(i)=becca.VsysRel(j);
    p.coVrotmax(i)=becca.VrotMaxRel(j); if (p.coVrotmax(i)==0), p.coVrotmax(i)=NaN; end
    p.coPA(i)=becca.PA(j);
    p.coInc(i)=becca.inc(j);
end
doc.coVsys='Systemic galaxy velocity (relativistic convention), derived from CO rotation curve or taken from LEDA if no rotation curve is possible.';
doc.coVrotmax='Maximum rotational velocity of the CO, determined from the outermost rings where the rotation velocity does not change by more than 1.8 km/s over 3 adjacent rings.';
doc.coPA='Galaxy PA determined from the CO rotation curve using the orientation that minimizes the radial component, or if not available from Glenn van de Ven table of orientation of outer isophotes, or otherwise from LEDA.';
doc.coInc='Galaxy PA determined from the CO rotation curve using the inclination that produces the best fit, or if not available from Glenn van de Ven table of orientation of outer isophotes, or otherwise from LEDA.';

prof=parseCsvTable('profiles.csv',1);
for i=1:length(p.Name),
    j=find(strcmp(p.Name{i},prof.Name));
    if isempty(j), fprintf(1,'Galaxy %s not found in profiles.csv\n',p.Name{i}); continue; end
    p.coScaleMol(i)=prof.ScaleMol(j);
    p.coScaleMolHi(i)=prof.ScaleMolHi(j);
    p.coScaleMolLo(i)=prof.ScaleMolLo(j);
    p.coNormMol(i)=prof.NormMol(j);
    p.coeNormMol(i)=prof.eNormMol(j);
    p.coScaleSt(i)=prof.ScaleSt(j);
    p.coeScaleSt(i)=prof.eScaleSt(j);
    p.coNormSt(i)=prof.NormSt(j);
    p.coeNormSt(i)=prof.eNormSt(j);
    p.coR50Mol(i)=prof.R50Mol(j);
    p.coeR50Mol(i)=prof.eR50Mol(j);
    p.coR50St(i)=prof.R50St(j);
    p.coeR50St(i)=prof.eR50St(j);
end
doc.coScaleMol='Exponential scale length for CO disk (kpc). Derived by filling in undetected values in annulii with 1-sigma.';
doc.coScaleMolHi='Upper limit to exponential scale length by filling in annulii with 2-sigma values (kpc).';
doc.coScaleMolLo='Lower limit to exponential scale length by filling in annulii with zeros (kpc).';
doc.coNormMol='Normalization of exponential disk profile, i.e. density at R=0 (Msun/pc^2).';
doc.coeNormMol='Error in normalization of exponential disk profile (Msun/pc^2).';
doc.coScaleSt='Exponential scale length for the mass of the stellar disk (kpc).';
doc.coeScaleSt='Formal error in stellar scale length fit (kpc).';
doc.coNormSt='Normalization in stellar exponential disk profile, i.e. density at R=0 (Msun/pc^2).';
doc.coeNormSt='Formal error in stellar disk normalization (Msun/pc^2).';
doc.coR50Mol='Radius enclosing 50% of the molecular mass (kpc).';
doc.coeR50Mol='Error in radius enclosing 50% of the molecular mass (kpc), including beam size.';
doc.coR50St='Radius enclosing 50% of the stellar mass (kpc).';
doc.coeR50St='Error in radius enclosing 50% of the stellar mass (kpc), including beam size.';

writeCsvTable(p,'DETableFinal.csv');
writeDoc(doc,'DETableFinal.doc.txt');


