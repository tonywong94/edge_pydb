edge_bbpars_natv.csv

Global fit parameters from Bbarolo on natv CO data

| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy name |
| bbRactr | deg | float64 | .4f | R.A. J2000 of center used by Bbarolo |
| bbDectr | deg | float64 | .4f | Dec. J2000 of center used by Bbarolo |
| bbVsys | km / s | float64 | .2f | Systemic velocity determined by Bbarolo (radio-LSR) |
| bbKinInc | deg | float64 | .2f | Inclination determined by Bbarolo |
| bbKinPA | arcsec | float64 | .2f | Position angle E from N determined by Bbarolo |
| bbMask |   | string |   | Identifier for cube mask |

edge_bbpars_smo7.csv

Global fit parameters from Bbarolo on smo7 CO data

| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy name |
| bbRactr | deg | float64 | .4f | R.A. J2000 of center used by Bbarolo |
| bbDectr | deg | float64 | .4f | Dec. J2000 of center used by Bbarolo |
| bbVsys | km / s | float64 | .2f | Systemic velocity determined by Bbarolo (radio-LSR) |
| bbKinInc | deg | float64 | .2f | Inclination determined by Bbarolo |
| bbKinPA | arcsec | float64 | .2f | Position angle E from N determined by Bbarolo |
| bbMask |   | string |   | Identifier for cube mask |

edge_hiflux.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| Refcode |   | string |   | Ref Code for spectrum |
| Vsys | km / s | float32 |   | 'Vsys from LEDA, optical barycentric' |
| Deltav | km / s | float32 |   | Velocity channel width |
| Robust_rms | mJy | float32 |   | rms noise from mad_std over noise region |
| RefInt | Jy km / s | float32 |   | integrated HI flux over fixed window of +/- 600.0 |
| RefUnc | Jy km / s | float32 |   | uncertainty in RefInt |
| SigInt | Jy km / s | float32 |   | integrated HI flux over get_signal_range window |
| SigUnc | Jy km / s | float32 |   | uncertainty in SigInt |
| SigVmin | km / s | float32 |   | signal start from get_signal_range |
| SigVmax | km / s | float32 |   | signal end from get_signal_range |
| BadFlag |   | bool |   | True for poor signal |

edge_rfpars.csv



date: '2019-10-13'

| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | standard EDGE-CALIFA galaxy name |
| rfVsys | km / s | float64 |   | CO relativistic helocentric vsys (km/s); see VsysFlag for how this was derived |
| rfVsysFlag |   | int64 |   | how is Vsys derived; 0=CO kinematics; 1=Ha kinematics; 3=LEDA |
| rfLSRK2helio | km / s | float64 |   | conversion from LSRK velocity to heliocentric velocity (km/s); V(helio)=V(LSRK)-LSRK2helio |
| rfRflat | arcsec | float64 |   | the radius (arcsec) at which the rotation curve flattens from the Lelli+16 algorithm |
| rfVrotMaxL | km / s | float64 |   | the maximum CO rotation velocity (relativistic km/s) from the Lelli+2016 algorithm |
| rfeVrotMaxL | km / s | float64 |   | the error in rfVrotMaxL from the Lelli+2016 algorithm |
| rfVrotMaxLFlag |   | int64 |   | 0 = either rfeVrotMaxL\rfVrotMaxL > 10% rfeVrotMaxL>rfVrotMaxL rfeVrotmaxL=0 so rfRflat rfVrotMaxL and rfeVrotMaxL values are unreliable; 1 = rfRflat rfVrotMaxL and rfeVrotMaxL values are reliable |
| rfVrotMax | km / s | float64 |   | the maximum CO rotation velocity (relativistic km/s) using median Vrot at radii larger than 12 arcsec as used in Levy+18 |
| rfeVrotMax | km / s | float64 |   | the error in rfVrotMax using standard deviation of Vrot at radii larger than 12 arcsec as used in Levy+18 |
| rfPA | deg | float64 |   | position angle (degrees) ranging from 0-365 degrees defined counterclockwise from North to approaching side of major axis; see PAFlag for how this was derived |
| rfPAFlag |   | int64 |   | this describes where the PA listed is from; 0=CO kinematics; 1=Ha kinematics; 2=outer optical isophotes; 3=LEDA |
| rfInc | deg | float64 |   | inclination (degrees); see IncFlag for how this was derived |
| rfIncFlag |   | int64 |   | this describes where the inclination listed is from; 0=CO kinematics; 1=Ha kinematics; 2=outer optical isophotes; 3=LEDA |
| rfKinXoff | arcsec | float64 |   | offset of CO kinematic center from ledaRA (arcsec) |
| rfKinYoff | arcsec | float64 |   | offset of CO kinematic center from ledaDec (arcsec) |
| rfKinRA | hourangle | float64 |   | RA of kinematic center including rfKinXoff (J2000 hours); rfKinRA=ledaRA+rfKinXoff |
| rfKinDecl | deg | float64 |   | Decl of kinematic center including rfKinYoff (J2000 degrees); rfKinDecl=ledaDE+rfKinYoff |

edge_califa.csv



date: '2019-08-31'

| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   |   |
| caMass | dex(solMass) | float64 |   | Stellar mass from col 149 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| caeMass | dex(solMass) | float64 |   | Error in stellar mass from col 150 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| caSFR | dex(solMass / yr) | float64 |   | SFR from col 151 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| caeSFR | dex(solMass / yr) | float64 |   | Error in SFR from col 152 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| caOH | dex | float64 |   | Oxygen abundance as 12+log(O/H) from col 153 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| caeOH | dex | float64 |   | Error in oxygen abundance from col 154 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| caAvgas | mag | float64 |   | Nebular extinction as Av from col 155 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| caeAvgas | mag | float64 |   | Error in nebular extinction from col 156 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| caAvstars | mag | float64 |   | Stellar extinction as Av from col 157 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| caeAvstars | mag | float64 |   | Error in stellar extinction from col 158 in Pipe3D_NSA_CALIFA-DR3_candidates.csv |
| Su | mag | float64 |   | SDSS u magnitude from col 4 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction |
| Sg | mag | float64 |   | SDSS g magnitude from col 8 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction |
| Sr | mag | float64 |   | SDSS r magnitude from col 12 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction |
| Si | mag | float64 |   | SDSS i magnitude from col 16 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction |
| caB | mag | float64 |   | B magnitude from col 20 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction |
| caV | mag | float64 |   | V magnitude from col 24 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction |
| caR | mag | float64 |   | R magnitude from col 28 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction |
| caRe | arcsec | float64 |   | Equivalent radius Re from col 34 in get_mag_cubes_v2.2.csv |
| caeRe | arcsec | float64 |   | Error in equivalent radius from col 35 in get_mag_cubes_v2.2.csv |
| caEllipticity |   | float64 |   | Ellipticity sqrt(1-b^2/a^2) from col 38 in get_mag_cubes_v2.2.csv |
| caPA | deg | float64 |   | PA from col 39 in get_mag_cubes_v2.2.csv |
| caR50 | arcsec | float64 |   | R50 from col 40 in get_mag_cubes_v2.2.csv |
| caeR50 | arcsec | float64 |   | Error in R50 from col 41 in get_mag_cubes_v2.2.csv |
| caR90 | arcsec | float64 |   | R90 from col 42 in get_mag_cubes_v2.2.csv |
| caeR90 | arcsec | float64 |   | Error in R90 from col 43 in get_mag_cubes_v2.2.csv |
| caOH_O3N2 | dex | float64 |   | O3N2-based metallicity from <OH_O3N2> col 4 in get_proc_elines_CALIFA.csv |
| caZgas |   | float64 |   | 'Redshift for gas lines, from <z_gas> col 14 in get_proc_elines_CALIFA.csv' |
| caZstars |   | float64 |   | 'Redshift for stars, from <z_stars> col 15 in get_proc_elines_CALIFA.csv' |
| caAge | dex(Gyr) | float64 |   | Mean stellar age from <log_age_mean_LW> col 37 in get_proc_elines_CALIFA.csv |
| caeAge | dex(Gyr) | float64 |   | Error in mean stellar age from <s_log_age_mean_LW> col 38 in get_proc_elines_CALIFA.csv |
| caFHa | dex(1e-16 erg / (cm2 s)) | float64 |   | 'Log of Halpha flux, from <log_F_Ha> column 145 in get_proc_elines_CALIFA.csv' |
| caFHacorr | dex(1e-16 erg / (cm2 s)) | float64 |   | 'Log of Halpha flux, extinction corrected, from <log_F_Ha_cor> column 146 in get_proc_elines_CALIFA.csv' |
| caLHacorr | dex(erg / s) | float64 |   | 'Log of Halpha luminosity, extinction corrected, from <log_L_Ha_cor> column 147 in get_proc_elines_CALIFA.csv' |
| caMstars | dex(solMass) | float64 |   | 'Log of stellar mass, from column 73 log_Mass in get_proc_elines_CALIFA.csv' |
| caDistMpc | Mpc | float64 |   | 'Luminosity distance in Mpc computed from caZgas assuming Ho=70, Om=0.27, Ol=0.73' |
| caDistP3d | Mpc | float64 |   | Luminosity distance in Mpc from get_proc_elines_CALIFA.csv |
| caFlgWav5 |   | int64 |   | Flag (-1/0/1/2=NA/good/minor/bad) for wavelength calibration V500 |
| caFlgReg5 |   | int64 |   | Flag (-1/0/1/2=NA/good/minor/bad) for 2D registration rel to SDSS V500 |
| caFlgImg5 |   | int64 |   | Flag (-1/0/1/2=NA/good/minor/bad) for reconstructed image quality V500 |
| caFlgWav12 |   | int64 |   | Flag (-1/0/1/2=NA/good/minor/bad) for wavelength calibration V1200 |
| caFlgReg12 |   | int64 |   | Flag (-1/0/1/2=NA/good/minor/bad) for 2D registration rel to SDSS V1200 |
| caFlgImg12 |   | int64 |   | Flag (-1/0/1/2=NA/good/minor/bad) for reconstructed image quality V1200 |

edge_leda.csv



date: '2019-10-13'

| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   |   |
| ledaRA | hourangle | float64 |   | RA J2000 from LEDA <celpos> |
| ledaDE | deg | float64 |   | DEC J2000 from LEDA <celpos> |
| ledaA_Bgal | mag | float64 |   | Galactic A_B from LEDA <ag> |
| ledaType |   | float64 |   | Morphological type from LEDA <t> |
| ledaMorph |   | string |   | Hubble type from LEDA <type> |
| ledaD25 | arcmin | float64 |   | Apparent B diameter from LEDA <logd25> linearized |
| ledaAxRatio |   | float64 |   | Maj/min axis ratio from LEDA <logr25> linearized |
| ledaPA | deg | float64 |   | 'PA from LEDA <pa>, N to E' |
| ledaIncl | deg | float64 |   | Morph inclination from LEDA <incl> |
| ledaVrad | km / s | float64 |   | cz from radio data from LEDA <vrad> |
| ledaVmaxg | km / s | float64 |   | HI max v_rot not corr for incl from LEDA <vmaxg> |
| ledaVrot | km / s | float64 |   | HI max v_rot corr for incl from LEDA <vrot> |
| ledaBar |   | string |   | B = bar present from LEDA <bar> |
| ledaRing |   | string |   | R = ring present from LEDA <ring> |
| ledaMultiple |   | string |   | M = multiple system from LEDA <multiple> |
| ledaBt | mag | float64 |   | Apparent B total magnitude from LEDA <bt> |
| ledaIt | mag | float64 |   | Apparent I total magnitude from LEDA <it> |
| ledaMfir | mag | float64 |   | FIR flux as magnitude from LEDA <mfir> |
| ledaM21 | mag | float64 |   | HI line flux as magnitude from LEDA <m21> |
| ledaVvir | km / s | float64 |   | Virgo infall corr cz from LEDA <vvir> |
| ledaModz | mag | float64 |   | Dist modulus from LEDA <modz> based on <vvir>; NGC2880 NGC4211 and UGC05498 substituted with NED scaled to same cosmology |
| ledaDistMpc | Mpc | float64 |   | Distance in Mpc corresponding to ledaModz |

edge_ned.csv



date: '2019-08-31'

| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   |   |
| nedRA | deg | float64 |   | J2000 RA from NED |
| nedDE | deg | float64 |   | J2000 DEC from NED |
| nedVopt | km / s | int64 |   | cz from NED |

edge_nsa.csv



date: '2019-08-31'

| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   |   |
| nsaZdist |   | float64 |   | Redshift from NASA-Sloan Atlas using pecular velocity model of Willick+ 97 |
| nsaAu | mag | float64 |   | Galactic extinction in u from Schlegel+ 97 |
| nsaAg | mag | float64 |   | Galactic extinction in g from Schlegel+ 97 |
| nsaAr | mag | float64 |   | Galactic extinction in r from Schlegel+ 97 |
| nsaAi | mag | float64 |   | Galactic extinction in i from Schlegel+ 97 |
| nsaAz | mag | float64 |   | Galactic extinction in z from Schlegel+ 97 |

edge_rdist.csv



date: '2019-08-31'

| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   |   |
| rdScaleMol | kpc | float64 |   | Exponential scale length for CO disk derived by filling in undetected values in annulii with 1-sigma |
| rdeScaleMol | kpc | float64 |   | Statistical error from fit to exponential scale length |
| rdScaleMolHi | kpc | float64 |   | Upper limit to exponential scale length by filling in annuli with 2-sigma values |
| rdScaleMolLo | kpc | float64 |   | Lower limit to exponential scale length by filling in annulii with zeros |
| rdNormMol | solMass / pc2 | float64 |   | 'Normalization of CO exponential disk profile, i.e. density at R=0' |
| rdeNormMol | solMass / pc2 | float64 |   | Error in normalization of CO exponential disk profile |
| rdScaleSt | kpc | float64 |   | Exponential scale length for the mass of the stellar disk |
| rdeScaleSt | kpc | float64 |   | Formal error in stellar scale length fit |
| rdNormSt | solMass / pc2 | float64 |   | 'Normalization of stellar exponential disk profile, i.e. density at R=0' |
| rdeNormSt | solMass / pc2 | float64 |   | Formal error in stellar disk normalization |
| rdR50Mol | kpc | float64 |   | Radius enclosing 50% of the molecular mass |
| rdeR50Mol | kpc | float64 |   | 'Error in radius enclosing 50% of the molecular mass, including beam size' |
| rdR50St | kpc | float64 |   | Radius enclosing 50% of the stellar mass |
| rdeR50St | kpc | float64 |   | 'Error in radius enclosing 50% of the stellar mass, including beam size' |
| rdScaleSFR | kpc | float64 |   | Exponential scale length for SFR from extinction corrected Ha |
| rdeScaleSFR | kpc | float64 |   | Error in exponential scale length for SFR from extinction corrected Ha |

edge_wise.csv



date: '2019-08-31'

| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   |   |
| W1 | mag | float64 |   | 3.4 um Vega magnitude from image photometry from table by Bitsakis. NGC0598 and NGC4676A are from w1mpro in allWISE catalog |
| eW1 | mag | float64 |   | Error in W1 in Vega magnitudes from table by Bitsakis. NGC0598 and NGC4676A are from w1mpro+2.5*log(w1snr) in allWISE catalog |
| W2 | mag | float64 |   | 4.6 um Vega magnitude from Bitsakis |
| eW2 | mag | float64 |   | Error in W2 in Vega magnitudes from Bitsakis |
| W3 | mag | float64 |   | 12 um Vega magnitude from Bitsakis |
| eW3 | mag | float64 |   | Error in W3 in Vega magnitudes from Bitsakis |
| W4 | mag | float64 |   | 22 um Vega magnitude from Bitsakis |
| eW4 | mag | float64 |   | Error in W4 in Vega magnitudes from Bitsakis |
| W1lum | dex(erg / s) | float64 |   | 'Luminosity in W1 from caDistMpc and W1 magnitude using the zero point, frequency, and bandwidth from Jarrett+11' |
| W2lum | dex(erg / s) | float64 |   | 'Luminosity in W2 from caDistMpc and W2 magnitude using the zero point, frequency, and bandwidth from Jarrett+11' |
| W3lum | dex(erg / s) | float64 |   | 'Luminosity in W3 from caDistMpc and W3 magnitude using the zero point, frequency, and bandwidth from Jarrett+11' |
| W4lum | dex(erg / s) | float64 |   | 'Luminosity in W4 from caDistMpc and W4 magnitude using the zero point, frequency, and bandwidth from Jarrett+11' |
| W4SFR | solMass / yr | float64 |   | SFR from W4lum using Catalan-Torrecilla+14 calibration |

edge_coflux_e20.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| coRactr_e20 | deg | float64 | .4f | Reference R.A. of e20 CARMA cube |
| coDectr_e20 | deg | float64 | .4f | Reference Dec. of e20 CARMA cube |
| coCtrint_e20 | K km / s | float64 | .3f | Unmasked CO intensity at reference pixel |
| coBmaj_e20 | arcsec | float64 | .3f | Beam major axis of e20 co cube |
| coBmin_e20 | arcsec | float64 | .3f | Beam minor axis of e20 co cube |
| coBpa_e20 | deg | float64 | .3f | Beam pos ang (deg E of N) of e20 co cube |
| coNomask_e20 | Jy km / s | float64 | .3f | co flux from unmasked e20 cube |
| coeNomask_e20 | Jy km / s | float64 | .3f | co flux uncertainty from unmasked e20 cube |
| coNomaskDv_e20 | km / s | float64 | .3f | co velocity width for unmasked e20 cube |
| coDilated_e20 | Jy km / s | float64 | .3f | co flux from dilated-masked e20 cube |
| coeDilated_e20 | Jy km / s | float64 | .3f | co flux uncertainty from dilated-masked e20 cube |
| coSmooth_e20 | Jy km / s | float64 | .3f | co flux from smoothed-masked e20 cube |
| coeSmooth_e20 | Jy km / s | float64 | .3f | co flux uncertainty from smoothed-masked e20 cube |
| coSmoothDv_e20 | km / s | float64 | .3f | co velocity width for smooth-masked e20 cube |
| coMask2d_e20 | Jy km / s | float64 | .3f | co flux from 2D masked e20 cube |
| coeMask2d_e20 | Jy km / s | float64 | .3f | co flux uncertainty from 2D masked e20 cube |
| coSNRmax_e20 |   | float64 | .3f | Maximum SNR of co in e20 cube |
| coSNR4pix_e20 |   | float64 | .3f | Number of XY pixels with co peak Tb > 4 sigma |
| coSNR5pix_e20 |   | float64 | .3f | Number of XY pixels with co peak Tb > 5 sigma |
| cottBmaj_e20 | arcsec | float64 | .3f | Beam major axis of e20 13co cube |
| cottBmin_e20 | arcsec | float64 | .3f | Beam minor axis of e20 13co cube |
| cottBpa_e20 | deg | float64 | .3f | Beam pos ang (deg E of N) of e20 13co cube |
| cottNomask_e20 | Jy km / s | float64 | .3f | 13co flux from unmasked e20 cube |
| cotteNomask_e20 | Jy km / s | float64 | .3f | 13co flux uncertainty from unmasked e20 cube |
| cottNomaskDv_e20 | km / s | float64 | .3f | 13co velocity width for unmasked e20 cube |
| cottDilated_e20 | Jy km / s | float64 | .3f | 13co flux from dilated-masked e20 cube |
| cotteDilated_e20 | Jy km / s | float64 | .3f | 13co flux uncertainty from dilated-masked e20 cube |
| cottSmooth_e20 | Jy km / s | float64 | .3f | 13co flux from smoothed-masked e20 cube |
| cotteSmooth_e20 | Jy km / s | float64 | .3f | 13co flux uncertainty from smoothed-masked e20 cube |
| cottSmoothDv_e20 | km / s | float64 | .3f | 13co velocity width for smooth-masked e20 cube |
| cottMask2d_e20 | Jy km / s | float64 | .3f | 13co flux from 2D masked e20 cube |
| cotteMask2d_e20 | Jy km / s | float64 | .3f | 13co flux uncertainty from 2D masked e20 cube |
| cottSNRmax_e20 |   | float64 | .3f | Maximum SNR of 13co in e20 cube |
| cottSNR4pix_e20 |   | float64 | .3f | Number of XY pixels with 13co peak Tb > 4 sigma |
| cottSNR5pix_e20 |   | float64 | .3f | Number of XY pixels with 13co peak Tb > 5 sigma |

edge_coflux_iram.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| coRactr_iram | deg | float64 | .4f | Reference R.A. of iram cube |
| coDectr_iram | deg | float64 | .4f | Reference Dec. of iram cube |
| coCtrint_iram | K km / s | float64 | .3f | Unmasked CO intensity at reference pixel |
| coNomask_iram | Jy km / s | float64 | .3f | 12co flux from unmasked iram cube |
| coeNomask_iram | Jy km / s | float64 | .3f | 12co flux uncertainty from unmasked iram cube |
| coNomaskDv_iram | km / s | float64 | .3f | 12co velocity width for unmasked iram cube |
| coDilated_iram | Jy km / s | float64 | .3f | 12co flux from dilated-masked iram cube |
| coeDilated_iram | Jy km / s | float64 | .3f | 12co flux uncertainty from dilated-masked iram cube |
| coSmooth_iram | Jy km / s | float64 | .3f | 12co flux from smoothed-masked iram cube |
| coeSmooth_iram | Jy km / s | float64 | .3f | 12co flux uncertainty from smoothed-masked iram cube |
| coSmoothDv_iram | km / s | float64 | .3f | 12co velocity width for smooth-masked iram cube |
| coMask2d_iram | Jy km / s | float64 | .3f | 12co flux from 2D masked iram cube |
| coeMask2d_iram | Jy km / s | float64 | .3f | 12co flux uncertainty from 2D masked iram cube |
| coSNRmax_iram |   | float64 | .3f | Maximum SNR of 12co in iram cube |
| coSNR4pix_iram |   | float64 | .3f | Number of XY pixels with 12co peak Tb > 4 sigma |
| coSNR5pix_iram |   | float64 | .3f | Number of XY pixels with 12co peak Tb > 5 sigma |

edge_coflux_natv.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| coRactr_natv | deg | float64 | .4f | Reference R.A. of natv CARMA cube |
| coDectr_natv | deg | float64 | .4f | Reference Dec. of natv CARMA cube |
| coCtrint_natv | K km / s | float64 | .3f | Unmasked CO intensity at reference pixel |
| coDvhel | km / s | float64 | .3f | Add this to go from LSR to Barycentric frame |
| coBmaj_natv | arcsec | float64 | .3f | Beam major axis of natv co cube |
| coBmin_natv | arcsec | float64 | .3f | Beam minor axis of natv co cube |
| coBpa_natv | deg | float64 | .3f | Beam pos ang (deg E of N) of natv co cube |
| coNomask_natv | Jy km / s | float64 | .3f | co flux from unmasked natv cube |
| coeNomask_natv | Jy km / s | float64 | .3f | co flux uncertainty from unmasked natv cube |
| coNomaskDv_natv | km / s | float64 | .3f | co velocity width for unmasked natv cube |
| coDilated_natv | Jy km / s | float64 | .3f | co flux from dilated-masked natv cube |
| coeDilated_natv | Jy km / s | float64 | .3f | co flux uncertainty from dilated-masked natv cube |
| coSmooth_natv | Jy km / s | float64 | .3f | co flux from smoothed-masked natv cube |
| coeSmooth_natv | Jy km / s | float64 | .3f | co flux uncertainty from smoothed-masked natv cube |
| coSmoothDv_natv | km / s | float64 | .3f | co velocity width for smooth-masked natv cube |
| coMask2d_natv | Jy km / s | float64 | .3f | co flux from 2D masked natv cube |
| coeMask2d_natv | Jy km / s | float64 | .3f | co flux uncertainty from 2D masked natv cube |
| coSNRmax_natv |   | float64 | .3f | Maximum SNR of co in natv cube |
| coSNR4pix_natv |   | float64 | .3f | Number of XY pixels with co peak Tb > 4 sigma |
| coSNR5pix_natv |   | float64 | .3f | Number of XY pixels with co peak Tb > 5 sigma |
| cottBmaj_natv | arcsec | float64 | .3f | Beam major axis of natv 13co cube |
| cottBmin_natv | arcsec | float64 | .3f | Beam minor axis of natv 13co cube |
| cottBpa_natv | deg | float64 | .3f | Beam pos ang (deg E of N) of natv 13co cube |
| cottNomask_natv | Jy km / s | float64 | .3f | 13co flux from unmasked natv cube |
| cotteNomask_natv | Jy km / s | float64 | .3f | 13co flux uncertainty from unmasked natv cube |
| cottNomaskDv_natv | km / s | float64 | .3f | 13co velocity width for unmasked natv cube |
| cottDilated_natv | Jy km / s | float64 | .3f | 13co flux from dilated-masked natv cube |
| cotteDilated_natv | Jy km / s | float64 | .3f | 13co flux uncertainty from dilated-masked natv cube |
| cottSmooth_natv | Jy km / s | float64 | .3f | 13co flux from smoothed-masked natv cube |
| cotteSmooth_natv | Jy km / s | float64 | .3f | 13co flux uncertainty from smoothed-masked natv cube |
| cottSmoothDv_natv | km / s | float64 | .3f | 13co velocity width for smooth-masked natv cube |
| cottMask2d_natv | Jy km / s | float64 | .3f | 13co flux from 2D masked natv cube |
| cotteMask2d_natv | Jy km / s | float64 | .3f | 13co flux uncertainty from 2D masked natv cube |
| cottSNRmax_natv |   | float64 | .3f | Maximum SNR of 13co in natv cube |
| cottSNR4pix_natv |   | float64 | .3f | Number of XY pixels with 13co peak Tb > 4 sigma |
| cottSNR5pix_natv |   | float64 | .3f | Number of XY pixels with 13co peak Tb > 5 sigma |

edge_coflux_smo6.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| coRactr_smo6 | deg | float64 | .4f | Reference R.A. of smo6 CARMA cube |
| coDectr_smo6 | deg | float64 | .4f | Reference Dec. of smo6 CARMA cube |
| coCtrint_smo6 | K km / s | float64 | .3f | Unmasked CO intensity at reference pixel |
| coNomask_smo6 | Jy km / s | float64 | .3f | co flux from unmasked smo6 cube |
| coeNomask_smo6 | Jy km / s | float64 | .3f | co flux uncertainty from unmasked smo6 cube |
| coNomaskDv_smo6 | km / s | float64 | .3f | co velocity width for unmasked smo6 cube |
| coDilated_smo6 | Jy km / s | float64 | .3f | co flux from dilated-masked smo6 cube |
| coeDilated_smo6 | Jy km / s | float64 | .3f | co flux uncertainty from dilated-masked smo6 cube |
| coSmooth_smo6 | Jy km / s | float64 | .3f | co flux from smoothed-masked smo6 cube |
| coeSmooth_smo6 | Jy km / s | float64 | .3f | co flux uncertainty from smoothed-masked smo6 cube |
| coSmoothDv_smo6 | km / s | float64 | .3f | co velocity width for smooth-masked smo6 cube |
| coMask2d_smo6 | Jy km / s | float64 | .3f | co flux from 2D masked smo6 cube |
| coeMask2d_smo6 | Jy km / s | float64 | .3f | co flux uncertainty from 2D masked smo6 cube |
| coSNRmax_smo6 |   | float64 | .3f | Maximum SNR of co in smo6 cube |
| coSNR4pix_smo6 |   | float64 | .3f | Number of XY pixels with co peak Tb > 4 sigma |
| coSNR5pix_smo6 |   | float64 | .3f | Number of XY pixels with co peak Tb > 5 sigma |
| cottNomask_smo6 | Jy km / s | float64 | .3f | 13co flux from unmasked smo6 cube |
| cotteNomask_smo6 | Jy km / s | float64 | .3f | 13co flux uncertainty from unmasked smo6 cube |
| cottNomaskDv_smo6 | km / s | float64 | .3f | 13co velocity width for unmasked smo6 cube |
| cottDilated_smo6 | Jy km / s | float64 | .3f | 13co flux from dilated-masked smo6 cube |
| cotteDilated_smo6 | Jy km / s | float64 | .3f | 13co flux uncertainty from dilated-masked smo6 cube |
| cottSmooth_smo6 | Jy km / s | float64 | .3f | 13co flux from smoothed-masked smo6 cube |
| cotteSmooth_smo6 | Jy km / s | float64 | .3f | 13co flux uncertainty from smoothed-masked smo6 cube |
| cottSmoothDv_smo6 | km / s | float64 | .3f | 13co velocity width for smooth-masked smo6 cube |
| cottMk12_di_smo6 | Jy km / s | float64 | .3f | 13co flux in 12co dilated mask smo6 cube |
| cotteMk12_di_smo6 | Jy km / s | float64 | .3f | 13co flux uncertainty in 12co dilated mask smo6 cube |
| cottMk12_sm_smo6 | Jy km / s | float64 | .3f | 13co flux in 12co smoothed mask smo6 cube |
| cotteMk12_sm_smo6 | Jy km / s | float64 | .3f | 13co flux uncertainty in 12co smoothed mask smo6 cube |
| cottMk12_2d_smo6 | Jy km / s | float64 | .3f | 13co flux in 12co projected mask smo6 cube |
| cotteMk12_2d_smo6 | Jy km / s | float64 | .3f | 13co flux uncertainty in 12co projected mask smo6 cube |
| cottSNRmax_smo6 |   | float64 | .3f | Maximum SNR of 13co in smo6 cube |
| cottSNR4pix_smo6 |   | float64 | .3f | Number of XY pixels with 13co peak Tb > 4 sigma |
| cottSNR5pix_smo6 |   | float64 | .3f | Number of XY pixels with 13co peak Tb > 5 sigma |

edge_coflux_smo7.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| coRactr_smo7 | deg | float64 | .4f | Reference R.A. of smo7 CARMA cube |
| coDectr_smo7 | deg | float64 | .4f | Reference Dec. of smo7 CARMA cube |
| coCtrint_smo7 | K km / s | float64 | .3f | Unmasked CO intensity at reference pixel |
| coNomask_smo7 | Jy km / s | float64 | .3f | co flux from unmasked smo7 cube |
| coeNomask_smo7 | Jy km / s | float64 | .3f | co flux uncertainty from unmasked smo7 cube |
| coNomaskDv_smo7 | km / s | float64 | .3f | co velocity width for unmasked smo7 cube |
| coDilated_smo7 | Jy km / s | float64 | .3f | co flux from dilated-masked smo7 cube |
| coeDilated_smo7 | Jy km / s | float64 | .3f | co flux uncertainty from dilated-masked smo7 cube |
| coSmooth_smo7 | Jy km / s | float64 | .3f | co flux from smoothed-masked smo7 cube |
| coeSmooth_smo7 | Jy km / s | float64 | .3f | co flux uncertainty from smoothed-masked smo7 cube |
| coSmoothDv_smo7 | km / s | float64 | .3f | co velocity width for smooth-masked smo7 cube |
| coMask2d_smo7 | Jy km / s | float64 | .3f | co flux from 2D masked smo7 cube |
| coeMask2d_smo7 | Jy km / s | float64 | .3f | co flux uncertainty from 2D masked smo7 cube |
| coSNRmax_smo7 |   | float64 | .3f | Maximum SNR of co in smo7 cube |
| coSNR4pix_smo7 |   | float64 | .3f | Number of XY pixels with co peak Tb > 4 sigma |
| coSNR5pix_smo7 |   | float64 | .3f | Number of XY pixels with co peak Tb > 5 sigma |
| cottNomask_smo7 | Jy km / s | float64 | .3f | 13co flux from unmasked smo7 cube |
| cotteNomask_smo7 | Jy km / s | float64 | .3f | 13co flux uncertainty from unmasked smo7 cube |
| cottNomaskDv_smo7 | km / s | float64 | .3f | 13co velocity width for unmasked smo7 cube |
| cottDilated_smo7 | Jy km / s | float64 | .3f | 13co flux from dilated-masked smo7 cube |
| cotteDilated_smo7 | Jy km / s | float64 | .3f | 13co flux uncertainty from dilated-masked smo7 cube |
| cottSmooth_smo7 | Jy km / s | float64 | .3f | 13co flux from smoothed-masked smo7 cube |
| cotteSmooth_smo7 | Jy km / s | float64 | .3f | 13co flux uncertainty from smoothed-masked smo7 cube |
| cottSmoothDv_smo7 | km / s | float64 | .3f | 13co velocity width for smooth-masked smo7 cube |
| cottMk12_di_smo7 | Jy km / s | float64 | .3f | 13co flux in 12co dilated mask smo7 cube |
| cotteMk12_di_smo7 | Jy km / s | float64 | .3f | 13co flux uncertainty in 12co dilated mask smo7 cube |
| cottMk12_sm_smo7 | Jy km / s | float64 | .3f | 13co flux in 12co smoothed mask smo7 cube |
| cotteMk12_sm_smo7 | Jy km / s | float64 | .3f | 13co flux uncertainty in 12co smoothed mask smo7 cube |
| cottMk12_2d_smo7 | Jy km / s | float64 | .3f | 13co flux in 12co projected mask smo7 cube |
| cotteMk12_2d_smo7 | Jy km / s | float64 | .3f | 13co flux uncertainty in 12co projected mask smo7 cube |
| cottSNRmax_smo7 |   | float64 | .3f | Maximum SNR of 13co in smo7 cube |
| cottSNR4pix_smo7 |   | float64 | .3f | Number of XY pixels with 13co peak Tb > 4 sigma |
| cottSNR5pix_smo7 |   | float64 | .3f | Number of XY pixels with 13co peak Tb > 5 sigma |

bb_natv_fitvd_bbmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_natv_fitvd_bbmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_natv_fitvd_dilmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_natv_fitvd_dilmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_natv_fixvd_bbmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_natv_fixvd_bbmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_natv_fixvd_dilmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_natv_fixvd_dilmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_smo5_fitvd_bbmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_smo5_fitvd_bbmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_smo5_fitvd_dilmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_smo5_fitvd_dilmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_smo5_fixvd_bbmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_smo5_fixvd_bbmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_smo5_fixvd_dilmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_smo5_fixvd_dilmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_smo7_fitvd_bbmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_smo7_fitvd_bbmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_smo7_fitvd_dilmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_smo7_fitvd_dilmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbVdisp_e1 | km / s | float64 |   |   |
| bbVdisp_e2 | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_smo7_fixvd_bbmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_smo7_fixvd_bbmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

bb_smo7_fixvd_dilmsk.csv

''

| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbPA | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |
| bbNpix |   | int64 |   | number of pixels in ring |
| bbIntens | (Jy*km/s)/arcsec**2 | float64 |   | 'average intensity in ring, not corrected for inclination' |
| bbIntensRMS | (Jy*km/s)/arcsec**2 | float64 |   | standard deviation of intensity in ring |

bb_smo7_fixvd_dilmsk_freepa.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| bbName |   | string |   | Galaxy Name |
| bbRad | arcsec | float64 |   | Galactocentric radius of ring |
| bbVrot | km / s | float64 |   |   |
| bbVrot_e1 | km / s | float64 |   |   |
| bbVrot_e2 | km / s | float64 |   |   |
| bbVdisp | km / s | float64 |   |   |
| bbXpos | pix | float64 |   |   |
| bbYpos | pix | float64 |   |   |
| bbVmean | km / s | float64 |   | 'Systemic velocity in radio defn, LSR frame' |
| bbVmean_e1 | km / s | float64 |   |   |
| bbVmean_e2 | km / s | float64 |   |   |
| bbPA | deg | float64 |   |   |
| bbPA_e1 | deg | float64 |   |   |
| bbPA_e2 | deg | float64 |   |   |
| bbInc | deg | float64 |   |   |
| bbZ0 | arcsec | float64 |   | Scale height of disk in arcsec |
| bbVrotfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVdispfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbXposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbYposfit | pix | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbVsysfit | km / s | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbPAfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbIncfit | deg | int64 |   | 0=fixed to input value; 1=fitted per ring; -1=fixed to mean value of previous fit |
| bbFtype |   | int64 |   |   |
| bbWfunc |   | int64 |   |   |
| bbFlag |   | int64 |   | '=1 when some ringlog exits, but Bbarolo fails to produce at least one of the plots' |

rprof_de20_smo.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| r_arcs | arcsec | float64 |   | outer radius of ring in arcsec |
| r_kpc | kpc | float64 | 10.3f | outer radius of ring in kpc |
| r_r25 | '' | float64 | 10.3f | outer radius of ring as frac of opt radius |
| npix |   | int64 |   | number of available pixels in the ring |
| ngood |   | int64 |   | number of unmasked pixels in the ring |
| ico_avg | K km / s | float64 | 10.3f | average face-on intensity in ring if ngood/npix > 0.1 and setting blanks to zero |
| ico_rms | K km / s | float32 |   | rms face-on intensity in ring for unmasked pixels |
| detlim | K km / s | float64 | 10.3f | 3 sigma detection limit based on emom0max |
| sigmol | solMass / pc2 | float64 | 10.3f | average face-on surface density including He with alphaco=4.3 |
| cumlum | K km pc2 / s | float32 | 10.4e | total CO luminosity within the given radius |
| cummass | solMass | float32 | 10.4e | total molecular gas mass within the given radius |

rprof_smo7_smo.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| r_arcs | arcsec | float64 |   | outer radius of ring in arcsec |
| r_kpc | kpc | float64 | 10.3f | outer radius of ring in kpc |
| r_r25 | '' | float64 | 10.3f | outer radius of ring as frac of opt radius |
| npix |   | int64 |   | number of available pixels in the ring |
| ngood |   | int64 |   | number of unmasked pixels in the ring |
| ico_avg | K km / s | float64 | 10.3f | average face-on intensity in ring if ngood/npix > 0.1 and setting blanks to zero |
| ico_rms | K km / s | float32 |   | rms face-on intensity in ring for unmasked pixels |
| detlim | K km / s | float64 | 10.3f | 3 sigma detection limit based on emom0max |
| sigmol | solMass / pc2 | float64 | 10.3f | average face-on surface density including He with alphaco=4.3 |
| cumlum | K km pc2 / s | float32 | 10.4e | total CO luminosity within the given radius |
| cummass | solMass | float32 | 10.4e | total molecular gas mass within the given radius |

cospec_iram.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| coVlsr_iram | km / s | float64 |   |   |
| coNomaskSpec_iram | Jy | float64 | .5f | channel-averaged flux |
| coNomaskUnc_iram | Jy | float64 | .5f | formal uncertainty |
| coDilatedSpec_iram | Jy | float64 | .5f | channel-averaged flux |
| coDilatedUnc_iram | Jy | float64 | .5f | formal uncertainty |
| coSmoothSpec_iram | Jy | float64 | .5f | channel-averaged flux |
| coSmoothUnc_iram | Jy | float64 | .5f | formal uncertainty |
| coMask2dSpec_iram | Jy | float64 | .5f | channel-averaged flux |
| coMask2dUnc_iram | Jy | float64 | .5f | formal uncertainty |

cospec_smo7.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| coVlsr_smo7 | km / s | float64 |   |   |
| coNomaskSpec_smo7 | Jy | float64 | .5f | channel-averaged flux |
| coNomaskUnc_smo7 | Jy | float64 | .5f | formal uncertainty |
| coDilatedSpec_smo7 | Jy | float64 | .5f | channel-averaged flux |
| coDilatedUnc_smo7 | Jy | float64 | .5f | formal uncertainty |
| coSmoothSpec_smo7 | Jy | float64 | .5f | channel-averaged flux |
| coSmoothUnc_smo7 | Jy | float64 | .5f | formal uncertainty |
| coMask2dSpec_smo7 | Jy | float64 | .5f | channel-averaged flux |
| coMask2dUnc_smo7 | Jy | float64 | .5f | formal uncertainty |

edge_hispectra.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| source |   | string |   | Origin of spectrum |
| vbary | km / s | float64 |   | Optical barycentric velocity |
| hiflux | mJy | float64 |   | HI flux |

edgealma_hispectra.csv



| name | unit | datatype | format | description |
|---|---|---|---|---|
| Name |   | string |   | Galaxy Name |
| source |   | string |   | Origin of spectrum |
| vbary | km / s | float64 |   | Optical barycentric velocity |
| hiflux | mJy | float64 |   | HI flux |
