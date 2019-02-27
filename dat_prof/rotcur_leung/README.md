Note: Rotation curves published in Leung et al., http://adsabs.harvard.edu/doi/10.1093/mnras/sty288


From: Gigi Leung  
Subject: Re: Revised manuscript and reply to referee (EDGE: comparison of stellar and CO V_circ)  
Date: May 18, 2018 at 4:46:10 PM CDT  
To: "Wong, Tony"

Hi Tony,

Attached please find the rotation curve values I used in my paper. There are in total 54 galaxies. Each file has three columns: R (arcsec), vrot (km/s), vrot error (km/s).

The files in the folder 'EDGE_CO_vrot_bsc' are beam-smearing-corrected, while the ones in EDGE_CO_vrot are not.

Let me know if anything is unclear or if I can help in anyways.

Cheers,
Gigi


> On Feb 7, 2018, at 8:54 PM, Wong, Tony wrote:
> 
> Gigi,
> 
> Thanks for your reply!  Text or CSV files are fine; again it would be useful to have the profiles both with and without the beam smearing correction.  Also a machine-readable version of your Table 1 would be useful.
> 
> Cheers,
> Tony
> 
>> On Feb 5, 2018, at 11:14 AM, Gigi Leung wrote:
>> 
>> Hi Tony,
>> 
>> Yes, the moment maps I used are the Gaussian fitted ones provided by Alberto and Rebecca. I then convert the values to the optical convention by myself.
>> 
>> For PA and inclinations, for most galaxies I simply adopt the values in the CALIFA master table and they are indeed derived from photometry. There are a number of galaxies as marked in table 1 in my paper which I fit the PA kinematically by myself. In the first step, each ring is allowed to have its own PA, and central position. Then I fix the kinematic centre as the average of the central positions of all rings and refit the kinematics, again allowing each ring to have its own PA. And during the process the location of the rings (meaning which pixels belong to which rings and the radii of the individual rings) can change. These first 2 iterations are simple fitting of sin θ around the ring. In the third iteration, PA is also fixed (as the average of all the rings from the 2nd iteration) and I now apply a harmonic decomposition.
>> 
>> In addition to the 20 pixel requirement I also for the latest version imposed an additional requirement: the thickness of the ring has to be at least 1/2 of the beam size. And the outer rings are cutoff sometimes because the data become patchy there or the signal to noise become low. 
>> 
>> The rotation curves and dispersion profiles are now just in text files, written in 3 columns: [radius (arcsec), Vrot, error]. Which format would you prefer to have it in for the data archive?
>> 
>> Cheers.
>> Gigi
>> 
>> 
>> 
>> On Thu, Jan 25, 2018 at 6:00 AM, Wong, Tony wrote:
>>
>> Gigi,
>> 
>> I am sorry it has taken so long to get back to you, and perhaps you have already resubmitted the paper.  It looks fine to me, and mostly I am keen to make sure that I understand how the CO data have been handled.  I assume you are using the Gaussian fit velocity maps from Alberto and Becca?  Did you then shifted the velocity values to the optical convention?
>> 
>> I also would like to understand a bit more the rotation curve fitting process.  When you say you fit ellipses to the mean velocity map, I guess you start by assuming a PA and INC from photometry, and then extract the mean velocities in a ring and use the sin θ component to determine an improved estimate of PA?  Is this then averaged across all rings to come up with a single value for PA_kin?  Also do you go through multiple iterations because when you change the PA you also change the location of the rings?  Is the requirement of 20 pixels per annulus adequate to decide how many rings per galaxy, or do you sometimes have to apply an outer radius cutoff because the fraction of a ring with good data becomes too small?
>> 
>> Let me know also what kind of format the rotation curves and velocity dispersion profiles are in, and how you suggest we integrate them into the data archive.  It would probably be useful to archive both the original and the beam smearing corrected curves.
>> 
>> Best,
>> Tony
>> 
