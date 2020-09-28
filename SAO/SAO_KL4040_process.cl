
noao
imred
ccdred

#SIMPLE  =                    T / C# FITS: 12/11/2019 10:39:23 PM 
#BITPIX  =                   16              
#NAXIS   =                    2 / Dimensionality  
#NAXIS1  =                 4096                    
#NAXIS2  =                 4096                    
#XBINNING=                    1 / Binning factor used on the X axis    
#EXPTIME =             11627906 / Duration of exposure in 10us periods 
#IMGSTPCL=                 4095 / Image Stop Column                        
#IMGSTPRW=                 4095 / Image Stop Row  
#FIRMVER =                    0 / Firmware Version
#ISSFTBIN=                    F / Is Software Binning flag
#ISHIFRM =                    F / Is High Frame flag 
#YBINNING=                    1 / Binning factor used on the Y axis 
#ISSTKFRM=                    F / Is Stacked Frame flag   
#USSFTAVE=                    F / Shifted Averaging used flag     
#GPS_LOCK=                    F / GPS Present at Capture Time
#BKSDILL =                    F / Back Side Illuminated      
#ISMRGFRM=                    T / Is Merged Frame flag           
#XORGSUBF=                    0 / Subframe origin on X axis  
#CCD-TEMP=                  -10 / CCD Temperature (C)            
#FPGATEMP=                  314 / FPGA Temperature (C)     
#SENSTEMP= '-1438.2 '           / Sensor Temperature (C)  
#CAM-MODE= 'Unknown(2)'         / Camera Mode  
#IMGTYPE =                    5 / The Image Type (ex. High Gain) 
#GAIN    =                 16.5 / Gain Value         
#SERIAL  = 'KL1914719'          / Camera serial number  
#DATE-OBS= '2019-12-11T22:39:21.873 NOGPS' / Date of observation (capture)
#DATE    = '2019-12-11T13:39:23.169000' / File create date 
#GEO_LAT =                    0 / [deg] Geocentric latitude  
#GEO_LONG=                    0 / [deg] Geocentric longitude  
#PLOTVER = '1.2.31.2'           / FliPilot Version   
#FLIBVER = '1.11.24 '           / LibFliPro Version 
#RECTARR =                    T / Rectangular Array 
#YORGSUBF=                    0 / Subframe origin on Y axis 
#BASETEMP=                   76 / Heatsink Temperature (C) 
#PREREFRW=                    0 / Number of Pre-Reference Rows
#PSTREFRW=                    0 / Number of Post-Reference Rows
#EXTEND  =                    T / Extensions are permitted  
#SIGNED  =                    F / False if unsigned array type
#METAVER =                    2 / Metadata Version 
#INSTRUME= 'KEPLER KL4040'      / Camera Model  
#SNRPXLDP=                   12 / Sensor Pixel Depth capability    
#RFPXLPRW=                    0 / Number of Post Reference pixels per row 
#PRFPXPRW=                    0 / Number of Pre Reference pixels per row 
#DEADPXCR=                    F / Dead Pixel Correction flag  
#PXLORDIM=                    T / Pixel Ordered Image flag  
#HDR-MODE=                    T / Image captured in HDR Mode   
#LOWNOISE=                    F / Low Noise flag  
#LOWDARK =                    T / Low Dark Current flag   
#HRZDIRIN=                    F / Horizontal Scan Direction Invert flag 
#VRTDIRIN=                    F / Vertical Scan Direction Invert flag 
#NONRWALL=                    F / Non-Row Alligned Image flag  
#BLKLVLAD=                15672 / Black Level Adjust value 
#BLKLVLSN=                   20 / Black Level Sun value  
#SHUTOPDL=                  150 / Shutter Open Delay 
#SHUTCLDL=                    0 / Shutter Close Delay  
#ILLSTRDL=                    0 / Illumination Start Delay  
#ILLSTPDL=                    0 / Illumination Stop Delay   
#TRKSTRRW=                    0 / Tracking Start Row   
#TRKSTRCL=                    0 / Tracking Start Column  
#TRKSTPRW=                    0 / Tracking Stop Row   
#TRKSTPCL=                    0 / Tracking Stop Column
#TRKFRMPR=                    0 / Tracking Frames Per Image Frame
#FRAMENUM=                    1 / Frame Number 
#STRTEXRW=                    0 / Starting Exposure Row 
#BSCALE  =                    1 / Linear factor in scaling equation 
#NUMDTCH =                    8 / Number Of Data Channels   
#BZERO   =                32768 / Zero point in scaling equation
#END    


# process for the one night

pwd
#!gethead *.fit FILTER OBJECT IMAGETYP EXPTIME UTDATE UTSTART RA DEC > fileinfo.txt

# bias combine


ls BIAS*.fit > zero.list
#!gethead BIAS*.fit FILTER OBJECT IMAGETYP EXPTIME

zerocom @zero.list output=zero.fits combine=median reject=crreject ccdtype='' scale=none gain=16.5 rdnoise=3.7

imstat BIAS*.fit
imstat zero.fits

# dark combine
ls DK1_*.fit > dark1.list
ls DK10*.fit > dark10.list
ls DK30*.fit > dark30.list
ls DK120*.fit > dark120.list
ls DK180*.fit > dark180.list
ls DARK120*.fit > dark120.list

# dark - bias
imarith @dark1.list - zero.fits z@dark1.list
imarith @dark10.list - zero.fits z@dark10.list
imarith @dark30.list - zero.fits z@dark30.list
imarith @dark120.list - zero.fits z@dark120.list
imarith @dark180.list - zero.fits z@dark180.list



darkcom z@dark1.list output=dark1.fits combine=median reject=crreject ccdtype='' scale=none gain=16.5 rdnoise=3.7
darkcom z@dark10.list output=dark10.fits combine=median reject=crreject ccdtype='' scale=none gain=16.5 rdnoise=3.7
darkcom z@dark30.list output=dark30.fits combine=median reject=crreject ccdtype='' scale=none gain=16.5 rdnoise=3.7
darkcom z@dark120.list output=dark120.fits combine=median reject=crreject ccdtype='' scale=none gain=16.5 rdnoise=3.7
darkcom z@dark180.list output=dark180.fits combine=median reject=crreject ccdtype='' scale=none gain=16.5 rdnoise=3.7
#darkcom @dark.list output=dark.fits combine=median reject=crreject ccdtype='' scale=none gain=16.5 rdnoise=3.7

#flat combine

ls DFlat*.fit > flat.list
ls DFlat_Merge*.fit > flat_Merge.list
#hselect @flat.list $I "filter= 'R(Ha)'" > raflat.list ; type raflat.list

imstat @flat_Merge.list
imarith @flat_Merge.list - zero.fits z@flat_Merge.list
imarith z@flat_Merge.list - dark1.fits dz@flat_Merge.list
imstat d@flat.list

flatcom dz@flat_Merge.list output=flat_Merge.fits combine=median reject=minmax ccdtype='' process- gain=16.5 rdnoise=3.7 scale=mode subset- ccdtype='' 

imstat flat_Merge.fits

imstat flat_Merge.fits field='mean' format- | scan(x)
=x
imarith flat_Merge.fits / (x) nflat.fits



