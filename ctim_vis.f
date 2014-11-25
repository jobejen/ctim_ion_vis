c.. this is an example program to read the MHD simulation data
c   and interpolate the to an arbitrary point within the simulation
c   box.
c     adpted from http://openggcm.sr.unh.edu/wiki/images/d/dc/Ioread.txt
c     
c      
c.. What you really need to know:
c
c  0.) The coordinate system is 'inverted GSE', i.e., -X_gse,-Y_gse,Z_gse.
c  1.) data are stored in a `compressed ascii' format, and thus are basically
c      machine independent (except for vintage mainframes). The 8'th bit
c      is used, thus `ftp' transfers should be done in `binary' mode since
c      some ftp clients strip the 8'th bit. However, there are no control
c      characters ( ascii 0 -- 31 ) in the files, except for line feed.
c  2.) open files as 'formatted'.
c  3.) read a 1d field with:
c      call getf11(iunit,a,n,l1,l2,it)
c      iunit (integer, input):            fortran file unit
c      a (real, output):                  field
c      n (integer, output):               size of a
c      l1 (character*80, input/output):   field identifier
c      l2 (character*80, output):         character tag
c      it (integer, output):              integer tag
c  4.) note:
c      a) you dont need to know the filed size, it is returned in `n',
c         however, `a' must be dimensioned sufficiently large.
c      b) `l1' and `l2' MUST be declared `character*80' in the calling 
c         program.
c      c) if `l1' is set to 'any', then the next field will be read.
c      d) if `l1' is set to a known identifier (like 'gridx'), the routine
c         will skip if necessary to that field and read it.
c      e) on `EOF' or any read error `n' returns less than zero.
c      f) `l2' usually contains some additional information, like the UT time.
c      g) `it' usually contains a value of seconds since the start of the run.
c  5.) read 2d and 3d fields with:
c      call getf21(iunit,a,nx,ny,l1,l2,it) and
c      call getf31(iunit,a,nx,ny,nz,l1,l2,it) , respectively.
c      Everything works the same way as with 1d fields, except that the multiple
c      field dimensions (nx,ny) or (nx,ny,nz) are returned.
c  6.) The returned values in `a' are of limited precision, about 3.e-4 in the
c      worst case.  That is good enough for virtually any purpose, because 
c      interpolation errors are generally larger.  `a' also contains no true
c      zeroes.  Zero values are returned as the smallest possible positive
c      value.
c  7.) The routines have been used on Intel-Linux boxes (pgf77 or fort77), SunOS,
c      Solaris, IBM RS6000 (AIX), CRAY-C90, DEC/Ultrix, CRAY-T3E, INTEL-Paragon,
c      and should run without problems on any UNIX workstation.
c      And, yes, GNU g77 works too.
c      I dont know about Microsoft, Apple, and what else there might exist.
c  8.) Corresponding write routines (datf11, datf21, datf31) are included.  Their
c      use should be straightforward.
c      
c.. What you may want to know:
c
c 1.)  Each field has a 5 line header that stores some pertinent
c      information (magic number, dimensions, etc.).
c 2.)  The field of values `a' is cropped to a dynamic range (ratio of
c      maximum to minimum abs value) of no more than 1.0e16, then is chopped 
c      into chunks of 64.  For each chunk the log of the abs value is taken and
c      linearly discretized in about 4400 steps (that number, with the dynamic range
c      determines the precision).  The resulting integer is factored into 2
c      integers and combined with the sign, resulting in 2 integers in the range
c      0-94. These are converted into ascii characters (33-126), which are then
c      run length encoded and written to the file (2 lines per 64 chunk).
c      The 8'th bit is used to indicate run length counts, thus all stored
c      characters lie in the range (33-126,172-254). In addition, a checksum is 
c      stored for each row and compared upun reading.  Thus, damaged files are 
c      easily detected.  Compression is usually
c      pretty good, `gzipping', i.e., nearly optimal compression of the files,
c      usually yields only 25% further compression.  In the worst case, 2 bytes
c      per value are used, typical is 1-1.5 bytes per value, and constant or nearly
c      constant fields compress to 0.1 bytes per value or better.
c 3.)  Using the software with C may be possible in at least 3 ways:
c      a) Compiling the routines with a f77 compiler, compiling the C routines
c      (including main()) with a C compiler and linking with the C compiler.
c      In this case, the proper calling sequence for the fortran routines must be
c      used, i.e., parameters passed by pointers, and character strings pass 2
c      pointers, possibly appending an underscore character to the fortran routine
c      name.  Also, the fortran libraries (something like libU77.a and libF77.a
c      must be provided in the link step.
c      b) Using f2c (from ftp.netlib.org) to convert the f77 code to C code and
c      compiling and linking everything with a C compiler.  In this case the
c      proper f2c libraries must be used because the f77 code uses IO and math
c      functions.
c      c) Rewrite the f77 routines in C.
c      I have not tried either and I would appreciate to hear about the results
c      if anyone tries out the above steps.
c 4.)  Things to do:
c      a) write a native C version
c      b) use different encoding (non-log) for fields with limited range, like
c      B-components and V-components.  That should improve compression and allow
c      for better precision.
c
c     -----------------------------------------------------
c     modified by Joseph Jensen July 2014 to read *ioc* files
c     some important information about the ioc files and the grid stsructure
c   21 grids are used longitudinally (from 0 to 360 degrees), which lead to longitudinal resolution
c   of 18 degrees. 91 grids are used between the two poles (from -90 to 90 degrees), which lead to  
c   latitudinal resolution of 2 degrees. Such grids are not necessarily at r=1.0 Re. Actually it    
c   can be at any altitudes, like the potential it follows the pressure gradients. also note that 
c   ctim is in geographic coordinates.  
c  
c..... read fields                                                                                  
c      10 2D parameters                                                                             
c      1  ctim_efdis  21 91                                                                         
c      2  ctim_medis  21 91                                                                         
c      3  ctim_efdif  21 91                                                                         
c      4  ctim_medif  21 91                                                                         
c      5  ctim_pot  21 91                                                                           
c      6  ctim_eden  21 91                                                                          
c      7  ctim_nmf2  21 91                                                                          
c      8  ctim_nmf2h  21 91                                                                         
c      9  ctim_hiped  21 91                                                                         
c      10  ctim_hihall  21 91                                                                       
c
c     this would be an example call strucure for 2d fields
cc      do i=1, 10                                                                                   
cc      l1='any'                                                                                     
cc      call ctim_getf21(10,val,nnx,nny,l1,l2,it)                                                    
cc      write(0,*)'c     ',i,'  ',l1(1:jlen(l1,80)),' ',nnx,nny                                      
cc            write(0,*)i, val(1), val(15)                                                           
cc      enddo                                                                                        
c                                                                                                   
c      19 3D parameters                                                                             
c 1  CTIM_DYN_3D_JX_3D_GEO  15 91 20 :                                                              
c 2  CTIM_DYN_3D_JY_3D_GEO  15 91 20 :                                                              
c 3  CTIM_NDENS_3D_GEO  15 91 20 : atmosphere density                                               
c 4  CTIM_VX_3D_GEO  15 91 20 : neutral wind velocity in theta direction                            
c 5  CTIM_VY_3D_GEO  15 91 20 : neutral wind velocity in phi direction                              
c 6  CTIM_WVZ_3D_GEO  15 91 20 : neutral wind vertical velocity                                     
c 7  CTIM_EPS_3D_GEO  15 91 20 : internal energy density                                            
c 8  CTIM_HT_3D_GEO  15 91 20 : height                                                              
c 9  CTIM_RMT_3D_GEO  15 91 20 : atmosphere molecular weight                                        
c 10  CTIM_OM1_3D_GEO  15 91 20 : dP/dt, the vertical wind in the pressure coordinate system        
c 11  CTIM_ALDEN_3D_GEO  15 91 20 : number density of electrons                                     
c 12  CTIM_TTS_3D_GEO  15 91 20 :                                                                   
c 13  CTIM_PSAO_3D_GEO  15 91 20 : mass percentage of atomic oxygen in the atmosphere               
c 14  CTIM_PSMO_3D_GEO  15 91 20 : mass percentage of molecular oxygen in the atmosphere            
c 15  CTIM_PSMN_3D_GEO  15 91 20 : mass percentage of molecular nitrogen in the atmosphere          
c 16  CTIM_D13D1_3D_GEO  90 91 20                                                                   
c 17  CTIM_D13D2_3D_GEO  90 91 20                                                                   
c 18  CTIM_D23D1_3D_GEO  90 91 20                                                                   
c 19  CTIM_D23D2_3D_GEO  90 91 20  
c
c     an example call structure for 3d files would be like this
cc     l1="any"
cc     call getf31(10,val,nnx,nny,nnz,l1,l2,it)
cc     write(0,*)' got ',l3(1:jlen(l3,80)),' ',nnx,nny,nnz
c
c--------------------------------------------------------
      program iocread
c--------------------------------------------------------
c     uses the script run_all_dens.scr to run multiple times and generate the files
c     This program takes the ioc file and reads density then interpolates to a 
c     satellites position and outputs the data in a form for gnuplot to plot the data.
c     Joseph B Jensen july 2014
c
      real den(15,91,20),ht(15,91,20),denb(15,91,20)
      real zkmheight(90,90,20),oprod_rate_nprec(90,90,20)
      real o2prod_rate_nprec(90,90,20),n2prod_rate_nprec(90,90,20)
      real oprod_rate_yprec(90,90,20),o2prod_rate_yprec(90,90,20)
      real n2prod_rate_yprec(90,90,20)
      real*8 ssecond,seconds,xgeo,ygeo,zgeo,xgse,ygse,zgse,lat,long
      real satx(10),saty(10),satz(10),satxgeo(10),satygeo(10)
      real satzgeo(10),interheight,mindis,val1,val2,avl3,val4
      integer altnum,syear,smonth,sday,shour,sminute,sattime(10)
      integer timenum,timenuminit,nyear,nmonth,nday,nhour,nminute,nsat
      integer pos1,pos2,pos3,pos4,pos5
      character*80 l1(10),fgrid,l2
      character*100 dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
      character*20 dummy,satname(10)
      character*200 f3d,backgroundfile
      character*10,altnumc,chtimenum

!!!!! read in the inputs!!!!!
c     get the file name of the ioc file to read in the data
      call getarg(1,f3d)
      open(10,file=f3d,status='old')
c     get the variables to look at (e.g. density or electron energy)
      do I=1,7
      call getarg(I+1,l1(I))
      enddo
c     get the altitude/pressure height layer (1-15)
      call getarg(9,chtimenum) 
c     convert to an integer
      read (chtimenum,'(I10)') timenum
      timenuminit=timenum
      print*,"timemun= ",timenum
c     read the l1 from the ioc file

      call getf31(10,zkmheight,nnx,nny,nnz,l1(1),l2,it)
      call getf31(10,oprod_rate_nprec,nnx,nny,nnz,l1(2),l2,it)
      call getf31(10,o2prod_rate_nprec,nnx,nny,nnz,l1(3),l2,it)
      call getf31(10,n2prod_rate_nprec,nnx,nny,nnz,l1(4),l2,it)
      call getf31(10,oprod_rate_yprec,nnx,nny,nnz,l1(5),l2,it)
      call getf31(10,o2prod_rate_yprec,nnx,nny,nnz,l1(6),l2,it)
      call getf31(10,n2prod_rate_yprec,nnx,nny,nnz,l1(7),l2,it)


!!!! get the imte right to make transformations
c     we need to initialize the coordinate transformation program by calling ctim_cotr_set
c     this !!!WILL HAVE PROBLEMS!! for ends of months and end of years, I'll get to it when I need it
      syear=2011
      smonth=8
      sday=4
      shour=12
      sminute=00
      ssecond=00.00
      nyear=2011
      nmonth=8
      nday=-1
      nhour=-1
      nminute=-1
      seconds=-1.0
c     take the timenum from the file and divide it up into year,month,day,hour,min,sec
      nday=int((timenum+sday*86400+shour*3600+sminute*60+ssecond)/86400)
      if (nday-sday .eq. 0) then
         nhour=INT((timenum+shour*3600+sminute*60+ssecond)/3600.0) 
         timenum=timenum-3600*(nhour-shour)
         nminute=INT((timenum+sminute*60+ssecond)/60.0)
         timenum=timenum - 60.0 * nminute 
         seconds=timenum
      else
         timenum=timenum-86400*(nday-sday)+3600*shour+60*sminute+ssecond
         nhour=INT((timenum)/3600.0) 
         timenum=timenum-3600*(nhour)
         nminute=INT((timenum)/60.0)
         timenum=timenum - 60.0 * nminute
         seconds=timenum
      endif
c     check the output and make sure it makes sense
      print*,"nday= ",nday
      print*,"nhour= ",nhour
      print*,"nminute= ",nminute
      print*,"second= ",seconds
      call ctim_cotr_set(nyear,nmonth,nday,nhour,nminute,seconds)

!!!!  output the data to plot
c     longitude slice of the ionization rates
      lslice=5
      open(unit=101, file="long_slice_000.txt")
      open(unit=106, file="long_slice_090.txt")
      open(unit=111, file="long_slice_180.txt")
      open(unit=116, file="long_slice_270.txt")
      open(unit=117, file="alt_slice_k-10.txt")
      do I=1,90           !long
         do J=1,20        !lat
            do k=1,90     !height

       if (K .eq. 15)then
c     create xyz geo coordinates
               the=(I-1)*2.*(3.14159/180)
               phi=(mod((((j-1)*18.)+180.0),360.0)-180)*(3.14159/180)
               R=6731+(zkmheight(K,I,J)/1000)
       xgeo=(6731+(zkmheight(K,I,J)/1000))*sin(the)*cos(phi)
       ygeo=(6731+(zkmheight(K,I,J)/1000))*sin(the)*sin(phi)
       zgeo=(6731+(zkmheight(K,I,J)/1000))*cos(the)
c     print out a polar profile

          call ctim_cotr('geo','gse',xgeo,ygeo,zgeo,xgse,ygse,zgse)
          write(117,*)xgeo,ygeo,zgeo,xgse,ygse,zgse,zkmheight(k,I,J),
     1oprod_rate_nprec(k,I,J),o2prod_rate_nprec(k,I,J),
     2n2prod_rate_nprec(k,I,J),oprod_rate_yprec(k,I,J),
     3o2prod_rate_yprec(k,I,J),n2prod_rate_yprec(k,I,J),I,J,K
          if (J .eq. 20) then 
               the=(I-1)*2.*(3.14159/180)
               phi=(mod((((1-1)*18.)+180.0),360.0)-180)*(3.14159/180)
               R=6731+(zkmheight(K,I,1)/1000)
       xgeo=(6731+(zkmheight(K,I,1)/1000))*sin(the)*cos(phi)
       ygeo=(6731+(zkmheight(K,I,1)/1000))*sin(the)*sin(phi)
       zgeo=(6731+(zkmheight(K,I,1)/1000))*cos(the)
             write(117,*)xgeo,ygeo,zgeo,xgse,ygse,zgse,zkmheight(k,I,1),
     1oprod_rate_nprec(k,I,1),o2prod_rate_nprec(k,I,1),
     2n2prod_rate_nprec(k,I,1),oprod_rate_yprec(k,I,1),
     3o2prod_rate_yprec(k,I,1),n2prod_rate_yprec(k,I,1),I,J,K
          endif
       endif
c     now print out lat slices
               if (j .eq. 1) then
      write(101,*)I,J,K,zkmheight(k,I,J),oprod_rate_nprec(k,I,J),
     1o2prod_rate_nprec(k,I,J),n2prod_rate_nprec(k,I,J),
     2oprod_rate_yprec(k,I,J),o2prod_rate_yprec(k,I,J),
     3n2prod_rate_yprec(k,I,J)
              endif
               if (j .eq. 6) then
      write(106,*)I,J,K,zkmheight(k,I,J),oprod_rate_nprec(k,I,J),
     1o2prod_rate_nprec(k,I,J),n2prod_rate_nprec(k,I,J),
     2oprod_rate_yprec(k,I,J),o2prod_rate_yprec(k,I,J),
     3n2prod_rate_yprec(k,I,J)
               endif
               if (j .eq. 11) then
      write(111,*)I,J,K,zkmheight(k,I,J),oprod_rate_nprec(k,I,J),
     1o2prod_rate_nprec(k,I,J),n2prod_rate_nprec(k,I,J),
     2oprod_rate_yprec(k,I,J),o2prod_rate_yprec(k,I,J),
     3n2prod_rate_yprec(k,I,J)
               endif
               if (j .eq. 16) then
      write(116,*)I,J,K,zkmheight(k,I,J),oprod_rate_nprec(k,I,J),
     1o2prod_rate_nprec(k,I,J),n2prod_rate_nprec(k,I,J),
     2oprod_rate_yprec(k,I,J),o2prod_rate_yprec(k,I,J),
     3n2prod_rate_yprec(k,I,J)
               endif
            enddo
            if (j .eq. 1) then
               write(101,*)" "
            endif
            if (j .eq. 6) then
               write(106,*)" "
            endif
            if (j .eq. 11) then
               write(111,*)" "
            endif
            if (j .eq. 16) then
               write(116,*)" "
            endif
         enddo
         write(117,*)" "
      enddo
      j=1
      do I=1,90                 !long                                                           
         do k=1,90              !height                                                          
            if (K .eq. 15)then
c     create xyz geo coordinates                                                                  
            the=(I-1)*2.*(3.14159/180)
            phi=(mod((((j-1)*18.)+180.0),360.0)-180)*(3.14159/180)
            R=6731+(zkmheight(K,I,J)/1000)
            xgeo=(6731+(zkmheight(K,I,J)/1000))*sin(the)*cos(phi)
            ygeo=(6731+(zkmheight(K,I,J)/1000))*sin(the)*sin(phi)
            zgeo=(6731+(zkmheight(K,I,J)/1000))*cos(the)
c     print out a polar profile                                                                     
            
            call ctim_cotr('geo','gse',xgeo,ygeo,zgeo,xgse,ygse,zgse)
c          write(117,*)xgeo,ygeo,zgeo,xgse,ygse,zgse,zkmheight(k,I,J),
c     1oprod_rate_nprec(k,I,J),o2prod_rate_nprec(k,I,J),
c     2n2prod_rate_nprec(k,I,J),oprod_rate_yprec(k,I,J),
c     3o2prod_rate_yprec(k,I,J),n2prod_rate_yprec(k,I,J),I,J,K
            endif
      enddo
      enddo
      



      close(unit=101)
      close(unit=106)
      close(unit=111)
      close(unit=116)

!!!!  now we need to read in the satellite files and get the orbit information !!!!
c      nsat=6   !!! number of satellites in the orbits.txt
c      open(unit=27,file="/home/joe/trillian/storm_2011_08_04_v4/storm-20
c     +14-AUG-05-jr12/orbits.txt", status="old") !input
c      open(unit=28, file="orbits_sattraj_gse.txt")  !output
c      open(unit=29, file="orbits_sattraj_geo.txt")  !output
c      do i=1,nsat
cc     really funny way to read in the orbits.txt file to deal with : seperation
c      do while (sattime(i) .ne. timenuminit)
c         read(27,*) dummy1
c         pos1 = index(dummy1, ":")
c         read(dummy1(pos1+1:), *) dummy2
c         pos2 = index(dummy2, ":")
c         satname(I) = dummy2(1:pos2-1)
c         read(dummy2(pos2+1:), *) dummy3
c         pos3 = index(dummy3, ":")
c         read (dummy3(1:pos3-1),'(I10)') sattime(I)
c         read(dummy3(pos3+1:), *) dummy4
c         pos4 = index(dummy4, ":")
c         read (dummy4(1:pos4-1),'(E16.7)') satx(I)
c         read(dummy4(pos4+1:), *) dummy5
c         pos5 = index(dummy5, ":")
c         read (dummy5(1:pos5-1),'(E16.7)') saty(I)
c         read(dummy5(pos5+1:), *) dummy6
c         read (dummy6,'(E16.7)') satz(I)
c      enddo
c      enddo
cc     convert to from gse to geo to make the interpolation calculations, but keep
cc     the gse info becuase it is easier to graph
c      do i=1,nsat
c         xgse=satx(I)
c         ygse=saty(I)
c         zgse=satz(I)
c         call ctim_cotr('gse','geo',xgse,ygse,zgse,xgeo,ygeo,zgeo) 
c         satxgeo(I)=xgeo
c         satygeo(I)=ygeo
c         satzgeo(I)=zgeo
c         write(29,*)satname(i),xgeo*6731,ygeo*6731,zgeo*6731,-0.0  !write sats position geo
c       write(28,*)satname(i),satx(i)*6731,saty(i)*6731,satz(i)*6731,-0.0 !write sats position gse
c      enddo
c      goto 104
c 103  print*,"end of sat file reading, continuing" !!! the end statement is needed if the sat dat isn't long enough
c      write(29,*)-9999,0.0,0.0,0.0  !!write dummy files for the satellites to remove from graph
c      write(28,*)-9999,0.0,0.0,0.0 
c 104  continue

!!!!  interpolating to a satellites coordinates !!!!
c     satellites position satname(i),satx(i)*6731,saty(i)*6731,satz(i)*6731
c     do the interpolation in GEO coordinates because it is easier (convert from gse to geoxyz)
c     !!!!WILL PROBABLY BREAK IF SATELLITE IS ABOVE ~500km!!!!
c      do I=1,nsat   
c      if (I .eq. 5 .or. I .eq. 6) then ! this is to look only at grace1 and 2 satellites, dmsp are too high
c      alt=sqrt(satxgeo(i)*satxgeo(i)+satygeo(i)*satygeo(i)+satzgeo(i)
c     1*satzgeo(i))
c      interheight=(alt-1.0)*6371000.0  !set interpolate height as altitude of sat
cc      interheight=410*1000.  !interpolate heigth as fixed value
c      print*,"alt of the sat: ",interheight
c      lat=acos(satzgeo(i)/alt)        !!theta in radians
c      long=acos(satxgeo(i)/(alt*sin(lat)))       !!phi in radians
c      lat=lat*(180/3.14159)
c      long=long*(180/3.14159)
c         j = (INT(lat) + 90)/2 + 1
c         k = INT(long)/18 + 1
c         print*,"j= ",j
c         print*,"k= ",k
c         print*,"lat= ",lat
c         print*,"long= ",long
c         iL=-10
c         mindis = 1000000.0
c         do iL1=1,15
c            if (mindis .gt. ABS(interheight-(ht(iL1,j,k)))) then
c               iL = iL1
c               mindis = ABS(interheight-(ht(iL1,j,k)))
c               if (mindis .lt. 5000.) then
c                  print*,"low mindis!!!"
c               endif
c            endif
c         enddo
c         print*,"here is the interheight ",interheight
c         print*,"here is the altitude slab ",iL
c         print*,"here is mindis ",mindis         
cc     now do the interpolation, a cubic spline I think is what the method is called
c          if (ht(iL,j,k) .lt. interheight) then
c            val1 = den(iL,j,k) +
c     *           ((interheight-ht(iL,j,k))/(ht(iL+1,j,k)-ht(iL,j,k)))*
c     *           (den(iL+1,j,k)-den(iL,j,k))
c          else
c            val1 = den(iL,j,k) +
c     *           ((interheight-ht(iL,j,k))/(ht(iL-1,j,k)-ht(iL,j,k)))*
c     *           (den(iL-1,j,k)-den(iL,j,k))
c          endif
c          val1 = val1*1.0E11
c          if (ht(iL,j,k+1) .lt. interheight) then
c            val2 = den(iL,j,k+1) +
c     *     ((interheight-ht(iL,j,k+1))/(ht(iL+1,j,k+1)-ht(iL,j,k+1)))*
c     *           (den(iL+1,j,k+1)-den(iL,j,k+1))
c          else
c            val2 = den(iL,j,k+1) +
c     *     ((interheight-ht(iL,j,k+1))/(ht(iL-1,j,k+1)-ht(iL,j,k+1)))*
c     *           (den(iL-1,j,k+1)-den(iL,j,k+1))
c          endif
c          val2 = val2*1.0E11
c          if (ht(iL,j+1,k) .lt. interheight) then
c            val3 = den(iL,j+1,k) +
c     *      ((interheight-ht(iL,j+1,k))/(ht(iL+1,j+1,k)-ht(iL,j+1,k)))*
c     *           (den(iL+1,j+1,k)-den(iL,j+1,k))
c          else
c            val3 = den(iL,j+1,k) +
c     *     ((interheight-ht(iL,j+1,k))/(ht(iL-1,j+1,k)-ht(iL,j+1,k)))*
c     *           (den(iL-1,j+1,k)-den(iL,j+1,k))
c          endif
c          val3 = val3*1.0E11
c          if (ht(iL,j+1,k+1) .lt. interheight) then
c            val4 = den(iL,j+1,k+1) +
c     * ((interheight-ht(iL,j+1,k+1))/(ht(iL+1,j+1,k+1)-ht(iL,j+1,k+1)))*
c     *           (den(iL+1,j+1,k+1)-den(iL,j+1,k+1))
c          else
c            val4 = den(iL,j+1,k+1) +
c     * ((interheight-ht(iL,j+1,k+1))/(ht(iL-1,j+1,k+1)-ht(iL,j+1,k+1)))*
c     *           (den(iL-1,j+1,k+1)-den(iL,j+1,k+1))
c          endif
c          val4 = val4*1.0E11
c          rate1=(long-18*int(long/18))/18.0
c          val12=val1+(val2-val1)*rate1
c          val34=val3+(val4-val3)*rate1
c          rate2=(lat-2*int(lat/2))/2.0
c          znden1 = val12+(val34-val12)*rate2
cc     error handeling for times it passes over the cap which ins't calculated
cc     or for times the satellite is at the edge of the grid and we can't interpolate
cc     these aren't very pretty and I hope that they are accurate
c      if (znden1 .lt. 0) then
c      print*,"density is negative interpolation error -999.999999"
c      znden1 = -999.999
c      goto 189
c      endif
c      if (znden1 .lt. 0.09) then
c      print*,"on edge of boundary at 89 deg so -999.999999"
c      znden1 = -999.999
c      goto 189
c      endif
c      if (den(il,J,K) .lt. 1E-20) then
c      print*,"Density is too small, probably at boundary -999.999999"
c      znden1 = -999.999
c      goto 189
c      endif
c      if(ht(iL,j,k) .lt. interheight .and. iL .eq.15)then
c      print*,"sat is outside of interpolatable grid, -999.999999"
c      znden1 = -999.999
c      endif
c 189  continue
c      print*,"here is the density ", den(iL,J,K)*1.0E11
c     1," at ",ht(iL,j,k)
c          print*,"here is the interpolated density ", znden1
c      print*,"here is the density ", den(iL+1,J,K)*1.0E11,
c     1" at  ", ht(iL+1,j,k)
c          if (I .eq. 5) then          
c          open(unit=30, file="Grace1_interpolated.txt",access="append")
c          write(30,*)timenuminit,lat,long,interheight/1000.0,znden1
c          close(unit=30)
c          endif
c          if (I .eq. 6) then          
c          open(unit=30, file="Grace2_interpolated.txt",access="append")
c          write(30,*)timenuminit,lat,long,interheight/1000.0,znden1
c          close(unit=30)
c          endif
c      write(29,*)satname(i),xgeo*6731,ygeo*6731,zgeo*6731,znden1 !write sats position geo
c      write(28,*)satname(i),satx(i)*6731,saty(i)*6731,satz(i)*6731,
c     1znden1 !write sats position gse
c      endif
      stop
      end

c-----------------------------------------------------------
      subroutine getf11(iu,a1,nx,l1,l2,it)
c-----------------------------------------------------------
      real a1(*)
      character*80 rec
      character*80 l1,l2
      m=len(l1)
      call end0(l1,m)
      isany=0
      if(l1(1:m).eq.'any') isany=1
100   continue
      read(iu,1000,end=190) rec
      if(rec(1:10).eq.'FIELD-1D-1') then
      read(iu,1000,end=190) rec
      if((isany.eq.0).and.l1(1:m).ne.rec(1:m)) goto 100
      l1=rec
      read(iu,1000,end=190) l2
      read(iu,*,end=190)it,nx
      call rdn2(iu,a1,nx,rec,it,rid)
      return
      endif
      goto 100
190   continue
      nx=-1
1000  format(a)
      return
      end
c-----------------------------------------------------------
      subroutine getf21(iu,a1,nx,ny,l1,l2,it)
c-----------------------------------------------------------
      real a1(*)
      character*80 rec
      character*80 l1,l2
      m=len(l1)
      call end0(l1,m)
      isany=0
      if(l1(1:m).eq.'any') isany=1
100   continue
      read(iu,1000,end=190) rec
      if(rec(1:10).eq.'FIELD-2D-1') then
      read(iu,1000,end=190) rec
      if((isany.eq.0).and.l1(1:m).ne.rec(1:m)) goto 100
      l1=rec
      read(iu,1000,end=190) l2
      read(iu,*,end=190)it,nx,ny
      call rdn2(iu,a1,nn,rec,it,rid)
      if(nn.lt.0)nx=nn
      return
      endif
      goto 100
190   continue
      nx=-1
1000  format(a)
      return
      end
c-----------------------------------------------------------
      subroutine end0(r,m)
c-----------------------------------------------------------
      character*80 r
      n=min0(m,len(r))
      do 100 i=1,n
      m=i-1
      if(r(i:i).eq.' ') return
100   continue
      m=n
      return
      end
c-----------------------------------------------------------
      subroutine datf11(iu,a1,nx,l1,l2,it)
c-----------------------------------------------------------
      character*80 l1,l2
      real a1(*)
      write(iu,'(a)')'FIELD-1D-1'
      write(iu,'(a)') l1
      write(iu,'(a)') l2
      write(iu,*) it,nx
      call wrn2(iu,a1,nx,'FUNC-1-1',it,float(it))
      return
      end
c-----------------------------------------------------------
      subroutine datf21(iu,a1,nx,ny,l1,l2,it)
c-----------------------------------------------------------
      character*80 l1,l2
      real a1(*)
      write(iu,'(a)')'FIELD-2D-1'
      write(iu,'(a)') l1
      write(iu,'(a)') l2
      write(iu,*)it,nx,ny
      call wrn2(iu,a1,nx*ny,'FUNC-2-1',it,float(ny))
      return
      end
c-----------------------------------------------------------
      subroutine datf31(iu,a1,nx,ny,nz,l1,l2,it)
c-----------------------------------------------------------
      character*80 l1,l2
      real a1(*)
      write(iu,'(a)')'FIELD-3D-1'
      write(iu,'(a)') l1
      write(iu,'(a)') l2
      write(iu,*)it,nx,ny,nz
      call wrn2(iu,a1,nx*ny*nz,'FUNC-3-1',it,float(ny))
      return
      end
c-----------------------------------------------------------
      subroutine getf31(iu,a1,nx,ny,nz,l1,l2,it)
c-----------------------------------------------------------
      real a1(*)
      character*80 rec
      character*80 l1,l2
      m=len(l1)
      call end0(l1,m)
      isany=0
      if(l1(1:m).eq.'any') isany=1
100   continue
      read(iu,1000,end=190) rec
      if(rec(1:10).eq.'FIELD-3D-1') then
      read(iu,1000,end=190) rec
      if((isany.eq.0).and.l1(1:m).ne.rec(1:m)) goto 100
      l1=rec
      read(iu,1000,end=190) l2
      read(iu,*,end=190) it,nx,ny,nz
      call rdn2(iu,a1,nn,rec,it,rid)
      if(nn.lt.0)nx=nn
      return
      endif
      goto 100
190   continue
      nx=-1
1000  format(a)
      return
      end
c.---------------------------------------
      subroutine wrn2(iu,a,n,cid,it,rid)
c.---------------------------------------
      real a(n)
      character*8 cid
      real a1(0:63)
      integer i1(0:63),i2(0:63),i3(0:63)
      zmin=1.e33
      zmax=1.e-33
      do 100 i=1,n
      b=abs(a(i))
      zmin=amin1(zmin,b)
      zmax=amax1(zmax,b)
100   continue
1000  format(a,i8,3e14.7,i8,a)
      zmin=amax1(1.e-33,zmin)
      zmax=amin1( 1.e33,zmax)
      z2=alog(zmax)
      z1=-76.
      if(zmin.ne.0.)z1=alog(zmin)
      z1=amax1(z1,z2-37.)
      if(abs(z2-z1).le.1.e-5) then
      z1=z1-1.0
      z2=z2+1.0
      endif
      z0=exp(z1)
      dz=float(4410)/(z2-z1)
      write(iu,1000)'WRN2',n,z1,z2,rid,it,cid
      do 200 k=1,n,64
      nk=min0(63,n-k)
      do 210 i=0,nk
      a1(i)=amin1(1.e33,abs(a(i+k)))
      a1(i)=amax1(z0,a1(i))
      a1(i)=dz*(alog(a1(i))-z1)+0.5
      i3(i)=a1(i)
      i3(i)=max0(0,i3(i))
      i3(i)=min0(4414,i3(i))
      i1(i)=i3(i)/94
      i2(i)=i3(i)-94*i1(i)
      if(a(i+k).lt.0.)i1(i)=i1(i)+47
      i1(i)=i1(i)+33
      i2(i)=i2(i)+33
210   continue
      call wrnenc(iu,i1,nk)
      call wrnenc(iu,i2,nk)
200   continue
      return
      end
c.---------------------------------------
      subroutine wrnenc(iu,i1,n)
c.---------------------------------------
      integer i1(0:63)
      integer i2(-1:63)
      character c*72
      ick=i1(0)
      ir=i1(0)
      ic=1
      k=-1
      do 100 i=1,n
      ick=ick+i1(i)
      if(i1(i).eq.ir) then
      ic=ic+1
      else
      if(ic.eq.1) then
      k=k+1
      i2(k)=ir
      else
      k=k+1
      i2(k)=ic+170
      k=k+1
      i2(k)=ir
      endif
      ic=1
      ir=i1(i)
      endif
100   continue
      if(ic.eq.1) then
      k=k+1
      i2(k)=ir
      else
      k=k+1
      i2(k)=ic+170
      k=k+1
      i2(k)=ir
      endif
      i2(-1)=33+mod(ick,92)
      j=0
      do 200 i=-1,k
      j=j+1
      c(j:j)=char(i2(i))
200   continue
      write(iu,'(a)')c(1:j)
      return
      end
c.---------------------------------------
      subroutine wrndec(iu,i1,n)
c.---------------------------------------
      integer i1(0:63)
      common /wrdqq9/i2(-2:67)
      character c*72
      c=' '
      read(iu,'(a)',end=900,err=900) c
      i=-2
      j=0
100   continue
      j=j+1
      if(j.gt.67) then
      write(0,*)'wrndec: cant find end of encoded record'
      n=-5
      return
      endif
      ir=ichar(c(j:j))
      if(ir.eq.32) goto 190
      if(ir.gt.127) then
      ic=ir-170
      j=j+1
      ir=ichar(c(j:j))
      do 110 l=1,ic
      i=i+1
      i2(i)=ir
110   continue
      else
      i=i+1
      i2(i)=ir
      endif
      goto 100
190   continue
      ick=0
      n=i
      if(n.gt.63) then
      write(0,*)'wrndec: n gt 63, n=',n
      n=-5
      write(0,'(a)')'rec:'
      write(0,'(a)')c
      return
      endif
      do 200 i=0,n
      i1(i)=i2(i)
      ick=ick+i1(i)
200   continue
      if(i2(-1).ne.33+mod(ick,92)) then
      write(0,*)'wrndec: checksum error '
      write(0,'(a,a)')'rec:',c
      write(0,*)i2(-1),33+mod(ick,92)
      n=-5
      return
      endif
      return
900   continue
      write(0,*)' wrndec eof/err '
      n=-5
      return
      end
c.---------------------------------------
      subroutine rdn2(iu,a,n,cid,it,rid)
c.---------------------------------------
c character decoding
      real a(*)
      character*8 cid
      character*4 did
      real a1(0:63)
      integer i1(0:63),i2(0:63),i3(0:63)
100   continue
      read(iu,1000,err=100,end=900)did,n,zmin,zmax,rid,it,cid
1000  format(a,i8,3e14.7,i8,a)
      if(did.ne.'WRN2') goto 100
      if(zmin.eq.zmax) then
      do 200 i=1,n
      a(i)=zmin
200   continue
      return
      endif
      l=0
      do 300 k=1,n,64
      nk=min0(63,n-k)
      call wrndec(iu,i1,nn)
      if(nn.ne.nk) then
      write(0,*)'rdn2: nn .ne. nk ',nn,nk,n,k
      n=-2
      return
      endif
      call wrndec(iu,i2,nn)
      if(nn.ne.nk) then
      write(0,*)'rdn2: nn .ne. nk ',nn,nk
      n=-2
      return
      endif
      dzi=(zmax-zmin)/float(4410)
      do 330 i=0,nk
      i1(i)=i1(i)-33
      i2(i)=i2(i)-33
      sig=1.
      if(i1(i).ge.47) then
      sig=-1.
      i1(i)=i1(i)-47
      endif
      i3(i)=i2(i)+94*i1(i)
      a1(i)=i3(i)
      a1(i)=dzi*a1(i)+zmin
      a1(i)=sig*exp(a1(i))
      l=l+1
      a(l)=a1(i)
330   continue
300   continue
      return
900   continue
      n=-1
      return
      end
c-------------------------------------
      integer function jlen(c,m)
c-------------------------------------
      character*80 c
      k=1
      do 100 i=1,m
      if(c(i:i).ne.' ')k=i
100   continue
      jlen=k
      return
      end
c----------------------------------------------------
      SUBROUTINE ctim_cotr_set(NYEAR,NMONTH,NDAY,NHOUR,NMINUTE,SEC)
c----------------------------------------------------
cdoc-
cdoc-------------------------------------------------------------------
cdoc-cotr_set
cdoc-------------------------------------------------------------------
cdoc-
cdoc-SYNOPSIS:
cdoc-  subroutine cotr_set(nyear,nmonth,nday,nhour,nminute,sec)
cdoc-
cdoc-FUNCTION:
cdoc-  computes coordinate transformation matrices
cdoc-  for cotr() for a given time
cdoc-
cdoc-CALL SEQUENCE:
cdoc-  must be called before cotr()
cdoc-
cdoc-PARAMETERS:
cdoc-  nyear (in,integer):  year (nn or nnnn notation) (1900<nyear<2000)
cdoc-  nmonth (in,integer):  month             (1<=nmonth<=12)
cdoc-  nday (in,integer):  day of month      (1<=nday<=31)
cdoc-  nhour (in,integer):  hour of day       (0<=nhour<24)
cdoc-  nminute (in,integer):  minute of hour    (0<=nminute<60)
cdoc-  sec (in,real):     seconds of minute (0.0<=sec<60.0)
cdoc-
cdoc-NOTES:
cdoc-  1.  no output, transformation matrices are stored internally
cdoc-      for use by cotr()
cdoc-  2.  dipole coordinates are computed from IGRF g10,g11,h11
cdoc-      coefficients between 1945 and 1990 with linear time
cdoc-      interpolation.  Before 1945 coefficients for 1945 are
cdoc-      used and after 1990 coefficients for 1990 are used.
cdoc-  3.  internal arithmetic is handled in double precision
cdoc-
cdoc-  4.  for nyear=1967 nmonth=1 nday=1 nhour=0 nminute=0 sec=any all transformations
cdoc-      are set to identity
cdoc-
cdoc-
cdoc-FILES:
cdoc-  none
cdoc-
cdoc-AUTHOR:
cdoc-  Joachim Raeder, IGPP/UCLA   June 1994
cdoc-
cdoc-REFERENCES:
cdoc-
cdoc-        Almanac for Computers 1991
cdoc-        Nautical Almanac Office
cdoc-        US Naval Observatory
cdoc-        Washington, DC
cdoc-
cdoc-        M A Hapgood
cdoc-        Space physics coordinate transformations: A user guide
cdoc-        Planet. Space Sci., 40, 711, 1992
cdoc-
      REAL*8 cmat,d1,eps,g10,g11,gang,gg10,gg11,gmstd, gmsth,h11,hh11
     *,one,pih,pmat,psi,qx,qxg,qy
      REAL*8 qyg,qz,qzg,rad,t,t0,tig,tmp,tmpco,tmpmat, tmpsi,tmpx,
     *tmpy,tmpz,uthour,w1,xjd, xjd0
      REAL*8 xl,xls,xm,xmjd,xmu,xxl,xxp,xy,zero
      INTEGER*4 i,i1,i2,i3,i4,ier,iscalled,j,jy,k,mjd
      COMMON /cotrd1/ t0,iscalled
      COMMON /cotrd2/ tmpmat(3,3,2),cmat(3,3),pmat(3,3,-7:7)
      COMMON /cotrd3/ ccfr,ccto
      CHARACTER*3 ccfr,ccto
      INTEGER*4 NYEAR,NMONTH,NDAY,NHOUR,NMINUTE
      REAL*8 SEC
      DIMENSION tig(12),g10(12),g11(12),h11(12)
      DATA tig(1),g10(1),g11(1),h11(1)/1000.00, - 30594.00, - 2285.00
     *,5810.00/
      DATA tig(2),g10(2),g11(2),h11(2)/1945.00, - 30594.00, - 2285.00
     *,5810.00/
      DATA tig(3),g10(3),g11(3),h11(3)/1950.00, - 30554.00, - 2250.00
     *,5815.00/
      DATA tig(4),g10(4),g11(4),h11(4)/1955.00, - 30500.00, - 2215.00
     *,5820.00/
      DATA tig(5),g10(5),g11(5),h11(5)/1960.00, - 30421.00, - 2169.00
     *,5791.00/
      DATA tig(6),g10(6),g11(6),h11(6)/1965.00, - 30334.00, - 2119.00
     *,5776.00/
      DATA tig(7),g10(7),g11(7),h11(7)/1970.00, - 30220.00, - 2068.00
     *,5737.00/
      DATA tig(8),g10(8),g11(8),h11(8)/1975.00, - 30100.00, - 2013.00
     *,5675.00/
      DATA tig(9),g10(9),g11(9),h11(9)/1980.00, - 29992.00, - 1956.00
     *,5604.00/
      DATA tig(10),g10(10),g11(10),h11(10)/1985.00, - 29873.00, - 
     *1905.00,5500.00/
      DATA tig(11),g10(11),g11(11),h11(11)/1990.00, - 29775.00, - 
     *1851.00,5411.00/
      DATA tig(12),g10(12),g11(12),h11(12)/3000.00, - 29775.00, - 
     *1851.00,5411.00/
      DATA pih/1.570796372D0/
      DATA rad/17.45329252D-3/
      DATA one/1.0D0/
      DATA zero/0.0D0/
c
c..... check input consistency
c
      jy = NYEAR
      IF ( jy.LT.149 ) jy = jy + 1900
      ier = 0
      IF ( jy.LT.1901 ) ier = 1
      IF ( jy.GT.2049 ) ier = 7
      IF ( NMONTH.LT.1 .OR. NMONTH.GT.12 ) ier = 2
      IF ( NDAY.LT.1 .OR. NDAY.GT.31 ) ier = 3
      IF ( NHOUR.LT.0 .OR. NHOUR.GT.23 ) ier = 4
      IF ( NMINUTE.LT.0 .OR. NMINUTE.GT.63 ) ier = 5
      IF ( SEC.LT.0.0 .OR. SEC.GT.60.0 ) ier = 6
      IF ( ier.NE.0 ) THEN
      WRITE (0,*) 'error cotr_set, ier= ',ier
      WRITE (0,*) NYEAR,jy,NMONTH,NDAY,NHOUR,NMINUTE,SEC
      STOP
      ENDIF
c
      iscalled = 1
c
c.....  universal time of day in hours (UT1)
c
      uthour = dble(float(NHOUR)) + (dble(float(NMINUTE))/60.0D0) + 
     *dble(SEC)/3600.0D0
c
c..... full julian day
c
      i1 = 367*jy
      i2 = (NMONTH+9)/12
      i2 = (i2+jy)*7
      i2 = i2/4
      i3 = (275*NMONTH)/9
      i4 = 100*jy + NMONTH
      d1 = one
      IF ( (dble(i4)-190002.5D0).LT.zero ) d1 = (-d1)
      xjd = dble(i1) - dble(i2) + dble(i3) + dble(NDAY) + 1721013.5D0
     * + uthour/24.0D0 - 0.5D0*d1 + 0.5
c
c..... julian day at 0 UT
c
      xjd0 = xjd - (uthour/24.0D0)
c
c.....  modified julian day
c
      xmjd = xjd - 2400000.5D0
      mjd = int(xmjd)
c
c.....  set to identity transformation if nyear=1967,....
c       except for the MHD transformation which still flips x- y- axes
c
      IF ( jy.LE.1967 ) THEN
      DO 302 k = -7,7
      DO 301 j = 1,3
      DO 300 i = 1,3
      pmat(i,j,k) = 0.0D0
      IF ( i.EQ.j ) pmat(i,j,k) = 1.0D0
300   CONTINUE
301   CONTINUE
302   CONTINUE
      DO 311 i = 1,3
      DO 310 j = 1,3
      pmat(j,i,6) = zero
310   CONTINUE
311   CONTINUE
      pmat(1,1,6) = (-1.0D0)
      pmat(2,2,6) = (-1.0D0)
      pmat(3,3,6) = (1.0D0)
      DO 321 i = 1,3
      DO 320 j = 1,3
      pmat(j,i,-6) = pmat(i,j,6)
320   CONTINUE
321   CONTINUE
      RETURN
c
      ENDIF
c
c
c..... solar coordinates
c
      t = (xjd-2451545.0D0)/36525.0D0
      t0 = (xjd0-2451545.0D0)/36525.0D0
      xm = 357.528D0 + 35999.050D0*t
      xl = 280.460D0 + 36000.772D0*t
      xm = dmod(xm+360.0D3,360.0D0)
      xl = dmod(xl+360.0D3,360.0D0)
      xls = xl + ((1.915D0-0.0048*t0)*dsin(rad*xm)) + (0.020D0*dsin(
     *2.0D0*rad*xm))
      eps = 23.439D0 - 0.013D0*t
      eps = dmod(eps+360.0D3,360.0D0)
c
c..... greenwich mean siderial time
c
      gmsth = 6.69737456D0 + 2400.051336D0*t0 + 0.0000258622D0*t0*t0 
     *+ 1.002737909D0*uthour
      gmsth = dmod(gmsth+24.0D3,24.0D0)
      gmstd = 15.0D0*gmsth
      gang = rad*gmstd
c
      eps = rad*eps
      xls = rad*xls
c
c.... gei to geo
c
      DO 94972 i = 1,3
      DO 94971 j = 1,3
      pmat(j,i,1) = zero
94971 CONTINUE
94972 CONTINUE
      tmpsi = dsin(gang)
      tmpco = dcos(gang)
      pmat(1,1,1) = tmpco
      pmat(2,2,1) = tmpco
      pmat(3,3,1) = one
      pmat(2,1,1) = -tmpsi
      pmat(1,2,1) = tmpsi
      DO 94914 i = 1,3
      DO 94913 j = 1,3
      pmat(j,i,-1) = pmat(i,j,1)
94913 CONTINUE
94914 CONTINUE
c
c.... gei to gse
c
      DO 94885 i = 1,3
      DO 94884 j = 1,3
      tmpmat(j,i,1) = zero
94884 CONTINUE
94885 CONTINUE
      tmpsi = dsin(xls)
      tmpco = dcos(xls)
      tmpmat(1,1,1) = tmpco
      tmpmat(2,2,1) = tmpco
      tmpmat(3,3,1) = one
      tmpmat(2,1,1) = -tmpsi
      tmpmat(1,2,1) = tmpsi
c
      DO 94827 i = 1,3
      DO 94826 j = 1,3
      tmpmat(j,i,2) = zero
94826 CONTINUE
94827 CONTINUE
      tmpsi = dsin(eps)
      tmpco = dcos(eps)
      tmpmat(1,1,2) = one
      tmpmat(2,2,2) = tmpco
      tmpmat(3,3,2) = tmpco
      tmpmat(3,2,2) = -tmpsi
      tmpmat(2,3,2) = tmpsi
c
      DO 94769 i = 1,3
      DO 94768 j = 1,3
      pmat(i,j,2) = tmpmat(i,1,1)*tmpmat(1,j,2) + tmpmat(i,2,1) *
     *tmpmat(2,j,2) + tmpmat(i,3,1)*tmpmat(3,j,2)
94768 CONTINUE
94769 CONTINUE
      DO 94740 i = 1,3
      DO 94739 j = 1,3
      pmat(j,i,-2) = pmat(i,j,2)
94739 CONTINUE
94740 CONTINUE
c
c
c .... gse to gsm
c
      xy = 100.0D0*(20.0D0+t)
      DO 100 k = 1,12 - 1
      IF ( xy.GE.tig(k) .AND. xy.LE.tig(k+1) ) THEN
      tmp = (xy-tig(k))/(tig(k+1)-tig(k))
      gg10 = (one-tmp)*g10(k) + tmp*g10(k+1)
      gg11 = (one-tmp)*g11(k) + tmp*g11(k+1)
      hh11 = (one-tmp)*h11(k) + tmp*h11(k+1)
      GOTO 190
      ENDIF
100   CONTINUE
190   xxl = datan(hh11/gg11)
      xxp = pih - dasin((gg11*dcos(xxl)+hh11*dsin(xxl))/gg10)
c
      qxg = dcos(xxp)*dcos(xxl)
      qyg = dcos(xxp)*dsin(xxl)
      qzg = dsin(xxp)
      tmpx = pmat(1,1,-1)*qxg + pmat(1,2,-1)*qyg + pmat(1,3,-1)*qzg
      tmpy = pmat(2,1,-1)*qxg + pmat(2,2,-1)*qyg + pmat(2,3,-1)*qzg
      tmpz = pmat(3,1,-1)*qxg + pmat(3,2,-1)*qyg + pmat(3,3,-1)*qzg
      qx = pmat(1,1,2)*tmpx + pmat(1,2,2)*tmpy + pmat(1,3,2)*tmpz
      qy = pmat(2,1,2)*tmpx + pmat(2,2,2)*tmpy + pmat(2,3,2)*tmpz
      qz = pmat(3,1,2)*tmpx + pmat(3,2,2)*tmpy + pmat(3,3,2)*tmpz
      psi = datan2(qy,qz)
      DO 94653 i = 1,3
      DO 94652 j = 1,3
      pmat(j,i,-3) = zero
94652 CONTINUE
94653 CONTINUE
      tmpsi = dsin(psi)
      tmpco = dcos(psi)
      pmat(1,1,-3) = one
      pmat(2,2,-3) = tmpco
      pmat(3,3,-3) = tmpco
      pmat(3,2,-3) = -tmpsi
      pmat(2,3,-3) = tmpsi
      DO 94595 i = 1,3
      DO 94594 j = 1,3
      pmat(j,i,3) = pmat(i,j,-3)
94594 CONTINUE
94595 CONTINUE
c
c.... gsm to sm
c
      xmu = datan(qx/sqrt(qy*qy+qz*qz))
      DO 94566 i = 1,3
      DO 94565 j = 1,3
      pmat(j,i,-4) = zero
94565 CONTINUE
94566 CONTINUE
      tmpsi = dsin(xmu)
      tmpco = dcos(xmu)
      pmat(1,1,-4) = tmpco
      pmat(2,2,-4) = one
      pmat(3,3,-4) = tmpco
      pmat(3,1,-4) = -tmpsi
      pmat(1,3,-4) = tmpsi
      DO 94508 i = 1,3
      DO 94507 j = 1,3
      pmat(j,i,4) = pmat(i,j,-4)
94507 CONTINUE
94508 CONTINUE
c
c.... geo to mag
c
      w1 = xxp - pih
      DO 94479 i = 1,3
      DO 94478 j = 1,3
      tmpmat(j,i,1) = zero
94478 CONTINUE
94479 CONTINUE
      DO 94450 i = 1,3
      DO 94449 j = 1,3
      tmpmat(j,i,2) = zero
94449 CONTINUE
94450 CONTINUE
      tmpsi = dsin(w1)
      tmpco = dcos(w1)
      tmpmat(1,1,1) = tmpco
      tmpmat(2,2,1) = one
      tmpmat(3,3,1) = tmpco
      tmpmat(3,1,1) = -tmpsi
      tmpmat(1,3,1) = tmpsi
      tmpsi = dsin(xxl)
      tmpco = dcos(xxl)
      tmpmat(1,1,2) = tmpco
      tmpmat(2,2,2) = tmpco
      tmpmat(3,3,2) = one
      tmpmat(2,1,2) = -tmpsi
      tmpmat(1,2,2) = tmpsi
      DO 94363 i = 1,3
      DO 94362 j = 1,3
      pmat(i,j,5) = tmpmat(i,1,1)*tmpmat(1,j,2) + tmpmat(i,2,1) *
     *tmpmat(2,j,2) + tmpmat(i,3,1)*tmpmat(3,j,2)
94362 CONTINUE
94363 CONTINUE
      DO 94334 i = 1,3
      DO 94333 j = 1,3
      pmat(j,i,-5) = pmat(i,j,5)
94333 CONTINUE
94334 CONTINUE
c
c.....  gse to mhd
c
      DO 94305 i = 1,3
      DO 94304 j = 1,3
      pmat(j,i,6) = zero
94304 CONTINUE
94305 CONTINUE
      pmat(1,1,6) = (-1.0D0)
      pmat(2,2,6) = (-1.0D0)
      pmat(3,3,6) = (1.0D0)
      DO 94276 i = 1,3
      DO 94275 j = 1,3
      pmat(j,i,-6) = pmat(i,j,6)
94275 CONTINUE
94276 CONTINUE
c
      RETURN
      END
c----------------------------------------------------
      SUBROUTINE ctim_cotr(CFR,CTO,X1,Y1,Z1,X2,Y2,Z2)
c----------------------------------------------------
cdoc-
cdoc-------------------------------------------------------------------
cdoc-cotr
cdoc-------------------------------------------------------------------
cdoc-
cdoc-SYNOPSIS:
cdoc-  subroutine cotr(cfr,cto,x1,y1,z1,x2,y2,z2)
cdoc-
cdoc-FUNCTION:
cdoc-  transform vector (x1,y1,z1) in coordinate system 'cfr'
cdoc-  to vector (x2,y2,z2) in coordinate system 'cto'
cdoc-  coordinate systems 'cfr' and 'cto' can be any of the following:
cdoc-  'gei'  :  geocentric equatorial inertial
cdoc-  'geo'  :  geographic
cdoc-  'gse'  :  geocentric solar ecliptic
cdoc-  'gsm'  :  geocentric solar magnetospheric
cdoc-  'sm '  :  solar magnetic
cdoc-  'mag'  :  geomagnetic
cdoc-  'mhd'  :  global mhd simulation system, like 'gse' with
cdoc-            x and y axes mirrored
cdoc-
cdoc-CALL SEQUENCE:
cdoc-  cotr_set() must be called before the first call to
cdoc-  cotr() or if the transformation is for a different time
cdoc-
cdoc-PARAMETERS:
cdoc-  cfr (in,character*3)   'from' coordinate system, for ex. 'gse'
cdoc-  cto (in,character*3)   'to' coordinate system, for ex. 'gsm'
cdoc-  x1,y1,z1 (in,real)     input vector
cdoc-  x2,y2,z2 (out,real)    transformed vector
cdoc-
cdoc-NOTES:
cdoc-  see NOTES and REFERENCES for cotr_set()
cdoc-
cdoc-FILES:
cdoc-  none
cdoc-
cdoc-AUTHOR:
cdoc-  Joachim Raeder, IGPP/UCLA   June 1994
cdoc-
c
c
      REAL*8 cmat,pmat,t0,tmpmat,xx1,xx2,yy1,yy2,zz1,zz2
      INTEGER*4 i,icalled,ifr,iscalled,ito,j,k,kmat,l
      CHARACTER*(*) CFR,CTO
      REAL*8 X1,Y1,Z1,X2,Y2,Z2
      COMMON /cotrd1/ t0,iscalled
      COMMON /cotrd2/ tmpmat(3,3,2),cmat(3,3),pmat(3,3,-7:7)
      COMMON /cotrd3/ ccfr,ccto
      CHARACTER*3 ccfr,ccto
      CHARACTER*3 clfr,clto
      INTEGER*4 imat(7,7,7)
      DATA icalled/0/
      DATA (imat(l,1,1),l=1,7)/0,0,0,0,0,0,0/
      DATA (imat(l,2,1),l=1,7)/ - 1,0,0,0,0,0,0/
      DATA (imat(l,3,1),l=1,7)/ - 2,0,0,0,0,0,0/
      DATA (imat(l,4,1),l=1,7)/ - 3, - 2,0,0,0,0,0/
      DATA (imat(l,5,1),l=1,7)/ - 4, - 3, - 2,0,0,0,0/
      DATA (imat(l,6,1),l=1,7)/ - 5, - 1,0,0,0,0,0/
      DATA (imat(l,7,1),l=1,7)/6, - 2,0,0,0,0,0/
      DATA (imat(l,1,2),l=1,7)/1,0,0,0,0,0,0/
      DATA (imat(l,2,2),l=1,7)/0,0,0,0,0,0,0/
      DATA (imat(l,3,2),l=1,7)/ - 2,1,0,0,0,0,0/
      DATA (imat(l,4,2),l=1,7)/ - 3, - 2,1,0,0,0,0/
      DATA (imat(l,5,2),l=1,7)/ - 4, - 3, - 2,1,0,0,0/
      DATA (imat(l,6,2),l=1,7)/ - 5,0,0,0,0,0,0/
      DATA (imat(l,7,2),l=1,7)/6, - 2,1,0,0,0,0/
      DATA (imat(l,1,3),l=1,7)/2,0,0,0,0,0,0/
      DATA (imat(l,2,3),l=1,7)/ - 1,2,0,0,0,0,0/
      DATA (imat(l,3,3),l=1,7)/0,0,0,0,0,0,0/
      DATA (imat(l,4,3),l=1,7)/ - 3,0,0,0,0,0,0/
      DATA (imat(l,5,3),l=1,7)/ - 4, - 3,0,0,0,0,0/
      DATA (imat(l,6,3),l=1,7)/ - 5, - 1,2,0,0,0,0/
      DATA (imat(l,7,3),l=1,7)/6,0,0,0,0,0,0/
      DATA (imat(l,1,4),l=1,7)/2,3,0,0,0,0,0/
      DATA (imat(l,2,4),l=1,7)/ - 1,2,3,0,0,0,0/
      DATA (imat(l,3,4),l=1,7)/3,0,0,0,0,0,0/
      DATA (imat(l,4,4),l=1,7)/0,0,0,0,0,0,0/
      DATA (imat(l,5,4),l=1,7)/ - 4,0,0,0,0,0,0/
      DATA (imat(l,6,4),l=1,7)/ - 5, - 1,2,3,0,0,0/
      DATA (imat(l,7,4),l=1,7)/6,3,0,0,0,0,0/
      DATA (imat(l,1,5),l=1,7)/2,3,4,0,0,0,0/
      DATA (imat(l,2,5),l=1,7)/ - 1,2,3,4,0,0,0/
      DATA (imat(l,3,5),l=1,7)/3,4,0,0,0,0,0/
      DATA (imat(l,4,5),l=1,7)/4,0,0,0,0,0,0/
      DATA (imat(l,5,5),l=1,7)/0,0,0,0,0,0,0/
      DATA (imat(l,6,5),l=1,7)/ - 5, - 1,2,3,4,0,0/
      DATA (imat(l,7,5),l=1,7)/6,3,4,0,0,0,0/
      DATA (imat(l,1,6),l=1,7)/1,5,0,0,0,0,0/
      DATA (imat(l,2,6),l=1,7)/5,0,0,0,0,0,0/
      DATA (imat(l,3,6),l=1,7)/ - 2,1,5,0,0,0,0/
      DATA (imat(l,4,6),l=1,7)/ - 3, - 2,1,5,0,0,0/
      DATA (imat(l,5,6),l=1,7)/ - 4, - 3, - 2,1,5,0,0/
      DATA (imat(l,6,6),l=1,7)/0,0,0,0,0,0,0/
      DATA (imat(l,7,6),l=1,7)/6, - 2,1,5,0,0,0/
      DATA (imat(l,1,7),l=1,7)/2, - 6,0,0,0,0,0/
      DATA (imat(l,2,7),l=1,7)/ - 1,2,6,0,0,0,0/
      DATA (imat(l,3,7),l=1,7)/ - 6,0,0,0,0,0,0/
      DATA (imat(l,4,7),l=1,7)/ - 3, - 6,0,0,0,0,0/
      DATA (imat(l,5,7),l=1,7)/ - 4, - 3, - 6,0,0,0,0/
      DATA (imat(l,6,7),l=1,7)/ - 5, - 1,2, - 6,0,0,0/
      DATA (imat(l,7,7),l=1,7)/0,0,0,0,0,0,0/
      IF ( icalled.EQ.0 .OR. iscalled.NE.0 ) THEN
      icalled = 1
      ccfr = 'XXX'
      ccto = 'xyz'
      ENDIF
      iscalled = 0
      xx1 = dble(X1)
      yy1 = dble(Y1)
      zz1 = dble(Z1)
      clfr = CFR
      clto = CTO
      IF ( clfr.NE.ccfr .OR. clto.NE.ccto ) THEN
      ifr = 0
      IF ( clfr.EQ.'GEI' .OR. clfr.EQ.'gei' ) ifr = 1
      IF ( clfr.EQ.'GEO' .OR. clfr.EQ.'geo' ) ifr = 2
      IF ( clfr.EQ.'GSE' .OR. clfr.EQ.'gse' ) ifr = 3
      IF ( clfr.EQ.'GSM' .OR. clfr.EQ.'gsm' ) ifr = 4
      IF ( clfr.EQ.'SM ' .OR. clfr.EQ.'sm ' ) ifr = 5
      IF ( clfr.EQ.'MAG' .OR. clfr.EQ.'mag' ) ifr = 6
      IF ( clfr.EQ.'MHD' .OR. clfr.EQ.'mhd' ) ifr = 7
      ito = 0
      IF ( clto.EQ.'GEI' .OR. clto.EQ.'gei' ) ito = 1
      IF ( clto.EQ.'GEO' .OR. clto.EQ.'geo' ) ito = 2
      IF ( clto.EQ.'GSE' .OR. clto.EQ.'gse' ) ito = 3
      IF ( clto.EQ.'GSM' .OR. clto.EQ.'gsm' ) ito = 4
      IF ( clto.EQ.'SM ' .OR. clto.EQ.'sm ' ) ito = 5
      IF ( clto.EQ.'MAG' .OR. clto.EQ.'mag' ) ito = 6
      IF ( clto.EQ.'MHD' .OR. clto.EQ.'mhd' ) ito = 7
      IF ( ifr.LT.1 .OR. ito.LT.1 ) THEN
      WRITE (0,*) ' error cotr: ifr,ito ',ifr,ito
      STOP
      ENDIF
c
      ccfr = clfr
      ccto = clto
      DO 101 i = 1,3
      DO 100 j = 1,3
      cmat(i,j) = 0.0D0
      IF ( i.EQ.j ) cmat(i,j) = 1.0D0
100   CONTINUE
101   CONTINUE
      DO 200 k = 1,7
      kmat = imat(k,ifr,ito)
      IF ( kmat.EQ.0 ) GOTO 290
      DO 211 i = 1,3
      DO 210 j = 1,3
      tmpmat(i,j,1) = pmat(i,1,kmat)*cmat(1,j) + pmat(i,2,kmat)*cmat(
     *2,j) + pmat(i,3,kmat)*cmat(3,j)
210   CONTINUE
211   CONTINUE
      DO 221 i = 1,3
      DO 220 j = 1,3
      cmat(i,j) = tmpmat(i,j,1)
220   CONTINUE
221   CONTINUE
200   CONTINUE
      ENDIF
290   xx2 = cmat(1,1)*xx1 + cmat(1,2)*yy1 + cmat(1,3)*zz1
      yy2 = cmat(2,1)*xx1 + cmat(2,2)*yy1 + cmat(2,3)*zz1
      zz2 = cmat(3,1)*xx1 + cmat(3,2)*yy1 + cmat(3,3)*zz1
      X2 = sngl(xx2)
      Y2 = sngl(yy2)
      Z2 = sngl(zz2)
      RETURN
      END
