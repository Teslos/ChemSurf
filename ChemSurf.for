        Program Chemsurf
        implicit double precision(a-h,o-z)
        double precision x2s(3), x1s(3), ac1s(3), 
     +   ac2s(3), sl(3), sr(3), ddx(3)
        character filetens*60, filether*60, phasename*24, name*24
        character nm1*24, nm2*24, fileplo*60
        character*24 names(5)
        integer   icom, iph, nct, ncp, ncn, ncmur, ncx, ncx1 nc1, nc2, j, k
        integer   ic1s, ic2s,ic1b,ic2b, iconst
	parameter (nwg=80000, nc=8, ncom=3)
	dimension stoi(ncom,ncom)
	dimension ncw(nc)
	dimension iwsg(nwg)
	logical SG1ERR, SG2ERR, TQGSP, pstat
c  ------------------------------------------
c       Initialize ThermoCalc read the surface tension data &
c       thermochemical data
c  ------------------------------------------
	open(13,file='debug',form='FORMATTED', 
     +  status='NEW', ERR=9)
        write(*,*) 'Surface tension for binary melts'
        call TQINI(nwg,iwsg)
        call TQSIO('OUTPUT',6)

        filetens = 'SURFTENS.DAT'
        open(12,file=filetens(:12), form='FORMATTED', 
     +  status='OLD', ERR=9)
   1    format(a)
   2    format(F12.6)
   3    format(F12.6,F12.6,F12.6)
   4    format(f12.6,f12.6)
        read(12,1) filether
        call TQRFIL(filether(:12), iwsg)
        write(*,*) 'Reading Thermo Data!'
c  ------------------------------------------
c       Read Phase Name
c  ------------------------------------------
        read(12,1) phasename
c  ------------------------------------------
c       Read BETA-MIX for equation
c  ------------------------------------------
        read(12,2) beta
c  ------------------------------------------
c       Read Name of constituent 1
c  ------------------------------------------
        read(12,1) nm1
c  ------------------------------------------
c       Read coefficients of surface tension equation
c       and Melting Temp for Constituent 1 
c  ------------------------------------------
        read(12,*) surf11, surf12, temp1
c  ------------------------------------------
c       Read coefficients of Molar volume equation for Constituent 1
c  ------------------------------------------
        read(12,*) v11, v12
c  ------------------------------------------
c       Read name of constituent 2
c  ------------------------------------------
        read(12,1) nm2
c  ------------------------------------------
c       Read coefficients of surface tension equation
c       and melting Temp. for constituent 2
c  ------------------------------------------
        read(12,*) surf21, surf22, temp2
c  ------------------------------------------
c       Read coefficients of molar volume equation for constituent 2
c  ------------------------------------------
        read(12,*) v21, v22
c  ------------------------------------------
c       Read Temperature and Composition step of Calculation
c  ------------------------------------------
        read(12,*) temp
        read(12,*) step        
c  ------------------------------------------
c       Preparation of Physical properties of components
c  ------------------------------------------
        R=8.314
        an0=6.02E+23**(1./3.)
c  ------------------------------------------
c       Surface tension
c  ------------------------------------------
        surf1 = surf11 - surf12*(temp - temp1)
        surf2 = surf21 - surf22*(temp - temp2)
c  ------------------------------------------
c       Molar volume
c  ------------------------------------------
        v1 = v11*(1. + v12*(temp - temp1))
        v2 = v21*(1. + v22*(temp - temp2))
c  ------------------------------------------
c       Molar surface area
c  ------------------------------------------
        s1 = 1.091*((v1)**(2./3.))*an0
        s2 = 1.091*((v2)**(2./3.))*an0
c  ------------------------------------------
c       Write header section for plot file including
c       values for pure component 1
c  ------------------------------------------
        write(*,*) 'Start calculation of surface tension'
        fileplo = 'surftens.plo'
        open(11, file = fileplo(:12), access='SEQUENTIAL', status='NEW'
     +  , err=11)
        write(11,*) -1
        write(11,*) 1
        write(11,9999) nm2(:2)
9999    format('"X(',A,')"','"SIGMA in N/m"')
        write(11,9998) nm1(:2),nm2(:2)
9998    format('2', 'Surface tension for ',A,'~',A,'"')
        write(11,9997) temp
9997    format('" for T=',F7.2,' K"')
        write(11,9996) 1,-1,0,2,2,0,0, '"Calc."'
9996    format(7(I2,1X),A3)
        write(11,*) 1.0, surf2
c  ------------------------------------------
c       Set new system components
c  ------------------------------------------
        names(1)='CEO2'
        names(2)='COO'
        names(3)='O'
        stoi(1,1)=1.0; stoi(1,2)=0.0; stoi(1,3)=2.0
        stoi(2,1)=0.0; stoi(2,2)=1.0; stoi(2,3)=1.0
        stoi(3,1)=0.0; stoi(3,2)=0.0; stoi(3,3)=1.0

        call TQSCOM(ncom, names, stoi, iwsg)
        
c  ------------------------------------------        
c       Get component names in the system
c  ------------------------------------------
        call TQGCOM(icom, names, iwsg)
        print *, 'This system has the following components:'
        print *, (names(i),i=1,icom)
        print *
c  ------------------------------------------
c       Get number of phases in the system
c  ------------------------------------------
        call TQGNP(iph, iwsg)
        print *, 'This system has', iph, ' phases:'
c  ------------------------------------------
c       get names and status of the phases in the system
c  ------------------------------------------
        do  i=1,iph
          call TQGPN(i, name, iwsg)
c          pstat = TQGSP(i, sname, an, iwsg)
          print *, i, '   ', name, '  ', sname, '  ', an
        end do
c  ------------------------------------------
c       Calculation of Activites in Bulk Phase
c  ------------------------------------------
c  ------------------------------------------
c       Find Phase Index for Phase Name
c  ------------------------------------------	
        call TQGPI( indexp, phasename(:9), iwsg )
        ipliq = indexp
        print *, 'IONIC_LIQ:', ipliq
        
	write(*,*) 'LIQ: ',phasename(:9),' phase num:',ipliq, 
     +	' COMPONENT 1: ',nm1(:2), ' COMPONENT 2: ',nm2(:2)
        name = nm1(:4)
       
        call TQGSCI( indexc, nm1, iwsg )
        ic1 = indexc
	name = nm2(:4)        
        call TQGSCI( indexc, nm2, iwsg )
        ic2 = indexc
        
        call TQGSCI( indexc, 'O', iwsg)
        ic3 = indexc
        
        write(*,*) 'ic1:', ic1,' ic2:', ic2, 'ic3:',ic3
        call TQGPI( indexp, 'O2GAS', iwsg )
	igas = indexp
	call TQCSP( indexp, 'DORMANT',1.0, iwsg )
c  ------------------------------------------
c       Change component of oxygen to special
c  ------------------------------------------
	call TQGSCI( indexc, 'O', iwsg )
	write(*,*) 'Index of oxygen:',indexc
	ic4 = indexc
c	call TQCSSC( 3, 'SPEC',iwsg )
c  ------------------------------------------
c       Set refernce state to liquid state
c  ------------------------------------------
	call TQSETR(ic4,igas,-1.0D00,200000.0D00,iwsg)
        call TQSETR(ic1, ipliq, -1.0D00, 100000D00,iwsg)
        call TQSETR(ic2, ipliq, -1.0D00, 100000D00,iwsg)
c  -----------------------------------------
c        List status
c  -----------------------------------------
        call TQLS(iwsg)
c  ------------------------------------------
c       Use component names to find indices of components
c       and define initial compositions for bulk phase
c  ------------------------------------------
	j = 1
	k = 1
        do 500 x2b = step, 1.- step, step
        x1b = 1. - x2b
c
	
        write(*,*) 'Indexc: ', ic1         
        write(*,*) 'Indexc: ', ic2
c  -----------------------------------------
c        Remove all conditions
c  -----------------------------------------
        call TQREMAC(iwsg)
        call TQSETC('P',-1,-1,1.0D5,ncp,iwsg)
	call TQSETC('MUR',-1,ic3,0.0D00,ncmur,iwsg)

        if ( j .eq. 1 ) then
c	   call TQSETC('IN', ipliq, ic1,1.0D+0-step, nc1, iwsg)
           j = j + 1
        else
c      	   call TQSETC('IN', ipliq, ic1, 0.0D00, nc1, iwsg)
       	endif
    
c         
        call TQSETC( 'X', -1, ic1, x2b, nc2, iwsg )
	
c        write(*,*) nc1,nc2
c  ------------------------------------------
c       Set temperature for calculation
c  ------------------------------------------
        call TQSETC('T',0,0,temp,nct,iwsg)
	call TQSETC('N', -1,-1,1.0D00,ncn,iwsg)
c	call TQSETC('IN', ipliq,2,step,ncx,iwsg)
c       call TQSETC('IN', ipliq,ic2,1.D00-step,ncx1,iwsg)
        
        write(*,*) '----------Cond. for bulk-----------'
        call TQLC(iwsg)
        write(*,*) '-----------------------------------'
c        read(5,*)
c  ------------------------------------------
c       Execute phase equilibrium calculation for bulk phase
c  ------------------------------------------
        call TQCE('',0,0,0D+0,iwsg)
c  ------------------------------------------
c       check for errors
c  ------------------------------------------
        if (sg2err(ierr)) then
           call RESERR
           call TQCE('',0,0,0D+0,iwsg)
           if (sg2err(ierr)) then
           	write(*,*) 'Error in calculation of equilibrium, 
     +      	check starting conditions'
           endif
        endif
           
        call TQLE(iwsg)
c  ------------------------------------------
c       Get Activity of constituent 1 in "Bulk Phase"
c  ------------------------------------------
        call TQGPCI(ipliq, ic1b, 'CE+4', iwsg)
	write(*,*) 'Index constituent CE for bulk:',ic1b
	call TQGPCI(ipliq, ic2b, 'CO+2', iwsg)
        write(*,*) 'Index constituent CO for bulk:',ic2b
        
        call TQGET1('ACR',-1, ic1b-1, val, iwsg)
        ac1b = val
        write(*,*) 'ac1b',ac1b
c  ------------------------------------------
c       Get Activity of constituent 2 in "Bulk Phase"
c  ------------------------------------------
        call TQGET1('ACR',-1, ic2b-1, val, iwsg)
	ac2b = val
	write(*,*) 'ac2b',ac2b
c  ------------------------------------------
c       Calculation of Surface activites using the Newton-Raphson method
c  ------------------------------------------
	xx2s = 0.0002
	dx2s = 0.0001
c  ------------------------------------------
c       Prepare variables for Newton-Raphson method
c  ------------------------------------------
 300    continue
          x2s(1) = xx2s-dx2s
          x1s(1) = 1.-x2s(1)
c
          x2s(2) = xx2s
          x1s(2) = 1.-x2s(2)
c
          x2s(3) = xx2s+dx2s
          x1s(3) = 1.-x2s(3)
c
	  
        do 400 i=1,3,1
c  ------------------------------------------
c       Input of component name & intial compositions of "surface phase"
c  ------------------------------------------
	
c  --------reinitialize poly module----------
c
        call TQREMAC(iwsg)

c      	call TQREMC(ncp, iwsg)
c      	call TQREMC(ncmur, iwsg)
c      	call TQREMC(nc1, iwsg)
c	call TQREMC(nc2, iwsg)
c	call TQREMC(nct, iwsg)
c	call TQREMC(ncn, iwsg)
c	call TQREMC(ncx, iwsg)
c	call TQREMC(ncx1, iwsg)
c
c  -------------------------------------------
c       Find constituents of liq phase.
c       In this case it is CE and CO constituent indexes. 
c  -------------------------------------------
        call TQGPI(ipliq,'ION',iwsg)
        if (SG2ERR(IERR)) then
           write(*,*) 'Error in getting liq phase'
	endif
c  --------------------------------------------
c       Get number of constituents in liq phase
c  --------------------------------------------
        call TQGNPC(ipliq, iconst, iwsg)
        do  j=1,iconst
	  call TQGPCN(ipliq, j, name, iwsg)
	  print *, j, '   ', name
        end do
        call TQGPCI(ipliq, ic1s, 'CE+4', iwsg)
        write(*,*) 'Index constituent CE:',ic1s
        call TQGPCI(ipliq, ic2s, 'CO+2', iwsg)
        write(*,*) 'Index constituent CO:',ic2s
        if (SG2ERR(IERR)) then
	   write(*,*) 'Error in constituent, 
     +	   check constituents of liquid phase'
	endif
        ic1s=3
        ic2s=2
	call TQSETC('X',ipliq,ic2s,x2s(i),nc3,iwsg)
c        call TQSETC('X',ipliq,ic1,x1s(i),nc3,iwsg)
c        call TQSETC('IN',ipliq,ic2,x2s(i),nc4,iwsg)
c        call TQSETC('IN',ipliq,ic1,x1s(i),nc3,iwsg)
        call TQSETC('T',0,0,temp,nct,iwsg)
        call TQSETC('P',-1,-1,1.0D5,ncp,iwsg)
        call TQSETC('N',-1,-1,1.0D00,ncn,iwsg)
        call TQSETC('MUR',-1,ic3,0.0D00,ncmur,iwsg)

        write(*,*) '---------Cond. for Surface Phase----------'
        call TQLC(iwsg)
        write(*,*) '------------------------------------------'
c        read(5,*)
c  ------------------------------------------
c       Execute phase equilibrium calculation for "Surface Phase"
c  ------------------------------------------
        call TQCE('',0,0,0D+0,iwsg)
        if (SG1ERR(ierr)) goto 900
        call TQLE(iwsg)
c  ------------------------------------------
c       Calc. activity of constituent 1 in "Surface phase"
c  ------------------------------------------
	call TQGPCN(ipliq,ic1s, name, iwsg)
	write(*,*) 'DEBUG:',ic1s
        call TQGET1('AC',ipliq,ic1s-1,val,iwsg)
c  ------------------------------------------
        ac1s(i) = (val/x1s(i))**beta*x1s(i)
        write(*,*) 'ac1s ',name,' ', i,':',val
	write(*,*) ac1s(i)
c  ------------------------------------------
c       Calc. activity of constituent 2 in "Surface phase"
c  ------------------------------------------
 	call TQGPCN(ipliq,ic2s, name, iwsg)
 	write(*,*) 'DEBUG:',ic2s
 	call TQGET1('AC',ipliq, ic2s-1, val, iwsg)
c  ------------------------------------------
        ac2s(i) = (val/x2s(i))**beta*x2s(i)
        write(*,*) 'ac2s ',name,' ',i,':',val
        write(*,*) ac2s(i)
c  ------------------------------------------
c       Calc. of Butler's equation
c  ------------------------------------------
        sl(i) = surf1 + (R*temp/s1)*LOG(ac1s(i)/ac1b)
        sr(i) = surf2 + (R*temp/s2)*LOG(ac2s(i)/ac2b)
 400    continue
c  ------------------------------------------
c       Calculation of Newton-Raphson step
c  ------------------------------------------
        dd = ((sl(3)-sr(3))-(sl(1)-sr(1)))/(2.*dx2s)
        d = (sl(2)-sr(2))
        write(*,*) 'difference of surface tens:',d
        if(d .lt. 0.001) goto 200
        xx2s=x2s(2)-(d/dd)
        write(*,*) 'New concentration:', xx2s
        ddx(1) = (d/dd)+dx2s
	ddx(2) = (d/dd)
	ddx(3) = (d/dd)-dx2s
	  
c  ------------------------------------------
c       Make sure only premited mole fractions are used
c       i.e 0=<xx2s=<1
c  -------------------------------------------
        if(xx2s .lt. 0) xx2s=.001
        if(xx2s .gt. 1.) xx2s=.999
        goto 300        
c
  200   write(13,*) sl(2), sr(2) 
  	surf = ( sl(2) + sr(2) )/2.
c  -------------------------------------------
c       Write present point of curve to plot file
c  -------------------------------------------
  20    format(e12.4,e12.4)
        write(11,20) 1-x2b, surf
c  -------------------------------------------
c       Reinitialize the poly-3 module
c  -------------------------------------------
  500   continue
c  -------------------------------------------
c    Insert value for pure component 2 as last point of curve
c  -------------------------------------------
        write(11,*) 0.0,surf1
        write(11,*)'ENDCURVE'
        stop
  10    write(*,*) 'Error while reading thermochemical data file!'
        stop
  9	write(*,*) 'Error while reading SURFTENS.DAT file !'
	stop    
  11    write(*,*) 'Error while writing SURFTENS.PLO file'
  900   continue
        end

