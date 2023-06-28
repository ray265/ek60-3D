!!!! version 0.24
 
module datagramformat
 implicit none
!configuration header datagram
	type :: headerdatagram
		character*4 deheader
 		integer*8 datetime  
		character*128 surveyname,transectname,soundername
		character*30 version
		character*98 spare
 		integer*4  transducercount
	end type headerdatagram
	type (headerdatagram):: header
!configuration transducer datagram
	type :: transducerdatagram
		character*128 channelid
 		integer*4  beamtype
		real*4 frequency,gain,equivalentbeamangle
		real*4 beam(2),angle(4),pos(3),dir(3)
		real*4 pulselengthtable(5)
		character*8 spare1
		real*4 gaintable(5)
		character*8 spare2
		real*4 sacorrectiontable(5)
		character*8 spare3
		character*16 gptsoftwareversion
		character*28 spare4
	end type transducerdatagram
!sample datagram
	type ::sampledatagram
 		integer*8 datetime  
		integer*2 channel,mode
		real*4 TD,FQ,TP,PL,BW,SI,SV,AC,heave,TxR,TxP,tempera
		integer*2 spare1,spare2
		real*4 RxR,RxP
		integer*4 offset,countnum
	end type sampledatagram
	type (sampledatagram) :: test
!NMEA
	type :: NMEA
		integer*8 datetime
		character*120 lenline
		integer*4 hhmmss,date
		real*4 lon,lat,vel,dir
	end type NMEA
	

end module datagramformat

!! module linearalgebra
module linearalgebra
contains

subroutine Gauss_Jordan(A,S,ANS,Row,Width)
 implicit none
 integer :: Row,Width,i
 real :: A(Row,Width),B(Row,Width),S(Row),ANS(Row)

 B=A
 ANS=S
 call Upper(B,ANS,Row,Width)
 call Lower(B,ANS,Row,Width)

  do i=1,Row
   ANS(i)=ANS(i)/B(i,2)
  end do
 return
end subroutine Gauss_Jordan

subroutine Upper(M,S,Row,Width)
 implicit none
 integer :: Row,Width,i,j
 real :: M(Row,Width),S(Row)
 real :: E
  do i=1,Row-1
     J=i+1
     E=M(j,1)/M(i,2)
     M(j,1:2)=M(j,1:2)-M(i,2:3)*E
     S(j)=S(j)-S(i)*E
  end do
 return
end subroutine Upper

subroutine Lower(M,S,Row,Width)
 implicit none
 integer :: Row,Width,i,j
 real :: M(Row,Width),S(Row)
 real :: E
 do i=Row,2,-1
     J=i-1
     E=M(j,3)/M(i,2)
     M(j,3)=M(j,3)-M(i,2)*E
     S(j)=S(j)-S(i)*E
 end do
 return
end subroutine Lower
end module


!!!!!!!! module spline
module spline
 use linearalgebra
 implicit none
 integer :: POINTS,SEGMENTS,envint
 real,allocatable  :: spln(:),yp(:)
 integer,allocatable  :: xp(:)
 real,allocatable :: a(:),b(:),c(:),d(:),h(:)
 real,allocatable :: matrix(:,:)
 real,allocatable :: f(:),ans(:)

contains

subroutine cublic_spline()
 implicit none
 integer :: i
  do i=1,SEGMENTS
   h(i)=real(xp(i+1))-real(xp(i))
  end do

  do i=2,POINTS-2
   matrix(i,1)=h(i)
   matrix(i,2)=2.0*( h(i+1)+h(i) )
   matrix(i,3)=h(i)
   f(i+1)=6.0*( (yp(i+2)-yp(i+1))/h(i+1) - (yp(i+1)-yp(i))/h(i) )
  end do

 f(1)=0.0
 f(POINTS)=0.0
 ans(1)=0.0
 ans(POINTS)=0.0

 call Gauss_Jordan(matrix,f(2),ans(2),POINTS-2,3)
  do i=1,SEGMENTS
   a(i)=( ans(i+1)-ans(i) )/( 6.0*h(i) )
   b(i)=ans(i)/2.0
   c(i)=( yp(i+1)-yp(i) )/h(i) -( ans(i+1)+2.0*ans(i) )*h(i)/6.0
   d(i)=yp(i)
  end do
 return
end subroutine cublic_spline

real function Spline_Evalue(x)
 implicit none
 integer :: i,x
 real :: diff
 if ( x <=xp(1) ) then
  i=1
 else if ( x>=xp(SEGMENTS) ) then
  i=SEGMENTS
 else
  do i=1,SEGMENTS
   if ( x>=xp(i) .and. x<=xp(i+1) ) exit
  end do
 end if
 diff=real(x)-real(xp(i))
 spline_Evalue=( (a(i)*diff+b(i))*diff+c(i) )*diff+d(i)
return
end function Spline_Evalue


subroutine SampleMax() 
 implicit none
 integer i,j

 xp(1)=1
 xp(POINTS)=size(spln)
 xp(SEGMENTS)=maxloc(spln(2+(SEGMENTS-2)*envint:size(spln)-1),1)+1+(SEGMENTS-2)*envint

 yp(1)=spln(1)
 yp(POINTS)=spln(POINTS)
 yp(SEGMENTS)=maxval(spln(2+(SEGMENTS-2)*envint:size(spln)-1) )

  do j=2,SEGMENTS-1,1
   xp(j)=maxloc(spln(2+(j-2)*envint:1+(j-1)*envint),1)+ 1+(j-2)*envint
   yp(j)=maxval(spln(2+(j-2)*envint:1+(j-1)*envint) )
  end do

 call Cublic_Spline()
return
end subroutine SampleMax

end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! main program
program read_ek60
!use typedef
use datagramformat
use spline
implicit none

		character (len=32) :: arg1,arg2,arg3,tempstr,tempstr1,tempstr2,tempstr3
		integer*4 i,j,k,l,chnum,pic,NMEAlength,maxlocation
 		integer*4  datelengthbegin,datelengthend
 		integer*8  datetime
		character*4 dgheader
		character*1,allocatable :: tracechar(:)
		integer*2,allocatable :: trace(:)
            integer*1,allocatable :: traceangle (:,:)
		real*4 pi,o,p,q
		character*120 tempchar,tempchar1
	type (transducerdatagram),allocatable :: khz(:)
	parameter (pi=4.0*atan(1.0))

!data
	type (sampledatagram),allocatable :: sample1(:),sample2(:)
		integer*4 :: RAWnum(2),RAWlength(2)
		integer*2,allocatable :: RAW1(:,:),mask1(:,:),RAW2(:,:),mask2(:,:),POSnum(:)
            integer*1,allocatable :: RAWangle1(:,:,:),RAWangle2(:,:,:)
		real*4,allocatable :: dB1(:,:),TS1(:,:),Sv1(:,:),depth1(:,:),TG1(:),Sa1(:),temp1(:,:),envelope1(:,:),duration1(:,:),mean1(:)
		real*4,allocatable :: dB2(:,:),TS2(:,:),Sv2(:,:),depth2(:,:),TG2(:),Sa2(:),temp2(:,:),envelope2(:,:),duration2(:,:),mean2(:)
		real*4,allocatable :: diff1(:,:),diff2(:,:)
		real*4,allocatable :: POS(:,:),fixPOS1(:,:,:),fixPOS2(:,:,:),tempPOS(:,:,:),mechanicalAngle1(:,:,:),mechanicalAngle2(:,:,:)
!NMEA
!		integer*8,allocatable :: GGAdt(:),VTGdt(:),GSAdt(:),ZDAdt(:),RMCdt(:),VLWdt(:),GLLdt(:)
!		character*1,allocatable :: GGA(:,:),VTG(:,:),GSA(:,:),ZDA(:,:),RMC(:,:),VLW(:,:),GLL(:,:)
		integer *4 GGAl,VTGl,GSAl,ZDAl,RMCl,VLWl,GLLl
		integer *4 GGAnum,VTGnum,GSAnum,ZDAnum,RMCnum,VLWnum,GLLnum
	type (NMEA),allocatable :: GGA(:),VTG(:),GSA(:),ZDA(:),RMC(:),VLW(:),GLL(:)



!!!!! test the trace number and gps ping number
! read configuration datagram
 ! read header
 call get_command_argument (1,arg1)
 call get_command_argument (2,arg2)
 call get_command_argument (3,arg3)
 arg1=trim(arg1)
 arg2=trim(arg2)
 arg3=trim(arg3)
 read(arg2,*) chnum
 read(arg3,*) pic
 GGAnum=0
 VTGnum=0
 GSAnum=0
 ZDAnum=0
 RMCnum=0
 VLWnum=0
 RAWnum=0
 GLLnum=0
 j=0
! open(unit=9910,file="zzzcx")
 open(unit=10,file=arg1,form='UNFORMATTED',access='stream')
	read(10) datelengthbegin,header
 	allocate(khz(header%transducercount))
	write(*,*) "file:",arg1
	write(*,*) "transducercount:",header%transducercount
 	read(10) khz,datelengthend
	write(*,*) "channelID:"
	write(*,'(G20.0)') khz%channelid

!read NMEA and sample
do while (.true.)
	read(10,end=100) datelengthbegin,dgheader
!	write (9910,*) dgheader,datelengthbegin
 if (dgheader == "NME0") then
	NMEAlength=datelengthbegin-12
 	allocate(tracechar(NMEAlength))
	read(10) datetime,tracechar,datelengthend
	 select case(tracechar(4)//tracechar(5)//tracechar(6))
	  case("GGA")
		if (GGAl < NMEAlength) then
		  GGAl=NMEAlength
		end if
		  GGAnum=GGAnum+1
	  case("VTG")
		if (VTGl < NMEAlength) then
		  VTGl=NMEAlength
		end if
		  VTGnum=VTGnum+1
	  case("GSA")
		if (GSAl < NMEAlength) then
		  GSAl=NMEAlength
		end if
		  GSAnum=GSAnum+1
	  case("ZDA")
		if (ZDAl < NMEAlength) then
		  ZDAl=NMEAlength
		end if
		  ZDAnum=ZDAnum+1
	  case("RMC")
		if (RMCl < NMEAlength) then
		  RMCl=NMEAlength
		end if
		  RMCnum=RMCnum+1
	  case("VLW")
		if (VLWl < NMEAlength) then
		  VLWl=NMEAlength
		end if
		  VLWnum=VLWnum+1	
	  case("GLL")
		if (GLLl < NMEAlength) then
		  GLLl=NMEAlength
		end if
		  GLLnum=GLLnum+1
	 end select
	deallocate(tracechar)
 else if (dgheader == "RAW0") then 
	read(10) test
    !! find max RAWlenght in each channel
  if ( test%countnum .GT. RAWlength(test%channel) ) then
	RAWlength(test%channel)=test%countnum
  end if

   if ( test%mode == 3 ) then
	allocate(trace(test%countnum),traceangle(2,test%countnum))
      read(10) trace,traceangle,datelengthend
      deallocate(trace,traceangle)
   else
	allocate(trace(test%countnum))
      read(10) trace,datelengthend
	deallocate(trace)
   end if
	write (9910,*) dgheader,datelengthbegin,datelengthend,test%countnum
 RAWnum(test%channel)=RAWnum(test%channel)+1
 end if 
!!! test datelength
 if (datelengthend /= datelengthbegin) then 
  write(*,*) "WARNING!!!!!!datalength not equivalent, stop program"
  stop
 end if

end do
100 continue
 close(10)
write(*,*) "(GPS number):GGA       VTG       GSA       ZDA       RMC       GLL"
write(*,*) GGAnum,VTGnum,GSAnum,ZDAnum,RMCnum,GLLnum
write(*,*) "(Length):GGA       VTG       GSA       ZDA       RMC       GLL"
write(*,*) GGAl,VTGl,GSAl,ZDAl,RMCl,GLLl
write(*,*) "VLW (heading number):",VLWnum,"RAWnum:",RAWnum
write(*,*) "VLWl:",VLWl,"RAWlength:",RAWlength
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! get  GGAnum,VLWnum,RAWlength,RAWnum end





!!!!!!!!!!!!!! read file
! allocate(GGA(GGAnum,GGAl),VTG(VTGnum,VTGl),GSA(GSAnum,GSAl),ZDA(ZDAnum,ZDAl),RMC(RMCnum,RMCl),VLW(VLWnum,VLWl),GLL(GLLnum,GLLl))
 allocate(GGA(GGAnum),VTG(VTGnum),GSA(GSAnum),ZDA(ZDAnum),RMC(RMCnum),VLW(VLWnum),GLL(GLLnum))
! allocate(GGAdt(GGAnum),VTGdt(VTGnum),GSAdt(GSAnum),ZDAdt(ZDAnum),RMCdt(RMCnum),VLWdt(VLWnum),GLLdt(GLLnum))
 allocate(RAW1(RAWnum(1),RAWlength(1)),TS1(RAWnum(1),RAWlength(1)),depth1(RAWnum(1),RAWlength(1)),fixPOS1(3,RAWnum(1),RAWlength(1)))
 allocate(RAW2(RAWnum(2),RAWlength(2)),TS2(RAWnum(2),RAWlength(2)),depth2(RAWnum(2),RAWlength(2)),fixPOS2(3,RAWnum(2),RAWlength(2)))
 allocate(Sv1(RAWnum(1),RAWlength(1)),dB1(RAWnum(1),RAWlength(1)),temp1(RAWnum(1),RAWlength(1)))
 allocate(Sv2(RAWnum(1),RAWlength(1)),dB2(RAWnum(1),RAWlength(1)),temp2(RAWnum(1),RAWlength(1)))
 allocate(RAWangle1(2,RAWnum(1),RAWlength(1)),RAWangle2(2,RAWnum(2),RAWlength(2)),tempPOS(2,RAWnum(1),RAWlength(1)))
 allocate(mechanicalAngle1(2,RAWnum(1),RAWlength(1)),mechanicalAngle2(2,RAWnum(2),RAWlength(2)))
 allocate(sample1(RAWnum(1)),TG1(RAWnum(1)),Sa1(RAWnum(1)),sample2(RAWnum(2)),TG2(RAWnum(2)),Sa2(RAWnum(2)))

 allocate(envelope1(RAWnum(1),RAWlength(1)))
 allocate(envelope2(RAWnum(2),RAWlength(1)))
 allocate(mean1(RAWlength(1)),mask1(RAWnum(1),RAWlength(1)),diff1(RAWnum(1),RAWlength(1)))
 allocate(mean2(RAWlength(2)),mask2(RAWnum(2),RAWlength(2)),diff2(RAWnum(2),RAWlength(2)))
 allocate(duration1(RAWnum(1),4))
 allocate(duration2(RAWnum(2),4))
 allocate(POSnum(RAWnum(1)),POS(RAWnum(1),4))

 duration1=huge(duration1)
 duration1=huge(duration2)
 mask1=-1
 mask2=0
 mean1=0
 mean2=0
 GGAnum=0
 VTGnum=0
 GSAnum=0
 ZDAnum=0
 RMCnum=0
 VLWnum=0
 RAWnum=0
 GLLnum=0
 RAW1=-9999
 RAW2=-9999
  !! read header , transducer header
 open(unit=10,file=arg1,form='UNFORMATTED',access='stream')
	read(10) datelengthbegin,header,khz,datelengthend
  !! read NMEA and sample
do while (.true.)
	read(10,end=200) datelengthbegin,dgheader
 if (dgheader == "NME0") then
	NMEAlength=datelengthbegin-12
 	allocate(tracechar(NMEAlength))
  !! read NMEA to array
	read(10) datetime,tracechar,datelengthend
           tempchar=''
   !!array to string
	  do i=1,NMEAlength-2,1
             tempchar=trim(tempchar)//tracechar(i)
	  end do
    !! replace camma to space
        do i=1,len_trim(tempchar),1
	   if(tempchar(i:i) .eq. ',') then
	    tempchar(i:i)=' '
	   end if 
	  end do

	 select case(tempchar(4:6))
!	  case("GGA")
!		GGAnum=GGAnum+1
!		GGAdt(GGAnum)=datetime
!		GGA(GGAnum,:)=tracechar(:)
!	  case("VTG")
!		VTGnum=VTGnum+1
!		VTGdt(VTGnum)=datetime
!		VTG(VTGnum,:)=tracechar(:)
!	  case("GSA")
!		GSAnum=GSAnum+1
!		VTGdt(GSAnum)=datetime
!		GSA(GSAnum,:)=tracechar(:)
!	  case("ZDA")
!		ZDAnum=ZDAnum+1
!		ZDAdt(ZDAnum)=datetime
!		ZDA(ZDAnum,:)=tracechar(:)
	  case("RMC")
		RMCnum=RMCnum+1
		RMC(RMCnum)%datetime=datetime
		RMC(RMCnum)%lenline=trim(tempchar)
		read(tempchar,*) tempchar1,RMC(RMCnum)%hhmmss,tempchar1,RMC(RMCnum)%lat,tempchar1,RMC(RMCnum)%lon&
                            &,tempchar1,RMC(RMCnum)%vel,RMC(RMCnum)%dir,RMC(RMCnum)%date
!	  case("VLW")
!		VLWnum=VLWnum+1
!		VLWdt(VLWnum)=datetime
!		VLW(VLWnum,:)=tracechar(:)
!	  case("GLL")
!		GLLnum=GLLnum+1
!		GLLdt(GLLnum)=datetime
!		GLL(GLLnum,:)=tracechar(:)
	 end select
	deallocate(tracechar)
 else if (dgheader == "RAW0") then 
   !! sample header select
	read(10) test
	RAWnum(test%channel)=RAWnum(test%channel)+1
	allocate(trace(test%countnum),traceangle(2,test%countnum))
   if ( test%mode == 3 ) then
      read(10) trace,traceangle,datelengthend
   else
      read(10) trace,datelengthend
   end if

	select case(test%channel)
	 case(1)
		sample1(RAWnum(1))=test
	 	 if(test%mode == 3 ) then 
		   RAW1(RAWnum(1),1:test%countnum) =trace
               RAWangle1(:,RAWnum(1),1:test%countnum)=traceangle
		    deallocate(trace,traceangle)
	 	 else
	         RAW1(RAWnum(1),1:test%countnum)=trace
		    deallocate(trace)
	 	 end if
	 case(2)
		sample2(RAWnum(2))=test
	 	 if(test%mode == 3 ) then 
	         RAW2(RAWnum(2),1:test%countnum)=trace
               RAWangle2(:,RAWnum(2),1:test%countnum)=traceangle
		    deallocate(trace,traceangle)
	 	 else
		   RAW2(RAWnum(2),1:test%countnum)=trace
		    deallocate(trace)
	 	 end if
	end select
 end if 
end do

200 continue
 close(10)
write(*,*) "read over"

!!!! calculate Position(Lon Lat vel direction)
k=1
do i=1,RAWnum(1),1
 do j=k,RMCnum,1
   if ( sample1(i)%datetime .LE. RMC(j)%datetime ) then
     POSnum(i)=j
       if ( j-1 .eq. 0   ) then
            POS(i,1)=RMC(j)%lon
		POS(i,2)=RMC(j)%lat 
	 end if
     k=j
     o=real(sample1(i)%datetime-RMC(j-1)%datetime)/real(RMC(j)%datetime-RMC(j-1)%datetime)
     p=floor(RMC(j-1)%lon/100)+mod(RMC(j-1)%lon,100.0)/60
     q=floor(RMC(j)%lon/100)+mod(RMC(j)%lon,100.0)/60
     POS(i,1)=p+o*(q-p)
     p=floor(RMC(j-1)%lat/100)+mod(RMC(j-1)%lat,100.0)/60
     q=floor(RMC(j)%lat/100)+mod(RMC(j)%lat,100.0)/60
     POS(i,2)=p+o*(q-p)
     POS(i,3)=RMC(j-1)%vel+o*(RMC(j)%vel-RMC(j-1)%vel)
!!!  may have bug for degree?
     POS(i,4)=RMC(j-1)%dir+o*(RMC(j)%dir-RMC(j-1)%dir)
     goto 777
    end if
 end do 
 777 continue



end do 


!!!! calculate fixPOS(x,y,z)  >ENZ 
!!!! mechanicalAngle(degree  along athwar)
mechanicalAngle1(1,:,:)=real(RAWangle1(1,:,:)*180.0/128.0)/khz(1)%angle(1)-khz(1)%angle(3)
mechanicalAngle1(2,:,:)=real(RAWangle1(2,:,:)*180.0/128.0)/khz(1)%angle(2)-khz(1)%angle(4)
mechanicalAngle2(1,:,:)=real(RAWangle2(1,:,:)*180.0/128.0)/khz(2)%angle(1)-khz(2)%angle(3)
mechanicalAngle2(2,:,:)=real(RAWangle2(2,:,:)*180.0/128.0)/khz(2)%angle(2)-khz(2)%angle(4)


!!!! calculate dB,TS,Sv
!!!find transducer peak gain
write(*,*) "calculate channel1"
!!! channel == 1
do i=1,RAWnum(1)
 do j=1,5
  if (sample1(i)%PL == khz(1)%pulselengthtable(j) )  then
	TG1(i)=khz(1)%gaintable(j)
	Sa1(i)=khz(1)%sacorrectiontable(j)
  end if
 end do
end do

!! calculate R (distance)

forall(i=1:RAWnum(1),j=1:RAWlength(1))
   depth1(i,j)=sample1(i)%SI*1500*j/2
end forall

!!!!  Z,y,x   tempPOS(1)-athwar-y    POS(2)-along-x 
 fixPOS1(3,:,:)=depth1(:,:)/sqrt(tan(mechanicalAngle1(1,:,:)*pi/180)**2+tan(mechanicalAngle1(2,:,:)*pi/180)**2+1) 
  !! y_athwar (alpha)
   tempPOS(2,:,:)=tan(mechanicalAngle1(1,:,:)*pi/180)*fixPOS1(3,:,:)
  !! x_along (beta)
   tempPOS(1,:,:)=tan(mechanicalAngle1(2,:,:)*pi/180)*fixPOS1(3,:,:)
do i=1,RAWnum(1)
 !!  N
  fixPOS1(2,i,:)=tempPOS(1,i,:)*cos(POS(i,4)*pi/180)-tempPOS(2,i,:)*sin(POS(i,4)*pi/180)
 !!  E
  fixPOS1(1,i,:)=tempPOS(1,i,:)*sin(POS(i,4)*pi/180)+tempPOS(2,i,:)*cos(POS(i,4)*pi/180)
end do

!!dB,TS,Sv
write(*,*) "dB,TS,SV1..."
!open(98765,file='zzxc_xxxx')
!open(98764,file='zzxc_xxxxx')
dB1=RAW1*10*log10(2.0)/256
do i=1,RAWnum(1)
 temp1(i,:)=dB1(i,:)+2*sample1(i)%AC*depth1(i,:)&
         &-10*log10(sample1(i)%TP*(10**(TG1(i)/ 10)*sample1(i)%SV/sample1(i)%FQ/4/pi)**2)
 TS1(i,:)=temp1(i,:)+40*log10(depth1(i,:))
 Sv1(i,:)=temp1(i,:)+20*log10(depth1(i,:))&
	   &-10*log10((sample1(i)%SV*sample1(i)%PL*10**(khz(1)%equivalentbeamangle/10))/2)-2*Sa1(i)

! do j=700,RAWlength(1),1
!   If (Sv1(i,j) .GE. -50) then
!      write(98765,*) i,depth1(i,j)
!        exit
!   end if 
! end do
! do j=RAWlength(1),700,-1
!   If (Sv1(i,j) .GE. -50) then
!      write(98764,*) i,depth1(i,j)
!        exit
!   end if 
! end do
end do

!!!envelope
 envint=20
 SEGMENTS=ceiling(real((RAWlength(1)-2))/real(envint))+1
 POINTS=SEGMENTS+1

 allocate(a(SEGMENTS),b(SEGMENTS),c(SEGMENTS),d(SEGMENTS),h(SEGMENTS))
 allocate(matrix(POINTS-2,3),spln(RAWlength(1)))
 allocate(xp(POINTS),yp(POINTS),f(POINTS),ans(POINTS))


 do i=1,RAWnum(1)
  spln=Sv1(i,:) 
  call SampleMax()
  do j=1,RAWlength(1)
   envelope1(i,j)= Spline_Evalue(j)
  end do
 end do 


!!  duration   50dB 1.maxloc(m) 2.maxvalue(dB) 3.upper limit(m) 4. lower limit(m)
write(*,*) "duration1......"
do i=1,RAWnum(1),1
maxlocation=maxloc(Sv1(i,50:RAWlength(1)-5),1)+49
 duration1(i,1)= depth1(i,maxlocation)
 duration1(i,2)= Sv1(i,maxlocation)

if ( duration1(i,2) .lt. -60) then
 mask1(i,:)=-999
 cycle
end if
 !!upper limit & mask value below the upper limit to -999
 j=1
 o=0.0
 p=0.0
  do while ( o .GE. -55.0 .or. p .ge. o )
   p=o 
   o=maxval(envelope1(i,maxlocation-j-8:maxlocation-j))
   j=j+1
  end do
  duration1(i,3)=depth1(i,maxlocation-j+1)
  mask1(i,maxlocation-j+1:)=-999
  mean1(1:maxlocation-j)=mean1(1:maxlocation-j)+Sv1(i,1:maxlocation-j)
 !!lower limit
 j=1
 o=0.0

 if ( maxval(envelope1(i,maxlocation:RAWlength(1))) .GE. -55.0 )  then
   do while ( o .GE. -55.0  ) 
      if ( (maxlocation+j+8) .GE. RAWlength(1) ) then
        o=maxval(envelope1(i,maxlocation+j:RAWlength(1)))
      else
        o=maxval(envelope1(i,maxlocation+j:maxlocation+j+8))
      end if 
     j=j+1 
   end do
   duration1(i,4)=depth1(i,maxlocation+j-1)
 else
   duration1(i,4)=depth1(i,RAWlength(1))
 end if
end do
  mean1=mean1/(RAWnum(1)-count(mask1==-999,1))

!!!

 do j=1,RAWlength(1),1
   if ( depth1(1,j) .ge. 1000.0 .and. mean1(j) .ge. -95.0) then 
     mean1(j)=-95.0
   end if
 end do


  !!set mask (dectect  plump  dB   mask=0)  

do i=1,RAWnum(1),1
 do j=1,RAWlength(1),1
   if (mask1(i,j) .ge. -1 ) then
    diff1(i,j)=Sv1(i,j)- mean1(j)
    if (Sv1(i,j) .GE. -90) then
     mask1(i,j)=0
      if ( diff1(i,j) > 12.0 ) then 
       mask1(i,j)=1
      end if
     end if
   end if
 end do
end do

!!!!! decrease noise point   mask value (=1 great than -90 and diff great than 10, =2  no neighbor grid great than 10)

forall (i=1:RAWnum(1),j=2:RAWlength(1)-1, (count(mask1(i-2:i+2,j-2:j+2) == 1) > 6 .AND. mask1(i,j) == 1) )  mask1(i,j)=2

!!!!! channel 1 end




write(*,*) "calculate channel2"
!!!! channel 2 
if ( header%transducercount .eq. chnum) then


do i=1,RAWnum(2)
 do j=1,5
  if (sample2(i)%PL == khz(1)%pulselengthtable(j) )  then
	TG2(i)=khz(2)%gaintable(j)
	Sa2(i)=khz(2)%sacorrectiontable(j)
  end if
 end do
end do


!! calculate R (distance)
do i=1,RAWnum(2)
 do j=1,RAWlength(2)
   depth2(i,j)=sample2(i)%SI*1500*j/2
 end do 
end do

!!dB,TS,Sv
write(*,*) "dB,TS,SV2..."
  dB2=RAW2*10*log10(2.0)/256
do i=1,RAWnum(2)
 temp2(i,:)=dB2(i,:)+2*sample2(i)%AC*depth2(i,:)&
         &-10*log10(sample2(i)%TP*(10**(TG2(i)/ 10)*sample2(i)%SV/sample2(i)%FQ/4/pi)**2)
 TS2(i,:)=temp2(i,:)+40*log10(depth2(i,:))
 Sv2(i,:)=temp2(i,:)+20*log10(depth2(i,:))&
	   &-10*log10((sample2(i)%SV*sample2(i)%PL*10**(khz(2)%equivalentbeamangle/10))/2)-2*Sa2(i)

!! upper envelope
  !!set mask (dectect  pump   _70dB)  
  do j=1,floor(RAWlength(2)/real(envint)),1
      maxlocation=maxloc(Sv2(i,1+(j-1)*envint:1+j*envint),1)+(j-1)*envint
!	envelope2(i,1,j)=depth2(i,maxlocation)
!	envelope2(i,2,j)=Sv2(i,maxlocation)
!	 if (envelope2(i,2,j) .GE. -100 ) then
!	   mask2(i,j)=1
!	 end if
  end do
end do


!!  duration   50dB 1.maxloc(m) 2.maxvalue(dB) 3.upper limit(m) 4. lower limit(m)
write(*,*) "duration2......"
do i=1,RAWnum(2),1
!maxlocation=maxloc(envelope2(i,2,7:),1)+6
! duration2(i,1)= envelope2(i,1,maxlocation)
! duration2(i,2)= envelope2(i,2,maxlocation)

 !!upper limit & mask value below the upper limit to -999
 j=1
 o=0.0
  do while ( o .GE. -60.0) 
   o=Sv2(i,maxlocation-j)
   j=j+1
  end do

  duration2(i,3)=Sv2(i,maxlocation-j+1)
  mask2(i,maxlocation-j+1:)=-999
  mean2(1:maxlocation-j)=mean2(1:maxlocation-j)+Sv2(i,1:maxlocation-j)


 !!lower limit
 j=1
 o=0.0
  do while ( o .GE. -60.0) 
!   o=envelope2(i,2,maxlocation+j)
   j=j+1
  end do
!  duration2(i,4)=envelope2(i,1,maxlocation+j-1)
end do

end if
!!!! channel 2 end
write(*,*) "calculate over"





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! output

open(901,file='pos')
open(1101,file='mask1_1')
open(1301,file='mean1')
open(1401,file='mask1_2')
open(1001,file='zz_testpos')




!! mean
write(*,*) "output mean1"

 do i=1,size(mean1),1
    write(1301,'(2G20.10)') depth1(1,i),mean1(i)
 end do

!!mask
write(*,*) "output mask"
l=0
do i=1,RAWnum(1),1
 k=0
 do j=1,RAWlength(1),1
  select case (mask1(i,j))
   case (1)
   !!mask1
      write(1101,*) i,POS(i,1:2),depth1(i,j),Sv1(i,j),sample1(i)%datetime
   case (2)
	k=k+1
   !!mask2
      write(1401,*) i,POS(i,1:2),depth1(i,j),Sv1(i,j),fixPOS1(1:3,i,j),duration1(i,3),sample1(i)%datetime
  end select
 end do
 !!pos: i,lon,lat,depth,max dB, upper depth,lower depth,sum
 write(901,'(I6,2F14.8,4(F12.3,1X),2I6,I20)') i,POS(i,1:2),duration1(i,:),k,k+l,sample1(i)%datetime
 l=k
end do


!! output envelope
 if ( pic .eq. 1 ) then
	write(*,*) "output trace &envelope "
	do i=1,RAWnum(1),1
	  write(tempstr,*) i
	  tempstr1='trace1_'//adjustl(tempstr)
	  tempstr2='envelope1_'//adjustl(tempstr)
	  open(701,file=tempstr1)
	  open(801,file=tempstr2)
	    do j=2,RAWlength(1),1
   	     write(701,*) sample1(i)%SI*j,depth1(i,j),RAW1(i,j),dB1(i,j),TS1(i,j),Sv1(i,j)
 	     write(801,'(2G15.5)') depth1(i,j),envelope1(i,j)
 	    end do
	end do
 else if ( pic .eq. 2 ) then
	write(*,*) "output xyz_file"
!! output xyz, angular (along athw)
	open(1201,file='xyz')
	open(2201,file='xyz_bi',FORM='UNFORMATTED',ACCESS='STREAM')
	open(2401,file='xyz_bi_num',FORM='UNFORMATTED',ACCESS='STREAM')
	open(2301,file='angle')
	do i=1,RAWnum(1),1
	    do j=2,RAWlength(1),1
           write(1201,*) POS(i,1:2),depth1(i,j),Sv1(i,j)
           write(2201) POS(i,1),depth1(i,j),Sv1(i,j)
           write(2401) i,depth1(i,j),Sv1(i,j)
           write(2301,*) POS(i,1:2),depth1(i,j),(real(RAWangle1(1,i,j))*180.0/128.0)/khz(1)%angle(1)-khz(1)%angle(3)&
                         &,(real(RAWangle1(2,i,j))*180.0/128.0)/khz(1)%angle(2)-khz(1)%angle(4)
 	    end do
	end do
 end if

!!NMEA

!if (GGAnum .ne. 0) then
! open(1501,file='zz_GGA')
! do i=1,GGAnum,1
!  write(1501,*) GGA(i,:)
! end do
! close(1501)
!end if

!if (VTGnum .ne. 0) then
! open(1601,file='zz_VTG')
! do i=1,VTGnum,1
!  write(1601,*) VTG(i,:)
! end do
! close(1601)
!end if

!if (GSAnum .ne. 0) then
! open(1701,file='zz_GSA')
! do i=1,GSAnum,1
!  write(1701,*) GSA(i,:)
! end do
! close(1701)
!end if

!if (ZDAnum .ne. 0) then
! open(1801,file='zz_ZDA')
! do i=1,ZDAnum,1
!  write(1801,*) ZDA(i,:)
! end do
! close(1801)
!end if

if (RMCnum .ne. 0) then
 open(1901,file='RMC')
 do i=1,RMCnum,1
   write(1901,102) RMC(i)%hhmmss,RMC(i)%date,RMC(i)%lat,RMC(i)%lon,RMC(i)%vel,RMC(i)%dir,RMC(i)%datetime
 end do
 close(1901)
end if
102 FORMAT(I10.6,I10.6,F15.6,F16.6,F8.2,F7.1,I24.18)



!if (GLLnum .ne. 0) then
! open(2001,file='zz_GLL')
! do i=1,GLLnum,1
!  write(2001,*) GLL(i,:)
! end do
! close(2001)
!end if

!if (VLWnum .ne. 0) then
! open(2101,file='zz_VLW')
! do i=1,VLWnum,1
!  write(2101,*) VLW(i,:)
! end do
! close(2101)
!end if


write(*,*) "RawNum_1",RAWnum(1)
write(*,*) "transducer_1 header"
write(*,*) khz(1)%beamtype
write(*,*) khz(1)%equivalentbeamangle,khz(1)%beam(1),khz(1)%beam(2)
write(*,*) khz(1)%angle(:)
write(*,*) khz(1)%pos(:)
write(*,*) khz(1)%dir(:)
write(*,*) khz(1)%pulselengthtable(:)
write(*,*) khz(1)%gaintable(:)
write(*,*) khz(1)%sacorrectiontable(:)
write(*,*) "sample_1 header"
write(*,*) sample1(1)%TD, sample1(1)%FQ,sample1(1)%TP,sample1(1)%PL,sample1(1)%BW
write(*,*) sample1(1)%SI, sample1(1)%SV,sample1(1)%AC,sample1(1)%heave
write(*,*) sample1(1)%TxR, sample1(1)%TxP,sample1(1)%tempera,sample1(1)%spare1,sample1(1)%spare2
write(*,*) sample1(1)%RxR, sample1(1)%RxP,sample1(1)%offset,sample1(1)%countnum



stop
end program


