module param
! this module is used to define the dimensions of large domain and the time duration of available data
  integer,parameter :: bigx=149          !number of grids in x direction
  integer,parameter :: bigy=129          !number of grids in y direction
  integer,parameter :: days=152          !number of days for input fields E, W, U, V and P
  real   ,parameter :: dx=25000.00       !unit(m), the grid-spacing of each grid in x direction
  real   ,parameter :: dy=25000.00       !unit(m), the grid-spacing of each grid in x direction
  integer,parameter :: day1=16           !starting date for drm simulation
  integer,parameter :: day2=152          !ending date for drm simulation
  integer,parameter :: domsize=950       !number of grids within the target region
  integer,parameter :: region_num=5      !number of subregions in map file
  integer,parameter :: dt=1800           !unit(second), time interval of each back-tracing step
  integer,parameter :: max_tracing=1000  ! maximum times of back trajactory(>48*15)
  integer,parameter :: max_iteration=100 ! maximum times of iteration
  integer,parameter :: warmup=15         !days that to be excluded at the beginning used for back-tracking (= day1-1)
  integer,parameter :: velocity_t=3      !unit(hr), temporal resolution for velocity U and V fields
  integer,parameter :: valve1=6          !times of tracing per unit velocity time interval (=3600*3/1800)
  integer,parameter :: valve2=48         !times of tracing per day (=3600*24/1800)
  integer,parameter :: NAN_value=-9.99E8
  integer,parameter :: ithyear=2004      !year to be simulated
end module param
!----------------------------------------------------
module typedef
  implicit none
  type start_end
    integer :: start_l,end_l,renum
  end type
end module

!====================================================
program Recycling

! load data first
! data includes daily data and several hourly data(in this case 3 hourly):
! Daily: precipitable water(PW),evaportranspiration(ET),precipitation(PP)
! 3-hourly: u and v
  use param
  use typedef
  implicit none
  real, dimension(bigx,bigy,days,3)                 :: grid_days    ! 1-PW,2-ET,3-PP
  real, dimension(bigx,bigy,days*24/velocity_t,4:5) :: grid_hours   ! 4-U3,5-V3
! thus grid_days(:,:,:,1)=PW(mm), grid_days(:,:,:,2)=ET(mm/3h), grid_days(:,:,:,3)=PP(mm/day)
! grid_hours(:,:,:,4)=U3(m/s),grid_hours(:,:,:,5)=V3(m/s)
! modeling date from day1 to day2, domain from x1,y1 to x2,y2
  real                                              :: TPW=0,TET=0,TPP=0,TU_1=0,TV_1=0,TU_2=0,TV_2=0
  real, dimension(max_tracing,3)                    :: EW_vector=0 !matrix used to save the E/W data of each grid on specific day and the region it belongs as well
  real                                              :: E_W=0  !variable save integrated EW_vector of each grid on specific day
  real, dimension(domsize,day1:day2)                :: rrcalc  !matrix to save recycling ratio of each grid on specific day
  real, dimension(max_tracing)                      :: ppx,ppy  !variable used to record the position of tracing
  real, dimension(day1:day2)                        :: rr
  integer                                           :: i,j,k,day,time,it,intx,inty,intx1,inty1,flag1,flag2,rn,t
  real                                              :: oldppx,oldppy,diffx,diffy
  integer                                           :: daycount,hourcount,count1,count2
  real                                              :: temprr,sumrr,sumpp,s

  type(start_end), dimension(max_tracing)           :: reg_count
  real           , dimension(domsize,region_num)    :: grid_rr     ! record contribution to one grid box
  real           , dimension(day1:day2,region_num)  :: daily_rr  !record daily recycling ratio of different regions
  integer        , dimension(bigx,bigy)             :: map
  integer        , dimension(domsize,2)             :: domainij
  real           , dimension(region_num)            :: sum_region
  real           , dimension(domsize)               :: precip_grid=0
  real           , allocatable                      :: various_reg(:,:)
  integer                                           :: reg_change,old_reg,new_reg,zeroi,ein,xx,ind,L1,L2
  real                                              :: s1,pral

! build files to output data from running the model
  call buildfile

! load data into array grid_days and grid_hours
  call openfile(grid_days,grid_hours)  !including fixing data <0
  call readmap(map)
  call monsoon(domainij)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!############################################################################################
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! start the daily loop during the chosen duration from day1 to day2
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'!!!  Daily Loop Starts  !!!'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  daily : do day=day1,day2

    print *,'day:', day
    write(500,*) day    ! file 500 is output to record information for each back-trajectory

    !initialized vars to be zero
    sumrr=0
    sumpp=0
   
    do xx=1,domsize 
      precip_grid(xx)=0
    end do    

    !start the domain loop
    domain : do k=1,domsize
        i=domainij(k,1)
        j=domainij(k,2)
        ppx(1)=(i-1)*dx+dx/2
        ppy(1)=(j-1)*dy+dy/2
        TPP=grid_days(i,j,day,1)
        E_W=0
        
        write(500,*) k,day,ppx(1),ppy(1)
        
        count2=0
        reg_change=0
        hourcount=day*24/velocity_t   !starting temporal position in velocity files
        daycount=day                  !starting temporal position in other files
     
        ! initialize vars to be zero
        do zeroi=1,max_tracing
          do ein=1,3
            EW_vector(zeroi,ein)=0
          end do
          reg_count(zeroi)%start_l=0
          reg_count(zeroi)%end_l=0
          reg_count(zeroi)%renum=0
        end do

!====================================================================================================
!****************************************************************************************************
!====================================================================================================
        ! start the time loop
        timeloop : do time=1,max_tracing

          ! tell ppx and ppy are within large domain or not-----------------------------------
          if ((ppx(time)<=0).or.(ppx(time)>=(dx*bigx))) then
           ! print *,'Beyond boundary during timeloop when time is:',time
            exit timeloop
          else if ((ppy(time)<=0).or.(ppy(time)>=(dy*bigy))) then
           ! print *,'Beyond boundary during timeloop when time is:',time
            exit timeloop
          end if
          !-----------------------------------------------------------------------------------

          !read the Evaporation,precipitable water and wind velocity by subroutine 'finduvEW'
          call finduvEW(ppx(time),ppy(time),hourcount,daycount,grid_days(:,:,daycount,1),grid_days(:,:,daycount,2),&
              grid_hours(:,:,hourcount,4),grid_hours(:,:,hourcount,5),TPW,TET,TU_1,TV_1,intx1,inty1,flag1)

          if (flag1==1) then
            exit timeloop
          else
            EW_vector(time,1)=(time-1)*dt
            EW_vector(time,2)=TET/TPW
            EW_vector(time,3)=map(intx1,inty1)
            !TET in unit of mm/3h thus TET*8 and EW_vector*8 has the unit of mm/day and day-1 individually
          end if
 
          if (time .eq. 1) then
            reg_count(1)%start_l=1
            reg_count(1)%renum=map(i,j)
            reg_change=1
          else
            old_reg=EW_vector(time-1,3)
            new_reg=EW_vector(time,3)
            if (old_reg .ne. new_reg) then
              reg_change=reg_change+1
              reg_count(reg_change-1)%end_l=time
              reg_count(reg_change)%start_l=time
              reg_count(reg_change)%renum=new_reg      
            end if
          end if
  
          ! find the next point using Iterative Technique--------------------------------------
          call finduv(ppx(time),ppy(time),hourcount,grid_hours(:,:,hourcount,4),& 
               grid_hours(:,:,hourcount,5),TU_2,TV_2)
          oldppx=ppx(time)-(TU_1+TU_2)*dt/2
          oldppy=ppy(time)-(TV_1+TV_2)*dt/2
          if ((oldppx<=0).or.(oldppx>=(dx*bigx))) exit timeloop
          if ((oldppy<=0).or.(oldppy>=(dy*bigy))) exit timeloop
!================================================================================================
          iteration : do it=1,max_iteration
            !tell if the point still within the large domain
            if ((oldppx<=0).or.(oldppx>=(dx*bigx))) then
            !  print *,'Beyond boundary during iteration when it is',it
              exit iteration
            else if ((oldppy<=0).or.(oldppy>=(dy*bigy))) then
            !  print *,'Beyond boundary during iteration when it is',it
              exit iteration
            end if           

            call finduv(oldppx,oldppy,hourcount,grid_hours(:,:,hourcount,4),&
                 grid_hours(:,:,hourcount,5),TU_2,TV_2,intx,inty,flag2)

            if (flag2==1) then
              exit iteration
            else
              ppx(time+1)=ppx(time)-(TU_1+TU_2)*dt/2
              ppy(time+1)=ppy(time)-(TV_1+TV_2)*dt/2
            
              diffx=abs(ppx(time+1)-oldppx)
              diffy=abs(ppy(time+1)-oldppy)
              if ((diffx<0.0001).and.(diffy<0.0001)) exit  iteration
              oldppx=ppx(time+1)
              oldppy=ppy(time+1)
            end if
          end do iteration
!================================================================================================
        ! change the temporal data to use based on back tracing
          count1=count1+1
          if (mod(count1,valve1)==0) then
            hourcount=hourcount-1
            if(hourcount<=0) exit timeloop
          end if
          count2=count2+1
          if (mod(count2,valve2)==0)  then
            daycount=daycount-1
            if(daycount<=0) exit timeloop
          end if
 
        end do timeloop

!===================================================================================================
!***************************************************************************************************
!===================================================================================================
        ! do the integral of the entire path
        allocate(various_reg(reg_change,3))
        pral=1
        if (count2>0) then
          reg_count(reg_change)%end_l=count2
          call SPTNQ(EW_vector(1:count2,1),EW_vector(1:count2,2),count2,s)
          E_W=s
          do ein=1,reg_change
            L1=reg_count(ein)%start_l
            L2=reg_count(ein)%end_l
            call SPTNQ(EW_vector(L1:L2,1),EW_vector(L1:L2,2),L2-L1+1,s1)
            various_reg(ein,1)=s1
            various_reg(ein,2)=exp(-s1)
            various_reg(ein,3)=pral*(1-various_reg(ein,2))
            pral=pral*various_reg(ein,2)
            write(500,'(I5,5X,I7,5X,I7,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,I2)')k,L1,L2,various_reg(ein,1:3),&
               pral,reg_count(ein)%renum
          end do
        else
          E_W=0
        end if
   
        ! compute the contribution of each region
        do ein=1,region_num
          sum_region(ein)=0;
        end do
        do ein=1,reg_change
          ind=reg_count(ein)%renum
          sum_region(ind)=sum_region(ind)+various_reg(ein,3)
        end do
        do ein=1,region_num
          grid_rr(k,ein)=sum_region(ein)
        end do
        write(200,'(F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3)') grid_rr(k,1:region_num)
   
        deallocate(various_reg)
        !calculate recycling ratio
        temprr=1-exp(-E_W)
        sumrr=sumrr+temprr*TPP
        sumpp=sumpp+TPP
        precip_grid(k)=TPP

    end do domain
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    print *,'tracing completed'
    rr(day)=sumrr/sumpp
    print *,rr(day)
   
    do ein=1,region_num
      sumrr=0
      sumpp=0
      do k=1,domsize
        sumrr=sumrr+grid_rr(k,ein)*precip_grid(k)
        sumpp=sumpp+precip_grid(k)
      end do
      daily_rr(day,ein)=sumrr/sumpp
    end do
    write(300,'(I3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3)') day,daily_rr(day,1:region_num)   !FILE 300 is output to save contribution from each subregion to target region
    
  end do daily

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!##################################################################################################
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   close(200)
   close(300)
   close(500)
end program Recycling


!===========================================================================
!subroutine to build files for output
!---------------------------------------------------------------------------
subroutine buildfile
  use param
  implicit none
  integer,parameter :: filen_out=3
  character(len=4)  :: items(filen_out)=(/'GRDR','DAYR','BCTJ'/)
  ! 'BT'-back trajactory position,'EW'-E/W,'RR'-recycling ratio
  !  The last three files are test files
  character(len=70) :: chfile1,chfile2,path,chfile3
  
  write(path,'(A29,I4,A1)') '/data/francina/b/hhu18/2LDRM/',ithyear,'/'
  write(chfile1,'(A34,A11,I4,A1,A4,A7)') path,'WRF_OUTPUT_',ithyear,'_',items(1),'_PW.txt'
  open(unit=200,file=chfile1)

  write(chfile2,'(A34,A11,I4,A1,A4,A7)') path,'WRF_OUTPUT_',ithyear,'_',items(2),'_PW.txt'
  open(unit=300,file=chfile2)

  write(chfile3,'(A34,A11,I4,A1,A4,A7)') path,'WRF_OUTPUT_',ithyear,'_',items(3),'_PW.txt'
  open(unit=500,file=chfile3)

end subroutine
!==========================================================================
! subroutine to read in the i and j indices for the target region in the domain
!--------------------------------------------------------------------------
subroutine monsoon(domain)
  use param
  implicit none
  integer,dimension(domsize,2),intent(out) :: domain
  integer :: ii,error
  character(len=70) :: file_name
  file_name='/data/francina/a/francina/Pileus/Tracer_NAM/WRF_Tracer_NAM_ij.txt'
  open(14,file=file_name,action='read',status='old',access='sequential')
  ii=1
  do while(.true.)
    read(14,'(I2,4X,I3)',iostat=error)domain(ii,1),domain(ii,2)
    if(error/=0) exit
    ii=ii+1
  end do
  write(*,*) 'Monsoon region has been loaded.'
  close(14)
end subroutine
!===========================================================================
! subroutine to read in the input fields E, W, P, U, V
!---------------------------------------------------------------------------
subroutine openfile(xdays,xhours)
  use param
  implicit none
  integer :: file_n=5, i,j,k,l,flag
  real, dimension(bigx,bigy,days,3),intent(out) :: xdays
  real, dimension(bigx,bigy,days*24/velocity_t,4:5),intent(out) :: xhours
  character :: n_c
  character(len=70) :: chfile, path
  character(len=90) :: whpath
  character(len=2) :: items(5)=(/'PW','ET','PP','U3','V3'/)
  write(path,'(A44,I4,A1)') '/data/francina/a/francina/Pileus/Tracer_NAM/',ithyear,'/'

  do l=1, file_n
    write(chfile,'(A3,I4,A7,I3,A1,A2,A4)') 'DRM',ithyear,'0601WRF',days,'_',items(l),'.txt'
    write(whpath,'(A49,A24)') path,chfile
    print *, chfile
    open(10,file=whpath)
    if(l<=3)then
      read(10,*) (((xdays(i,j,k,l),i=1,bigx),j=1,bigy),k=1,days)
    else
      read(10,*) (((xhours(i,j,k,l),i=1,bigx),j=1,bigy),k=1,days*24/velocity_t)
    end if
    write(*,'(A30,A20)') chfile, 'has been opened.'
    close(10)
  end do
end subroutine
!==========================================================================
! subroutine to read in masks for each subregion
!__________________________________________________________________________
subroutine readmap(regions)
  use param
  implicit none
  integer,dimension(bigx,bigy),intent(out) :: regions
  character(len=70) :: file_name
  integer :: i,j
    
  file_name='/data/francina/a/francina/Pileus/Tracer_NAM/WRF_Tracer_5reg.txt'
  open(12,file=file_name)
  read(12,*) ((regions(i,j),i=1,bigx),j=1,bigy)
  write(*,*) 'Map has been loaded.'
  close(12)
end subroutine
!==========================================================================
! subroutine to find U, V, E and W for a certain grid at a certain time
!__________________________________________________________________________
subroutine finduvEW(ppx,ppy,hourcount,daycount,W,E,U,V,outW,outE,outU,outV,xint,yint,flag)
  use param
  implicit none
  real,intent(in) :: ppx,ppy
  integer,intent(in) :: daycount,hourcount
  real,dimension(bigx,bigy),intent(in) ::U,V,E,W
  real,intent(out) :: outW,outE,outU,outV
  integer,intent(out) :: xint,yint,flag
  flag=0
  xint=ceiling(ppx/dx)
  yint=ceiling(ppy/dy)
  outW=W(xint,yint)
  outE=E(xint,yint)
  outU=U(xint,yint)
  outV=V(xint,yint)
  if (outW==NAN_value) then
    flag=1
  elseif (outE==NAN_value) then
    flag=1
  elseif (outU==NAN_value) then
    flag=1
  elseif (outV==NAN_value) then
    flag=1
  end if
end subroutine
!==========================================================================
! subroutine to find U and V for a certain grid at a certain time
! called during back-tracking
!__________________________________________________________________________
subroutine finduv(xp,yp,hourcount,u,v,outu,outv,xint,yint,flag)
  use param
  implicit none
  real,intent(in) ::xp,yp
  integer,intent(in) ::hourcount
  real,dimension(bigx,bigy),intent(in) ::u,v
  real,intent(out)::outu,outv
  integer,intent(out)::xint,yint,flag
  flag=0
  xint=ceiling(xp/dx)
  yint=ceiling(yp/dy)
  outu=u(xint,yint)
  outv=v(xint,yint)
  if (outu==NAN_value) then
    flag=1
  else if (outv==NAN_value) then
    flag=1
  end if
end subroutine
!============================================================================
! subroutine to do the integral of E/W along the back-trajectory
!____________________________________________________________________________
subroutine SPTNQ(x,y,n,s)
  use param
  implicit none
  integer,intent(in):: n
  real,dimension(n),intent(in)::x,y
  real,intent(out) :: s
  integer:: i,j,k
  real :: a,b1,b2
  
  s=0
  do i=1,n-1
    a=x(i+1)-x(i)
    b1=y(i)
    b2=y(i+1)
    s=s+a*(b1+b2)/(2*3600*velocity_t)
  end do
end subroutine
!_____________________________________________________________
