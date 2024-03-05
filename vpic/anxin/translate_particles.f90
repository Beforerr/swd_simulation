!---------------------------------------
! parallel conversion; Using Bill's data type at LANL
! 
! this code convert VPIC output into gda files,
! which are "bricks" of data
! 
!---------------------------------------

module MPI
  include "mpif.h"                                                                                         
  integer myid,numprocs,ierr                                                                               
  integer master  

  ! MPI IO stuff
  integer nfiles, nbands
  parameter(nfiles=36)!66)
  parameter(nbands=0)!0)

  integer sizes(3), subsizes(3), starts(3)
  integer fileinfo, ierror, fh, filetype, status(MPI_STATUS_SIZE), output_format, continuous, file_per_slice
  integer(kind=MPI_OFFSET_KIND) :: disp, offset
  integer err_length,ierror2
  character*(256) err_msg
  character*(40), dimension (nfiles+2*nbands) :: fnames
  CHARACTER*(40) cfname

  parameter(master=0)                                                                                      
  parameter(continuous=1)
  parameter(file_per_slice=2)
  
  integer append_to_files

end module MPI

module topology
implicit none
  integer topology_x,topology_y,topology_z
  real(kind=8) tx,ty,tz
  integer, allocatable :: domain_size(:,:,:,:), idxstart(:,:), idxstop(:,:) 

  type :: ht_type
     
     integer(kind=4) :: tx,ty,tz         ! number of processes in x and y
     integer(kind=4) :: nx,ny,nz         ! number of cells in each direction that belong to this process
     integer(kind=4) :: start_x, stop_x, start_z, stop_z, start_y,stop_y ! where to start/stop in x/y/z
     integer(kind=4) :: ix,iy,iz 

  end type ht_type

  type(ht_type) :: ht

end module topology


module data_struct
  implicit none

  ! define structure for V0 header
  type :: v0header
    integer(kind=4) :: step, nx, ny, nz
    real(kind=4) :: dt, dx, dy, dz, x0, y0, z0
    real(kind=4) :: cvac, eps0, damp
    integer(kind=4) :: mp_rank, mp_num_proc
    integer(kind=4) :: species_id
    real(kind=4) :: q_m
  end type v0header

  ! define structure for array header
  type :: array_header
    integer(kind=4) :: size, ndim, dim
  end type array_header

  ! define structure for input particle data
  type :: ipdata
    real(kind=4) :: dx, dy, dz
    integer(kind=4) :: i
    real(kind=4) :: ux, uy, uz
    real(kind=4) :: q
  end type ipdata

  ! define structure for output particle data
  type :: opdata
    real(kind=4) :: x, y, z
    real(kind=4) :: ux, uy, uz
    real(kind=4) :: q
  end type opdata
end module data_struct

program pppp
  use data_struct
  use topology
  use MPI
  implicit none

  ! declare v0 header
  type(v0header) :: v0
  type(array_header) :: ah
  type(ipdata) :: pi
  type(opdata) :: po

  integer(kind=4) :: ndomains, n, tindex, ip
  integer(kind=4) :: tindex_new, tindex_start, tindex_stop
  integer(kind=4) :: output_record, nskip, nout
  character(len=60) :: fname, ofname2, ofname3, ofname4, ofname5
  character(len=1) :: species = "H"
  logical :: check, dfile
  integer(kind=4) :: version, dump_type
  integer(kind=4) :: nxp2, nyp2, nzp2
  integer(kind=4) :: ix, iy, iz, i
  real(kind=4) :: x, y, z

  real(kind=8) :: nx_d, ny_d, nz_d, dt
  real(kind=8) :: xmax, ymax, zmax
  integer(kind=4) :: httx, htty, httz
  integer(kind=4) :: nx, ny, nz

  ! phase space data
  real(kind=4), dimension(:,:), allocatable :: fylim, fxlim, fypara, fyperp

  ! magnetic field data
  integer(kind=4) :: noutB
  integer(kind=8) :: pos
  real(kind=4), dimension(:,:,:), allocatable :: Bx, By, Bz

  ! parameters for particle diagnostics
  real(kind=4) :: xmax_diag, xmin_diag, yup_diag, ylow_diag, umax_diag, umin_diag, dx, du
  integer(kind=4) :: num_x, num_u
  real(kind=4) :: ymax_diag, ymin_diag, xup_diag, xlow_diag, dy
  integer(kind=4) :: num_y

  ! variables for timing
  ! integer :: c1,c2,cr,cm
  real(kind=8) :: mp_elapsed

  ! init MPI 
  call MPI_INIT(ierr)                                                                                      
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)                                                           
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)                                                       

  namelist /datum/ httx,htty,httz,tindex_start,tindex_stop,&
    output_format,append_to_files,&
    xmax_diag,xmin_diag,yup_diag,ylow_diag,umax_diag,umin_diag,num_x,num_u,&
    ymax_diag,ymin_diag,xup_diag,xlow_diag,num_y, noutB

  ! read the configuration file
  if (myid == master) then
    open(unit=10, file='conf_particles.dat', form='formatted', status='old')
    read(10, datum)
    close(10)
  endif

  call MPI_BCAST(httx,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(htty,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(httz,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(tindex_start,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tindex_stop,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(output_format,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(append_to_files,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(xmax_diag,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xmin_diag,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(yup_diag, 1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ylow_diag,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(umax_diag,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(umin_diag,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_x,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_u,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(ymax_diag,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ymin_diag,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xup_diag, 1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xlow_diag,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_y,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(noutB,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  
  
  ht%tx = httx
  ht%ty = htty
  ht%tz = httz

 ! check the topology for consistency

  if ( (ht%tx*ht%ty*ht%tz /= numprocs).or.(topology_x/ht%tx*ht%tx /= topology_x).or.&
       (topology_z/ht%tz*ht%tz /= topology_z).or.(topology_y/ht%ty*ht%ty /= topology_y) ) then

     if (myid == master) print *, "invalid converter topology"
     call MPI_FINALIZE(ierr)
     stop

  endif
  

  ! read info.bin
  ! read the info file
  if (myid == master) then
    open(unit=10, file='info.bin', status='old',&
      form='unformatted', access='stream')
    read(10) tx
    read(10) ty
    read(10) tz

    read(10) xmax
    read(10) ymax
    read(10) zmax
    
    read(10) nx_d
    read(10) ny_d
    read(10) nz_d

    read(10) dt

    close(10)
  endif

  call MPI_BCAST(tx, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ty, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tz, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(xmax, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ymax, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(zmax, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nx_d, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ny_d, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nz_d, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(dt, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)

  ! convert to integers
  topology_x = floor(tx + 0.5)
  topology_y = floor(ty + 0.5)
  topology_z = floor(tz + 0.5)

  nx = floor(nx_d + 0.5)
  ny = floor(ny_d + 0.5)
  nz = floor(nz_d + 0.5)

  ! numner of cells for each process

  ht%nx = nx/ht%tx
  ht%ny = ny/ht%ty
  ht%nz = nz/ht%tz

  ! number of domains = number of processors in simulations
  ndomains = topology_x*topology_y*topology_z

  ! specify parameters for diagnostics
  if (num_x == 0) then
    num_x = nx_d
  endif
  dx = (xmax_diag - xmin_diag) / num_x
  du = (umax_diag - umin_diag) / num_u
  dy = (ymax_diag - ymin_diag) / num_y

  if (myid == master) then
    write(*,*) "-----------------------------------------------"
    write(*,*) " Topology: ", topology_x, topology_y, topology_z
    write(*,*) " nx,nz,nz: ", nx,ny,nz
    write(*,*) " ht: nx,ny,nz: ", ht%nx,ht%ny,ht%nz
    write(*,*) "xmin_diag, xmax_diag, num_x: ", xmin_diag, xmax_diag, num_x
    write(*,*) "umin_diag, umax_diag, num_u: ", umin_diag, umax_diag, num_u
    write(*,*) "-----------------------------------------------"
  endif

  ! Determine number of iterations between output files
  if (myid == master) then
  
    dfile=.false.
    tindex= tindex_start
    do while(.not.dfile)
       tindex=tindex+1
       write(fname,"(A,I0,A,I0,A)")"particle/T.",tindex,&
         "/"//species//"particle.",tindex,".0"
       if (tindex .ne. 1) inquire(file=trim(fname),exist=dfile)
    enddo
    nskip = 1
    nout = (tindex-tindex_start)*nskip
    !nout = 1 !1time part L.O.  
  
  ! Total size of domain
  
    print *,"---------------------------------------------------"
    print *
    print *,"xmax=",xmax,"   ymax=",ymax,"   zmax=",zmax
    print *
    print *,"Iterations between output=",nout
    print *,"---------------------------------------------------"
  
  endif

  call MPI_BCAST(nout,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  
  ! Need to determine the last record written,so we know which time slice to process next

  if (append_to_files==1) then
   output_record = (tindex_start/nout) + 1
  else
    output_record = 1
  endif

  tindex = tindex_start

  ! allocate memory for phase space data
  allocate(fylim(num_x,num_u))
  allocate(fxlim(num_y,num_u))
  allocate(fypara(num_y,num_u))
  allocate(fyperp(num_y,num_u))

  ! allocate memory for magnetic field data
  allocate(Bx(nx,ny,nz))
  allocate(By(nx,ny,nz))
  allocate(Bz(nx,ny,nz))

  ! initialize timer
  ! call system_clock(count_rate=cr)
  ! call system_clock(count_max=cm)

  !loop over time slices
  dfile = .true.
  do while(dfile)

    ! initialize phase space data
    fylim = 0.0
    fxlim = 0.0
    fypara = 0.0
    fyperp = 0.0
    if (myid == master) then
      print *, " processing particles; time slice: ", tindex
    endif

    ! read magnetic field data
    if (myid == master) then
      ! Bx
      fname = 'data/bx.gda'
      open(unit=1, file=trim(fname), access='stream', status='unknown', form='unformatted', action='read')
      pos = tindex/noutB*sizeof(Bx) + 1
      read(1, pos=pos) Bx
      close(1)
      ! By
      fname = 'data/by.gda'
      open(unit=1, file=trim(fname), access='stream', status='unknown', form='unformatted', action='read')
      read(1, pos=pos) By
      close(1)
      ! Bz
      fname = 'data/bz.gda'
      open(unit=1, file=trim(fname), access='stream', status='unknown', form='unformatted', action='read')
      read(1, pos=pos) Bz
      close(1)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Bx,size(Bx),MPI_REAL,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(By,size(By),MPI_REAL,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Bz,size(Bz),MPI_REAL,master,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! time for start of loop
    ! call system_clock(c1)
    mp_elapsed = MPI_WTIME()
    
    ! loop over domains/processors
    do n = myid+1, ndomains, numprocs
      write(fname,"(A,I0,A,I0,A,I0)") "particle/T.",tindex,"/"//species//"particle.",&
        tindex,".",n-1
      inquire(file=trim(fname), exist=check)

      if (check) then
        open(unit=10, file=trim(fname), status='unknown', form='binary', access='stream')
      else
        write(*,*) "Can't find file:", fname
        write(*,*) 
        write(*,*) " ***  Terminating ***"
        stop
      endif

      ! Binary compatibility information
      call read_boilerplate(10)
      ! Dump type and header format version
      read(10) version
      read(10) dump_type
      ! High level information
      read(10) v0
      ! Particle array information
      read(10) ah

      ! Now loop over particles
      nxp2 = v0%nx + 2
      nyp2 = v0%ny + 2
      nzp2 = v0%nz + 2
      do ip = 1, ah%dim
        read(10) pi

        i = pi%i
        iz = i / (nxp2*nyp2)
        iy = (i - iz*nxp2*nyp2) / nxp2
        ix = i - nxp2*(iy+nyp2*iz)

        ! Compute real particle position from relative coords and grid data
        ! Store data first in tmp particle struct po so we can further process
        po%x = v0%x0 + ((ix-1) + (pi%dx+1)*0.5)*v0%dx
        po%y = v0%y0 + ((iy-1) + (pi%dy+1)*0.5)*v0%dy
        po%z = v0%z0 + ((iz-1) + (pi%dz+1)*0.5)*v0%dz
        po%ux = pi%ux
        po%uy = pi%uy
        po%uz = pi%uz
        po%q  = pi%q

        ! deposit particle distribution
        call deposit_xux(po, fylim, xmin_diag, dx, num_x, umin_diag, du, num_u, yup_diag, ylow_diag)
        call deposit_yux(po, fxlim, ymin_diag, dy, num_y, umin_diag, du, num_u, xup_diag, xlow_diag)
        call deposit_yucl_2d(po, fypara, fyperp, ymin_diag, dy, num_y, umin_diag, du, num_u,&
          xup_diag, xlow_diag, Bx, By, Bz, -xmax/2.0, -ymax/2.0, v0%dx, v0%dy, nx, ny)
      enddo ! ip; particle loop
    enddo ! n; ndomains loop

    if (myid == master) then
      call MPI_REDUCE(MPI_IN_PLACE, fylim, size(fylim), MPI_REAL, MPI_SUM, master, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(MPI_IN_PLACE, fxlim, size(fxlim), MPI_REAL, MPI_SUM, master, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(MPI_IN_PLACE, fypara, size(fypara), MPI_REAL, MPI_SUM, master, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(MPI_IN_PLACE, fyperp, size(fyperp), MPI_REAL, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    else
      call MPI_REDUCE(fylim, fylim, size(fylim), MPI_REAL, MPI_SUM, master, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(fxlim, fxlim, size(fxlim), MPI_REAL, MPI_SUM, master, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(fypara, fypara, size(fypara), MPI_REAL, MPI_SUM, master, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(fyperp, fyperp, size(fyperp), MPI_REAL, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    endif

    if (myid == master) then
      ! open output file
      write(ofname2,"(A,I0,A)") "data/xux.",tindex,'.bin'
      write(ofname3,"(A,I0,A)") "data/yux.",tindex,'.bin'
      write(ofname4,"(A,I0,A)") "data/yupara.",tindex,'.bin'
      write(ofname5,"(A,I0,A)") "data/yuperp.",tindex,'.bin'
      open(unit=2, file=trim(ofname2), action='write', form='unformatted',&
        access='stream', status='replace')
      open(unit=3, file=trim(ofname3), action='write', form='unformatted',&
        access='stream', status='replace')
      open(unit=4, file=trim(ofname4), action='write', form='unformatted',&
        access='stream', status='replace')
      open(unit=5, file=trim(ofname5), action='write', form='unformatted',&
        access='stream', status='replace')
      ! write phase space data
      write(2) fylim
      write(3) fxlim
      write(4) fypara
      write(5) fyperp
      ! close output file
      close(2)
      close(3) 
      close(4) 
      close(5) 
    endif

    ! time for one time slice
    ! call system_clock(c2)
    ! write(*,*) "system_clock: ", (c2-c1)/real(cr)
    mp_elapsed = MPI_WTIME() - mp_elapsed
    if (myid==master) write(*,'(A, F5.1)') " => time(s):",mp_elapsed
    
    ! check if there is another time slice to read
    dfile = .false.
    tindex_new = tindex + nout
    if ((tindex_new .gt. 1) .and. (tindex_new<=tindex_stop)) then
      write(fname,"(A,I0,A,I0,A,I0)") "particle/T.",tindex_new,"/"//species//&
        "particle.",tindex_new,".0"
      inquire(file=trim(fname), exist=dfile)
    endif

    tindex = tindex_new
    if (dfile) output_record = output_record + 1

  enddo ! time loop

  contains

  subroutine read_boilerplate(iunit)
    implicit none
    integer(kind=1)sizearr(5)
    integer(kind=2)cafevar 
    integer(kind=4)deadbeefvar
    real(kind=4)realone
    real(kind=8)doubleone
    integer iunit
    read(iunit)sizearr
    read(iunit)cafevar
    read(iunit)deadbeefvar
    read(iunit)realone
    read(iunit)doubleone
  end subroutine read_boilerplate

end program pppp


subroutine deposit_xux(po, f, xmin, dx, nx, vmin, dv, nv, yup, ylow)
  use data_struct
  implicit none

  integer, intent(in) :: nx, nv
  type(opdata), intent(in) :: po
  real, dimension(nx,nv), intent(inout) :: f
  real, intent(in) :: xmin, dx, vmin, dv, yup, ylow

  ! local data
  integer :: ix, iv
  real :: rx, rv

  rx = (po%x  - xmin) / dx + 1
  ix = rx
  rv = (po%ux - vmin) / dv + 1
  iv = rv

  if ((ix>=1) .and. (ix<=nx) .and. (iv>=1) .and. (iv<=nv)) then
    rx = rx - ix
    rv = rv - iv
    if ((po%y >= ylow) .and. (po%y < yup)) then
      f(ix,  iv  ) = f(ix,  iv  ) + (1.0-rx)*(1.0-rv)*po%q
      f(ix+1,iv  ) = f(ix+1,iv  ) +      rx *(1.0-rv)*po%q
      f(ix,  iv+1) = f(ix,  iv+1) + (1.0-rx)*     rv *po%q
      f(ix+1,iv+1) = f(ix+1,iv+1) +      rx *     rv *po%q
    endif
  endif

  return
end subroutine deposit_xux


subroutine deposit_yux(po, f, ymin, dy, numy, vmin, dv, numv, xup, xlow)
  use data_struct
  implicit none

  integer, intent(in) :: numy, numv
  type(opdata), intent(in) :: po
  real, dimension(numy,numv), intent(inout) :: f
  real, intent(in) :: ymin, dy, vmin, dv, xup, xlow

  ! local data
  integer :: iy, iv
  real :: ry, rv

  ry = (po%y  - ymin) / dy + 1
  iy = ry
  rv = (po%ux - vmin) / dv + 1
  iv = rv

  if ((iy>=1) .and. (iy<=numy) .and. (iv>=1) .and. (iv<=numv)) then
    ry = ry - iy
    rv = rv - iv
    if ((po%x >= xlow) .and. (po%x < xup)) then
      f(iy,  iv  ) = f(iy,  iv  ) + (1.0-ry)*(1.0-rv)*po%q
      f(iy+1,iv  ) = f(iy+1,iv  ) +      ry *(1.0-rv)*po%q
      f(iy,  iv+1) = f(iy,  iv+1) + (1.0-ry)*     rv *po%q
      f(iy+1,iv+1) = f(iy+1,iv+1) +      ry *     rv *po%q
    endif
  endif

  return
end subroutine deposit_yux


subroutine deposit_yucl_2d(po, fpara, fperp, ymind, dyd, numy, vmind, dvd, numv, xup, xlow,&
    Bx, By, Bz, xmin, ymin, dx, dy, nx, ny)
  ! 2d means the configuration space is 2d.
  use data_struct
  use MPI
  implicit none

  integer, intent(in) :: numy, numv, nx, ny
  type(opdata), intent(in) :: po
  real, dimension(numy,numv), intent(inout) :: fpara, fperp
  real, intent(in) :: ymind, dyd, vmind, dvd, xup, xlow
  real, intent(in) :: dx, dy
  real(kind=8), intent(in) :: xmin, ymin
  real, dimension(nx,ny), intent(in) :: Bx, By, Bz

  ! local data
  integer :: ix, iy, ivpara, ivperp, ix1, iy1
  real :: rx, ry, rvpara, rvperp, vperp, vpara, bxp, byp, bzp, btot

  if ((po%x >= xlow) .and. (po%x < xup)) then
    ! interpolate magnetic field to particle's position
    rx = (po%x - xmin) / dx + 1
    ix = rx
    rx = rx - ix
    ix1 = ix + 1
    ry = (po%y - ymin) / dy + 1
    iy = ry
    ry = ry - iy
    iy1 = iy + 1
    bxp = Bx(ix,iy)*(1.0-rx)*(1.0-ry) + Bx(ix1,iy)*rx*(1.0-ry) + Bx(ix,iy1)*(1.0-rx)*ry + Bx(ix1,iy1)*rx*ry
    byp = By(ix,iy)*(1.0-rx)*(1.0-ry) + By(ix1,iy)*rx*(1.0-ry) + By(ix,iy1)*(1.0-rx)*ry + By(ix1,iy1)*rx*ry
    bzp = Bz(ix,iy)*(1.0-rx)*(1.0-ry) + Bz(ix1,iy)*rx*(1.0-ry) + Bz(ix,iy1)*(1.0-rx)*ry + Bz(ix1,iy1)*rx*ry

    ry = (po%y  - ymind) / dyd + 1
    iy = ry
    ry = ry - iy
    
    btot = sqrt(bxp*bxp + byp*byp + bzp*bzp)
    if (btot /= 0) then
      vpara = (po%ux*bxp + po%uy*byp + po%uz*bzp) / btot
    else
      vpara = 0.0
    endif
    vperp = sqrt(po%ux*po%ux + po%uy*po%uy + po%uz*po%uz - vpara*vpara)
    rvpara = (vpara - vmind) / dvd + 1
    ivpara = rvpara
    rvpara = rvpara - ivpara
    rvperp = (vperp - vmind) / dvd + 1
    ivperp = rvperp
    rvperp = rvperp - ivperp

    if ((iy>=1) .and. (iy<numy)) then

      if ((ivpara>=1) .and. (ivpara<numv)) then
        fpara(iy,  ivpara  ) = fpara(iy,  ivpara  ) + (1.0-ry)*(1.0-rvpara)*po%q
        fpara(iy+1,ivpara  ) = fpara(iy+1,ivpara  ) +      ry *(1.0-rvpara)*po%q
        fpara(iy,  ivpara+1) = fpara(iy,  ivpara+1) + (1.0-ry)*     rvpara *po%q
        fpara(iy+1,ivpara+1) = fpara(iy+1,ivpara+1) +      ry *     rvpara *po%q
      endif

      if ((ivperp>=1) .and. (ivperp<numv)) then
        fperp(iy,  ivperp  ) = fperp(iy,  ivperp  ) + (1.0-ry)*(1.0-rvperp)*po%q
        fperp(iy+1,ivperp  ) = fperp(iy+1,ivperp  ) +      ry *(1.0-rvperp)*po%q
        fperp(iy,  ivperp+1) = fperp(iy,  ivperp+1) + (1.0-ry)*     rvperp *po%q
        fperp(iy+1,ivperp+1) = fperp(iy+1,ivperp+1) +      ry *     rvperp *po%q
      endif
    
    endif

  endif  ! xlow <= po%x < xup

  return
end subroutine deposit_yucl_2d
