!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module printings.f90: ...description...
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Yvelin Giret, Thomas Pope, Julien Eng
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module printings
    use constants
    use param
    implicit none
    ! equilibration parameters
    double precision, allocatable, dimension(:)   :: time_save
    double precision                              :: diff_time
    double precision, allocatable, dimension(:,:) :: flux_save
    double precision                              :: diff_flux
    integer                               :: i_equilibration
    ! timing parameters
    double precision                              :: t1,t2
    integer                               :: iteration_threshold = 10
    ! converter parameters 
    double precision                              :: flux_to_cur
    double precision                              :: flux_to_mob
    contains
!=====================================================================!
!                           !
!=====================================================================!
    subroutine print_logo()
    implicit none
    write(screen,fmdlm2) delim0, blank1
    write(screen,fmchr1) "                ____ _               _____ _____ ____        _  ____  __  ____       "
    write(screen,fmchr1) "               / ___| |__   __ _ _ _|_   _| ____|  _ \      | |/ /  \/  |/ ___|      "
    write(screen,fmchr1) "              | |   | '_ \ / _` | '__|| | |  _| | | | |_____| ' /| |\/| | |          "
    write(screen,fmchr1) "              | |___| | | | (_| | |   | | | |___| |_| |_____| . \| |  | | |___       "
    write(screen,fmchr1) "               \____|_| |_|\__,_|_|   |_| |_____|____/      |_|\_\_|  |_|\____|      "
    write(screen,fmdlm1) blank1               
    write(screen,fmchr1) "                Charge Transfer and Exciton Dynamics with Kinetic Monte-Carlo       "
    write(screen,fmdlm1) blank1               
    write(screen,fmchr1) "                      *****    ******.                                                      "
    write(screen,fmchr1) "              , *  *..**.*************        **,                    ,                      "
    write(screen,fmchr1) "           **    ,,.        .********,                   ,       .,,,,,        .            "
    write(screen,fmchr1) "./////*,,**** *....* * .*     ******.          ,,,,,,      ,.,.,,,,,,,,,,,,,,,,,,,,,.,,,,.  "
    write(screen,fmchr1) ",/////*************,    .*    **             ***.**,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ."
    write(screen,fmchr1) "        /***********   ****              .,  ** ,*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,    ,,       "
    write(screen,fmchr1) "          ****************, *              /,/***,,,,,,.,,,,,,,,,,,.,,.,,,//,,,.            "
    write(screen,fmchr1) "           /////////..//*                /// // ** .  *...,,,***/////.//////.               "
    write(screen,fmchr1) "            ///////////                  **,,*      *,,////,,.*,//////////  .,              "
    write(screen,fmchr1) "              .,,,                     *,.,,,,***,,,.,,, ...,*****.*/////,                  "
    write(screen,fmchr1) "                   ,.                  *,,**,,,,,,,,* **      **   .**.                     "
    write(screen,fmchr1) "                       ,,,,*            ***,,,****,,,,,,            .                       "
    write(screen,fmchr1) "                      ,,,******,*             *,,,.,,                , .,,  ,.,.            "
    write(screen,fmchr1) "                       ,,,.******             .,,,,.,  .                    .* *            "
    write(screen,fmchr1) "                         .,*****.              ,**.,  ,                 *********    .      "
    write(screen,fmchr1) "                        ,,,,,*                 .,,,                     **********          "
    write(screen,fmchr1) "                        ,,,                                                     ,       ,   "
    write(screen,fmchr1) "                        .,                                                                  "
    write(screen,fmchr1) "                        .,       ,                                                          "
    write(screen,fmdlm3) blank1, delim1, blank1
    write(screen,fmchr1) "version 1.4"
    write(screen,fmchr2) "Authors:", "Thomas Pope and Yvelin Giret"
    write(screen,fmdlm1) blank1
    endsubroutine print_logo
!=====================================================================!
!                           !
!=====================================================================!   
    subroutine initialize_print_timer()
    implicit none
    call cpu_time(t1)
    endsubroutine initialize_print_timer
!=====================================================================!
!                           !
!=====================================================================!   
    subroutine initialize_print_converters(c1)
    use type_charge
    implicit none
    type(charge)                          :: c1
    flux_to_cur = e_to_C / (dble(c1%Nx * c1%Ny * c1%Nz) * (lattice * a0_to_nm * 1e-9)**2)
    flux_to_mob = mb_to_SI * lattice / field
    endsubroutine initialize_print_converters
!=====================================================================!
!                                                                     !
!=====================================================================!
    subroutine open_file(u,substr,ip,exts)
    implicit none
    character(*), intent(in)           :: substr
    character(*), intent(in), optional :: exts
    integer,      intent(in)           :: u, ip
    character(100)                     :: fname
    if(PRESENT(exts)) then
      write(fname,'("OUTPUT_",a,"_",i0,a)') substr, ip, exts
    else
      write(fname,'("OUTPUT_",a,"_",i0,".dat")') substr, ip
    endif
    open(u,file=fname,access='append')
    endsubroutine open_file

    subroutine allocate_print_arrays(Nchrg)
    implicit none
    integer :: Nchrg
    allocate(time_save(N_equilibration))
    allocate(flux_save(Nchrg,N_equilibration))
    time_save         = 0.0d0
    flux_save         = 0.0d0
    i_equilibration   = 1
    endsubroutine allocate_print_arrays
!=====================================================================!
!                           !
!=====================================================================!
    subroutine print_step_output(time_AU,iteration,c1,c2,e1,e2)
    use my_mpi
    use type_charge
    implicit none
    type(charge),     intent(in)           :: c1
    type(charge),     intent(in), optional :: c2
    type(charge),     intent(in), optional :: e1
    type(charge),     intent(in), optional :: e2
    integer,          intent(in)           :: iteration
    double precision, intent(in)           :: time_AU
    double precision                       :: cpu_timer
    double precision                       :: time, flux_over_time
    double precision                       :: mob1, mob2
    double precision                       :: cur1, cur2
    character(1)                           :: creturn=achar(13)
    time         = time_AU * hEh_to_s
    if(iteration * 100 / Nsteps .ge. iteration_threshold .and. Nsteps .gt. 0 .and. i_am_root) then
      call cpu_time(t2)
      cpu_timer           = t2 - t1
      if(print_verb.eq.1) then
        if(cpu_timer.lt.60) then
          write(screen,101,advance='NO') creturn, iteration_threshold, cpu_timer
        else
          t1                  = t2
          write(screen,111,advance='NO') creturn, iteration_threshold, int(cpu_timer/60), mod(cpu_timer,60.0d0)
        endif
        if(iteration_threshold.eq.100) write(screen,*) 
      else
        t1                  = t2
        if(cpu_timer.lt.60) then
          write(screen,102) iteration_threshold, cpu_timer
        else
          t1                  = t2
          write(screen,112) iteration_threshold, int(cpu_timer/60), mod(cpu_timer,60.0d0)
        endif
      endif
      iteration_threshold = iteration_threshold + 10
      flush(screen)
    endif 
    if (PRESENT(c2)) then 
      flux_over_time = dble(c1%flux) / time_AU
      cur1           = flux_over_time * flux_to_cur
      mob1           = flux_over_time * flux_to_mob / dble(c1%N)
      flux_over_time = - dble(c2%flux) / time_AU
      cur2           = flux_over_time * flux_to_cur
      mob2           = flux_over_time * flux_to_mob / dble(c2%N)   
      write(111,*) time, c1%flux, c2%flux, c1%relax_count, c2%relax_count
      write(114,*) time, mob1, mob2
      write(117,*) time, e1%N, e2%N
      time_save(i_equilibration)   = time_Au
      flux_save(1,i_equilibration) = c1%flux
      flux_save(2,i_equilibration) = c2%flux
      if(i_equilibration.lt.N_equilibration) then
        i_equilibration = i_equilibration + 1
      else
        diff_time                       = time_save(N_equilibration) - time_save(1)
        time_save(:N_equilibration-1)   = time_save(2:)
        diff_flux                       = flux_save(1,N_equilibration) - flux_save(1,1)
        flux_save(1,:N_equilibration-1) = flux_save(1,2:)
        flux_over_time                  = dble(diff_flux) / diff_time
        cur1                            = flux_over_time * flux_to_cur
        mob1                            = flux_over_time * flux_to_mob / dble(c1%N)
        diff_flux                       = flux_save(2,N_equilibration) - flux_save(2,1)
        flux_save(2,:N_equilibration-1) = flux_save(2,2:)
        flux_over_time                  = - dble(diff_flux) / diff_time
        cur2                            = flux_over_time * flux_to_cur
        mob2                            = flux_over_time * flux_to_mob / dble(c2%N)
      endif 
      write(112,*) time, cur1, cur2
      write(115,*) time, mob1, mob2
      flush(117)
    else
      flux_over_time = dble(c1%flux) / time_AU
      cur1           = flux_over_time * flux_to_cur
      mob1           = flux_over_time * flux_to_mob / dble(c1%N)
      write(111,*) time, c1%flux, c1%relax_count
      write(114,*) time, mob1
      time_save(i_equilibration)   = time_Au
      flux_save(1,i_equilibration) = c1%flux
      if(i_equilibration.lt.N_equilibration) then
        i_equilibration = i_equilibration + 1
      else
        diff_time                       = time_save(N_equilibration) - time_save(1)
        time_save(:N_equilibration-1)   = time_save(2:)
        diff_flux                       = flux_save(1,N_equilibration) - flux_save(1,1)
        flux_save(1,:N_equilibration-1) = flux_save(1,2:)
        flux_over_time                  = dble(diff_flux) / diff_time
        cur1                            = flux_over_time * flux_to_cur
        mob1                            = flux_over_time * flux_to_mob / dble(c1%N)
      endif 
      write(112,*) time, cur1
      write(115,*) time, mob1
    endif
    flush(111); flush(112); flush(114); flush(115)
!    call print_restart("CHRG.RESTART",c1,time_AU)
101 format(a,t3,"|",t6,i3,"% complete. Time elapsed =",f10.5," seconds",t101,"|")
111 format(a,t3,"|",t6,i3,"% complete. Time elapsed =",i0," minutes",f10.5," seconds",t101,"|")
102 format(t2,"|",t5,i3,"% complete. Time elapsed since last print=",f10.5," seconds",t100,"|")
112 format(t2,"|",t5,i3,"% complete. Time elapsed since last print =",i0," minutes",f10.5," seconds",t100,"|")
    endsubroutine print_step_output
!=====================================================================!
!                           !
!=====================================================================!
    subroutine print_mobility_analysis()
    use param
    implicit none
    integer                                          :: npoint
    integer                                          :: i,j,k,istep
    double precision                                 :: time
    character(100)                                   :: fname
    double precision, dimension(Nsteps/modulo_print) :: mob_inp
    double precision, dimension(9,0:MPI_N_calc-1)    :: mob_ave, mob_var, flux
    double precision, dimension(9)                   :: bmob_ave, bmob_std
    double precision                                 :: rwsd, var_diff, flux_ave, flux_std, intersec
    double precision, dimension(9)                   :: RWO
    npoint = Nsteps / modulo_print
    istep=int(npoint/10)
    write(screen,fmdlm3) blank1, delim1, blank1
    do j=0,MPI_N_calc-1
     write(fname,'("OUTPUT_mobility_b_",i0,".dat")') j
     open(115,file=fname)
     do i=1,npoint
      read(115,*) time, mob_inp(i)
     enddo
     close(115)
     do i=1,9
      mob_ave(i,j)=sum(mob_inp(npoint - (istep * i):)) / (istep * i)
      mob_var(i,j)=sum(mob_inp(npoint - (istep * i):) - mob_ave(i,j))**2 / (istep * i - 1)
     enddo
    enddo
    do i=1,9
     bmob_ave(i)=sum(mob_ave(i,:)) / MPI_N_calc
     bmob_std(i)=sqrt( (sum((mob_ave(i,:) - bmob_ave(i))**2) + sum(mob_var(i,:))) / (MPI_N_calc-1))
    enddo   
    do j=0,MPI_N_calc-1
     write(fname,'("OUTPUT_flux_",i0,".dat")') j
     open(115,file=fname)
     do k=1,9; do i=1,istep
       read(115,*) time, flux(10-k,j)
     enddo; enddo
    enddo
    do k=1,9
     flux_ave = sum(flux(k,:)) / MPI_N_calc
     flux_std = sqrt( sum((flux(k,:) - flux_ave)**2) / (MPI_N_calc-1))
     rwsd     = 0.5d0 * sqrt(dble(nsteps) / dble(k))
     var_diff = rwsd**2-flux_std**2
     intersec = (flux_ave * rwsd**2 - flux_std * rwsd * sqrt(flux_ave**2 + 2 * var_diff * log(rwsd / flux_std))) / var_diff
     RWO(k)   = 50.d0 * abs(erf(intersec / (sqrt(2.d0) * rwsd)) + erf((intersec - flux_ave) / (sqrt(2.d0) * flux_std)))
    enddo
    do i=1,9
     write(screen,100) i*10,bmob_ave(i),bmob_std(i),RWO(i)
    enddo  
    write(screen,fmdlm1) blank1  
100 format(t2,"|",t5,"Mobility of the Final ",i0,"% of Data:",t40,es13.6," +/- ",es12.6,t76,"RWO:",2x,f6.2,"%",t100,"|")
    endsubroutine print_mobility_analysis
!=====================================================================!
!       collect and print the singlet or triplet types                !
!=====================================================================!
    subroutine print_exciton_analysis(sings,trips)
    use my_mpi
    use type_charge
    implicit none
    type(charge), intent(in)                               :: sings, trips
    integer,          dimension(0:sings%NE)                :: shist
    integer,          dimension(0:sings%NE,0:MPI_N_calc-1) :: shmpi
    double precision, dimension(0:sings%NE)                :: shout, shstd
    integer,          dimension(0:trips%NE)                :: thist
    integer,          dimension(0:trips%NE,0:MPI_N_calc-1) :: thmpi
    double precision, dimension(0:trips%NE)                :: thout, thstd
    double precision                                       :: hsum
    integer                                                :: i
    shist = 0; shmpi = 0
    thist = 0; thmpi = 0
    do i=1,sings%N; shist(sings%iE(i)) = shist(sings%iE(i)) + 1; enddo
    do i=1,trips%N; thist(trips%iE(i)) = thist(trips%iE(i)) + 1; enddo
    if(i_am_root) then
     shmpi(:,0) = shist
     thmpi(:,0) = thist
     do i=1,we_np-1
      call mpi_recv(shmpi(:,i),sings%NE,mpi_integer,i,0,we_comm,istat,ierr)
     enddo
     do i=1,we_np-1
      call mpi_recv(thmpi(:,i),trips%NE,mpi_integer,i,0,we_comm,istat,ierr)
     enddo
     do i=0,sings%NE
      shout(i) = dble(sum(shmpi(i,:))) / dble(we_np)
      shstd(i) = sqrt(dble(sum((shmpi(i,:)-shout(i))**2)) / dble(we_np-1))
     enddo
     do i=0,trips%NE
      thout(i) = dble(sum(thmpi(i,:))) / dble(we_np)
      thstd(i) = sqrt(dble(sum((thmpi(i,:)-thout(i))**2)) / dble(we_np-1))
     enddo
     write(screen,fmdlm2) delim1, blank1
     write(screen,fmchr1) "Singlet Excited State Distribution:"
     write(screen,fmdlm1) blank1
     hsum  = sum(shout)
     write(screen,fmchfl) "Average number of singlets:", hsum
     write(screen,fmdlm1) blank1
     do i=0,sings%NE
      write(screen,100) "S",i+1,sings%iH(i),sings%iL(i),sings%E0(i)*Eh_to_eV,shout(i),shstd(i)
     enddo
     write(screen,fmdlm3) blank1, delim2, blank1
     write(screen,fmchr1) "Triplet Excited State Distribution:"
     write(screen,fmdlm1) blank1
     hsum  = sum(thout)
     write(screen,fmchfl) "Average number of triplets:", hsum
     write(screen,fmdlm1) blank1
     do i=0,trips%NE
      write(screen,100) "T",i+1,trips%iH(i),trips%iL(i),trips%E0(i)*Eh_to_eV,thout(i),thstd(i)
     enddo
    else
     call mpi_ssend(shist,sings%NE,mpi_integer,0,0,we_comm,ierr)
     call mpi_ssend(thist,trips%NE,mpi_integer,0,0,we_comm,ierr)
    endif
100 format(t2,"|",t5,a1,"_",i0,2x,"H",i0,"->L",i0,t20,"Energy = ",f7.2,"eV", &
        & t45,"Population = ",f9.2" +/- ",f9.2,t100,"|")
    endsubroutine print_exciton_analysis
!=====================================================================!
!       collect and print the singlet or triplet types                !
!=====================================================================!
        subroutine print_exciton_checks(extn,namefile)
        use my_mpi
        use type_charge
        implicit none
        type(charge), intent(in)      :: extn
        character(*), intent(in)      :: namefile
        integer, dimension(0:extn%NE) :: hist, hmpi
        double precision, dimension(0:extn%NE) :: hout
        integer                       :: i
        hist = 0
        hmpi = 0
        do i=1,extn%N
         hist(extn%iE(i)) = hist(extn%iE(i)) + 1
        enddo
        if(i_am_root) then
         do i=1,we_np-1
          call mpi_recv(hmpi,extn%NE,mpi_integer,i,0,we_comm,istat,ierr)
          hist = hmpi + hist
         enddo
         hout = dble(hist) / dble(we_np)
         open(333, file=namefile)
         do i=0,extn%NE
          write(333,*) i,extn%iH(i),extn%iL(i),extn%E0(i)*Eh_to_eV,hout(i)
         enddo
         close(333)
        else
         call mpi_ssend(hist,extn%NE,mpi_integer,0,0,we_comm,ierr)
        endif
        endsubroutine print_exciton_checks
!=====================================================================!
!                           !
!=====================================================================!
        subroutine print_3D_matrix(Nx,Ny,Nz,E_site,namefile)
        implicit none
        integer,          intent(in) :: Nx, Ny, Nz
        double precision, intent(in) :: E_site(Nx,Ny,Nz)
        character(*),     intent(in) :: namefile
        ! Local variables:
        integer :: ix, iy, iz
        open(333, file=namefile)
        do iz=1, Nz
         do iy=1, Ny
          do ix=1, Nx
           write(333,*) ix, iy, iz, E_site(ix,iy,iz) * Eh_to_eV
          enddo
         enddo
        enddo
        flush(333)
        close(333)
        write(screen,fmchr3) 'The file', namefile, 'has been printed.'
        flush(screen)
        endsubroutine print_3D_matrix
!=====================================================================!
!                           !
!=====================================================================!
        subroutine print_trajectory1(u,chrg,iter,time_AU)
        use type_charge
        implicit none
        type(charge),     intent(in) :: chrg
        integer,          intent(in) :: u, iter
        double precision, intent(in) :: time_AU
        double precision             :: time_SI
        integer                      :: i
        double precision             :: a
        a = 2.0d0
        time_SI = time_AU * hEh_to_s
        write(u,*) chrg%N
        write(u,'(" iteration:",i10,"; time: ",e20.12)') iter, time_SI
        do i=1,chrg%N
         write(u,'(1x,"O",3(3x,f12.1))') chrg%R(i)%x*a, chrg%R(i)%y*a, chrg%R(i)%z*a
        enddo
        flush(u)
        end subroutine print_trajectory1
!=====================================================================!
!                           !
!=====================================================================!
        subroutine print_trajectory2(u,holes,elecs,sings,trips,iter,time_AU)
        use type_charge
        implicit none
        type(charge),     intent(in) :: holes, elecs, sings, trips
        integer,          intent(in) :: u, iter
        double precision, intent(in) :: time_AU
        double precision             :: time_SI
        integer                      :: i, N
        double precision             :: a
        a       = 2.0d0
        time_SI = time_AU * hEh_to_s
        N       = 0
        do i=1,holes%N; if(holes%exists(i)) N = N + 1; enddo        
        do i=1,elecs%N; if(elecs%exists(i)) N = N + 1; enddo        
        N       = N + sings%N + trips%N
        write(u,*) N
        write(u,'(" iteration:",i10,"; time: ",e20.12)') iter, time_SI
        do i=1,holes%N
         if(holes%exists(i)) write(u,10) "O", holes%R(i)%x*a, holes%R(i)%y*a, holes%R(i)%z*a
        enddo
        do i=1,elecs%N
         if(elecs%exists(i)) write(u,10) "N", elecs%R(i)%x*a, elecs%R(i)%y*a, elecs%R(i)%z*a
        enddo
        do i=1,sings%N
         write(u,10) "Au", sings%R(i)%x*a, sings%R(i)%y*a, sings%R(i)%z*a
        enddo
        do i=1,trips%N
         write(u,10) "Cu", trips%R(i)%x*a, trips%R(i)%y*a, trips%R(i)%z*a
        enddo
        flush(u)
10      format(1x,a,3(3x,f12.1))
        endsubroutine print_trajectory2
!--------------------------------------------------------------------!
! Cube file writer adapted from code by J. Eng
!--------------------------------------------------------------------!
     SUBROUTINE print_site_flux(u,site_flux,iter)
     use kinds
     use param
     implicit none 
     integer, intent(in)                      :: u
     integer, intent(in)                      :: iter
     integer, dimension(Nx,Ny,Nz), intent(in) :: site_flux
     integer                                  :: ix,iy,iz
     double precision                                 :: x,y,z
     write(u,*) 'Site Flux X10^6 in .cube format '
     write(u,*) 'Generated by Thomas Pope, Yvelin Giret and Julien Eng'
     write(u,100) Nx*Ny*Nz, 0.0d0, 0.0d0, 0.0d0
     write(u,100) Nx, lattice, 0.0d0, 0.0d0 
     write(u,100) Ny, 0.0d0, lattice, 0.0d0
     write(u,100) Nz, 0.0d0, 0.0d0, lattice 
     do ix=1,Nx; x=ix*lattice
       do iy=1,Ny; y=iy*lattice
         do iz=1,Nz; z=iz*lattice
       write(u,'(I5,4F13.6)') 1, 1.0d0,x,y,z
     enddo; enddo; enddo
     write(u,'(i5)') 1
     do ix=1,Nx; do iy=1,Ny
       write(u,'(6E19.9)') 1e6*dble(site_flux(ix,iy,:))/dble(iter)
     enddo; enddo
     flush(u)
100  format(I10,3F13.6)
     end subroutine print_site_flux
!--------------------------------------------------------------------!
! Cube file writer adapted from code by J. Eng
!--------------------------------------------------------------------!
     SUBROUTINE print_site_fluxes(u,site_flux1,site_flux2,iter)
     use kinds
     use param
     implicit none 
     integer, intent(in)                      :: u
     integer, intent(in)                      :: iter
     integer, dimension(Nx,Ny,Nz), intent(in) :: site_flux1,site_flux2
     integer                                  :: ix,iy,iz
     double precision                                 :: x,y,z
     write(u,*) 'Site Fluxes in .cube format '
     write(u,*) 'Generated by Thomas Pope, Yvelin Giret and Julien Eng'
     write(u,100) Nx*Ny*Nz, 0.0d0, 0.0d0, 0.0d0
     write(u,100) Nx, lattice, 0.0d0, 0.0d0 
     write(u,100) Ny, 0.0d0, lattice, 0.0d0
     write(u,100) Nz, 0.0d0, 0.0d0, lattice 
     do ix=1,Nx; x=ix*lattice
       do iy=1,Ny; y=iy*lattice
         do iz=1,Nz; z=iz*lattice
       write(u,'(I5,4F13.6)') 1, 1.0d0,x,y,z
     enddo; enddo; enddo
     write(u,'(i5)') 2
     do ix=1,Nx; do iy=1,Ny
       write(u,'(6E19.9)') (dble(site_flux1(ix,iy,iz))/dble(iter),dble(site_flux2(ix,iy,iz))/dble(iter),iz=1,Nz)
     enddo; enddo
     flush(u)
100  format(I10,3F13.6)
     end subroutine print_site_fluxes
end module printings
