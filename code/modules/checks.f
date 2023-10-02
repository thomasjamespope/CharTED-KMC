!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module checks.f90: ...description...
!%
!% Copyright (C) 2020 Yvelin Giret
!% \author Yvelin Giret
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module checks
   use kinds
   use constants
   implicit none

   contains

   subroutine check_events(max_events,events,iteration,open_path,closed_path)
       implicit none
       integer, intent(in) :: max_events, events, iteration
       integer, intent(in) :: open_path, closed_path

       if(max_events .ne. events) then
           write(screen,*) 'WARNING: max_events .NE. events at &
                  & iteration:', iteration
           flush(screen)
        end if
        if(events .ne. open_path + closed_path) then
           write(screen,*) 'WARNING: events .NE. open_path + closed_path &
                  & at iteration:', iteration, events, open_path, closed_path
           flush(screen)
        end if
   end subroutine check_events

   subroutine check_site_number(Nx,Ny,Nz,correlated_type)
        implicit none
        integer,      intent(inout) :: Nx, Ny, Nz
        character(*), intent(in)    :: correlated_type

        ! local variables:
        integer :: Nx_keep, Ny_keep, Nz_keep

        Nx_keep = Nx
        Ny_keep = Ny
        Nz_keep = Nz

        if(correlated_type == 'ANTIFERRO') then
           if(mod(Nx,2) .ne. 0) Nx = Nx + 1
           if(mod(Ny,2) .ne. 0) Ny = Ny + 1
           if(mod(Nz,2) .ne. 0) Nz = Nz + 1
           if(Nx.ne.Nx_keep .or. Ny.ne.Ny_keep .or. Nz.ne.Nz_keep) then
            write(screen,fmdlm1) blank1
            write(screen,fmchr1) 'For antiferro-like correlated disorder the number &
                            &of sites in each direction must be even:'
            if(Nx.ne.Nx_keep) write(screen,10) 'X: ', Nx_keep, ' ===> ', Nx
            if(Ny.ne.Ny_keep) write(screen,10) 'Y: ', Ny_keep, ' ===> ', Ny
            if(Nz.ne.Nz_keep) write(screen,10) 'Z: ', Nz_keep, ' ===> ', Nz
            write(screen,fmdlm2) blank1, delim1             
           endif
        end if
10      format(t2,"|",t5,a,2x,i0,2x,a,2x,i0,t100,"|")
   end subroutine check_site_number
end module checks
