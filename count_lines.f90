       subroutine count_lines(nunit,nlines)
       implicit none
       integer, intent(in)            :: nunit
       integer, intent(out)           :: nlines
        integer                        :: check
        nlines = 0
        do
        read(nunit,*,iostat=check)
        if(check /= 0)exit
         nlines = nlines + 1
        end do
        rewind(nunit)
        end subroutine count_lines
