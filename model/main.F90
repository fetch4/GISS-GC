program main
  implicit none
  include 'netcdf.inc'

  integer :: rc, fid
  character(len=256) :: fname

  fname = 'TOPO'
  rc = nf_open(trim(fname), nf_nowrite, fid)
  if (rc /= nf_noerr) then
    print *, 'Error opening file: ', trim(fname)
    stop 1
  end if
  print *, "Opened ", trim(fname)

  rc = nf_close(fid)
  if (rc /= nf_noerr) then
    print *, 'Error closing file: ', fid
    stop 1
  end if
  print *, "Closed ", trim(fname)
  call modelE_mainDriver()
end program main
