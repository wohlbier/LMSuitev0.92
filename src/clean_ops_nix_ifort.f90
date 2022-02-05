subroutine clean_outputs
  use iflport !! this gives the systemqq function
  use parameters, only : max_movies
  use parameters, only : run_parameters
  use parameters, only : output_parameters

  logical :: file_exists, shell
  integer :: i
  character(10) :: dir_name
  character(25) :: dir_path

  !! remove directories movie_i for all i up to max_movies
  clean:do i = 1, max_movies
     if (i <= 9) then
        write(dir_name,fmt="(A6,I1)")'movie_',i
        write(dir_path,fmt="(A10,A7,A1)") './outputs/',dir_name,'/'
     elseif (i >= 10 .and. i <= 99) then
        write(dir_name,fmt="(A6,I2)")'movie_',i
        write(dir_path,fmt="(A10,A8,A1)") './outputs/',dir_name,'/'
     endif
     
     !! if the directory exists, remove it
     inquire(file=dir_path, exist=file_exists)
     if (file_exists) then
        shell = systemqq("cd outputs ; rm -rf "//dir_name)
     endif
  end do clean

  !! Delete all files in outputs directory
  !! put this between removal of directories and creation so *nix version
  !! doesn't report that it can't remove directories
  if (output_parameters % clean_outputs) then
     print*, 'CLEAN_OUTPUTS=TRUE: REMOVING ALL FILES IN ./outputs/'
     print*, ''
     shell = systemqq("cd ./outputs ; rm -f *")
  end if


  !! check if movie_i exists. no: make it, yes: clean it out
  make_dirs:do i=1,run_parameters % num_movie_namelists 
     if (i <= 9) then
        write(dir_name,fmt="(A6,I1)")'movie_',i
        write(dir_path,fmt="(A10,A7,A1)") './outputs/',dir_name,'/'
     elseif (i >= 10 .and. i <= 99) then
        write(dir_name,fmt="(A6,I2)")'movie_',i
        write(dir_path,fmt="(A10,A8,A1)") './outputs/',dir_name,'/'
     endif

     !! make the new directory
     shell = systemqq("cd outputs ; mkdir "//dir_name//" ; cd ..")
  end do make_dirs

end subroutine clean_outputs
