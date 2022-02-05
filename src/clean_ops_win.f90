subroutine clean_outputs
  use parameters, only : max_movies
  use parameters, only : run_parameters
  use parameters, only : output_parameters

  logical :: file_exists
  integer :: i
  character(10) :: dir_name
  character(25) :: dir_path

  !! Delete all files in outputs directory
  if (output_parameters % clean_outputs) then
     print*, 'CLEAN_OUTPUTS=TRUE: REMOVING ALL FILES IN ./outputs/'
     print*, ''
     call system ("cd outputs & echo y|del *.* >nul & cd ..")
  end if

  !! remove directories movie_i for i larger than num_movie_namelists
  clean:do i = run_parameters % num_movie_namelists + 1, max_movies
     if (i <= 9) then
        write(dir_name,fmt="(A6,I1)")'movie_',i
        write(dir_path,fmt="(A10,A7,A4)")'.\outputs\',dir_name,'\NUL'
     elseif (i >= 10 .and. i <= 99) then
        write(dir_name,fmt="(A6,I2)")'movie_',i
        write(dir_path,fmt="(A10,A8,A4)")'.\outputs\',dir_name,'\NUL'
     endif
     
     !! if the directory exists, remove it
     inquire(file=dir_path, exist=file_exists)
     if (file_exists) then
        call system("cd outputs & cd "//dir_name//"& &
             &echo y|del *.* >nul & cd .. & rmdir "//dir_name//" & cd ..")
     endif
  end do clean

  make_dirs:do i=1,run_parameters % num_movie_namelists 
     if (i <= 9) then
        write(dir_name,fmt="(A6,I1)")'movie_',i
        write(dir_path,fmt="(A10,A7,A4)")'.\outputs\',dir_name,'\NUL'
     elseif (i >= 10 .and. i <= 99) then
        write(dir_name,fmt="(A6,I2)")'movie_',i
        write(dir_path,fmt="(A10,A8,A4)")'.\outputs\',dir_name,'\NUL'
     endif

     inquire(file=dir_path, exist=file_exists)
     if (.not. file_exists) then
        !! create the new directory
        call system("cd outputs & mkdir "//dir_name//" & cd ..")
     else
        !! clean up the existing directories
        call system("cd outputs & cd "//dir_name//"& &
             &echo y|del *.* >nul & cd .. & cd ..")
     end if
  end do make_dirs


end subroutine clean_outputs
