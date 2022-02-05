##########################################################
##		Main makefile for LMSuite		##
##########################################################

# choose platform
PLATFORM = linux
#PLATFORM = unix
#PLATFORM = windows

# edit file ./make_files/makefile_$(PLATFORM) to select
# your compiler

.PHONY : all
all:
	cd make_files && make -f makefile_$(PLATFORM)

.PHONY : realclean clean

realclean : clean
	cd make_files && make -f makefile_$(PLATFORM) realclean

clean :
	cd make_files && make -f makefile_$(PLATFORM) clean
	$(RM) lmsuite
