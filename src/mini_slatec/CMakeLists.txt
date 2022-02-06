set (LIBSRC
     d1mach.f
     dbesi0.f
     dbesi1.f
     dbesk0.f
     dbesk1.f
     dbsi0e.f
     dbsi1e.f
     dbsk0e.f
     dbsk1e.f
     dcsevl.f
     fdump.f
     i1mach.f
     initds.f
     j4save.f
     xercnt.f
     xerhlt.f
     xermsg.f
     xerprn.f
     xersve.f
     xgetua.f)

add_library (mini_slatec_s ${LIBSRC})
