C  Test the consistency of the CRD FORTRAN read/write routines.
c
      character*256
     .          argv(10),
     .          in_name,
     .          out_name
      character*512
     .          line_in, line_out
       integer  argc,
     .          argl(12),
     .          trimlen,
     .          lengtht


      argc= iargc()
      do 10 i=1,argc
      call getarg(i,argv(i))
 10   argl(i)= trimlen(argv(i))

      in_name= argv(1)(1:argl(1))
      out_name= argv(2)(1:argl(2))
      lengtht= trimlen (in_name)
      write(*,*) "in: ",in_name(1:lengtht)
      lengtht= trimlen (out_name)
      write(*,*) "out: ",out_name(1:lengtht)

      open(1,file=in_name,status='old',err=98,iostat=ioerr)
      open(2,file=out_name,status='unknown',err=98,iostat=ioerr)

 20   read(1,'(a)',end=999,err=99,iostat=ioerr) line_in
        if (line_in(1:2).eq.'h1') then
           call read_h1 (line_in)
           call write_h1 (line_out)
        endif
        if (line_in(1:2).eq.'h2') then
           call read_h2 (line_in)
           call write_h2 (line_out)
        endif
        if (line_in(1:2).eq.'h3') then
           call read_h3 (line_in)
           call write_h3 (line_out)
        endif
        if (line_in(1:2).eq.'h4') then
           call read_h4 (line_in)
           call write_h4 (line_out)
        endif
        if (line_in(1:2).eq.'h8') then
           call read_h8 (line_in)
           call write_h8 (line_out)
        endif
        if (line_in(1:2).eq.'h9') then
           call read_h9 (line_in)
           call write_h9 (line_out)
        endif
        if (line_in(1:2).eq.'c0') then
           call read_c0 (line_in)
           call write_c0 (line_out)
        endif
        if (line_in(1:2).eq.'c1') then
           call read_c1 (line_in)
           call write_c1 (line_out)
        endif
        if (line_in(1:2).eq.'c2') then
           call read_c2 (line_in)
           call write_c2 (line_out)
        endif
        if (line_in(1:2).eq.'c3') then
           call read_c3 (line_in)
           call write_c3 (line_out)
        endif
        if (line_in(1:2).eq.'c4') then
           call read_c4 (line_in)
           call write_c4 (line_out)
        endif
        if (line_in(1:2).eq.'00') then
           call read_00 (line_in)
           call write_00 (line_out)
        endif
        if (line_in(1:2).eq.'10') then
           call read_10 (line_in)
           call write_10 (line_out)
        endif
        if (line_in(1:2).eq.'11') then
           call read_11 (line_in)
           call write_11 (line_out)
        endif
        if (line_in(1:2).eq.'12') then
           call read_12 (line_in)
           call write_12 (line_out)
        endif
        if (line_in(1:2).eq.'20') then
           call read_20 (line_in)
           call write_20 (line_out)
        endif
        if (line_in(1:2).eq.'21') then
           call read_21 (line_in)
           call write_21 (line_out)
        endif
        if (line_in(1:2).eq.'30') then
           call read_30 (line_in)
           call write_30 (line_out)
        endif
        if (line_in(1:2).eq.'40') then
           call read_40 (line_in)
           call write_40 (line_out)
        endif
        if (line_in(1:2).eq.'50') then
           call read_50 (line_in)
           call write_50 (line_out)
        endif
        if (line_in(1:2).eq.'60') then
           call read_60 (line_in)
           call write_60 (line_out)
        endif
      lengtht= trimlen(line_out)
      write(2,'(a)') line_out(1:lengtht)
      go to 20

 98   write(*,*) "Open Error"
      stop
 99   write(*,*) "Read Error"
      stop
 999  write(*,*) "Done!"
      stop
      end
