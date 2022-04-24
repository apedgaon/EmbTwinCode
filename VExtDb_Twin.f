      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
        implicit none
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       Variables passed from interface
        integer lOp, niArray, nrArray
        integer i_Array(niArray)
        real*8 r_Array(nrArray)

C       Contents of i_Array
        integer i_int_nTotalNodes, i_int_nTotalElements,
     1          i_int_kStep, i_int_kInc, i_int_iStatus,
     2          i_int_lWriteRestart
        parameter( i_int_nTotalNodes     = 1,
     *            i_int_nTotalElements  = 2,
     *            i_int_kStep           = 3,
     *            i_int_kInc            = 4,
     *            i_int_iStatus         = 5,
     *            i_int_lWriteRestart   = 6  )

C       Possible values for the lOp argument
        integer j_int_StartAnalysis, j_int_StartStep,
     1          j_int_SetupIncrement, j_int_StartIncrement,
     2          j_int_EndIncrement, j_int_EndStep,
     3          j_int_EndAnalysis
        parameter(j_int_StartAnalysis    = 0,
     *            j_int_StartStep        = 1,
     *            j_int_SetupIncrement   = 2,
     *            j_int_StartIncrement   = 3,
     *            j_int_EndIncrement     = 4,
     *            j_int_EndStep          = 5,
     *            j_int_EndAnalysis      = 6 )

C       Possible values for i_Array(i_int_iStatus)
        integer j_int_Continue, j_int_TerminateStep,
     1          j_int_TerminateAnalysis
        parameter(j_int_Continue          = 0,
     *            j_int_TerminateStep     = 1,
     *            j_int_TerminateAnalysis = 2)

C       Contents of r_Array
        integer i_flt_TotalTime, i_flt_StepTime,
     1          i_flt_dTime
        parameter(i_flt_TotalTime   = 1,
     1            i_flt_StepTime    = 2,
     2            i_flt_dTime       = 3 )

        integer kStep
        integer kInc

        kStep = i_Array(i_int_kStep)
        kInc  = i_Array(i_int_kInc)

        if (lop.eq.j_int_StartAnalysis) then
          call sleep(20)
          call read_gb_elements()
        end if

      end subroutine


      subroutine read_gb_elements()
        implicit none
        character(len = 50) :: gbfilename
        character (len = 1) :: pathseperator = '/'
        character (len = 256) :: outdir, gbfilepath
        integer :: lenoutdir, i, io, num_gb_elems

        gbfilename = 'gb_elems_1000.txt'
        call vgetoutdir(outdir, lenoutdir)
        gbfilepath = trim(outdir)//pathseperator//trim(gbfilename)

        num_gb_elems = 0
        open(120, file = gbfilepath, status = 'old', iostat = io)
        do
          read(120, *, iostat = io)
          if (io /= 0) exit
          num_gb_elems = num_gb_elems + 1
        end do
        close(120)

      end subroutine