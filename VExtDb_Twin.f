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

      end subroutine