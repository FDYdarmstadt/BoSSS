      program sparsesolverdemo
      implicit none
      include 'mpif.h'
      integer SolverRef, SolverRef2, ierr, mtxref
      integer mpiRank, mpiSize
      integer i, i0, iE, l, NN, icol, irow
      real*8 d
      real*8 X(100), rhs(100)

      call ilPSPinit()

c     ---------------------------------------
c     Define sparse solver
c     ---------------------------------------

c     Option a): directly in the code:
C     (a bit unhandy in fixed fromat FORTRAN)
      call ISparseSolver_FromXMLBegin('\\',ierr)
      call ISparseSolver_FromXMLSubmit('<sparsesolver \\',ierr)
      call ISparseSolver_FromXMLSubmit(' name="g">  \\ ',ierr)
      call ISparseSolver_FromXMLSubmit('<type>PCG</type>\\',ierr)
      call ISparseSolver_FromXMLSubmit('<library>monkey</library> \\ '
     * ,ierr)
      call ISparseSolver_FromXMLSubmit('<specific>\\',ierr)
      call ISparseSolver_FromXMLSubmit('<MaxIterations\\',ierr)
      call ISparseSolver_FromXMLSubmit('>20000</MaxIterations>\\ ',ierr)
      call ISparseSolver_FromXMLSubmit('<ConvergenceType>\\ ',ierr)
      call ISparseSolver_FromXMLSubmit('Relative \\ ',ierr)
      call ISparseSolver_FromXMLSubmit('</ConvergenceType> \\ ',ierr)
      call ISparseSolver_FromXMLSubmit('<Tolerance>\\',ierr)
      call ISparseSolver_FromXMLSubmit('1.0e-9\\',ierr)
      call ISparseSolver_FromXMLSubmit('</Tolerance>\\',ierr)
      call ISparseSolver_FromXMLSubmit('<DevType>CPU</DevType>\\',ierr)
      call ISparseSolver_FromXMLSubmit('</specific>\\',ierr)
      call ISparseSolver_FromXMLSubmit('</sparsesolver> \\ ',ierr)
      call ISparseSolver_FromXMLEnd(SolverRef2,ierr)

c     Option b): load from an XML file:
      call ISparseSolver_FromXMLFileF(SolverRef,'\\',
     *   '../slvcfgex.xml\\','hypre-pcg\\',ierr)

c     ------------------------
c     Create and define matrix
c     ------------------------
      call mpi_comm_rank(mpi_comm_world,mpirank,ierr)
      call mpi_comm_size(mpi_comm_world,mpisize,ierr)

      NN = 100                       ! matrix size
      i0 = NN*MPIrank/MPIsize;       ! first local index
      iE = NN*(MPIrank + 1)/MPIsize; ! first row on next proc
      L=iE-i0;                       ! number of rows on this processor

c     create matrix
      call MsrMatrix_New(MtxRef,L,NN,1,1,ierr) 


c     set matrix entries: remember to use C-indices, starting at 0
c     (we set an identity matrix)
      do i = 0,L-1
        iRow = i + i0;
        iCol = iRow;
        d = 1
        call IMutableMatrix_SetEntry(MtxRef,iRow,iCol,1.0d0,ierr)
      enddo

c     save matrix to textfile (e.g. to watch in MATLAB)
      call IMutableMatrixEx_SaveToTextFileSparseF(mtxref,
     * '\\','mtx.txt\\', ierr) ! dont forget to spec and USE termination char.


c     pass matrix to solver
      call ISparseSolver_DefineMatrix(SolverRef,MtxRef, ierr)

c     release matrix object
c     the MsrMatrix is for assambly -> mutual
c     the solver creates its own, fixed but more efficient format
      call Common_releaseObject(mtxref,ierr)

c     ------------
c     solve system
c     ------------

c     define rhs:
      do i = 1,L
        X(i) = 0
        rhs(i) = 1
      enddo

c     call solver
      call ISparseSolver_Solve(Solverref,L,X,rhs,ierr)


c     ------------
c     shutdown
c     ------------

      call ilPSPshutdown()


      end


