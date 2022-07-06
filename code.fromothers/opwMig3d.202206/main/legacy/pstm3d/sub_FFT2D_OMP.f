      SUBROUTINE FFT2D_OMP(A2D, AT2D, N1, N2, INV)
C     Perform 2D FFT with flag: INV. =0, forward; =1, backward.
C     A2D(N1,N2)  is the input 2d buffer, overwritten by output.
C     AT2D(N1,N2) is a temporal buffer.
C     Author: FengBo, 2010-07-17.
      IMPLICIT NONE
      INTEGER N1, N2, INV
      COMPLEX A2D(N1,N2), AT2D(N2,N1)

      INTEGER I1, I2

C     1D FFT to the fast-dimension.
C$OMP PARALLEL
C$OMP DO
      DO I2 = 1, N2
         CALL FFT1D(A2D(1,I2), N1, INV)
      END DO
C$OMP END DO
C$OMP END PARALLEL

C     Transpose A2D to  AT2D
      DO I1 = 1, N1
      DO I2 = 1, N2
         AT2D(I2,I1) = A2D(I1,I2)
      END DO
      END DO

C     1D FFT to the slow-dimension, using the transposed buffer AT2D
C$OMP PARALLEL
C$OMP DO
      DO I1 = 1, N1
         CALL FFT1D(AT2D(1,I1), N2, INV)
      END DO
C$OMP END DO
C$OMP END PARALLEL

C     Transpose AT2D to  A2D
      DO I1 = 1, N1
      DO I2 = 1, N2
         A2D(I1,I2) = AT2D(I2,I1)
      END DO
      END DO

      RETURN
      END SUBROUTINE FFT2D_OMP
*======================================================================*

        SUBROUTINE FFT2D(AA, Ahx, Amx, KhxNum, KmxNum, INV)

        INTEGER KmxNum, KhxNum, INV
        COMPLEX AA(KhxNum, KmxNum), Amx(KmxNum), Ahx(KhxNum)

        DO 567 Imx=1, KmxNum

         DO 667 Ihx=1, KhxNum
          Ahx(Ihx)=AA(Ihx, Imx)
667      CONTINUE

         CALL FFT1D(Ahx, KhxNum, INV)

         DO 767 Ihx=1, KhxNum
          AA(Ihx, Imx)=Ahx(Ihx)
767      CONTINUE

567     CONTINUE

        DO 568 Ihx=1, KhxNum

         DO 668 Imx=1, KmxNum
          Amx(Imx)=AA(Ihx, Imx)
668      CONTINUE

         CALL FFT1D(Amx, KmxNum, INV)

         DO 868 Imx=1, KmxNum
          AA(Ihx, Imx)=Amx(Imx)
868      CONTINUE

568     CONTINUE

        RETURN
        END SUBROUTINE FFT2D
*======================================================================*
