MODULE AnalysisGLI
PUBLIC::InvGLI
PUBLIC::ScanningGLI

CONTAINS
!*=======================================================================================!
!*                          Convert GLI into Internal File                               !
!*=======================================================================================!
SUBROUTINE InvGLI(fn2,fn_tmp)
	USE global
	USE gli

IMPLICIT NONE
CHARACTER(LEN=256)::na,fn2,fn_tmp
INTEGER(KIND=4)::id,x,y,z,istat,i,nu
REAL(KIND=4)::fb,unass1,unass2

i=1

fn_tmp=TRIM(ADJUSTL(fn2))//'_tmp.dat'
INQUIRE(IOLENGTH=nu) geo
OPEN(11,FILE=fn2,ACCESS='SEQUENTIAL',STATUS='OLD')
OPEN(12,FILE=fn_tmp,ACCESS='DIRECT',STATUS='REPLACE',RECL=nu*length)
	DO  WHILE(.TRUE.)
		READ(11,*,IOSTAT=istat) na
		IF(istat/=0) EXIT
		IF(na=='LINE') THEN
			BACKSPACE(11)
			READ(11,*) geo%na,geo%id
			WRITE(12,REC=i) geo
			i=i+1
		ELSE IF(na=='SHOT') THEN
			BACKSPACE(11)
			READ(11,*) geo%na,geo%id,geo%x,geo%y,geo%z,geo%fb,geo%unass1,geo%unass2
			WRITE(12,REC=i) geo
			i=i+1
		ELSE IF(na=='TRACE') THEN
			BACKSPACE(11)
			READ(11,*) geo%na,geo%id,geo%x,geo%y,geo%z,geo%fb
			WRITE(12,REC=i) geo
			i=i+1
		END IF
		
	END DO
CLOSE(12)
CLOSE(11)
END SUBROUTINE InvGLI
!========================================================================================!

SUBROUTINE ScanningGLI(line_start,lines,line_interval,&
				 shot_start,shots,shot_interval,&
				 fn_tmp,table,ntr)
	USE global
	USE gli

IMPLICIT NONE
CHARACTER(LEN=256),INTENT(IN)::fn_tmp
INTEGER::line_start,lines,line_interval,&
				 shot_start,shots,shot_interval,&
				 l,s,nu,istat,m,n,i,ntr
INTEGER::table(lines,shots,3),table1(lines)

i=1
m=0
n=0
ntr=0
INQUIRE(IOLENGTH=nu) geo
OPEN(11,FILE=fn_tmp,ACCESS='DIRECT',STATUS='OLD',RECL=nu*length)
DO l=1,lines
	DO WHILE(.TRUE.)
		READ(11,REC=i,IOSTAT=istat) geo
		IF(istat/=0) EXIT
		i=i+1
		IF(geo%na=='LINE')THEN
			m=m+1
			IF(m==(line_start+(l-1)*line_interval))THEN
				table1(l)=i-1
				EXIT
			END IF
		END IF
	END DO
END DO

DO l=1,lines
	i=table1(l)+1
	n=0
	DO s=1,shots
		DO WHILE(.TRUE.)
			READ(11,REC=i,IOSTAT=istat) geo
			IF(istat/=0) EXIT
			i=i+1
			IF(geo%na=='SHOT')THEN
				n=n+1
				IF(n==(shot_start+(s-1)*shot_interval))THEN
					table(l,s,1)=i-1
					EXIT
				END IF
			END IF
		END DO
	END DO
END DO

DO l=1,lines
	DO s=1,shots
		i=table(l,s,1)+1
		DO WHILE(.TRUE.)
			READ(11,REC=i,IOSTAT=istat) geo
			IF(istat/=0) THEN
				table(l,s,2)=i-table(l,s,1)-1
				EXIT
			END IF
			i=i+1
			IF(geo%na=='SHOT'.OR.geo%na=='LINE')THEN
				table(l,s,2)=i-table(l,s,1)-2
				EXIT
			END IF
		END DO
	END DO
END DO

DO l=1,lines
	DO s=1,shots
		table(l,s,3)=1+ntr
		ntr=ntr+table(l,s,2)
	END DO
END DO

CLOSE(11)
END SUBROUTINE ScanningGLI
END MODULE AnalysisGLI
