        parameter(ns=201)
        character a,b,c,d

        open(10,file='shot_cnooc_fault_id.txt',status='unknown')        

        do is=1,ns
        a=char((is-is/10000*10000)/1000+48)
        b=char((is-is/1000*1000)/100+48)
        c=char((is-is/100*100)/10+48)
        d=char((is-is/10*10)+48)
        write(*,*)  a,b,c,d
        write(10,*)  a,b,c,d
        enddo

        close(10)
 
        END
