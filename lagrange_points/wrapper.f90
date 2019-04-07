module wrapper
    
    
    use DOPRI
    use Diff_ope
    
   
    

    contains

    !SUBRUTINA PROPAGATOR: Subrutina que adapta los distintos propagadores para simplificar la interfaz.
    !ARGUMENTOS:
    !PROP: character ("ODEX", "DOP853", "ODE113") que define el propagador a utilizar
    !SYS: character ("non_linear", "linear") que define si se utiliza el problema completo o linealizado
    !GL: integer, número de ecuaciones del sistema, o grados de libertad
    !t: real, valor inicial de la variable independiente (tiempo)
    !x: vector of reals, de tamaño GL, contiene las condiciones iniciales del problema, y devuelve los valores finales en tf
    !tf: real, valor final de la variable independiente.
    
    subroutine propagator( prop, sys, GL, t, x, tf ) 
    
    character(len=*) :: prop, sys
    
    integer, intent (in) :: GL
    real, intent(in) :: t, tf
    real,intent(inout) ::  x(:)
    
    
    
    real :: WORK(226),RPAR(1),RTOL(1),ATOL(1), H
    real :: relerr,abserr
    integer:: i, N, LWORK, IDID, IWORK(39), LIWORK,ITOL,IPAR(1), IOUT, iflag!LWORK87!LIWORK21
    
    
    WORK=0.
    RPAR=0.
    RTOL=1d-11
    ATOL=1d-11
    H = 0.
    LWORK=226
    IWORK=0
    LIWORK=39
    ITOL=0
    IPAR=0
    IOUT=0
    
    relerr = 1d-11
    abserr = 1d-11
    iflag = 1
    
    if (prop == "ODEX") then
        
        if (sys == "non_linear")then
        
        call ODEX(GL,F_prop,t,x ,tf ,H,& !DOP853
            RTOL,ATOL,ITOL,&
            SOLOUT,IOUT,&
            WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
        
        else if (sys == "linear")then
            
            
        call ODEX(GL,F_lin,t,x ,tf ,H,& !DOP853
            RTOL,ATOL,ITOL,&
            SOLOUT,IOUT,&
            WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
            
            
            
        end if
        
    
    
    else if (prop == "DOP853") then
        
        
        if (sys == "non_linear")then
        
        call DOP853(GL,F_prop,t,x ,tf ,&
            RTOL,ATOL,ITOL,&
            SOLOUT,IOUT,&
            WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
        
        
        else if (sys == "linear")then
            
            
        call DOP853(GL,F_lin,t,x ,tf ,&
            RTOL,ATOL,ITOL,&
            SOLOUT,IOUT,&
            WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
            
            
        end if
        
        
        
    else if (prop == "ODE113")then
        
        
        if (sys == "non_linear")then
        
            call ode(F_prop_ODE,GL,x,t,tf,relerr,abserr,iflag,work,iwork)
        
        
        else if (sys == "linear")then
            
            call ode(F_lin_ODE,GL,x,t,tf,relerr,abserr,iflag,work,iwork)
            
        end if
        
    end if
    
    end subroutine
    
    
    
    
    
    
     SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,&
                                      RPAR,IPAR,IRTRN,XOUT)
                    real RPAR(*),IPAR(*),XOUT,X,CON,Y,XOLD
                    integer NR,ND,N,ICOMP,IRTRN
                    
                    
     end subroutine



end module
    