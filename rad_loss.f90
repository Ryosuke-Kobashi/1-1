! Calculate the radiative losses using Townsend et al. (2009) exact integration scheme20190718
! Uses piecewise powerlaw cooling curve
subroutine rad_loss

use sweepsize , only : maxsweep
use global, only : i_CD, i_FS, dt, smallp, xm_p, xm_e, k_b, year!gamm is comment-out
use zone , only : imax, jmax, kmax, zro, zpr_tot, zpr_gas, zpr_CR, zpr_B, zT_p, zT_e, zpr_gas_ratio, zpr_CR_ratio, zpr_PB_ratio, zgam, zgamm, zn_pro, zn_ele, zn_ion !zgam=gas+CR
use cool

implicit none
integer :: i, j, k, i_min, i_max
integer :: k_broken
double precision :: TminCoolCurveLog, TmaxCoolCurveLog, scaleCoolCurve, scaleCoolCurveI, ytef, ytef2, yinv
double precision, dimension(nCoolCurveBins) :: y_tef, a_broken, temp_ref, cool_ref
double precision, dimension(imax) :: temp_hydro, coolfunc, t_cool, n_a, n_b, n_c !mu_e, mu_H
double precision, dimension(imax) :: temp_hydronew, dtemp, dpres, P_eff

double precision, parameter :: kBoltzmann = 1.3806d-16      ! Boltz in cgs
double precision, parameter :: mSol       = 1.989d+33       ! Msun, cgs
double precision, parameter :: yr         = 3.156d+07       ! sec per year
double precision, parameter :: pc         = 3.08d+18        ! cm per parsec
double precision, parameter :: mu         = 1.035d-24       ! mean particle mass in g
double precision, parameter :: cmpkm      = 1.0d+05         ! 10^5 cm per km
double precision, parameter :: gamm       = 2.0d0/3.0d0     ! gamma-1, where gamma = 5/3


!-----------------------------------------------------------------------

!Implicit cooling code w/ CIE cooling curve
!!i=401
!!write(*,*) "zpr_ratios"; write(*,*) zpr_PB_ratio(i,1,1),zpr_CR_ratio(i,1,1),zpr_gas_ratio(i,1,1); write(*,*) zpr_B(i,1,1),zpr_CR(i,1,1),zpr_gas(i,1,1),zpr_tot(i,1,1)
! Calculate cooling-curve-dependent parameters
!!write(*,*) "a"
TminCoolCurveLog = log10(TminCoolCurve)
TmaxCoolCurveLog = log10(TmaxCoolCurve)
scaleCoolCurve   = (dble(nCoolCurveBins)-1.d0) / (TmaxCoolCurveLog-TminCoolCurveLog)
scaleCoolCurveI  = 1.d0 / scaleCoolCurve !Twidth of one bin

! Set up necessary parameters to calculate cooling
! use a double timestep so radiate need only be called once

!(assumes solar abundances, fully ionized H and He)
!(this happens at logT~6, see Table 6 in S&D))
!muOverk = mu/kBoltzmann
!kOvermu = kBoltzmann/mu
!convFactor = -0.25d0 * gamm * (2.d0*dt) / (mu*kBoltzmann)
!the units of this are K*s*erg^-1*g^-1

!smallTmpLoss = 1.0d-30
!nIterationsSecant = 8

!coolCurve(:) = 1.d-22

i_min = i_CD+1
i_max = i_FS-1
!do k = 1, kmax, 1  !no need?-no need
!    do j = 1, jmax, 1   !jmax=kmax=1 by zonemod.f90


!logT_k, logLambda_k, alpha_k
do k = 1, nCoolCurveBins
   temp_ref(k) = TminCoolCurveLog + scaleCoolCurveI*(k-1) !log
   cool_ref(k) = log10(coolCurve(k)) !log
enddo
!!write(*,*) "b"
do k = 1, nCoolCurveBins-1
   a_broken(k) = (cool_ref(k+1)-cool_ref(k))/(temp_ref(k+1)-temp_ref(k))
enddo

!!write(*,*) "c"
    ! broken power-law
    y_tef(nCoolCurveBins) = 0.d0
    do k = nCoolCurveBins-1, 1, -1
        if(a_broken(k) /= 1.d0) then
          y_tef(k)=y_tef(k+1) - 10**cool_ref(nCoolCurveBins)*10**temp_ref(k)*(1.d0-(10**temp_ref(k)/10**temp_ref(k+1))**(a_broken(k)-1.d0))&
            /((1.d0-a_broken(k))*10**cool_ref(k)*10**temp_ref(nCoolCurveBins))
        else
          y_tef(k)=y_tef(k+1) - 10**cool_ref(nCoolCurveBins)*10**temp_ref(k)*log(10**temp_ref(k)/10**temp_ref(k+1))&
               /(10**cool_ref(k)*10**temp_ref(nCoolCurveBins))
       endif
!!write(*,*) "d"
    enddo

!!write(*,*) "e"
    do i = i_min, i_max
!write(*,*) "f,i=",i
!        mu_e(i) = zro(i,1,1)/zn_ele(i,1,1)   !soon after eq.(5)
!        mu_H(i) = zro(i,1,1)/zn_pro(i,1,1)   !zn_ele&zn_pro is full-ionized and CIE also mu
       n_a(i) = zn_ele(i,1,1)
       n_b(i) = zn_pro(i,1,1)
       n_c(i) = sum(zn_ion(:,i,1,1))+zn_ele(i,1,1)
!       n_b(i) = sum(zn_ion(:,i,1,1))
!       n_b(i) = sum(zn_ion(:,i,1,1)*z_Aatom(:,i,1,1))
!       n_c(i) = sum(zn_ion(:,i,1,1))+zn_ele

!!!!!!!!temp_hydro(i) = mu*zpr_gas(i,1,1)/(zro(i,1,1)*kBoltzmann)
          
       !determination of ytef and yinv
       !write(*,*) "zT_e",zT_e(i,1,1)
!!write(*,*) "zpr_gas,tot", zpr_gas(i,1,1),zpr_tot(i,1,1)
       if(zT_e(i,1,1) <= 10**TminCoolCurveLog) then
!write(*,*) "g,zT_e=",zT_e(i,1,1),zT_p(i,1,1)
cycle
!if(zpr_gas(i,1,1)<1.d-50)stop
         else
!write(*,*) "h,zT_e=",zT_e(i,1,1),zT_p(i,1,1)
            k=1   !want k to run
            do j=1,nCoolCurveBins-1
               if(zT_e(i,1,1) >= 10**temp_ref(j)) then
                  k = j !eq.(A5)
!write(*,*) "running k", k
               else
                   exit
               endif
            enddo
            !write(*,*) "i,k"   !for test
            !write(*,*) "zT_e", zT_e(i,1,1)
            !write(*,*) "n_a,n_b,n_c"; write(*,*) n_a(i),n_b(i),n_c(i)
!            write(*,*) "dt", dt
            
                if(a_broken(k) /= 1.d0) then
                    ytef = y_tef(k)+ 10**cool_ref(nCoolCurveBins)*10**temp_ref(k)*(1.d0-(10**temp_ref(k)/zT_e(i,1,1))**(a_broken(k)-1.d0))&
                        /((1.d0-a_broken(k))*10**cool_ref(k)*10**temp_ref(nCoolCurveBins))!!!k can be used????
                else
                    ytef = y_tef(k)+ 10**cool_ref(nCoolCurveBins)*10**temp_ref(k)*log(10**temp_ref(k)/zT_e(i,1,1))&
                        /(10**cool_ref(k)*10**temp_ref(nCoolCurveBins))
                 endif
!write(*,*) "cool_ref,cool_ref,ytef,ytef(k)"; write(*,*) cool_ref(nCoolCurveBins),cool_ref(k),ytef,y_tef(k)

                coolfunc(i) = 10**cool_ref(k)*(zT_e(i,1,1)/10**temp_ref(k))**(a_broken(k))   !eq.(A4)+(13)
                t_cool(i) = k_b*n_c(i)*zT_e(i,1,1)&
                    /(gamm*n_a(i)*n_b(i)*coolfunc(i))   !gamm=5./3.

                ytef2 = ytef+zT_e(i,1,1)*10**cool_ref(nCoolCurveBins)*(2.d0*dt)&
                     /(10**temp_ref(nCoolCurveBins)*coolfunc(i)*t_cool(i))                                 ! 2.d0*dt
!write(*,*) "coolfunc,t_cool,ytef2;YO"; write(*,*) coolfunc(i),t_cool(i)/year,ytef2
!!!!!!!!!!!!!!!!write(*,*) "test",1.d0/10**temp_ref(k),10**temp_ref(k)*10**temp_ref(k),1.d0/10**temp_ref(k)*10**temp_ref(k) !! =OK!!1/Tk,Tk^2,1
                
                if(a_broken(k) /= 1.d0) then
                    yinv = 10**temp_ref(k)*(1.d0-(1.d0-a_broken(k))*10**cool_ref(k)*10**temp_ref(nCoolCurveBins)*(ytef2-y_tef(k))&
                            /(10**cool_ref(nCoolCurveBins)*10**temp_ref(k)))**(1.d0/(1.d0-a_broken(k)))
                else
                    yinv = 10**temp_ref(k)*exp(-10**cool_ref(k)*10**temp_ref(nCoolCurveBins)*(ytef2-y_tef(k))&
                            /(10**cool_ref(nCoolCurveBins)*10**temp_ref(k)))
                 endif
!write(*,*) "yinv=", yinv
!write(*,*) "j"
                temp_hydronew(i) = yinv
!write(*,*) "temp_new", temp_hydronew(i)
                
                dtemp(i) = temp_hydronew(i)-zT_e(i,1,1)   !negative
                dpres(i) = dtemp(i)*zn_ele(i,1,1)*k_b
                
!write(*,*) "dT,dP", dtemp(i),dpres(i)
                  ! zT_e(i,1,1) = zT_e(i,1,1)+dtemp(i)   !zT_i??no zT??
                    zT_e(i,1,1) = max(zT_e(i,1,1)+dtemp(i),TminCoolCurve)
                 !zT_p(i,1,1) = max(zT_p(i,1,1)+dtemp(i),TminCoolCurve)
                 !zT_i(:,i,1,1) = zT_i(:,i,1,1)+dtemp(i) ! other elements are not included now
                    zpr_gas(i,1,1) = max(zpr_gas(i,1,1)+dpres(i),smallp)
                    zpr_tot(i,1,1) = max(zpr_tot(i,1,1)+dpres(i),smallp)
                    zpr_gas_ratio(i,1,1) = zpr_gas(i,1,1)/zpr_tot(i,1,1)
                    zpr_CR_ratio(i,1,1) = zpr_CR(i,1,1)/zpr_tot(i,1,1)
                    zpr_PB_ratio(i,1,1) = zpr_B(i,1,1)/zpr_tot(i,1,1)
!write(*,*) "zprnew_g,t", zpr_gas(i,1,1), zpr_tot(i,1,1)
                    
                    P_eff(i) = 2.5e0*zpr_gas_ratio(i,1,1)+4.e0*zpr_CR_ratio(i,1,1)+2.e0*zpr_PB_ratio(i,1,1)   !g_b=2!solving equation of pressure using gamma
                    zgam(i,1,1) = P_eff(i)/(P_eff(i)-1.e0)
                    zgamm(i,1,1) = zgam(i,1,1)-1.e0
!write(*,*) "zgamm?>1.33-1",zgamm(i,1,1)
!write(*,*) "i,zT_e,dT,zne,dP"; write(*,*) i, zT_e(i,1,1), dtemp(i), zn_ele(i,1,1), dpres(i)
!write(*,*) "zpr_ratios"; write(*,*) zpr_PB_ratio(i,1,1),zpr_CR_ratio(i,1,1),zpr_gas_ratio(i,1,1)!!; write(*,*) zpr_B(i,1,1),zpr_CR(i,1,1),zpr_gas(i,1,1),zpr_tot(i,1,1)
            !enddo
        endif
     enddo
!write(*,*) "l"
return
end subroutine

