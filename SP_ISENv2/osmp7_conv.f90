      subroutine osmp7_conv(dir, sigma, u_n, C_mp, flux)

      use mod_common

      IMPLICIT NONE

! +++++Variables arguments

      INTEGER, INTENT(IN) :: dir

      REAL(kind=8), INTENT(IN) :: sigma
      REAL(kind=8), INTENT(IN), DIMENSION(nvar, -nordre+1:nordre) :: u_n
      REAL(kind=8), INTENT(OUT), DIMENSION(nvar) ::  C_mp,flux

      INTEGER :: l, lm, m
      INTEGER :: is

! +++++Variables locales

      REAL(kind=8), DIMENSION(-nordre+1:nordre-1) :: cp12, hp12, ecp12
      REAL(kind=8), DIMENSION(3, -nordre+1:nordre-1) :: up12

      REAL(kind=8), DIMENSION(-1:1) :: d4i, d4ip1, d4min, dmm, dm4
      REAL(kind=8), DIMENSION(-1:1) :: dip12
      REAL(kind=8), DIMENSION(-nordre+1:nordre) :: d2
      REAL(kind=8), DIMENSION(-nordre+2:nordre-2) :: ri, rm, rp
      REAL(kind=8) :: signfi, signfip1
      REAL(kind=8) :: phi_o, phi_tvd, phi_tvdmin, phi_tvdmax
      REAL(kind=8) :: gmdp, glcm
      REAL(kind=8) :: phi_min, phi_max, phi_md, phi_lc
      REAL(kind=8), DIMENSION(nvar) ::  phi_lim
      REAL(kind=8), DIMENSION(-nordre+1:nordre) :: co7_1, co7_2, &
                                                   co7_3, co7_4, co7_5

      REAL(kind=8), DIMENSION(nvar) ::  delu, deluf
      REAL(kind=8), DIMENSION(-nordre+1:nordre,nvar) ::  delw, delwnu
      REAL(kind=8), DIMENSION(-nordre+1:nordre,nvar) :: vp, xnu


!       print*,"convection", u_n(:,0)


! ================== Definition des moyennes aux interfaces =================

		do l = -nordre+1,nordre-1
        do m = 2, nvar-1
           up12(m-1,l) = 0.5 * (u_n(m,l+1)/u_n(1,l+1) + u_n(m,l)/u_n(1,l))
        enddo
      enddo

! +++++valeurs propres aux demi-mailles + CFL local

     do l = -nordre+1,nordre-1
			do m = 1,nvar
			vp(l,m) =  up12(dir,l)
			xnu(l,m) = sigma * abs(up12(dir,l))
			end do
     end do

#if OS7 .OR. TVDMP
! +++++Invariants de Riemann

      do l = -nordre+1,nordre-1
					do m = 1, nvar
          delw(l,m) = u_n(m,l+1) - u_n(m,l)
          end do
		end do

! .....delwnu = (1.-xnu) * delw

		do l = -nordre+1,nordre-1
			do m = 1,nvar
			delwnu(l,m) = (1.-xnu(l,m))*delw(l,m)
			end do
		enddo
#endif

! +++++Terme correctif par onde

      do m = 1, nvar

#if OS7 .OR. TVDMP

! .....coefficients du schema O7
        do l = -nordre+1,nordre-1
          co7_1(l) = (1.+xnu(l,m))/3.
          co7_2(l) = co7_1(l) * (xnu(l,m)-2.)/4.
          co7_3(l) = co7_2(l) * (xnu(l,m)-3.)/5.
          co7_4(l) = co7_3(l) * (xnu(l,m)+2.)/6.
          co7_5(l) = co7_4(l) * (xnu(l,m)+3.)/7.
        enddo

! +++++rapport des invariants

        do l = -nordre+2, nordre-2

! .....calcul de r+
          signfi = sign(1.,delwnu(l-1,m))

          signfip1 = sign(1.,delwnu(l,m))

          rp(l) = signfi*signfip1* &
                  (abs(delwnu(l-1,m))+zero)/(abs(delwnu(l,m))+zero)

! .....calcul de r-
          signfi = sign(1.,delwnu(l+1,m))

          signfip1 = sign(1.,delwnu(l,m))

          rm(l) = signfi*signfip1* &
                  (abs(delwnu(l+1,m))+zero)/(abs(delwnu(l,m))+zero)

        enddo
#endif

#ifdef TVDMP
! .....mesure de la courbure
        do l = -nordre+2, nordre-1
          d2(l) = vp(l,m)*delw(l,m) - vp(l-1,m)*delw(l-1,m)
        enddo
        if( nordre == 2 ) then
          d2(-1) = d2(0)
          d2(2) = d2(1)
        endif

! .....minmod entre i et i+1
        do l = -1, 1
          if( d2(l)*d2(l+1) < 0. ) then
            dmm(l) = 0.
          else
            dmm(l) = min(abs(d2(l)) , abs(d2(l+1)))
            if(d2(l) < 0.) dmm(l) = -dmm(l)
          endif
        enddo

! .....mesure plus severe
        do l = -1, 1
          d4i(l) = 4.*d2(l) - d2(l+1)
          d4ip1(l) = 4.*d2(l+1) - d2(l)
          if( d4i(l)*d4ip1(l)  < 0. ) then
            d4min(l) = 0.
          else
            d4min(l) = min(abs(d4i(l)) , abs(d4ip1(l)))
            if(d4i(l) < 0.) d4min(l) = -d4min(l)
          endif
          if( d4min(l)*dmm(l) < 0.) then
            dm4(l) = 0.
          else
            dm4(l) = min(abs(d4min(l)) , abs(dmm(l)))
            if(d4min(l) < 0.) dm4(l) = -dm4(l)
          endif

          dip12(l) = dm4(l)
        enddo
#endif

#if OS7 .OR. TVDMP
! .....Cas vp > 0
        if(vp(0,m) > 0.) then
          is = 1

          do l = -nordre+2, nordre-2
            ri(l) = rp(l)
          enddo

! .....Cas vp <= 0
        else
          is = -1

          do l = -nordre+2, nordre-2
            ri(l) = rm(is*l)
          enddo

        endif

! +++++fonction de precision

        phi_o=1.
! .....ordre 3
        if(nordre >= 2) phi_o=phi_o-(co7_1(0)-ri(0)*co7_1(-is))

        if(nordre >= 3) then
! .....ordre 4
          phi_o=phi_o+(co7_2(0)-2.*ri(0)*co7_2(-is) &
                      +ri(0)*ri(-1)*co7_2(-2*is))
! .....ordre 5
          phi_o=phi_o-(co7_3(is)/ri(+1)-3.*co7_3(0) &
                      +3.*ri(0)*co7_3(-is)-ri(0)*ri(-1)*co7_3(-2*is))
        endif

        if(nordre == 4) then
! .....ordre 6
          phi_o=phi_o+(co7_4(2*is)/(ri(+1)*ri(+2)) &
                      -4.*co7_4(is)/ri(+1)+6.*co7_4(0) &
                      -4.*ri(0)*co7_4(-is)+ri(0)*ri(-1)*co7_4(-2*is))
! .....ordre 7
          phi_o=phi_o-(co7_5(2*is)/(ri(+1)*ri(+2)) &
                      -5.*co7_5(is)/ri(+1)+10.*co7_5(0) &
                      -10.*ri(0)*co7_5(-is) &
                      +5.*ri(0)*ri(-1)*co7_5(-2*is) &
                      -ri(0)*ri(-1)*ri(-2)*co7_5(-3*is))
        endif
#endif

#ifdef TVDMP

!+++++limiteur de flux TVD

! !        phi_tvdmin = 2.*ri(0)*(1.-xnu(-is,m))/ &
! !                     (xnu(0,m)*(1.-xnu(0,m))+zero)
        phi_tvdmin = 2.*ri(0)/(xnu(-is,m)+zero)
        phi_tvdmax = 2./(1.-xnu(0,m)+zero)
        phi_tvd = max(0.,min(phi_tvdmin,phi_o,phi_tvdmax))

! +++++limiteur de flux MP

! !        gmdp = 0.5 * (1. - dip12(0)/(delwnu(0,m)+zero))
! !        glcm = 0.5 * (1. + dip12(-is)/(delwnu(0,m)+zero))
        gmdp = 0.5 * (1. - dip12(0)/(vp(0,m)*delw(0,m)+zero))
        glcm = 0.5 * (1. + dip12(-is)/(vp(0,m)*delw(0,m)+zero))

        phi_md = 2.*gmdp/(1.-xnu(0,m)+zero)
        phi_lc = 2.*ri(0)*glcm*(1.-xnu(-is,m))/ &
                 (xnu(0,m)*(1.-xnu(0,m))+zero)

        phi_min = max(min(0.,phi_md), min(0., phi_tvdmin, phi_lc))
        phi_max = min(max(phi_tvdmax,phi_md), max(0., phi_tvdmin, phi_lc))

        phi_lim(m) = max(phi_min, min(phi_o, phi_max))
#endif

#ifdef ROE
! .....Schema de Roe
       phi_lim(m) = 0.
#endif
!!! .....Schema de Lax-Wendroff
!!        phi_lim(m) = 1.
!!! .....Schema de Beam-Warming
!!        phi_lim(m) = ri(0)
!!! .....Schema TVD
!!        phi_lim(m) = phi_tvd

#ifdef OS7
!.....Schema OS7
       phi_lim(m) = phi_o
#endif

      enddo
!
! +++++Calcul de la correction C(i+1/2) + F(i+1/2) convection

		do m = 1, nvar

			delu(m) = u_n(m,1) - u_n(m,0)

			C_mp(m) =  - 0.5 * (1.-phi_lim(m)*(1.-xnu(0,m))) * abs(vp(0,m)) * delu(m)

			flux(m) = 0.5 * vp(0,m) * delu(m)

		enddo


	end subroutine osmp7_conv