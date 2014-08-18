      recursive subroutine encoding_det(maille)

      use mod_common
      use mod_structure
      use mod_recherche

      IMPLICIT NONE

      TYPE(structure_maille), POINTER :: maille

      INTEGER :: i, j, k, l, lm, m, n
      INTEGER :: indice, increment
      INTEGER :: niv_sup, sj, sk
      INTEGER, DIMENSION(ndim) :: ind_loc, ind_sup, ind_cor
      INTEGER, DIMENSION(3) :: ind

      REAL(kind=8) :: det_loc
      REAL(kind=8), DIMENSION(nvar) :: qx, qy, qz,  qxy, qyz, qxz, qxyz
      REAL(kind=8), DIMENSION(nvar) :: u_tilde

      REAL(kind=8),DIMENSION(-s:s,-s:s,-s:s,nvar) :: f_support

      TYPE(structure_maille), pointer :: courant

! .....Allocation du support
      if(ndim == 3) then
        sj = s
        sk = s
      elseif(ndim == 2) then
        sj = s
        sk = 0
      elseif(ndim == 1) then
        sj = 0
        sk = 0
      endif

! ========== traitement des details sur le niveau grossier ===========

! .....seulement si la feuille n'est pas fictive

      if( associated(maille%fils(1)%ptr) ) then

      if( maille%fils(1)%ptr%feuille ) then

! +++++Initialisation

        do l = 1, 2**ndim

! .....initialisation du detail sur la maille courante

          do m = 1, nvar
            maille%det(m, l) = 0.
          enddo

! .....initialisation du detail sur les fils de la maille

          do lm = 1, 2**ndim
            do m = 1, nvar
              maille%fils(l)%ptr%det(m, lm) = 0.
            enddo
          enddo

        enddo

! +++++Interpolation de la valeur du fils gauche sur le niveau grossier
!
! .....Valorisation du support pour un arbre graduel (-s:+s)
!
! .....niveau et indice de la maille courante

        niv_sup = maille%maille_niv
        ind_sup(:) = maille%maille_i(:)

        do m = 1,nvar
          f_support(0,0,0,m) = maille%u(m)
        enddo

! *****Recherche de tous les points du support et valorisation

! .....polynome de degre >= 1
        if(s > 0 ) call recherche_support(maille, s, f_support)

! *****Calcul de utilde par le polynome d'interpolation

! .....calcul des differents termes du polynome

        if(s > 0) then
          do m = 1,nvar
            qx(m) = 0.
            do i = 1, s
              qx(m) = qx(m) + coef(i)* &
                      (f_support(i,0,0,m)-f_support(-i,0,0,m))
            enddo
          enddo
        else
          do m = 1,nvar
            qx(m) = 0.
          enddo
        endif

        if(sj > 0 ) then
          do m = 1,nvar
            qy(m) = 0.
            qxy(m) = 0.
            do j = 1, sj
              qy(m) = qy(m) + coef(j)* &
                      (f_support(0,j,0,m)-f_support(0,-j,0,m))

              do i = 1,s
                qxy(m) = qxy(m) + coef(j)*coef(i)* &
                         (f_support(i,j,0,m)-f_support(i,-j,0,m) - &
                          f_support(-i,j,0,m)+f_support(-i,-j,0,m))
              enddo
            enddo
          enddo
        else
          do m = 1,nvar
            qy(m) = 0.
            qxy(m) = 0.
          enddo
        endif

        if(sk > 0 ) then
          do m = 1,nvar
            qz(m) = 0.
            qxz(m) = 0.
            qyz(m) = 0.
            qxyz(m) = 0.
            do k = 1, sk
              qz(m) = qz(m) + coef(k)* &
                      (f_support(0,0,k,m)-f_support(0,0,k,m))

              do i = 1,s
                qxz(m) = qxz(m) + coef(k)*coef(i)* &
                         (f_support(i,0,k,m)-f_support(i,0,-k,m) - &
                          f_support(-i,0,k,m)+f_support(-i,0,-k,m))
              enddo

              do j = 1,sj
                qyz(m) = qyz(m) + coef(j)*coef(k)* &
                         (f_support(0,j,k,m)-f_support(0,-j,k,m) - &
                          f_support(0,j,-k,m)+f_support(0,-j,-k,m))

                do i = 1,s
                  qxyz(m) = qxyz(m) + coef(k)*coef(j)*coef(i)* &
                           (f_support(i,j,k,m)-f_support(i,j,-k,m) - &
                            f_support(i,-j,k,m)-f_support(-i,j,k,m) + &
                            f_support(i,-j,-k,m)+f_support(-i,j,-k,m) + &
                            f_support(-i,-j,k,m)-f_support(-i,-j,-k,m))
                enddo
              enddo
            enddo
          enddo
        else
          do m = 1,nvar
            qz(m) = 0.
            qxz(m) = 0.
            qyz(m) = 0.
            qxyz(m) = 0.
          enddo
        endif

! .....calcul de la valeur estimee sur le maillage grossier

        do l = 1, 2**ndim

          ind(:) = 0
          k = 0
          n = 0
          do m = ndim, 1, -1
            k = k+1
            i = ((l-1) - n)/2**(ndim-k)
            n = n + i*2**(ndim-k)

            ind(m) = i
          enddo

          do m = 1,nvar
            u_tilde(m) = f_support(0, 0, 0, m)
          enddo

          do m = 1,nvar

            u_tilde(m) = u_tilde(m) - ( &
             - (-1)**ind(1)*qx(m) - (-1)**ind(2)*qy(m) - (-1)**ind(3)*qz(m) &
             + (-1)**(ind(1)+ind(2))*qxy(m) + (-1)**(ind(2)+ind(3))*qyz(m) &
             + (-1)**(ind(2)+ind(3))*qxz(m) &
             - (-1)**(ind(1)+ind(2)+ind(3))*qxyz(m) )

! +++++calcul du detail => erreur d'approximation

            maille%det(m, l) = maille%fils(l)%ptr%u(m) - u_tilde(m)

          enddo

        enddo
!
! .....detail total pour seuillage et elagage

        det_loc = 0.
        do l = 1, 2**ndim

         do m = 1,nvar
! .....norme L2
!!!!            det_loc = det_loc + maille%det(m, l)**2
! .....norme L1
           det_loc = det_loc + abs( maille%det(m, l) )
         enddo

!          det_loc = det_loc + maille%det(1, l)**2
!           det_loc = det_loc + abs(maille%det(1, l))
        enddo

!!        maille%detail = sqrt(det_loc / float(2**ndim))
        maille%detail = det_loc / float(2**ndim)
!!        maille%detail = det_loc / float(nvar)

! .....valeur extreme du detail
        det_max = max( det_max, maille%detail )

! .....On passe au niveau plus fin ==> fils gauche et droite

        do l = 1, 2**ndim
          call encoding_det(maille%fils(l)%ptr)
        enddo

! +++++Sinon : on annule le detail

      else

! .....initialisation du detail sur la maille courante

        do l = 1, 2**ndim
          do m = 1, nvar
            maille%det(m, l) = 0.
          enddo
        enddo

        maille%detail = 0.

      endif

      endif

      end subroutine encoding_det
