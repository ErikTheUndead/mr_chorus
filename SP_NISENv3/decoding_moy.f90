      recursive subroutine decoding_moy(maille)

      use mod_common
      use mod_structure
      use mod_recherche

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      INTEGER :: i, j, k, l, m, n
      INTEGER :: indice, increment
      INTEGER :: niv_sup, sj, sk
      INTEGER, DIMENSION(ndim) :: ind_loc, ind_sup, ind_cor
      INTEGER, DIMENSION(3) :: ind

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

! ========== Parcours de l'arbre a partir de la racine =============

      if( associated(maille%fils(1)%ptr) ) then

! +++++Interpolation de la valeur du fils sur le niveau grossier

! *****Valorisation du support pour un arbre graduel (-s:+s)

! .....niveau et indice de la maille courante
        niv_sup = maille%maille_niv
        ind_sup(:) = maille%maille_i(:)

        do m = 1,nvar
          f_support(0,0,0,m) = maille%u(m)
        enddo
          
! .....Recherche de tous les points du support

! .....polynome de degre >= 1
        if(s > 0 ) call recherche_support(maille, s, f_support)

! *****calcul des polynomes d'interpolation

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
            
! *****calcul de utilde par polynome d'interpolation

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

! .....prediction sur le maillage grossier

            u_tilde(m) = u_tilde(m) - ( &
             - (-1)**ind(1)*qx(m) - (-1)**ind(2)*qy(m) - (-1)**ind(3)*qz(m) &
             + (-1)**(ind(1)+ind(2))*qxy(m) + (-1)**(ind(2)+ind(3))*qyz(m) &
             + (-1)**(ind(2)+ind(3))*qxz(m) &
             - (-1)**(ind(1)+ind(2)+ind(3))*qxyz(m) )

! .....valeurs sur les fils (au niveau plus fin)

            maille%fils(l)%ptr%u(m) = u_tilde(m) + maille%det(m, l)

          enddo

        enddo
!
! +++++On passe aux fils
!
        do l = 1, 2**ndim
          call decoding_moy(maille%fils(l)%ptr)
        enddo

      endif

      end subroutine decoding_moy
