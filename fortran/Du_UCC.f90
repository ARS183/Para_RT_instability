!
!============================Parallel==============================
!

subroutine Du2Dx_PADE4(u,d2f,pxflag)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
	integer :: pxflag

    
    if (pxflag==0) then
        if (npx==0) then
            call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,1)
            call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,1)
            call D2F_PADE4(u,d2f,nx,hx)

        elseif (npx==npx0-1) then
            call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,2)
            call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,2)
            call D2F_PADE4(u,d2f,nx,hx)

        else
            call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,0)
            call D2F_PADE4(u,d2f,nx,hx)

        endif
    elseif (pxflag==1) then
        call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,0)
        call D2F_PADE4(u,d2f,nx,hx)
    endif

!    deallocate(d2f,df)

end subroutine Du2Dx_PADE4

!------------------------------------------------------

subroutine Du2Dy_PADE4(u,d2f,pyflag)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:ny+LAP),d2f(1-LAP:ny+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:ny+LAP)
	integer :: pyflag

    if (pyflag==0) then
        if (npy==0) then
            call OCFD_D2F_BOUND_PADE4(u,d2f,ny,hy,1)
            call OCFD_D2F_SB_PADE4(u,d2f,ny,hy,1)
            call D2F_PADE4(u,d2f,ny,hy)

        elseif (npy==npy0-1) then
            call OCFD_D2F_BOUND_PADE4(u,d2f,ny,hy,2)
            call OCFD_D2F_SB_PADE4(u,d2f,ny,hy,2)
            call D2F_PADE4(u,d2f,ny,hy)

        else
            call OCFD_D2F_SB_PADE4(u,d2f,ny,hy,0)
            call D2F_PADE4(u,d2f,ny,hy)

        endif
    elseif (pyflag==1) then
        call OCFD_D2F_SB_PADE4(u,d2f,ny,hy,0)
        call D2F_PADE4(u,d2f,ny,hy)
    endif
!    deallocate(d2f,df)

end subroutine Du2Dy_PADE4


!
!========================================================-
!




!
!===============´®ÐÐ===================
!
subroutine Du2Dx_PADE4_serial(u,d2f)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
!	integer :: nx

        call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,1)

        call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,2)

        call D2F_PADE4(u,d2f,nx,hx)

!    deallocate(d2f,df)

end subroutine Du2Dx_PADE4_serial


!--------------------------------------------------------

subroutine Du2Dy_PADE4_serial(u,d2f)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:ny+LAP),d2f(1-LAP:ny+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:ny+LAP)
!	integer :: nx

        call OCFD_D2F_BOUND_PADE4(u,d2f,ny,hy,1)

        call OCFD_D2F_BOUND_PADE4(u,d2f,ny,hy,2)

        call D2F_PADE4(u,d2f,ny,hy)

!    deallocate(d2f,df)

end subroutine Du2Dy_PADE4_serial


!
!========================================
!




!====================================================================
!-------------------------------Parallel----------------------------
!======================================================================

subroutine DuDx_UCC_UpWind(u,d2f,df,pxflag)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
	integer :: pxflag

    if (pxflag==0) then
        if (npx==0) then
            call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)
            call DF_UCC45_P(u,d2f,df,nx,hx,1)

        elseif (npx==npx0-1) then
            call OCFD_DF_SB2_UCC45_P(u,df,d2f,nx,hx)
            call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)
            call DF_UCC45_P(u,d2f,df,nx,hx,2)

        else
            call OCFD_DF_SB2_UCC45_P(u,df,d2f,nx,hx)
            call DF_UCC45_P(u,d2f,df,nx,hx,0)
        endif
    elseif (pxflag==1) then
        call OCFD_DF_SB2_UCC45_P(u,df,d2f,nx,hx)
        call DF_UCC45_P(u,d2f,df,nx,hx,0)
    endif

!    deallocate(d2f,df)

end subroutine DuDx_UCC_UpWind

!------------------------------------------------------------------------------

subroutine DuDy_UCC_UpWind(u,d2f,df,pyflag)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:ny+LAP),d2f(1-LAP:ny+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:ny+LAP)
	integer :: pyflag

    if (pyflag==0) then
        if (npy==0) then
            call OCFD_DF_BOUND_UCC45(u,df,ny,hy,1)
            call DF_UCC45_P(u,d2f,df,ny,hy,1)

        elseif (npy==npy0-1) then
            call OCFD_DF_SB2_UCC45_P(u,df,d2f,ny,hy)
            call OCFD_DF_BOUND_UCC45(u,df,ny,hy,2)
            call DF_UCC45_P(u,d2f,df,ny,hy,2)

        else
            call OCFD_DF_SB2_UCC45_P(u,df,d2f,ny,hy)
            call DF_UCC45_P(u,d2f,df,ny,hy,0)
        endif
    elseif (pyflag==1) then
        call OCFD_DF_SB2_UCC45_P(u,df,d2f,ny,hy)
        call DF_UCC45_P(u,d2f,df,ny,hy,0)
    endif
!    deallocate(d2f,df)

end subroutine DuDy_UCC_UpWind



!------------------------------------------------------------------------------

subroutine DuDx_UCC_DownWind(u,d2f,df,pxflag)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
	integer :: pxflag

    if (pxflag==0) then
        if (npx==0) then

            call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)
            call OCFD_DF_SB2_UCC45_M(u,df,d2f,nx,hx)
            call DF_UCC45_M(u,d2f,df,nx,hx,1)

        elseif (npx==npx0-1) then

            call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)
            call DF_UCC45_M(u,d2f,df,nx,hx,2)

        else
            call OCFD_DF_SB2_UCC45_M(u,df,d2f,nx,hx)
            call DF_UCC45_M(u,d2f,df,nx,hx,0)
        
        endif
    elseif (pxflag==1) then
        call OCFD_DF_SB2_UCC45_M(u,df,d2f,nx,hx)
        call DF_UCC45_M(u,d2f,df,nx,hx,0)
    endif

!    deallocate(d2f,df)

end subroutine DuDx_UCC_DownWind




!---------------------------------------------------------------------

subroutine DuDy_UCC_DownWind(u,d2f,df,pyflag)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:ny+LAP),d2f(1-LAP:ny+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:ny+LAP)
	integer :: pyflag

    if (pyflag==0) then

        if (npy==0) then

            call OCFD_DF_BOUND_UCC45(u,df,ny,hy,1)
            call OCFD_DF_SB2_UCC45_M(u,df,d2f,ny,hy)
            call DF_UCC45_M(u,d2f,df,ny,hy,1)

        elseif (npy==npy0-1) then

            call OCFD_DF_BOUND_UCC45(u,df,ny,hy,2)
            call DF_UCC45_M(u,d2f,df,ny,hy,2)

        else
            call OCFD_DF_SB2_UCC45_M(u,df,d2f,ny,hy)
            call DF_UCC45_M(u,d2f,df,ny,hy,0)
        
        endif

    elseif (pyflag==1) then
        call OCFD_DF_SB2_UCC45_M(u,df,d2f,ny,hy)
        call DF_UCC45_M(u,d2f,df,ny,hy,0)
    endif

!    deallocate(d2f,df)

end subroutine DuDy_UCC_DownWind




!===============================================================
!------------------------------------------------------------
!===============================================================




!
!============================Serial==========================
!
subroutine DuDx_UCC_UpWind_serial(u,d2f,df)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
!	integer :: nx

 
        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)

        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)

        call DF_UCC45_P(u,d2f,df,nx,hx,12)

!    deallocate(d2f,df)

end subroutine DuDx_UCC_UpWind_serial

!-------------------------------------------------------
subroutine DuDy_UCC_UpWind_serial(u,d2f,df)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:ny+LAP),d2f(1-LAP:ny+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:ny+LAP)
!	integer :: nx

 
        call OCFD_DF_BOUND_UCC45(u,df,ny,hy,1)

        call OCFD_DF_BOUND_UCC45(u,df,ny,hy,2)

        call DF_UCC45_P(u,d2f,df,ny,hy,12)

!    deallocate(d2f,df)

end subroutine DuDy_UCC_UpWind_serial

!-----------------------------------------------------------

subroutine DuDx_UCC_DownWind_serial(u,d2f,df)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
!	integer :: nx

        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)

        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)
  
        call DF_UCC45_M(u,d2f,df,nx,hx,12)


!    deallocate(d2f,df)

end subroutine DuDx_UCC_DownWind_serial

!----------------------------------------------------------

subroutine DuDy_UCC_DownWind_serial(u,d2f,df)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:ny+LAP),d2f(1-LAP:ny+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:ny+LAP)
!	integer :: nx

        call OCFD_DF_BOUND_UCC45(u,df,ny,hy,1)

        call OCFD_DF_BOUND_UCC45(u,df,ny,hy,2)
  
        call DF_UCC45_M(u,d2f,df,ny,hy,12)


!    deallocate(d2f,df)

end subroutine DuDy_UCC_DownWind_serial

!
!============================================================
!