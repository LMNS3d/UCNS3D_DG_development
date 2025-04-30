MODULE PRESTORE
USE DECLARATION
USE DERIVATIVES
USE LIBRARY
USE BASIS
USE TRANSFORM
USE LOCAL
USE LAPCK
USE DG_FUNCTIONS
IMPLICIT NONE


 CONTAINS
 
 
 SUBROUTINE PRESTORE_1(N)
!> @brief
!> This subroutine calls other subroutines to prestore the pseudoinverse reconstruction least square matrices
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,i
KMAXE=XMPIELRANK(N)

if (dimensiona.eq.3)then
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO schedule(STATIC)
DO I=1,KMAXE

CALL LOCALISE_STENCIL(N,I)

CALL LOCALISE_STEN2(N,I)



call CHECKGRADS(N,I)



CALL FIND_ROT_ANGLES(N,I)

CALL PRESTORE_RECONSTRUCTION3(N,I)
iconsidered=i
 if (dg.eq.1)then
 
 CALL PRESTORE_DG1
    end if
 

END DO
!$OMP END DO
!$OMP END PARALLEL



else

!$OMP PARALLEL DEFAULT(SHARED) 



!$OMP DO
DO I=1,KMAXE


CALL LOCALISE_STENCIL2d(N,I)



CALL LOCALISE_STEN2d(N,I)


 call CHECKGRADS2d(N,I)
 


CALL FIND_ROT_ANGLES2D(N,I)


 CALL PRESTORE_RECONSTRUCTION2(N,I)

 
 ICONSIDERED=I
 if (dg.eq.1)then
 
 CALL PRESTORE_DG1
 
 
    end if
END DO
!$OMP END DO
!$OMP END PARALLEL





!takis this needs to be within a omp parallel region!




end if





END SUBROUTINE PRESTORE_1






SUBROUTINE PRESTORE_RECONSTRUCTION3(N,iconsi)
!> @brief
!> This subroutine prestores the pseudoinverse reconstruction least square matrices
IMPLICIT NONE
INTEGER,INTENT(IN)::N,iconsi
INTEGER::I,J,K,llco,ll,ii,igf,IGF2,IFD2,idum,idum2,iq,jq,lq,IHGT,IHGJ,iqp,iqp2,NND,k0,g0,lcou,lcc,iqqq,ICOND1,ICOND2
real::ssss,gggg,UPTEMP,LOTEMP,X_STENCIL,Y_STENCIL,Z_STENCIL,DIST_STEN,DIST_STEN2
real, dimension(IELEM(N,Iconsi)%inumneighbours-1, ielem(n,iconsi)%idegfree):: stencil
real, dimension(20,IELEM(N,Iconsi)%inumneighbours-1, ielem(n,iconsi)%idegfree):: stencilS
real, dimension(ielem(n,iconsi)%idegfree, ielem(n,iconsi)%idegfree):: invmat
REAL,DIMENSION(7,IELEM(N,Iconsi)%inumneighbours-1)::WLSQR

i=iconsi


IDUM=0;
                if (ielem(n,i)%interior.eq.1)then
                        DO j=1,IELEM(N,I)%IFCA
                        if (ielem(n,i)%ibounds(J).gt.0)then
                            if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
                                IDUM=1
                            end if
                        END IF
                        END DO
                end if
                





        ICONSIDERED=I
		
		INTBS=zero
		JXX=1;IXX=i;LXX1=1
		number_of_dog=ielem(n,i)%idegfree
		kxx=ielem(n,i)%iorder
		ELTYPE=ielem(n,i)%ishape
		compwrt=0
		INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1)
		
                INTEG_BASIS(I)%VALUE(1:ielem(n,i)%IDEGFREE)=INTBS(1:ielem(n,i)%IDEGFREE)
                
               
		IF (IWENO.EQ.1)THEN
		CALL INDICATORMATRIX(N,I)
		
		END IF
		
		if (ees.eq.5)then
		INTBS=zero
		JXX=1;IXX=i;LXX1=1
		number_of_dog=idegfree2
		kxx=IORDER2
		ELTYPE=ielem(n,i)%ishape
                compwrt=1
		INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1)
		
                INTEG_BASIS(I)%VALUEc(1:number_of_dog)=INTBS(1:number_of_dog)
                
                
                
                
            ICONSIDERED=I
                
		IF (IWENO.EQ.1)THEN
		CALL INDICATORMATRIX2(N,I)
		
		compwrt=0
		
		END IF
		
		end if

		
		LLCO=IELEM(N,I)%ADMIS


        IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	 INUMO=ielem(n,i)%iorder
	 imax2=numneighbours2-1
	 inum2=numneighbours2
	 ideg2=iDEGFREE2
	 inumo2=iorder2

                DIST_STEN=ZERO
				
				DIST_STEN=-tolbig; DIST_STEN2=tolbig
	       DO LL=1,ielem(n,i)%iNUMNEIGHBOURS
if (ilocal_elem(1)%VOLUME(1,LL).lt.DIST_STEN2)then
DIST_STEN2=ilocal_elem(1)%VOLUME(1,LL)
end if
if (ilocal_elem(1)%VOLUME(1,LL).gt.DIST_STEN)then
DIST_STEN=ilocal_elem(1)%VOLUME(1,LL)
end if
			end do
			
IELEM(N,I)%STENCIL_DIST=MAX((DIST_STEN/DIST_STEN2),(DIST_STEN2/DIST_STEN))
				
				
		DO LL=1,LLCO	!ADMIS

			
! !-------------------FOR DEBUGGING ONLY -----------------------------------------!
! ! 			
! !-------------------FOR DEBUGGING ONLY -----------------------------------------!


if((ees.ne.5).or.(ll.eq.1))then
IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	 INUMO=ielem(n,i)%iorder
 compwrt=0
 number_of_dog=IDEG

else
imax=numneighbours2-1
	 inum=numneighbours2
	 ideg=iDEGFREE2
	 inumo=IORDER2
	 number_of_dog=IDEG
 compwrt=1
end if
DO K=1,imax
ixx=i
kxx=INUMO


X_STENCIL=(ilocal_elem(1)%XXC(ll,k+1)-ilocal_elem(1)%XXC(ll,1))**2
Y_STENCIL=(ilocal_elem(1)%YYC(ll,k+1)-ilocal_elem(1)%YYC(ll,1))**2
Z_STENCIL=(ilocal_elem(1)%ZZC(ll,k+1)-ilocal_elem(1)%ZZC(ll,1))**2

DIST_STEN2=SQRT(X_STENCIL+Y_STENCIL+Z_STENCIL)

IF (WEIGHT_LSQR.EQ.1)THEN
	    WLSQR(ll,K)=1.0D0/((DIST_STEN2))
	    ELSE
	    WLSQR(ll,K)=1.0D0
	    END IF




	    !ILOCAL_RECON3(I)%WEIGHTL(LL,K)=WLSQR(K)

!DIST_STEN=MAX(DIST_STEN,DIST_STEN2)

!IELEM(N,I)%STENCIL_DIST=WLSQR(K)

			
			




if (fastest.eq.1)then
x1 = ilocal_elem(1)%XXC(ll,k+1)-ilocal_elem(1)%XXC(ll,1)
y1 = ilocal_elem(1)%YYC(ll,k+1)-ilocal_elem(1)%YYC(ll,1)
z1 = ilocal_elem(1)%ZZC(ll,k+1)-ilocal_elem(1)%ZZC(ll,1)


IF (GREENGO.EQ.0)THEN

if (idum.eq.1)then
    if((ees.ne.5).or.(ll.eq.1))then
    ILOCAL_RECON3(I)%STENCILS(LL,K,1:ielem(n,i)%idegfree)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,ielem(n,i)%iorder,IXX,ielem(n,i)%idegfree)
    ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
    else
    compwrt=1
    ILOCAL_RECON3(I)%STENCILSc(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG)
    ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
    compwrt=0
    end if
else
if((ees.ne.5).or.(ll.eq.1))then

 compwrt=0
STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG)
else
 compwrt=1
STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG)
 compwrt=0
end if
end if

ELSE
 if((ees.ne.5).or.(ll.eq.1))then
 compwrt=0
STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG)
else
 compwrt=1
STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG)
 compwrt=0

end if
END IF
else


IXX=i;jxx=k+1;lxx1=ll
ELTYPE=ilocal_elem(1)%ishape(ll,k+1)
IF (GREENGO.EQ.0)THEN
if (idum.eq.1)then
		if((ees.ne.5).or.(ll.eq.1))then
		 compwrt=0
		ILOCAL_RECON3(I)%STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IDEG)
		ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
		else
		 compwrt=1
		ILOCAL_RECON3(I)%STENCILSc(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IDEG)
		ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
		 compwrt=0
		end if
else
    if((ees.ne.5).or.(ll.eq.1))then
   compwrt=0
    STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IDEG)
    else
    compwrt=1
    STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IDEG)
    compwrt=0
    end if
end if
ELSE
if((ees.ne.5).or.(ll.eq.1))then
   compwrt=0
    STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IDEG)
    else
    compwrt=1
    STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IDEG)
    compwrt=0
    end if

END IF


end if




END DO	!IMAXEDEGFREE








IF ((GREENGO.EQ.0))THEN
     if (idum.eq.1)then

! DO IQ=1,ideg
! DO JQ=1,iDEG
invmat=zero
LSCQM=zero
! END DO	!IDEGFREE
! END DO	!IMAXDEGFREE

   

    DO IQ=1,IDEG
    DO JQ=1,IDEG
    DO LQ=1,IMAX
    if((ees.ne.5).or.(ll.eq.1))then
    LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
    +((ILOCAL_RECON3(I)%STENCILS(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)))
    else
    LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
    +((ILOCAL_RECON3(I)%STENCILSc(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILSc(LL,LQ,IQ)))
    end if
    END DO
    END DO
    END DO  
    else

    ! DO IQ=1,ideg
    ! DO JQ=1,iDEG
    invmat=zero
    LSCQM=zero
    ! END DO	!IDEGFREE
    ! END DO	!IMAXDEGFREE

    DO IQ=1,IDEG
    DO JQ=1,IDEG
    DO LQ=1,IMAX
    LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
    +((STENCILS(LL,LQ,JQ)*STENCILS(LL,LQ,IQ)))
    END DO
    END DO
    END DO 



    end if


ELSE

! DO IQ=1,ideg
! DO JQ=1,iDEG
invmat=zero
LSCQM=zero
! END DO	!IDEGFREE
! END DO	!IMAXDEGFREE

DO IQ=1,IDEG
DO JQ=1,IDEG
DO LQ=1,IMAX
LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
+((STENCILS(LL,LQ,JQ)*STENCILS(LL,LQ,IQ)))
END DO
END DO
END DO 

END IF




QFF(:,:)=zero; RFF(:,:)=zero; QTFF(:,:)=zero; RFF(:,:)=zero;  INVRFF(:,:)=zero 



CALL QRDECOMPOSITION(LSCQM,QFF,RFF,IDEG)
CALL TRANSPOSEMATRIX(QFF,QTFF,IDEG)
IVGT=IDEG+1
CALL INVERT(RFF,INVRFF,IVGT)
invmat(1:IDEg,1:ideg)=MATMUL(INVRFF(1:ideg,1:IDEG),QTFF(1:IDEG,1:IDEG))






 




! 
! 
! 
! 
IF (GREENGO.EQ.0)THEN
    if (idum.eq.1)then
        if((ees.ne.5).or.(ll.eq.1))then
        stencil(1:imax,1:ideg)=ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:ideg)

        else
        stencil(1:imax,1:ideg)=ILOCAL_RECON3(I)%STENCILSc(LL,1:imax,1:ideg)

        end if
    else

    stencil(1:imax,1:ideg)=STENCILS(LL,1:imax,1:ideg)
    end if
ELSE

    stencil(1:imax,1:ideg)=STENCILS(LL,1:imax,1:ideg)

END IF

if((ees.ne.5).or.(ll.eq.1))then
! call gemm(                                               &
!                invmat,                                               &
!                stencil,                                              &
!                ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),             &
!                'N',                                                  & ! transposition flag for invmat
!                'T'                                                   & ! transposition flag for stencil
!             )
            
            
call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL),IDEG)
			
			do iq=1,imax
			ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)&
			*ilocal_elem(1)%VOLUME(ll,iq+1)*WLSQR(ll,iq)
			end do
			
			
else
! call gemm(                                               &
!                invmat(1:ideg,1:ideg),                                               &
!                stencil(1:imax,1:ideg),                                              &
!                ILOCAL_RECON3(I)%invmat_stenciltc(1:IDEG,1:IMAX,LL),             &
!                'N',                                                  & ! transposition flag for invmat
!                'T'                                                   & ! transposition flag for stencil
!             )

            call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stenciltc(1:IDEG,1:IMAX,LL),IDEG)


do iq=1,imax
			ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)&
			*ilocal_elem(1)%VOLUME(ll,iq+1)*WLSQR(ll,iq)
			end do

end if



 





if (ielem(n,i)%ggs.ne.1)then	!ggs
    
IF (ITESTCASE.EQ.4)THEN	!TEST

	 

	 

IF (LL.EQ.1)THEN		!stencils
! 

  
  if (idum.eq.1)then		!for wall only
    
	 
	
!     	
    DO IHGT=1,3; DO IHGJ=1,3
    AINVJT(IHGT,IHGJ)=ILOCAL_RECON3(I)%INVCCJAC(IHGJ,IHGT)
   END DO; END DO

	 idum2=0
	 
	 
				BASEFACEVAL=zero
				BASEFACGVAL=zero
				PERMUTATION=zero
				PERMUTATIONg=zero
			
	 
	 
	 
	 DO j=1,IELEM(N,I)%IFCA		!for all faces
	    if (ielem(n,i)%ibounds(J).gt.0)then		!for bounded only
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then		!!for bounded only 2                         
					facex=J;iconsidered=i
								  CALL coordinates_face_inner(N,Iconsidered,facex)
								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
							    
								    AY=cords(2)
								    AX=cords(1)
								    AZ=cords(3)
				
					
					    VEXT(1,1)=AX;VEXT(1,2)=AY;VEXT(1,3)=AZ
					    VEXT(1,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(1,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  
					    AX=VEXT(1,1);AY=VEXT(1,2);AZ=VEXT(1,3)
				
				
				
				
				
				
				ANGLE1=IELEM(N,I)%FACEANGLEX(j)
				ANGLE2=IELEM(N,I)%FACEANGLEY(j)
				
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
				NNX=(NX*AINVJT(1,1))+(NY*AINVJT(2,1))+(NZ*AINVJT(3,1))
				NNY=(NX*AINVJT(1,2))+(NY*AINVJT(2,2))+(NZ*AINVJT(3,2))
				NNZ=(NX*AINVJT(1,3))+(NY*AINVJT(2,3))+(NZ*AINVJT(3,3))
				
				DO IQ=1, IDEG
				IF (POLY.EQ.1) THEN
				XDER(IQ)=DFX(AX,AY,AZ,IQ);  YDER(IQ)=DFY(AX,AY,AZ,IQ);  ZDER(IQ)=DFZ(AX,AY,AZ,IQ)
				END IF
				IF (POLY.EQ.2) THEN
				XDER(IQ)=DLX(AX,AY,AZ,IQ);  YDER(IQ)=DLY(AX,AY,AZ,IQ);  ZDER(IQ)=DLZ(AX,AY,AZ,IQ)
				END IF
				IF (POLY.EQ.4) THEN
				XDER(IQ)=TL3DX(AX,AY,AZ,IQ);  YDER(IQ)=TL3DY(AX,AY,AZ,IQ);  ZDER(IQ)=TL3DZ(AX,AY,AZ,IQ)
				END IF
				
				END DO

				
				BASEFACEVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,ielem(n,i)%Iorder,I,ielem(n,i)%IDEGFREE)

				if (thermal.eq.0)then
				BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=((NNX*XDER(1:ielem(n,i)%IDEGFREE))+(NNY*YDER(1:ielem(n,i)%IDEGFREE))+(NNZ*ZDER(1:ielem(n,i)%IDEGFREE)))
				ELSE
				BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,ielem(n,i)%Iorder,I,ielem(n,i)%IDEGFREE)

				end if


				DO IQ=1,IDEG
				ILOCAL_RECON3(I)%WALLCOEFF(IQ)=BASEFACEVAL(IQ)
				ILOCAL_RECON3(I)%WALLCOEFG(IQ)=BASEFACGVAL(IQ)
				PERMUTATION(IQ)=IQ
				PERMUTATIONG(IQ)=IQ
				END DO
		
				
		
		
		
		
		
		end if!for bounded only 2        
		end if!for bounded only
		end do!for all faces
		
		GGGG=-TOLBIG
G0=0
DO IQ=1,IDEG
IF (ABS(BASEFACGVAL(IQ)) >GGGG)THEN
GGGG=ABS(BASEFACGVAL(IQ))
G0=IQ
END IF
END DO




SSSS=-TOLBIG
K0=0
DO IQ=1,IDEG
IF (ABS(BASEFACEVAL(IQ)) >SSSS)THEN
SSSS=ABS(BASEFACEVAL(IQ))
K0=IQ
END IF
END DO

!  

! now swap basis functions and thus coefficients
PERMUTATION(1)=K0; PERMUTATION(K0)=1; SSSS=BASEFACEVAL(1)
BASEFACEVAL(1)=BASEFACEVAL(K0); BASEFACEVAL(K0)=SSSS

PERMUTATIONG(1)=G0; PERMUTATION(G0)=1; GGGG=BASEFACGVAL(1)
BASEFACGVAL(1)=BASEFACGVAL(G0); BASEFACGVAL(G0)=GGGG
! ILOCAL_RECON3(II)%K0=K0
! ILOCAL_RECON3(II)%G0=G0
ILOCAL_RECON3(I)%K0=K0
ILOCAL_RECON3(I)%G0=G0
		
		
		
		
LSQM = ZERO
DO LQ=1,IMAX
LCOU=0
    DO IQ=1,IDEG
	IF (IQ.EQ.G0) CYCLE
	    LCOU=LCOU+1
LSQM(LQ,LCOU)=ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)&
-ILOCAL_RECON3(I)%STENCILS(LL,LQ,G0)*ILOCAL_RECON3(I)%WALLCOEFG(IQ)/ILOCAL_RECON3(I)%WALLCOEFG(G0)
    END DO
END DO
ILOCAL_RECON3(I)%TEMPSQ(1:IMAX,1:IDEG-1)=LSQM(1:IMAX,1:IDEG-1)





VELLSQMAT=ZERO
DO IQ=1,IDEG-1; DO JQ=1,IDEG-1;	DO LCC=1,IMAX
!now store the least square matrix
VELLSQMAT(JQ,IQ)= VELLSQMAT(JQ,IQ)+(LSQM(LCC,JQ)*LSQM(LCC,IQ))
END DO;	END DO;	END DO

LSCQM=ZERO
LSCQM(1:IDEG-1,1:IDEG-1)=VELLSQMAT(1:IDEG-1,1:IDEG-1)
Qff=ZERO; Rff=ZERO; QTff=ZERO; INVRff=ZERO
	  IVGT=IDEG
	  
    CALL QRDECOMPOSITION(LSCQM,Qff,Rff,IVGT-1)
    CALL TRANSPOSEMATRIX(Qff,QTff,IVGT-1)
    
    CALL INVERT(Rff,INVRff,IVGT)
!   final inverted R^(-1)*Q^(-1)
ILOCAL_RECON3(I)%TEMPSQMAT(1:IDEG-1,1:IDEG-1) =&
 MATMUL(INVRff(1:IDEG-1,1:IDEG-1),QTff(1:IDEG-1,1:IDEG-1))

! definy the matrix defining the least-square reconstruction in this case for velocity
LSQM = ZERO
DO LQ=1,IMAX
LCOU=0
    DO IQ=1,IDEG
	IF (IQ.EQ.K0) CYCLE
	    LCOU=LCOU+1
	     LSQM(LQ,LCOU)=ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)&
 -ILOCAL_RECON3(I)%STENCILS(LL,LQ,K0)*ILOCAL_RECON3(I)%WALLCOEFF(IQ)/ILOCAL_RECON3(I)%WALLCOEFF(K0)
    END DO
END DO
ILOCAL_RECON3(I)%VELLSQ(1:IMAX,1:IDEG-1)=LSQM(1:IMAX,1:IDEG-1)





VELLSQMAT=ZERO
DO IQ=1,IDEG-1; DO JQ=1,IDEG-1;	DO LCC=1,IMAX
!now store the least square matrix
VELLSQMAT(JQ,IQ)= VELLSQMAT(JQ,IQ)+(LSQM(LCC,JQ)*LSQM(LCC,IQ))
END DO;	END DO;	END DO

LSCQM=ZERO
LSCQM(1:IDEG-1,1:IDEG-1)=VELLSQMAT(1:IDEG-1,1:IDEG-1)

Qff=ZERO; Rff=ZERO; QTff=ZERO; INVRff=ZERO
	  IVGT=IDEG
    CALL QRDECOMPOSITION(LSCQM,Qff,Rff,IVGT-1)
    CALL TRANSPOSEMATRIX(Qff,QTff,IVGT-1)
    
    CALL INVERT(Rff,INVRff,IVGT)
!   final inverted R^(-1)*Q^(-1)
ILOCAL_RECON3(I)%VELINVLSQMAT(1:IDEG-1,1:IDEG-1) = &
MATMUL(INVRff(1:IDEG-1,1:IDEG-1),QTff(1:IDEG-1,1:IDEG-1))
			
		
		
		
		
		
		
		end if!for wall only
		end if!stencils
		end if!for test
		end if!ggs
		end do
		
		



END SUBROUTINE PRESTORE_RECONSTRUCTION3


subroutine walls_higher(n)
!> @brief
!> This subroutine allocates the memory for wall bounded cells for their constrained least squares reconstruction
implicit none
integer,intent(in)::n
integer::i,j,k,imax,ideg,inumo,inum,kmaxe,idum,idum2

KMAXE=XMPIELRANK(N)

DO I=1,KMAXE

  ielem(n,i)%walls=0
 IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	 INUMO=ielem(n,i)%iorder

if (ielem(n,i)%ggs.ne.1)then

 IDUM=0;idum2=0
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
	        IDUM=1
		if (dimensiona.eq.3)then
		if (ielem(n,i)%types_faces(j).eq.5)then
		  idum2=idum2+qp_quad_n
		else
		  idum2=idum2+QP_TRIANGLE_n
		end if
		else
		  idum2=idum2+qp_line_n
		end if
	      end if
	  END IF
	END DO
  end if
  if (idum.eq.1)then
    allocate(ielem(n,i)%num_of_wall_gqp(1))
    ielem(n,i)%num_of_wall_gqp(1)=idum2
    ielem(n,i)%walls=1
	    
	    IF (FASTEST.NE.1)THEN
	    ALLOCATE (ILOCAL_RECON3(I)%VELINVLSQMAT(IDEG-1,IDEG-1))
	    ALLOCATE (ILOCAL_RECON3(I)%WALLCOEFF(ideg))
	    ALLOCATE (ILOCAL_RECON3(I)%VELLSQ(IMAX,IDEG-1))
	    ALLOCATE (ILOCAL_RECON3(I)%TEMPSQMAT(IDEG-1,IDEG-1))
	    ALLOCATE (ILOCAL_RECON3(I)%WALLCOEFG(ideg))
	    ALLOCATE (ILOCAL_RECON3(I)%TEMPSQ(IMAX,IDEg-1))
	  
	    ILOCAL_RECON3(I)%VELINVLSQMAT=zero
	    ILOCAL_RECON3(I)%WALLCOEFF=zero
	    ILOCAL_RECON3(I)%WALLCOEFG=zero
	    ILOCAL_RECON3(I)%TEMPSQMAT=zero
	   END IF
	    
  end if
  
  end if

END DO



end subroutine walls_higher


SUBROUTINE PRESTORE_RECONSTRUCTION2(N,iconsi)
!> @brief
!> This subroutine prestores the pseudoinverse reconstruction least square matrices in 2d
IMPLICIT NONE
INTEGER,INTENT(IN)::N,iconsi
REAL::distancex,distancey
INTEGER::I,J,K,llco,ll,ii,igf,IGF2,IFD2,idum,idum2,iq,jq,lq,IHGT,IHGJ,iqp,iqp2,NND,k0,g0,lcou,lcc,ICOND1,ICOND2,ai,aj
real::ssss,gggg,UPTEMP,LOTEMP,DIST_STEN2,X_STENCIL,Y_STENCIL,DIST_STEN,maxai,minai
real, dimension(IELEM(N,Iconsi)%inumneighbours-1, ielem(n,iconsi)%idegfree):: stencil
real, dimension(7,IELEM(N,Iconsi)%inumneighbours-1, ielem(n,iconsi)%idegfree):: stencilS
real, dimension(ielem(n,iconsi)%idegfree, ielem(n,iconsi)%idegfree):: invmat
REAL,DIMENSION(7,IELEM(N,Iconsi)%inumneighbours-1)::WLSQR
REAL,DIMENSION(9)::integrals_basis
real, dimension(9, 9)::testmatrix,testmatrix1
REAL::SCALLING_HUANG
integer::IGQRULES_HUANG,QP_TRIANGLE_huang,QP_QUAD_HUANG,QP_LINE_HUANG
real, dimension(3, 2):: loworder_matrix
real, dimension(2, 3):: loworder_matrix_t,loworder_matrix_final
real, dimension(2, 2):: loworder_matrix_inverse



i=iconsi

        
		IDUM=0;
                if (ielem(n,i)%interior.eq.1)then
                        DO j=1,IELEM(N,I)%IFCA
                        if (ielem(n,i)%ibounds(J).gt.0)then
                            if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
                                IDUM=1
                            end if
                        END IF
                        END DO
                end if

		
		ICONSIDERED=I
		
		
		
		
		
		INTBS=zero
		JXX=1;IXX=i;LXX1=1
		number_of_dog=ielem(n,i)%idegfree
		kxx=ielem(n,i)%iorder
		ELTYPE=ielem(n,i)%ishape
		compwrt=0
		INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1)
		
		
                INTEG_BASIS(I)%VALUE(1:ielem(n,i)%IDEGFREE)=INTBS(1:ielem(n,i)%IDEGFREE)
                
                
                 ICONSIDERED=I
		IF (IWENO.EQ.1)THEN
		CALL INDICATORMATRIX(N,I)
		END IF

		if (ees.eq.5)then
		INTBS=zero
		JXX=1;IXX=i;LXX1=1
		number_of_dog=idegfree2
		kxx=IORDER2
		ELTYPE=ielem(n,i)%ishape
		compwrt=1
		INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1)
		
                INTEG_BASIS(I)%VALUEc(1:number_of_dog)=INTBS(1:number_of_dog)
		IF (IWENO.EQ.1)THEN
		CALL INDICATORMATRIX2(N,I)
		compwrt=0
		END IF
		
		end if
		LLCO=IELEM(N,I)%ADMIS

                
IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	 INUMO=ielem(n,i)%iorder
         imax2=numneighbours2-1
	 inum2=numneighbours2
	 ideg2=iDEGFREE2
	 inumo2=iorder2
                DIST_STEN=ZERO
		DO LL=1,LLCO	!ADMIS

			
! !-------------------FOR DEBUGGING ONLY -----------------------------------------!

! !-------------------FOR DEBUGGING ONLY -----------------------------------------!
if((ees.ne.5).OR.(ll.eq.1))then
IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	 INUMO=ielem(n,i)%iorder
 compwrt=0

else
imax=numneighbours2-1
	 inum=numneighbours2
	 ideg=iDEGFREE2
	 inumo=IORDER2
 compwrt=1
end if




	    DO K=1,imax
	    ixx=i
	    kxx=INUMO

! 	    WLSQR(K)=ilocal_elem(1)%XXC(ll,k+1)
	    IF (WEIGHT_LSQR.EQ.1)THEN
	    WLSQR(LL,K)=1.0D0/((SQRT(((ilocal_elem(1)%XXC(ll,k+1)-ilocal_elem(1)%XXC(ll,1))**2)+((ilocal_elem(1)%YYC(ll,k+1)-ilocal_elem(1)%YYC(ll,1))**2))))
	    ELSE
	    WLSQR(LL,K)=1.0D0
	    END IF
	    X_STENCIL=(ilocal_elem(1)%XXC(ll,k+1)-ilocal_elem(1)%XXC(ll,1))**2
	    Y_STENCIL=(ilocal_elem(1)%YYC(ll,k+1)-ilocal_elem(1)%YYC(ll,1))**2


	    DIST_STEN2=SQRT(X_STENCIL+Y_STENCIL)

	    DIST_STEN=MAX(DIST_STEN,DIST_STEN2)
	    
	    
	    IF (WEIGHT_LSQR.EQ.1)THEN
	    WLSQR(ll,K)=1.0D0/SQRT(X_STENCIL+Y_STENCIL)
	    ELSE
	    WLSQR(ll,K)=1.0D0
	    END IF

	    IELEM(N,I)%STENCIL_DIST=DIST_STEN/(ilocal_elem(1)%VOLUME(1,1)**(1/2))

		      if (fastest.eq.1)then
			      x1 = ilocal_elem(1)%XXC(ll,k+1)-ilocal_elem(1)%XXC(ll,1)
			      y1 = ilocal_elem(1)%YYC(ll,k+1)-ilocal_elem(1)%YYC(ll,1)

					if((ees.ne.5).OR.(ll.eq.1))then
					compwrt=0
					ILOCAL_RECON3(I)%STENCILS(LL,K,1:ielem(n,i)%idegfree)=WLSQR(ll,K)*basis_rec2d(N,x1,y1,ielem(n,i)%iorder,IXX,ielem(n,i)%idegfree)
					ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
					else
					compwrt=1
					ILOCAL_RECON3(I)%STENCILSc(LL,K,1:ideg)=WLSQR(ll,K)*basis_rec2d(N,x1,y1,inumo,IXX,ideg)
					ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
					compwrt=0
					end if

		      else
			    IXX=i;jxx=k+1;lxx1=ll
			    ELTYPE=ilocal_elem(1)%ishape(ll,k+1)

			      IF (GREENGO.EQ.0)THEN

					  if (idum.eq.1)then
						  if((ees.ne.5).OR.(ll.eq.1))then
						  compwrt=0
						  ILOCAL_RECON3(I)%STENCILS(LL,K,1:ielem(n,i)%idegfree)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,ielem(n,i)%idegfree)
						  ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
						  else
						  compwrt=1
						  ILOCAL_RECON3(I)%STENCILSc(LL,K,1:ideg)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,ideg)
						  ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
						  compwrt=0
						  end if
					  else
						  if((ees.ne.5).OR.(ll.eq.1))then
						  compwrt=0
						  else
						  compwrt=1
						  end if
					  STENCILS(LL,K,1:ideg)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,ideg)
					  compwrt=0
					  end if
	    
			      else
				  if((ees.ne.5).OR.(ll.eq.1))then
				  compwrt=0
				  else
				  compwrt=1
				  end if
				  STENCILS(LL,K,1:ideg)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,ideg)
				  compwrt=0
			      end if
		      end if		!fastest
	    END DO	!IMAXEDEGFREE




     
      IF ((GREENGO.EQ.0))THEN
		if (idum.eq.1)then

		      invmat=zero
		      LSCQM(:,:)=zero
		      DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
		      if((ees.ne.5).OR.(ll.eq.1))then
		      LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
		      +((ILOCAL_RECON3(I)%STENCILS(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)))
		      else
		      LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
		      +((ILOCAL_RECON3(I)%STENCILSc(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILSc(LL,LQ,IQ)))
		      end if
		      END DO;END DO;END DO  
		else

		      
		      invmat=zero
		      LSCQM(:,:)=zero
		      DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
		      LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
		      +((STENCILS(LL,LQ,JQ)*STENCILS(LL,LQ,IQ)))
		      END DO;END DO;END DO 
		end if


      ELSE

		
		invmat=zero
		LSCQM(:,:)=zero
		DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
		LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
		+((STENCILS(LL,LQ,JQ)*STENCILS(LL,LQ,IQ)))
		END DO;END DO;END DO 

	END IF 








QFF(:,:)=zero; RFF(:,:)=zero; QTFF(:,:)=zero; RFF(:,:)=zero;  INVRFF(:,:)=zero 
CALL QRDECOMPOSITION(LSCQM,QFF,RFF,IDEG)
CALL TRANSPOSEMATRIX(QFF,QTFF,IDEG)
IVGT=IDEG+1
CALL INVERT(RFF,INVRFF,IVGT)
invmat(1:IDEg,1:ideg)=MATMUL(INVRFF(1:ideg,1:IDEG),QTFF(1:IDEG,1:IDEG))




! 
! 
! 

IF (GREENGO.EQ.0)THEN
    if (idum.eq.1)then
            if((ees.ne.5).OR.(ll.eq.1))then
            
            stencil(1:imax,1:ideg)=ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:ideg)
            else
           
            stencil(1:imax,1:ideg)=ILOCAL_RECON3(I)%STENCILSc(LL,1:imax,1:ideg)

            end if
    else
stencil(1:imax,1:ideg)=STENCILS(LL,1:imax,1:ideg)
    end if
ELSE

stencil(1:imax,1:ideg)=STENCILS(LL,1:imax,1:ideg)

END IF

if((ees.ne.5).OR.(ll.eq.1))then
! call gemm(                                               &
!                invmat,                                               &
!                stencil,                                              &
!                ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),             &
!                'N',                                                  & ! transposition flag for invmat
!                'T'                                                   & ! transposition flag for stencil
!             )
            
            
call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL),IDEG)

            
           do iq=1,imax
			ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)&
			*ilocal_elem(1)%VOLUME(ll,iq+1)*WLSQR(ll,iq)
			end do

else
! call gemm(                                               &
!                invmat(1:IDEG,1:IDEG),                                               &
!                stencil(1:imax,1:ideg),                                              &
!                ILOCAL_RECON3(I)%invmat_stenciltc(1:ideg,1:imax,LL),             &
!                'N',                                                  & ! transposition flag for invmat
!                'T'                                                   & ! transposition flag for stencil
!             )
            
call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stenciltc(1:IDEG,1:IMAX,LL),IDEG)

do iq=1,imax
			ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)&
			*ilocal_elem(1)%VOLUME(ll,iq+1)*WLSQR(ll,iq)
			end do


end if

	     if (initcond.eq.0)then
            maxai=zero
            minai=tolbig
            
            
            do ai=1,imax
	      do aj=1,ideg
		    if (abs(ILOCAL_RECON3(I)%invmat_stencilt(aj,ai,LL)).ne.zero)then
		    maxai=max(maxai,abs(ILOCAL_RECON3(I)%invmat_stencilt(aj,ai,LL)))
		    minai=min(minai,abs(ILOCAL_RECON3(I)%invmat_stencilt(aj,ai,LL)))
		    
		    end if
	      end do
	    end do
		ilocal_recon3(i)%cond(ll)=maxai/minai
            
            
            
            end if



if (ielem(n,i)%ggs.ne.1)then	!ggs
    
IF (ITESTCASE.EQ.4)THEN	!TEST

	

	 

IF (LL.EQ.1)THEN		!stencils
! 

  if (idum.eq.1)then		!for wall only
    
	 
	 

    	
    DO IHGT=1,2; DO IHGJ=1,2
    AINVJT(IHGT,IHGJ)=ILOCAL_RECON3(I)%INVCCJAC(IHGJ,IHGT)
   END DO; END DO

	 idum2=0
	 
	 
				BASEFACEVAL=zero
				BASEFACGVAL=zero
				PERMUTATION=zero
				PERMUTATIONg=zero
			
	 
	 
	 
	 DO j=1,IELEM(N,I)%IFCA		!for all faces
	    if (ielem(n,i)%ibounds(J).gt.0)then		!for bounded only
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then		!!for bounded only 2                         
					facex=J;iconsidered=i
								  CALL coordinates_face_inner2D(N,Iconsidered,facex)
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    AY=cords(2)
								    AX=cords(1)

					    VEXT(1,1)=AX;VEXT(1,2)=AY;
					    VEXT(1,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(1:2,1:2),VEXT(1,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  
					    AX=VEXT(1,1);AY=VEXT(1,2)
				



				
				
				
				
				ANGLE1=IELEM(N,I)%FACEANGLEX(j)
				ANGLE2=IELEM(N,I)%FACEANGLEY(j)
				
				NX=ANGLE1
				NY=ANGLE2
				
				NNX=(NX*AINVJT(1,1))+(NY*AINVJT(2,1))
				NNY=(NX*AINVJT(1,2))+(NY*AINVJT(2,2))
				
				IF (POLY.EQ.4)THEN


				DO IQ=1, IDEG

				XDER(IQ)=TL2dX(AX,AY,IQ);  YDER(IQ)=TL2dY(AX,AY,IQ);


				END DO


				ELSE



				DO IQ=1, IDEG
				
				XDER(IQ)=DF2dX(AX,AY,IQ);  YDER(IQ)=DF2dY(AX,AY,IQ);
				
				
				END DO

				END IF



				ICONSIDERED=I
				compwrt=0
				BASEFACEVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC2D(N,AX,AY,ielem(n,i)%Iorder,Iconsidered,ielem(n,i)%IDEGFREE)

				if (thermal.eq.0)then


				BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=((NNX*XDER(1:ielem(n,i)%IDEGFREE))+(NNY*YDER(1:ielem(n,i)%IDEGFREE)))
				ELSE
				BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC2D(N,AX,AY,ielem(n,i)%Iorder,Iconsidered,ielem(n,i)%IDEGFREE)

				end if


				DO IQ=1,IDEG
				ILOCAL_RECON3(I)%WALLCOEFF(IQ)=BASEFACEVAL(IQ)
				ILOCAL_RECON3(I)%WALLCOEFG(IQ)=BASEFACGVAL(IQ)
				PERMUTATION(IQ)=IQ
				PERMUTATIONG(IQ)=IQ
				END DO
		
				
		
		
		
		
		
		end if!for bounded only 2        
		end if!for bounded only
		end do!for all faces
		
		GGGG=-TOLBIG
G0=0
DO IQ=1,IDEG
IF (ABS(BASEFACGVAL(IQ)) >GGGG)THEN
GGGG=ABS(BASEFACGVAL(IQ))
G0=IQ
END IF
END DO




SSSS=-TOLBIG
K0=0
DO IQ=1,IDEG
IF (ABS(BASEFACEVAL(IQ)) >SSSS)THEN
SSSS=ABS(BASEFACEVAL(IQ))
K0=IQ
END IF
END DO



! now swap basis functions and thus coefficients
PERMUTATION(1)=K0; PERMUTATION(K0)=1; SSSS=BASEFACEVAL(1)
BASEFACEVAL(1)=BASEFACEVAL(K0); BASEFACEVAL(K0)=SSSS

PERMUTATIONG(1)=G0; PERMUTATION(G0)=1; GGGG=BASEFACGVAL(1)
BASEFACGVAL(1)=BASEFACGVAL(G0); BASEFACGVAL(G0)=GGGG
! ILOCAL_RECON3(II)%K0=K0
! ILOCAL_RECON3(II)%G0=G0
ILOCAL_RECON3(I)%K0=K0
ILOCAL_RECON3(I)%G0=G0
		
		
		
		
LSQM = ZERO
DO LQ=1,IMAX
LCOU=0
    DO IQ=1,IDEG
	IF (IQ.EQ.G0) CYCLE
	    LCOU=LCOU+1
LSQM(LQ,LCOU)=ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)&
-ILOCAL_RECON3(I)%STENCILS(LL,LQ,G0)*ILOCAL_RECON3(I)%WALLCOEFG(IQ)/ILOCAL_RECON3(I)%WALLCOEFG(G0)
    END DO
END DO
ILOCAL_RECON3(I)%TEMPSQ(1:IMAX,1:IDEG-1)=LSQM(1:IMAX,1:IDEG-1)






VELLSQMAT=ZERO
DO IQ=1,IDEG-1; DO JQ=1,IDEG-1;	DO LCC=1,IMAX
!now store the least square matrix
VELLSQMAT(JQ,IQ)= VELLSQMAT(JQ,IQ)+(LSQM(LCC,JQ)*LSQM(LCC,IQ))
END DO;	END DO;	END DO
LSCQM=ZERO
LSCQM(1:IDEG-1,1:IDEG-1)=VELLSQMAT(1:IDEG-1,1:IDEG-1)

QFF=ZERO; RFF=ZERO; QTFF=ZERO; INVRFF=ZERO
	  IVGT=IDEG
    CALL QRDECOMPOSITION(LSCQM,QFF,RFF,IVGT-1)
    CALL TRANSPOSEMATRIX(QFF,QTFF,IVGT-1)
    
    CALL INVERT(RFF,INVRFF,IVGT)
!   final inverted R^(-1)*Q^(-1)
ILOCAL_RECON3(I)%TEMPSQMAT(1:IDEG-1,1:IDEG-1) =&
 MATMUL(INVRFF(1:IDEG-1,1:IDEG-1),QTFF(1:IDEG-1,1:IDEG-1))

! definy the matrix defining the least-square reconstruction in this case for velocity
LSQM = ZERO
DO LQ=1,IMAX
LCOU=0
    DO IQ=1,IDEG
	IF (IQ.EQ.K0) CYCLE
	    LCOU=LCOU+1
	     LSQM(LQ,LCOU)=ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)&
 -ILOCAL_RECON3(I)%STENCILS(LL,LQ,K0)*ILOCAL_RECON3(I)%WALLCOEFF(IQ)/ILOCAL_RECON3(I)%WALLCOEFF(K0)
    END DO
END DO
ILOCAL_RECON3(I)%VELLSQ(1:IMAX,1:IDEG-1)=LSQM(1:IMAX,1:IDEG-1)
VELLSQMAT=ZERO




DO IQ=1,IDEG-1; DO JQ=1,IDEG-1;	DO LCC=1,IMAX
!now store the least square matrix
VELLSQMAT(JQ,IQ)= VELLSQMAT(JQ,IQ)+(LSQM(LCC,JQ)*LSQM(LCC,IQ))
END DO;	END DO;	END DO

LSCQM=ZERO
LSCQM(1:IDEG-1,1:IDEG-1)=VELLSQMAT(1:IDEG-1,1:IDEG-1)


QFF=ZERO; RFF=ZERO; QTFF=ZERO; INVRFF=ZERO
	  IVGT=IDEG
    CALL QRDECOMPOSITION(LSCQM,QFF,RFF,IVGT-1)
    CALL TRANSPOSEMATRIX(QFF,QTFF,IVGT-1)
    
    CALL INVERT(RFF,INVRFF,IVGT)
!   final inverted R^(-1)*Q^(-1)
ILOCAL_RECON3(I)%VELINVLSQMAT(1:IDEG-1,1:IDEG-1) = &
MATMUL(INVRFF(1:IDEG-1,1:IDEG-1),QTFF(1:IDEG-1,1:IDEG-1))
			
		
		
		
		
		
		
		end if!for wall only
		end if!stencils
		end if!for test
		end if!ggs
		end do
		
		


	!written by huang apr 29 to reconstruct high-order poly for DG
		SCALLING_HUANG = 1000.0d0
		IGQRULES_HUANG = IGQRULES
		QP_TRIANGLE_huang = QP_TRIANGLE
		QP_QUAD_HUANG = QP_QUAD
		QP_LINE_HUANG = QP_LINE
		IGQRULES = 7
		QP_TRIANGLE = 36
		QP_QUAD = 36
		QP_LINE = 9

		VEXT = zero
		DO IQ = 1,3
			NODES_LIST(IQ,1) = ILOCAL_NODE(1)%X(1, 1, IQ) * SCALLING_HUANG
			NODES_LIST(IQ,2) = ILOCAL_NODE(1)%Y(1, 1, IQ) * SCALLING_HUANG
			VEXT(IQ, 1:2) = NODES_LIST(IQ,1:2)
		end do
		integrals_basis(:) = integrals_basis_huang(N)
		VEXT = zero
		DO K = 1,3
			DO IQ = 1,3
				NODES_LIST(IQ,1) = ILOCAL_NODE(1)%X(K+1, 2, IQ) * SCALLING_HUANG
				NODES_LIST(IQ,2) = ILOCAL_NODE(1)%Y(K+1, 2, IQ) * SCALLING_HUANG
				VEXT(IQ, 1:2) = NODES_LIST(IQ,1:2)
	
			end do
			! OUTPUT: ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_HIHGORDER(1:IELEM(N,I)%IDEGFREE,((K-1) * 3 +1):(3*K),1)
			distancex = (ilocal_elem(1)%xxc(K+1,2) - ilocal_elem(1)%xxc(K+1,1)) * SCALLING_HUANG
			distancey = (ilocal_elem(1)%yyc(K+1,2) - ilocal_elem(1)%yyc(K+1,1)) * SCALLING_HUANG

			KXX = 3

			ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_HIHGORDER(((K-1) * 3 +1),1:9,1) = COMPBASTRI_HUANG(N,9,integrals_basis)
			ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_HIHGORDER(((K-1) * 3 +2),1:9,1) = COMPBASTRI_V1_HUANG(N,9,distancex,integrals_basis)
			ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_HIHGORDER(((K-1) * 3 +3),1:9,1) = COMPBASTRI_V2_HUANG(N,9,distancey,integrals_basis)

			ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_rightpoly(1:6,((K-1) * 2 +1),1) = COMPBASTRI_V1_HUANG_rightpoly(N,6,distancex,distancey)
			ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_rightpoly(1:6,((K-1) * 2 +2),1) = COMPBASTRI_V2_HUANG_rightpoly(N,6,distancex,distancey)

		end do

		loworder_matrix(:,:) = zero
		loworder_matrix(1,:) = ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_HIHGORDER(1,1:2,1)
		loworder_matrix(2,:) = ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_HIHGORDER(4,1:2,1)
		loworder_matrix(3,:) = ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_HIHGORDER(7,1:2,1)
		DO K = 1,3
			DO IQ = 1,2
				loworder_matrix_t(iq,k) = loworder_matrix(k,iq)
			END DO
		END DO

		loworder_matrix_inverse(:,:) = MATMUL(loworder_matrix_t(:,:),loworder_matrix(:,:)) 
		loworder_matrix_inverse(:,:) = loworder_matrix_inverse(:,:) / (SCALLING_HUANG **3)
		loworder_matrix_inverse(:,:)=inverse_loworder(loworder_matrix_inverse(:,:))
		loworder_matrix_final(:,:) = MATMUL(loworder_matrix_inverse(:,:),loworder_matrix_t(:,:)) 

		ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_lowORDER(:,:,1)=loworder_matrix_final(:,:)



		VEXT = zero
		integrals_basis(:) = zero
		DO IQ = 1,3
			NODES_LIST(IQ,1) = ILOCAL_NODE(1)%X(1, 1, IQ)
			NODES_LIST(IQ,2) = ILOCAL_NODE(1)%Y(1, 1, IQ) 
			VEXT(IQ, 1:2) = NODES_LIST(IQ,1:2)
		end do
		integrals_basis(:) = integrals_basis_huang(N)
		VEXT = zero		
		DO K = 1,3
			DO IQ = 1,3
				NODES_LIST(IQ,1) = ILOCAL_NODE(1)%X(K+1, 2, IQ) 
				NODES_LIST(IQ,2) = ILOCAL_NODE(1)%Y(K+1, 2, IQ)
				VEXT(IQ, 1:2) = NODES_LIST(IQ,1:2)
	
			end do
			! OUTPUT: ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_HIHGORDER(1:IELEM(N,I)%IDEGFREE,((K-1) * 3 +1):(3*K),1)
			distancex = (ilocal_elem(1)%xxc(K+1,2) - ilocal_elem(1)%xxc(K+1,1)) 
			distancey = (ilocal_elem(1)%yyc(K+1,2) - ilocal_elem(1)%yyc(K+1,1)) 

			KXX = 3

			testmatrix(((K-1) * 3 +1),1:9) = COMPBASTRI_HUANG(N,9,integrals_basis)
			testmatrix(((K-1) * 3 +2),1:9) = COMPBASTRI_V1_HUANG(N,9,distancex,integrals_basis)
			testmatrix(((K-1) * 3 +3),1:9) = COMPBASTRI_V2_HUANG(N,9,distancey,integrals_basis)
			
			! ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_rightpoly(1:6,((K-1) * 2 +1),1) = COMPBASTRI_V1_HUANG_rightpoly(N,6,distancex,distancey)
			! ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_rightpoly(1:6,((K-1) * 2 +2),1) = COMPBASTRI_V2_HUANG_rightpoly(N,6,distancex,distancey)

		end do


		ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_HIHGORDER(:,:,1) = COMPMASSINV_huang(N,I)

		integrals_basis(:) = zero
		VEXT = zero
		DO IQ = 1,3
			NODES_LIST(IQ,1) = ILOCAL_NODE(1)%X(1, 1, IQ) * SCALLING_HUANG
			NODES_LIST(IQ,2) = ILOCAL_NODE(1)%Y(1, 1, IQ) * SCALLING_HUANG
			VEXT(IQ, 1:2) = NODES_LIST(IQ,1:2)
		end do
		integrals_basis(:) = integrals_basis_huang(N)


		DO IQ = 1,3
			NODES_LIST(IQ,1) = ILOCAL_NODE(1)%X(1, 1, IQ) * SCALLING_HUANG
			NODES_LIST(IQ,2) = ILOCAL_NODE(1)%Y(1, 1, IQ) * SCALLING_HUANG
			VEXT(IQ, 1:2) = NODES_LIST(IQ,1:2)
		END DO

		ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_project(1:5,1:9,1) = COMPBASTRI_project_HUANG(N,9,integrals_basis)
		ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_projectinv(1:5,1:5,1) = ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_project(1:5,1:5,1)
		ILOCAL_RECON3(I)%invmat_stencilt_DG_HUANG_projectinv(1:5,1:5,1) = COMPMASSINV_huang1(N,I)
		IGQRULES = IGQRULES_HUANG
		QP_TRIANGLE = QP_TRIANGLE_huang
		QP_QUAD = QP_QUAD_HUANG
		QP_LINE = QP_LINE_HUANG


END SUBROUTINE PRESTORE_RECONSTRUCTION2


FUNCTION integrals_basis_huang(N)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL::VOL
	REAL,dimension(9)::INTEG,integrals_basis_huang
	integer::lc,iq

	REAL,dimension(2,qp_quad)::QPOINTS_huang
	REAL,dimension(qp_quad)::WEQUA3D_huang


	
	WEQUA3D_huang(:) = QUADRATURETRIANGLE_huang_weight(N,IGQRULES)
	QPOINTS_huang(:,:) = QUADRATURETRIANGLE_huang_POINT(N,IGQRULES)


	VOL=TRIANGLEVOLUME(N)

	INTEG=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(1)=INTEG(1)+((x1)*WEQUA3D_huang(Lc)*VOL)
		INTEG(2)=INTEG(2)+((Y1)*WEQUA3D_huang(Lc)*VOL)
		INTEG(3)=INTEG(3)+((x1)*X1*WEQUA3D_huang(Lc)*VOL)
		INTEG(4)=INTEG(4)+((x1)*Y1*WEQUA3D_huang(Lc)*VOL)
		INTEG(5)=INTEG(5)+((Y1)*Y1*WEQUA3D_huang(Lc)*VOL)
		INTEG(6)=INTEG(6)+((x1)*X1*X1*WEQUA3D_huang(Lc)*VOL)
		INTEG(7)=INTEG(7)+((x1)*X1*Y1*WEQUA3D_huang(Lc)*VOL)
		INTEG(8)=INTEG(8)+((x1)*Y1*Y1*WEQUA3D_huang(Lc)*VOL)
		INTEG(9)=INTEG(9)+((Y1)*Y1*Y1*WEQUA3D_huang(Lc)*VOL)
	END DO
	INTEG = INTEG / vol
	integrals_basis_huang = INTEG
	! test part
	INTEG=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(1)=INTEG(1)+((x1 - integrals_basis_huang(1))*WEQUA3D_huang(Lc)*VOL)
		INTEG(2)=INTEG(2)+((Y1 - integrals_basis_huang(2))*WEQUA3D_huang(Lc)*VOL)
		INTEG(3)=INTEG(3)+(((x1)*X1  - integrals_basis_huang(3))*WEQUA3D_huang(Lc)*VOL)
		INTEG(4)=INTEG(4)+(((x1)*Y1  - integrals_basis_huang(4))*WEQUA3D_huang(Lc)*VOL)
		INTEG(5)=INTEG(5)+(((Y1)*Y1- integrals_basis_huang(5))*WEQUA3D_huang(Lc)*VOL)
		INTEG(6)=INTEG(6)+(((x1)*X1*X1- integrals_basis_huang(6))*WEQUA3D_huang(Lc)*VOL)
		INTEG(7)=INTEG(7)+(((x1)*X1*Y1- integrals_basis_huang(7))*WEQUA3D_huang(Lc)*VOL)
		INTEG(8)=INTEG(8)+(((x1)*Y1*Y1- integrals_basis_huang(8))*WEQUA3D_huang(Lc)*VOL)
		INTEG(9)=INTEG(9)+(((Y1)*Y1*Y1- integrals_basis_huang(9))*WEQUA3D_huang(Lc)*VOL)
	END DO
	do lc = 1,9
		if(ABS(INTEG(LC)).GT.1.0E-10)then
			! write(*,*)"integrals_basis_huang = "
			! write(*,*)INTEG(lc)
		end if
	end do


END FUNCTION

FUNCTION COMPBASTRI_HUANG(N,number_of_dog,integrals_basis)
	!> @brief
	!> This subroutine computes basis functions for each triangle
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,number_of_dog
	REal,DIMENSION(9),INTENT(IN)::integrals_basis
	REAL::VOL
	REAL,dimension(number_of_dog)::INTEG,COMPBASTRI_HUANG
	integer::lc,iq
	REAL,dimension(2,qp_quad)::QPOINTS_huang
	REAL,dimension(qp_quad)::WEQUA3D_huang


	
	WEQUA3D_huang(:) = QUADRATURETRIANGLE_huang_weight(N,IGQRULES)
	QPOINTS_huang(:,:) = QUADRATURETRIANGLE_huang_POINT(N,IGQRULES)
	
	VOL=TRIANGLEVOLUME(N)




	INTEG=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(1)=INTEG(1)+((X1 - integrals_basis(1))*WEQUA3D_huang(Lc)*VOL)
		INTEG(2)=INTEG(2)+((Y1 - integrals_basis(2))*WEQUA3D_huang(Lc)*VOL)
		INTEG(3)=INTEG(3)+((X1*X1 - integrals_basis(3))*WEQUA3D_huang(Lc)*VOL)

		INTEG(4)=INTEG(4)+((X1*Y1 - integrals_basis(4))*WEQUA3D_huang(Lc)*VOL)
		INTEG(5)=INTEG(5)+((Y1*Y1 - integrals_basis(5))*WEQUA3D_huang(Lc)*VOL)
		INTEG(6)=INTEG(6)+((X1*X1*X1 - integrals_basis(6))*WEQUA3D_huang(Lc)*VOL)

		INTEG(7)=INTEG(7)+((X1*X1*Y1 - integrals_basis(7))*WEQUA3D_huang(Lc)*VOL)
		INTEG(8)=INTEG(8)+((X1*Y1*Y1 - integrals_basis(8))*WEQUA3D_huang(Lc)*VOL)
		INTEG(9)=INTEG(9)+((Y1*Y1*Y1 - integrals_basis(9))*WEQUA3D_huang(Lc)*VOL)
	END DO



	COMPBASTRI_HUANG= INTEG

END FUNCTION

FUNCTION COMPBASTRI_V1_HUANG(N,number_of_dog,distancex,integrals_basis)
	!> @brief
	!> This subroutine computes basis functions for each triangle
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,number_of_dog
	REal,INTENT(IN)::distancex
	REal,DIMENSION(9),INTENT(IN)::integrals_basis
	REAL::VOL,INTEG1
	REAL,dimension(number_of_dog)::INTEG,COMPBASTRI_V1_HUANG
	integer::lc,iq
	REAL,dimension(2,qp_quad)::QPOINTS_huang
	REAL,dimension(qp_quad)::WEQUA3D_huang


	
	WEQUA3D_huang(:) = QUADRATURETRIANGLE_huang_weight(N,IGQRULES)
	QPOINTS_huang(:,:) = QUADRATURETRIANGLE_huang_POINT(N,IGQRULES)


	VOL=TRIANGLEVOLUME(N)

	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		
		INTEG1=INTEG1+((x1 - distancex)*WEQUA3D_huang(Lc)*VOL)
	
	END DO
	INTEG1 = INTEG1 / vol


	INTEG=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(1)=INTEG(1)+((X1 - integrals_basis(1))*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG1)))
		INTEG(2)=INTEG(2)+((Y1 - integrals_basis(2))*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG1)))
		INTEG(3)=INTEG(3)+((X1*X1 - integrals_basis(3))*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG1)))

		INTEG(4)=INTEG(4)+((X1*Y1 - integrals_basis(4))*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG1)))
		INTEG(5)=INTEG(5)+((Y1*Y1 - integrals_basis(5))*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG1)))
		INTEG(6)=INTEG(6)+((X1*X1*X1 - integrals_basis(6))*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG1)))

		INTEG(7)=INTEG(7)+((X1*X1*Y1 - integrals_basis(7))*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG1)))
		INTEG(8)=INTEG(8)+((X1*Y1*Y1 - integrals_basis(8))*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG1)))
		INTEG(9)=INTEG(9)+((Y1*Y1*Y1 - integrals_basis(9))*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG1)))
	END DO



		COMPBASTRI_V1_HUANG= INTEG

END FUNCTION


FUNCTION COMPBASTRI_V2_HUANG(N,number_of_dog,distancey,integrals_basis)
	!> @brief
	!> This subroutine computes basis functions for each triangle
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,number_of_dog
	REal,INTENT(IN)::distancey
	REal,DIMENSION(9),INTENT(IN)::integrals_basis
	REAL::VOL,INTEG1
	REAL,dimension(number_of_dog)::INTEG,COMPBASTRI_V2_HUANG
	integer::lc
	REAL,dimension(2,qp_quad)::QPOINTS_huang
	REAL,dimension(qp_quad)::WEQUA3D_huang

	WEQUA3D_huang(:) = QUADRATURETRIANGLE_huang_weight(N,IGQRULES)
	QPOINTS_huang(:,:) = QUADRATURETRIANGLE_huang_POINT(N,IGQRULES)

	VOL=TRIANGLEVOLUME(N)


	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG1=INTEG1+((Y1 - distanceY)*WEQUA3D_huang(Lc)*VOL)
	END DO
	INTEG1 = INTEG1 / vol


	INTEG=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(1)=INTEG(1)+((X1 - integrals_basis(1))*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG1)))
		INTEG(2)=INTEG(2)+((Y1 - integrals_basis(2))*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG1)))
		INTEG(3)=INTEG(3)+((X1*X1 - integrals_basis(3))*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG1)))
		INTEG(4)=INTEG(4)+((X1*Y1 - integrals_basis(4))*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG1)))
		INTEG(5)=INTEG(5)+((Y1*Y1 - integrals_basis(5))*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG1)))
		INTEG(6)=INTEG(6)+((X1*X1*X1 - integrals_basis(6))*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG1)))
		INTEG(7)=INTEG(7)+((X1*X1*Y1 - integrals_basis(7))*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG1)))
		INTEG(8)=INTEG(8)+((X1*Y1*Y1 - integrals_basis(8))*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG1)))
		INTEG(9)=INTEG(9)+((Y1*Y1*Y1 - integrals_basis(9))*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG1)))
	END DO


	COMPBASTRI_V2_HUANG= INTEG

END FUNCTION

FUNCTION COMPBASTRI_V1_HUANG_rightpoly(N,number_of_dog,distancex,distancey)
	!> @brief
	!> This subroutine computes basis functions for each triangle
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,number_of_dog
	REal,INTENT(IN)::distancex,distancey
	REAL::VOL,INTEG1
	REAL::INTEG_HUANG

	REAL,dimension(number_of_dog)::INTEG,COMPBASTRI_V1_HUANG_rightpoly
	integer::lc,IGQRULES_HUANG,QP_TRIANGLE_huang
	REAL,dimension(2,qp_quad)::QPOINTS_huang
	REAL,dimension(qp_quad)::WEQUA3D_huang

	WEQUA3D_huang(:) = QUADRATURETRIANGLE_huang_weight(N,IGQRULES)
	QPOINTS_huang(:,:) = QUADRATURETRIANGLE_huang_POINT(N,IGQRULES)

	
	VOL=TRIANGLEVOLUME(N)
	INTEG=zero
	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((x1 - distancex)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG_HUANG = INTEG1
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(1)=INTEG(1)+(WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG_HUANG)))
	END DO
	if(abs(INTEG(1)).gt.1.0E-03)then
		! write(*,*)"INTEG(1)="
		! write(*,*)INTEG(1)
	end if

	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((x1 - distancex)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(2)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(2)=INTEG(2)+((x1 - distancex - INTEG1)*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG_HUANG)))
	END DO


	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((y1 - distancey)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(3)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(3)=INTEG(3)+((y1 - distancey - INTEG1)*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG_HUANG)))
	END DO
	

	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((X1 - distancex)*(X1 - distancex)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(4)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(4)=INTEG(4)+(((X1 - distancex)*(X1 - distancex) - INTEG1)*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG_HUANG)))
	END DO

	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((X1 - distancex)*(Y1 - distancey)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(5)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(5)=INTEG(5)+(((X1 - distancex)*(Y1 - distancey) - INTEG1)*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG_HUANG)))
	END DO

	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((Y1 - distancey)*(Y1 - distancey)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(6)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(6)=INTEG(6)+(((Y1 - distancey)*(Y1 - distancey) - INTEG1)*WEQUA3D_huang(Lc)*VOL*((x1 - distancex - INTEG_HUANG)))
	END DO
	
	
	COMPBASTRI_V1_HUANG_rightpoly= INTEG
	
END FUNCTION
	
	
FUNCTION COMPBASTRI_V2_HUANG_rightpoly(N,number_of_dog,distancex,distancey)
	!> @brief
	!> This subroutine computes basis functions for each triangle
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,number_of_dog
	REal,INTENT(IN)::distancex,distancey
	REAL::VOL,INTEG1,INTEG_HUANG
	REAL,dimension(number_of_dog)::INTEG,COMPBASTRI_V2_HUANG_rightpoly
	integer::lc
	REAL,dimension(2,qp_quad)::QPOINTS_huang
	REAL,dimension(qp_quad)::WEQUA3D_huang

	WEQUA3D_huang(:) = QUADRATURETRIANGLE_huang_weight(N,IGQRULES)
	QPOINTS_huang(:,:) = QUADRATURETRIANGLE_huang_POINT(N,IGQRULES)

	
	VOL=TRIANGLEVOLUME(N)
	INTEG=zero
	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((Y1 - distancey)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG_HUANG = INTEG1
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(1)=INTEG(1)+(WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG_HUANG)))
	END DO

	if(abs(INTEG(1)).gt.1.0E-03)then
		! write(*,*)"INTEG(1)="
		! write(*,*)INTEG(1)
	end if
	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((x1 - distancex)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(2)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(2)=INTEG(2)+((x1 - distancex - INTEG1)*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG_HUANG)))
	END DO


	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((y1 - distancey)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(3)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(3)=INTEG(3)+((y1 - distancey - INTEG1)*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG_HUANG)))
	END DO
	

	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((x1 - distancex)*(x1 - distancex)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(4)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(4)=INTEG(4)+(((x1 - distancex)*(x1 - distancex) - INTEG1)*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG_HUANG)))
	END DO

	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((x1 - distancex)*(Y1 - distancey)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(5)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(5)=INTEG(5)+(((x1 - distancex)*(Y1 - distancey) - INTEG1)*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG_HUANG)))
	END DO

	!
	INTEG1=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc) 
		INTEG1=INTEG1+((Y1 - distancey)*(Y1 - distancey)*WEQUA3D_huang(Lc)*VOL)
	 
	END DO
	INTEG1 = INTEG1 / vol
	INTEG(6)=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(6)=INTEG(6)+(((Y1 - distancey)*(Y1 - distancey) - INTEG1)*WEQUA3D_huang(Lc)*VOL*((y1 - distancey - INTEG_HUANG)))
	END DO
	COMPBASTRI_V2_HUANG_rightpoly = INTEG
END FUNCTION

FUNCTION COMPBASTRI_project_HUANG(N,number_of_dog,integrals_basis)
	!> @brief
	!> This subroutine computes basis functions for each triangle
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,number_of_dog
	REal,DIMENSION(9),INTENT(IN)::integrals_basis
	REAL::VOL,INTEG1
	REAL,dimension(5,9)::INTEG,COMPBASTRI_project_HUANG
	REAL,dimension(2,qp_quad)::QPOINTS_huang
	REAL,dimension(qp_quad)::WEQUA3D_huang
	integer::lc

	WEQUA3D_huang(:) = QUADRATURETRIANGLE_huang_weight(N,IGQRULES)
	QPOINTS_huang(:,:) = QUADRATURETRIANGLE_huang_POINT(N,IGQRULES)
	
	VOL=TRIANGLEVOLUME(N)




	INTEG=zero
	DO Lc=1,qp_triangle
		x1=QPOINTS_huang(1,Lc);y1=QPOINTS_huang(2,Lc)
		INTEG(1,1)=INTEG(1,1)+((X1 - integrals_basis(1))*WEQUA3D_huang(Lc)*VOL*(X1 - integrals_basis(1)))
		INTEG(1,2)=INTEG(1,2)+((Y1 - integrals_basis(2))*WEQUA3D_huang(Lc)*VOL*(X1 - integrals_basis(1)))
		
		INTEG(1,3)=INTEG(1,3)+((X1*X1 - integrals_basis(3))*WEQUA3D_huang(Lc)*VOL*(X1 - integrals_basis(1)))
		INTEG(1,4)=INTEG(1,4)+((X1*Y1 - integrals_basis(4))*WEQUA3D_huang(Lc)*VOL*(X1 - integrals_basis(1)))
		INTEG(1,5)=INTEG(1,5)+((Y1*Y1 - integrals_basis(5))*WEQUA3D_huang(Lc)*VOL*(X1 - integrals_basis(1)))
		
		INTEG(1,6)=INTEG(1,6)+((X1*X1*X1 - integrals_basis(6))*WEQUA3D_huang(Lc)*VOL*(X1 - integrals_basis(1)))
		INTEG(1,7)=INTEG(1,7)+((X1*X1*Y1 - integrals_basis(7))*WEQUA3D_huang(Lc)*VOL*(X1 - integrals_basis(1)))
		INTEG(1,8)=INTEG(1,8)+((X1*Y1*Y1 - integrals_basis(8))*WEQUA3D_huang(Lc)*VOL*(X1 - integrals_basis(1)))
		INTEG(1,9)=INTEG(1,9)+((Y1*Y1*Y1 - integrals_basis(9))*WEQUA3D_huang(Lc)*VOL*(X1 - integrals_basis(1)))

		
		INTEG(2,1)=INTEG(2,1)+((X1 - integrals_basis(1))*WEQUA3D_huang(Lc)*VOL*(Y1 - integrals_basis(2)))
		INTEG(2,2)=INTEG(2,2)+((Y1 - integrals_basis(2))*WEQUA3D_huang(Lc)*VOL*(Y1 - integrals_basis(2)))
		
		INTEG(2,3)=INTEG(2,3)+((X1*X1 - integrals_basis(3))*WEQUA3D_huang(Lc)*VOL*(Y1 - integrals_basis(2)))
		INTEG(2,4)=INTEG(2,4)+((X1*Y1 - integrals_basis(4))*WEQUA3D_huang(Lc)*VOL*(Y1 - integrals_basis(2)))
		INTEG(2,5)=INTEG(2,5)+((Y1*Y1 - integrals_basis(5))*WEQUA3D_huang(Lc)*VOL*(Y1 - integrals_basis(2)))
		
		INTEG(2,6)=INTEG(2,6)+((X1*X1*X1 - integrals_basis(6))*WEQUA3D_huang(Lc)*VOL*(Y1 - integrals_basis(2)))
		INTEG(2,7)=INTEG(2,7)+((X1*X1*Y1 - integrals_basis(7))*WEQUA3D_huang(Lc)*VOL*(Y1 - integrals_basis(2)))
		INTEG(2,8)=INTEG(2,8)+((X1*Y1*Y1 - integrals_basis(8))*WEQUA3D_huang(Lc)*VOL*(Y1 - integrals_basis(2)))
		INTEG(2,9)=INTEG(2,9)+((Y1*Y1*Y1 - integrals_basis(9))*WEQUA3D_huang(Lc)*VOL*(Y1 - integrals_basis(2)))

		
		INTEG(3,1)=INTEG(3,1)+((X1 - integrals_basis(1))*WEQUA3D_huang(Lc)*VOL*(X1*X1 - integrals_basis(3)))
		INTEG(3,2)=INTEG(3,2)+((Y1 - integrals_basis(2))*WEQUA3D_huang(Lc)*VOL*(X1*X1 - integrals_basis(3)))
		
		INTEG(3,3)=INTEG(3,3)+((X1*X1 - integrals_basis(3))*WEQUA3D_huang(Lc)*VOL*(X1*X1 - integrals_basis(3)))
		INTEG(3,4)=INTEG(3,4)+((X1*Y1 - integrals_basis(4))*WEQUA3D_huang(Lc)*VOL*(X1*X1 - integrals_basis(3)))
		INTEG(3,5)=INTEG(3,5)+((Y1*Y1 - integrals_basis(5))*WEQUA3D_huang(Lc)*VOL*(X1*X1 - integrals_basis(3)))
		
		INTEG(3,6)=INTEG(3,6)+((X1*X1*X1 - integrals_basis(6))*WEQUA3D_huang(Lc)*VOL*(X1*X1 - integrals_basis(3)))
		INTEG(3,7)=INTEG(3,7)+((X1*X1*Y1 - integrals_basis(7))*WEQUA3D_huang(Lc)*VOL*(X1*X1 - integrals_basis(3)))
		INTEG(3,8)=INTEG(3,8)+((X1*Y1*Y1 - integrals_basis(8))*WEQUA3D_huang(Lc)*VOL*(X1*X1 - integrals_basis(3)))
		INTEG(3,9)=INTEG(3,9)+((Y1*Y1*Y1 - integrals_basis(9))*WEQUA3D_huang(Lc)*VOL*(X1*X1 - integrals_basis(3)))		


		INTEG(4,1)=INTEG(4,1)+((X1 - integrals_basis(1))*WEQUA3D_huang(Lc)*VOL*(X1*Y1 - integrals_basis(4)))
		INTEG(4,2)=INTEG(4,2)+((Y1 - integrals_basis(2))*WEQUA3D_huang(Lc)*VOL*(X1*Y1 - integrals_basis(4)))

		INTEG(4,3)=INTEG(4,3)+((X1*X1 - integrals_basis(3))*WEQUA3D_huang(Lc)*VOL*(X1*Y1 - integrals_basis(4)))
		INTEG(4,4)=INTEG(4,4)+((X1*Y1 - integrals_basis(4))*WEQUA3D_huang(Lc)*VOL*(X1*Y1 - integrals_basis(4)))
		INTEG(4,5)=INTEG(4,5)+((Y1*Y1 - integrals_basis(5))*WEQUA3D_huang(Lc)*VOL*(X1*Y1 - integrals_basis(4)))

		INTEG(4,6)=INTEG(4,6)+((X1*X1*X1 - integrals_basis(6))*WEQUA3D_huang(Lc)*VOL*(X1*Y1 - integrals_basis(4)))
		INTEG(4,7)=INTEG(4,7)+((X1*X1*Y1 - integrals_basis(7))*WEQUA3D_huang(Lc)*VOL*(X1*Y1 - integrals_basis(4)))
		INTEG(4,8)=INTEG(4,8)+((X1*Y1*Y1 - integrals_basis(8))*WEQUA3D_huang(Lc)*VOL*(X1*Y1 - integrals_basis(4)))
		INTEG(4,9)=INTEG(4,9)+((Y1*Y1*Y1 - integrals_basis(9))*WEQUA3D_huang(Lc)*VOL*(X1*Y1 - integrals_basis(4)))		


		INTEG(5,1)=INTEG(5,1)+((X1 - integrals_basis(1))*WEQUA3D_huang(Lc)*VOL*(Y1*Y1 - integrals_basis(5)))
		INTEG(5,2)=INTEG(5,2)+((Y1 - integrals_basis(2))*WEQUA3D_huang(Lc)*VOL*(Y1*Y1 - integrals_basis(5)))

		INTEG(5,3)=INTEG(5,3)+((X1*X1 - integrals_basis(3))*WEQUA3D_huang(Lc)*VOL*(Y1*Y1 - integrals_basis(5)))
		INTEG(5,4)=INTEG(5,4)+((X1*Y1 - integrals_basis(4))*WEQUA3D_huang(Lc)*VOL*(Y1*Y1 - integrals_basis(5)))
		INTEG(5,5)=INTEG(5,5)+((Y1*Y1 - integrals_basis(5))*WEQUA3D_huang(Lc)*VOL*(Y1*Y1 - integrals_basis(5)))

		INTEG(5,6)=INTEG(5,6)+((X1*X1*X1 - integrals_basis(6))*WEQUA3D_huang(Lc)*VOL*(Y1*Y1 - integrals_basis(5)))
		INTEG(5,7)=INTEG(5,7)+((X1*X1*Y1 - integrals_basis(7))*WEQUA3D_huang(Lc)*VOL*(Y1*Y1 - integrals_basis(5)))
		INTEG(5,8)=INTEG(5,8)+((X1*Y1*Y1 - integrals_basis(8))*WEQUA3D_huang(Lc)*VOL*(Y1*Y1 - integrals_basis(5)))
		INTEG(5,9)=INTEG(5,9)+((Y1*Y1*Y1 - integrals_basis(9))*WEQUA3D_huang(Lc)*VOL*(Y1*Y1 - integrals_basis(5)))		


	END DO



	COMPBASTRI_project_HUANG= INTEG

END FUNCTION


FUNCTION COMPMASSINV_huang(N,ICONSIDERED)
	!Calculate the inverse of the input matrix with Gauss-Jordan Elimination
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,ICONSIDERED
	integer :: I,j,k,l,m,irow,num_dofs
	real:: big,dum
	real,DIMENSION(9,9)::a,b,c,d
	real,DIMENSION(9,9)::COMPMASSINV_huang
	num_dofs=9
	


	do i  = 1,9
		do j = 1,9
			a(i,j) = ILOCAL_RECON3(ICONSIDERED)%invmat_stencilt_DG_HUANG_HIHGORDER(i,j,1)
		end do 
	end do
	! WRITE(*,*)"A="
	! WRITE(*,*)a
	d=a
	b(:,:)=zero
	
	
	
	
	do i = 1,num_dofs
		do j = 1,num_dofs
			b(i,j) = 0.0d0
		end do
		b(i,i) = 1.0d0
	end do 
	
	do i = 1,num_dofs 
	   big = a(i,i)
	   do j = i,num_dofs
		 if (a(j,i).gt.big) then
		   big = a(j,i)
		   irow = j
		 end if
	   end do
	   ! interchange lines i with irow for both a() and b() matrices
	   if (big.gt.a(i,i)) then
		 do k = 1,num_dofs
		   dum = a(i,k)                      ! matrix a()
		   a(i,k) = a(irow,k)
		   a(irow,k) = dum
		   dum = b(i,k)                 ! matrix b()
		   b(i,k) = b(irow,k)
		   b(irow,k) = dum
		 end do
	   end if
	   ! divide all entries in line i from a(i,j) by the value a(i,i); 
	   ! same operation for the identity matrix
	   dum = a(i,i)
	   do j = 1,num_dofs
		 a(i,j) = a(i,j)/dum
		 b(i,j) = b(i,j)/dum
	   end do
	   ! make zero all entries in the column a(j,i); same operation for indent()
	   do j = i+1,num_dofs
		 dum = a(j,i)
		 do k = 1,num_dofs
		   a(j,k) = a(j,k) - dum*a(i,k)
		   b(j,k) = b(j,k) - dum*b(i,k)               
				
		 end do
	   end do
	end do
	  
	 do i = 1,num_dofs-1
	   do j = i+1,num_dofs
		 dum = a(i,j)
		 do l = 1,num_dofs
		   a(i,l) = a(i,l)-dum*a(j,l)
		   b(i,l) = b(i,l)-dum*b(j,l)
		 end do
	   end do
	 end do




	!  write(*,*)"b1="
	!  write(*,*)b(:,:)
	!  c=MATMUL(d(:,:),b(:,:))
	!  write(*,*)"c1="
	!  write(*,*)c(:,:)
	 
	 COMPMASSINV_huang(:,:)=b(:,:)
	 
	 
END FUNCTION COMPMASSINV_huang

FUNCTION COMPMASSINV_huang1(N,ICONSIDERED)
	!Calculate the inverse of the input matrix with Gauss-Jordan Elimination
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,ICONSIDERED
	integer :: I,j,k,l,m,irow,num_dofs
	real:: big,dum
	real,DIMENSION(5,5)::a,b,c,d
	real,DIMENSION(5,5)::COMPMASSINV_huang1
	num_dofs=5
	


	do i  = 1,5
		do j = 1,5
			a(i,j) = ILOCAL_RECON3(ICONSIDERED)%invmat_stencilt_DG_HUANG_project(i,j,1)
		end do 
	end do
	! WRITE(*,*)"A="
	! WRITE(*,*)a
	d=a
	b(:,:)=zero
	
	
	
	
	do i = 1,num_dofs
		do j = 1,num_dofs
			b(i,j) = 0.0d0
		end do
		b(i,i) = 1.0d0
	end do 
	
	do i = 1,num_dofs 
	   big = a(i,i)
	   do j = i,num_dofs
		 if (a(j,i).gt.big) then
		   big = a(j,i)
		   irow = j
		 end if
	   end do
	   ! interchange lines i with irow for both a() and b() matrices
	   if (big.gt.a(i,i)) then
		 do k = 1,num_dofs
		   dum = a(i,k)                      ! matrix a()
		   a(i,k) = a(irow,k)
		   a(irow,k) = dum
		   dum = b(i,k)                 ! matrix b()
		   b(i,k) = b(irow,k)
		   b(irow,k) = dum
		 end do
	   end if
	   ! divide all entries in line i from a(i,j) by the value a(i,i); 
	   ! same operation for the identity matrix
	   dum = a(i,i)
	   do j = 1,num_dofs
		 a(i,j) = a(i,j)/dum
		 b(i,j) = b(i,j)/dum
	   end do
	   ! make zero all entries in the column a(j,i); same operation for indent()
	   do j = i+1,num_dofs
		 dum = a(j,i)
		 do k = 1,num_dofs
		   a(j,k) = a(j,k) - dum*a(i,k)
		   b(j,k) = b(j,k) - dum*b(i,k)               
				
		 end do
	   end do
	end do
	  
	 do i = 1,num_dofs-1
	   do j = i+1,num_dofs
		 dum = a(i,j)
		 do l = 1,num_dofs
		   a(i,l) = a(i,l)-dum*a(j,l)
		   b(i,l) = b(i,l)-dum*b(j,l)
		 end do
	   end do
	 end do




	!  write(*,*)"b1="
	!  write(*,*)b(:,:)
	!  c=MATMUL(d(:,:),b(:,:))
	!  write(*,*)"c1="
	!  write(*,*)c(:,:)
	 
	 COMPMASSINV_huang1(:,:)=b(:,:)
	 
	 
END FUNCTION COMPMASSINV_huang1

FUNCTION inverse_loworder(matrix)
	!Calculate the inverse of the input matrix with Gauss-Jordan Elimination
	IMPLICIT NONE
	real,DIMENSION(2,2),INTENT(IN)::matrix
	integer :: I,j,k,l,m,irow,num_dofs
	real:: big,dum
	real,DIMENSION(2,2)::a,b,c,d
	real,DIMENSION(2,2)::inverse_loworder
	num_dofs=2
	


	do i  = 1,2
		do j = 1,2
			a(i,j) = matrix(i,j)
		end do 
	end do
	! WRITE(*,*)"A="
	! WRITE(*,*)a
	d=a
	b(:,:)=zero
	
	
	
	
	do i = 1,num_dofs
		do j = 1,num_dofs
			b(i,j) = 0.0d0
		end do
		b(i,i) = 1.0d0
	end do 
	
	do i = 1,num_dofs 
	   big = a(i,i)
	   do j = i,num_dofs
		 if (a(j,i).gt.big) then
		   big = a(j,i)
		   irow = j
		 end if
	   end do
	   ! interchange lines i with irow for both a() and b() matrices
	   if (big.gt.a(i,i)) then
		 do k = 1,num_dofs
		   dum = a(i,k)                      ! matrix a()
		   a(i,k) = a(irow,k)
		   a(irow,k) = dum
		   dum = b(i,k)                 ! matrix b()
		   b(i,k) = b(irow,k)
		   b(irow,k) = dum
		 end do
	   end if
	   ! divide all entries in line i from a(i,j) by the value a(i,i); 
	   ! same operation for the identity matrix
	   dum = a(i,i)
	   do j = 1,num_dofs
		 a(i,j) = a(i,j)/dum
		 b(i,j) = b(i,j)/dum
	   end do
	   ! make zero all entries in the column a(j,i); same operation for indent()
	   do j = i+1,num_dofs
		 dum = a(j,i)
		 do k = 1,num_dofs
		   a(j,k) = a(j,k) - dum*a(i,k)
		   b(j,k) = b(j,k) - dum*b(i,k)               
				
		 end do
	   end do
	end do
	  
	 do i = 1,num_dofs-1
	   do j = i+1,num_dofs
		 dum = a(i,j)
		 do l = 1,num_dofs
		   a(i,l) = a(i,l)-dum*a(j,l)
		   b(i,l) = b(i,l)-dum*b(j,l)
		 end do
	   end do
	 end do




	!  write(*,*)"b1="
	!  write(*,*)b(:,:)
	!  c=MATMUL(d(:,:),b(:,:))
	!  write(*,*)"c1="
	!  write(*,*)c(:,:)
	 
	 inverse_loworder(:,:)=b(:,:)
	 
	 
END FUNCTION inverse_loworder


FUNCTION matrix_inverse(n1,a)
	!> @brief
	!> This subroutine computes basis functions for each triangle
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N1
	REal,DIMENSION(n1,n1),INTENT(IN)::A
    real :: matrix(n1, n1), inverse(n1, n1),matrix_inverse(n1, n1)
    integer :: i, j, k
    real :: pivot, factor
	matrix = a


    inverse = 0.0
    do i = 1, n1
        inverse(i, i) = 1.0
    end do


    do k = 1, n1
        pivot = matrix(k, k)
        matrix(k, :) = matrix(k, :) / pivot
        inverse(k, :) = inverse(k, :) / pivot
        do i = 1, n1
            if (i /= k) then
                factor = matrix(i, k)
                matrix(i, :) = matrix(i, :) - factor * matrix(k, :)
                inverse(i, :) = inverse(i, :) - factor * inverse(k, :)
            end if
        end do
    end do

	matrix_inverse = inverse

	

END FUNCTION




SUBROUTINE INDICATORMATRIX(N,Iconsi)
!> @brief
!> This subroutine computes the indicator matrices for weno reconstructions
IMPLICIT NONE
INTEGER,INTENT(IN)::N,iconsi
INTEGER::I,J,K,L,M,jx,jx2
i=iconsi

	IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	 INUMO=ielem(n,i)%iorder
	 
	iconsidered=i
	 

	VEXT=ZERO
    NODES_LIST=ZERO
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
    ELEM_LISTD=ZERO
          ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)=ZERO
      jx=IELEM(N,I)%NONODES
     
	  if (dimensiona.eq.3)then
	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILOCAL_NODE(1)%X(1,1,k)
	    NODES_LIST(k,2)=ILOCAL_NODE(1)%y(1,1,k)
	    NODES_LIST(k,3)=ILOCAL_NODE(1)%z(1,1,k)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE3!(N,ELTYPE,ELEM_DEC)
	  
	  else

	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILOCAL_NODE(1)%X(1,1,k)
	    NODES_LIST(k,2)=ILOCAL_NODE(1)%y(1,1,k)
	    VEXT(K,1:2)=NODES_LIST(k,1:2)
	  END DO
	  CALL DECOMPOSE2
	  	  end if

	 SELECT CASE(ielem(n,i)%ishape)

      CASE(1)
      CALL QUADRATUREHEXA(N,IGQRULES)
      voltemp=HEXAVOLUME(N)
      inump=qp_hexa
      CASE(2)
      CALL QUADRATURETETRA(N,IGQRULES)
      voltemp=tetraVOLUME(N)
      inump=qp_tetra
      CASE(3)
      CALL QUADRATUREPYRA(N,IGQRULES)
      voltemp=pyraVOLUME(N)
      inump=qp_pyra
      CASE(4)
      CALL QUADRATUREPRISM(N,IGQRULES)
      voltemp=prismVOLUME(N)
      inump=qp_prism

      CASE(5)
      CALL QUADRATUREquad(N,IGQRULES)
      voltemp=quadVOLUME(N)
      inump=qp_quad


      CASE(6)
      CALL QUADRATUREtriangle(N,IGQRULES)
      voltemp=TRIANGLEVOLUME(N)
      inump=qp_triangle

      END SELECT      
      
	if (dimensiona.eq.3)then
	IF (ielem(n,i)%mode.EQ.1)THEN
	  do K=1,ELEM_DEC
	      VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
	      CALL QUADRATURETETRA(N,IGQRULES)
	      voltemp=tetraVOLUME(N)
	      inump=qp_tetra
	      CALL WENOTET(N,WEFF)
	      ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)
	
	  END DO
	ELSE
	  

	    CALL WENOTET(N,WEFF)
	      ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)






	END IF
	else

	    do K=1,ELEM_DEC
            VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
            voltemp=TRIANGLEVOLUME(N)
            inump=qp_triangle
           
	    CALL WENOTET2d(N,WEFF)
             
	      ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)
            END DO




	end if




END SUBROUTINE INDICATORMATRIX


SUBROUTINE INDICATORMATRIX2(N,Iconsi)
!> @brief
!> This subroutine computes the indicator matrices for cweno reconstructions of the lower-order polynomials
IMPLICIT NONE
INTEGER,INTENT(IN)::N,iconsi
INTEGER::I,J,K,L,M,jx,jx2
i=iconsi

	IMAX=numneighbours2-1
	INUM=numneighbours2
	IDEG=iDEGFREE2
	 INUMO=IORDER2
	 
	
	 iconsidered=i

	VEXT=ZERO
    NODES_LIST=ZERO
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
    ELEM_LISTD=ZERO
          ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)=ZERO
      jx=IELEM(N,I)%NONODES
     
	  if (dimensiona.eq.3)then
	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILOCAL_NODE(1)%X(1,1,k)
	    NODES_LIST(k,2)=ILOCAL_NODE(1)%y(1,1,k)
	    NODES_LIST(k,3)=ILOCAL_NODE(1)%z(1,1,k)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE3!(N,ELTYPE,ELEM_DEC)
	  
	  else

	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILOCAL_NODE(1)%X(1,1,k)
	    NODES_LIST(k,2)=ILOCAL_NODE(1)%y(1,1,k)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE2
	  	  end if

	 SELECT CASE(ielem(n,i)%ishape)

      CASE(1)
      CALL QUADRATUREHEXA(N,IGQRULES)
      voltemp=HEXAVOLUME(N)
      inump=qp_hexa
      CASE(2)
      CALL QUADRATURETETRA(N,IGQRULES)
      voltemp=tetraVOLUME(N)
      inump=qp_tetra
      CASE(3)
      CALL QUADRATUREPYRA(N,IGQRULES)
      voltemp=pyraVOLUME(N)
      inump=qp_pyra
      CASE(4)
      CALL QUADRATUREPRISM(N,IGQRULES)
      voltemp=prismVOLUME(N)
      inump=qp_prism

      CASE(5)
      CALL QUADRATUREquad(N,IGQRULES)
      voltemp=quadVOLUME(N)
      inump=qp_quad


      CASE(6)
      CALL QUADRATUREtriangle(N,IGQRULES)
      voltemp=TRIANGLEVOLUME(N)
      inump=qp_triangle

      END SELECT      
      
	if (dimensiona.eq.3)then
	IF (ielem(n,i)%mode.EQ.1)THEN
	  do K=1,ELEM_DEC
	      VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
	      CALL QUADRATURETETRA(N,IGQRULES)
	      voltemp=tetraVOLUME(N)
	      inump=qp_tetra
	      CALL WENOTET(N,WEFF)
	      ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)
	
	  END DO
	ELSE
	  

	    CALL WENOTET(N,WEFF)
	      ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)






	END IF
	else
            
            do K=1,ELEM_DEC
            VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
            voltemp=TRIANGLEVOLUME(N)
            inump=qp_triangle

	    CALL WENOTET2d(N,WEFF)

	      ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)
            END DO
            
	end if

	
	



END SUBROUTINE INDICATORMATRIX2

	



SUBROUTINE WENOTET(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
WEFF=zero


		IF (POLY.EQ.1)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL PH1(N,WEFF)

		CASE(2)
		CALL PH2(N,WEFF)

		CASE(3)
		CALL PH3(N,WEFF)

		CASE(4)
		CALL PH4(N,WEFF)

		CASE(5)
		CALL PH5(N,WEFF)

		CASE(6)
		CALL PH6(N,WEFF)
		
	END SELECT
	END IF
	IF (POLY.EQ.2)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL PL1(N,WEFF)

		CASE(2)
		CALL PL2(N,WEFF)

		CASE(3)
		CALL PL3(N,WEFF)

		CASE(4)
		CALL PL4(N,WEFF)

		CASE(5)
		CALL PL5(N,WEFF)

		CASE(6)
		CALL PL6(N,WEFF)
		
	END SELECT
	END IF
	
	IF (POLY.EQ.4)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL TL3D1(N,WEFF)

		CASE(2)
		CALL TL3D2(N,WEFF)

		CASE(3)
		CALL TL3D3(N,WEFF)

		CASE(4)
		CALL TL3D4(N,WEFF)

		CASE(5)
		CALL TL3D5(N,WEFF)

		CASE(6)
		CALL TL3D6(N,WEFF)
		
	END SELECT
	END IF
	
	
	

END SUBROUTINE WENOTET

SUBROUTINE PH1(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=ZERO
		
 	DO I=1,IDEG
        	DO J=1,IDEG
			INTEG =ZERO
                scalerx=1.0d0
                
                
       			 DO K=1,inump
	
	    PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

	    INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH1
SUBROUTINE PH2(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=ZERO
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =ZERO
  DO K=1,iNUMp
	PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 
	
	INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
  ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH2
SUBROUTINE PH3(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
  DO K=1,inump
	  PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
	DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 
	

	INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
  ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH3
SUBROUTINE PH4(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
DO K=1,inump
    PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
    DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

    INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH4
SUBROUTINE PH5(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
DO I=1,IDEG
DO J=1,IDEG
scalerx=1.0d0
                
	INTEG =zero
	  DO K=1,inump
PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
  DFZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

		INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
	  ENDDO
	WEFF(I,J)=WEFF(I,J)+INTEG
ENDDO
ENDDO

END SUBROUTINE PH5
SUBROUTINE PH6(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
			PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		    DFX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		    DFX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DFZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DFZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH6

SUBROUTINE PL1(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
	
				PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				   DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				   DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL1
SUBROUTINE PL2(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 
				
				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL2
SUBROUTINE PL3(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				 PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 
				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL3
SUBROUTINE PL4(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL4
SUBROUTINE PL5(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
		PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		  DLX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL5
SUBROUTINE PL6(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
			     PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		    DLX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		    DLX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  DLZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*DLZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL6


SUBROUTINE TL3D1(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
	
				PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				   TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				   TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D1
SUBROUTINE TL3D2(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
				PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 
				
				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D2
SUBROUTINE TL3D3(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
				 PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 
				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D3
SUBROUTINE TL3D4(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
				PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
				TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D4
SUBROUTINE TL3D5(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
		PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		  TL3DX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D5
SUBROUTINE TL3D6(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
			     PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		    TL3DX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
		    TL3DX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) + &
                  TL3DZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I)*TL3DZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J) 

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D6





SUBROUTINE WENOTET2D(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
WEFF=zero

 
		IF (POLY.EQ.1)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL P2DH1(N,WEFF)

		CASE(2)
		CALL P2DH2(N,WEFF)

		CASE(3)
		CALL P2DH3(N,WEFF)

		CASE(4)
		CALL P2DH4(N,WEFF)

		CASE(5)
		CALL P2DH5(N,WEFF)

		CASE(6)
		CALL P2DH6(N,WEFF)
		
	END SELECT
	
	END IF
	
	
	
	IF (POLY.EQ.4)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL TL2DH1(N,WEFF)

		CASE(2)
		CALL TL2DH2(N,WEFF)

		CASE(3)
		CALL TL2DH3(N,WEFF)

		CASE(4)
		CALL TL2DH4(N,WEFF)

		CASE(5)
		CALL TL2DH5(N,WEFF)

		CASE(6)
		CALL TL2DH6(N,WEFF)
		
	END SELECT
	
	END IF
	
	
	

END SUBROUTINE WENOTET2D
! 
SUBROUTINE P2DH1(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
	
				PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
				   DF2DY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J) 
				   

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
!           		
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH1
SUBROUTINE P2DH2(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
			
			INTEG =zero
       			 DO K=1,inump
				PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
				(DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J)) * voltemp 	!huang
				
				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH2
SUBROUTINE P2DH3(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
        	
        	
        	
        	
			INTEG =zero
       			 DO K=1,inump
				 PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
				(DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J)) * voltemp + &
				(DF2DX3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DY3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY3(QPOINTS(1,K),QPOINTS(2,K),J))*voltemp * voltemp!huang
				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH3
SUBROUTINE P2DH4(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DX3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DY3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DX4(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX4(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
				DF2DY4(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY4(QPOINTS(1,K),QPOINTS(2,K),J) 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH4
SUBROUTINE P2DH5(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
		PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX4(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY4(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY4(QPOINTS(1,K),QPOINTS(2,K),J) + &
		  DF2DX5(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX5(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX4Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX4Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DXY4(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY5(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY5(QPOINTS(1,K),QPOINTS(2,K),J) 
                 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH5
SUBROUTINE P2DH6(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
			PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX4(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY4(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY4(QPOINTS(1,K),QPOINTS(2,K),J) + &
		    DF2DX5(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX5(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX4Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX4Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DXY4(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY5(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY5(QPOINTS(1,K),QPOINTS(2,K),J) + &
		    DF2DX6(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX6(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX5Y(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX5Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX4Y2(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX4Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX3Y3(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX3Y3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DX2Y4(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DX2Y4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DXY5(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DXY5(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  DF2DY6(QPOINTS(1,K),QPOINTS(2,K),I)*DF2DY6(QPOINTS(1,K),QPOINTS(2,K),J) 

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH6

SUBROUTINE TL2DH1(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
	
				PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
				   TL2DY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J) 
				   

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
!           		
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH1
SUBROUTINE TL2DH2(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
				PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J) 
				
				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH2
SUBROUTINE TL2DH3(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
        	
        	
        	
        	
			INTEG =zero
       			 DO K=1,inump
				 PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DY3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY3(QPOINTS(1,K),QPOINTS(2,K),J) 
				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH3
SUBROUTINE TL2DH4(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
				PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DY3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX4(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX4(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
				TL2DY4(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY4(QPOINTS(1,K),QPOINTS(2,K),J) 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH4
SUBROUTINE TL2DH5(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
		PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX4(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY4(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY4(QPOINTS(1,K),QPOINTS(2,K),J) + &
		  TL2DX5(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX5(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX4Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX4Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DXY4(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY5(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY5(QPOINTS(1,K),QPOINTS(2,K),J) 
                 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH5
SUBROUTINE TL2DH6(N,WEFF)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=ielem(n,iconsidered)%totvolume
                
			INTEG =zero
       			 DO K=1,inump
			PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX4(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY4(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY4(QPOINTS(1,K),QPOINTS(2,K),J) + &
		    TL2DX5(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX5(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX4Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX4Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DXY4(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY5(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY5(QPOINTS(1,K),QPOINTS(2,K),J) + &
		    TL2DX6(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX6(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX5Y(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX5Y(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX4Y2(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX4Y2(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX3Y3(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX3Y3(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DX2Y4(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DX2Y4(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DXY5(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DXY5(QPOINTS(1,K),QPOINTS(2,K),J) + &
                  TL2DY6(QPOINTS(1,K),QPOINTS(2,K),I)*TL2DY6(QPOINTS(1,K),QPOINTS(2,K),J) 

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH6



FUNCTION CALINTBASIS(N,IXX,JXX,KXX,LXX1)
!> @brief
!> This subroutine computes basis functions for each element
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,INTENT(INOUT):: IXX,JXX,KXX,LXX1
INTEGER::K
REAL,DIMENSION(1:number_of_dog)::CALINTBASIS


! 	kxx=ielem(n,ixx)%iorder
	CALINTBASIS = ZERO
     	
	CALINTBASIS(1:number_of_dog)=COMPBASEL(N,eltype,number_of_dog)
	
   END FUNCTION


FUNCTION COMPBASEL(N,ELTYPE,number_of_dog)
!> @brief
!> This subroutine computes basis functions for each element
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,INTENT(INout)::ELTYPE,number_of_dog
INTEGER::JX,K
real,dimension(1:number_of_dog)::compbasel,s1
S1=ZERO

 SELECT CASE(ELTYPE)
      CASE(1)
      JX=8;   ELEM_DEC=6
		ELTYPE=1
		 if ((JXX.eq.1))then
		  
		      do K=1,jx
			  NODES_LIST(k,1)=ILOCAL_NODE(1)%X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILOCAL_NODE(1)%Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILOCAL_NODE(1)%Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
		      END DO
		  if (ielem(n,ixx)%mode.eq.0)then
		  
		  COMPBASEL=COMPBASHEX(N,number_of_dog)

		    else
		  CALL DECOMPOSE3
		  do K=1,ELEM_DEC
! 		    VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,number_of_dog)
		  end do
		    COMPBASEL=S1
		    end if
		  
		  ELSE

		   do K=1,jx
			  NODES_LIST(k,1)=ILOCAL_NODE(1)%X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILOCAL_NODE(1)%Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILOCAL_NODE(1)%Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
			
		      END DO
		      
		      
		  CALL DECOMPOSE3
		  do K=1,ELEM_DEC
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,number_of_dog)
		    
		  end do
		    COMPBASEL=S1
		 END IF

      CASE(2)	
		  ELTYPE=2
		JX=4;   ELEM_DEC=1
		do K=1,jx
			  NODES_LIST(k,1)=ILOCAL_NODE(1)%X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILOCAL_NODE(1)%Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILOCAL_NODE(1)%Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
		      END DO
		  CALL DECOMPOSE3
		  do K=1,ELEM_DEC
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,number_of_dog)
		  end do
		    COMPBASEL=S1
		
      CASE(3)
		ELTYPE=3
		JX=5;   ELEM_DEC=2
		do K=1,jx
			  NODES_LIST(k,1)=ILOCAL_NODE(1)%X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILOCAL_NODE(1)%Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILOCAL_NODE(1)%Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
		      END DO
		  CALL DECOMPOSE3
		  do K=1,ELEM_DEC
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,number_of_dog)
		  end do
		    COMPBASEL=S1
	CASE(4)
		JX=6;   ELEM_DEC=3
	      ELTYPE=4
	      if ((JXX.eq.1))then
		  
		      do K=1,jx
			  NODES_LIST(k,1)=ILOCAL_NODE(1)%X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILOCAL_NODE(1)%Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILOCAL_NODE(1)%Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
		      END DO
		      CALL DECOMPOSE3
		  if (ielem(n,ixx)%mode.eq.0)then
		  
		  COMPBASEL=COMPBASPR(N,number_of_dog)

		    else
		  
		  do K=1,ELEM_DEC
		   VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,number_of_dog)
! 		   
		  end do
		    COMPBASEL=S1
		    end if
		  
		  ELSE
			    
		   do K=1,jx
			  NODES_LIST(k,1)=ILOCAL_NODE(1)%X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILOCAL_NODE(1)%Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILOCAL_NODE(1)%Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
			
		      END DO
		  CALL DECOMPOSE3
		  
		  do K=1,ELEM_DEC
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		   
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,number_of_dog)
		    
		  end do
		    COMPBASEL=S1
		 END IF

	      CASE(5)
		JX=4;   ELTYPE=5;ELEM_DEC=2
		     if ((JXX.eq.1))then
		     
		     
		   do K=1,JX
			  NODES_LIST(k,1)=ILOCAL_NODE(1)%X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILOCAL_NODE(1)%Y(LXX1,JXX,K)
		    vext(k,1:2)=NODES_LIST(k,1:2)
		    
		    END DO
		    CALL DECOMPOSE2
		    if (ielem(n,ixx)%mode.eq.0)then
		  COMPBASEL=COMPBASquad(N,number_of_dog)
		    else
		    
		  do K=1,ELEM_DEC
		    VEXT(1:3,1)=ELEM_LISTD(k,1:3,1)
		    VEXT(1:3,2)=ELEM_LISTD(k,1:3,2)
		    
		    
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTRi(N,number_of_dog)
		    
		  end do
		    COMPBASEL=S1
		    
		    
		    end if
		    
		    
		    
		else
		  do K=1,jx
			  NODES_LIST(k,1)=ILOCAL_NODE(1)%X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILOCAL_NODE(1)%Y(LXX1,JXX,K)
			 
			VEXT(K,:)=NODES_LIST(k,:)
			
		      END DO
		      
		      
		  CALL DECOMPOSE2
		  do K=1,ELEM_DEC
		    VEXT(1:3,1)=ELEM_LISTD(k,1:3,1)
		    VEXT(1:3,2)=ELEM_LISTD(k,1:3,2)
		    
		    
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTRi(N,number_of_dog)
		    
		  end do
		    COMPBASEL=S1	
		
		end if

	      CASE(6)

		JX=3;   ELTYPE=6
		
		
		   do K=1,JX
			  NODES_LIST(k,1)=ILOCAL_NODE(1)%X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILOCAL_NODE(1)%Y(LXX1,JXX,K)
		    vext(k,1:2)=NODES_LIST(k,1:2)
		    
		    END DO
		    
		  COMPBASEL=COMPBAStri(N,number_of_dog)


     END SELECT

END FUNCTION

FUNCTION COMPBASTR(N,number_of_dog)
!> @brief
!> This subroutine computes basis functions for each triangle
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog
integer::lc
REAL::VOL
REAL,dimension(1:number_of_dog)::INTEG,COMPBASTR


QPOINTS(1:3,1:QP_TETRA)=ZERO
WEQUA3D(1:QP_TETRA)=ZERO

VOL=TETRAVOLUME(N)

CALL QUADRATURETETRA(N,IGQRULES)

INTEG=ZERO
DO Lc=1,QP_TETRA
	  x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc);z1=QPOINTS(3,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC(N,x1,y1,z1,kxx,IXX,number_of_dog)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASTR = INTEG

END FUNCTION


FUNCTION COMPBASHEX(N,number_of_dog)
!> @brief
!> This subroutine computes basis functions for hexahedral
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog
REAL::VOL
REAL,dimension(number_of_dog)::INTEG,COMPBASHEX
integer::lc
QPOINTS(1:3,1:QP_HEXA)=zero
WEQUA3D(1:QP_HEXA)=zero

CALL QUADRATUREHEXA(N,IGQRULES)


 VOL=hexaVOLUME(N)

INTEG=zero
DO Lc=1,QP_HEXA
	x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc);z1=QPOINTS(3,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC(N,x1,y1,z1,kxx,IXX,number_of_dog)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASHEX = INTEG

END FUNCTION

FUNCTION COMPBASPR(N,number_of_dog)
!> @brief
!> This subroutine computes basis functions for prisms
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog
REAL::VOL
REAL,dimension(number_of_dog)::INTEG,COMPBASPR
integer::lc
QPOINTS(1:3,1:QP_PRISM)=zero
WEQUA3D(1:QP_PRISM)=zero

CALL QUADRATUREPRISM(N,IGQRULES)


 VOL=PRISMVOLUME(N)

INTEG=zero
DO Lc=1,QP_PRISM
	x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc);z1=QPOINTS(3,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC(N,x1,y1,z1,kxx,IXX,number_of_dog)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASPR= INTEG

END FUNCTION

FUNCTION COMPBASQUAD(N,number_of_dog)
!> @brief
!> This subroutine computes basis functions for quadrilateral
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog
REAL::VOL
REAL,dimension(number_of_dog)::INTEG,COMPBASQUAD
integer::lc
QPOINTS(1:2,1:QP_QUAD)=zero
WEQUA3D(1:QP_QUAD)=zero

CALL QUADRATUREQUAD(N,IGQRULES)


 VOL=QUADVOLUME(N)

INTEG=zero
DO Lc=1,qp_quad
	x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC2D(N,x1,y1,kxx,IXX,number_of_dog)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASQUAD= INTEG

END FUNCTION

FUNCTION COMPBASTRI(N,number_of_dog)
!> @brief
!> This subroutine computes basis functions for each triangle
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog
REAL::VOL
REAL,dimension(number_of_dog)::INTEG,COMPBAStri
integer::lc
QPOINTS(1:2,1:QP_QUAD)=zero
WEQUA3D(1:QP_QUAD)=zero

CALL QUADRATURETRIANGLE(N,IGQRULES)


 VOL=TRIANGLEVOLUME(N)

INTEG=zero
DO Lc=1,qp_triangle
	x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC2D(N,x1,y1,kxx,IXX,number_of_dog)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASTRI= INTEG

END FUNCTION


SUBROUTINE INVERT(RFF,INVRFF,IVGT)
!> @brief
!> This subroutine inverts a special matrix
  IMPLICIT NONE
  INTEGER,INTENT(IN)::IVGT
  REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(IN) ::RFF
  REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(inOUT)::INVRFF
  INTEGER I,J,K,GT
 
  INVRff(1:IVGT-1,1:IVGT-1) = zero
 GT=IVGT-1
  DO I=GT,1,-1
    INVRFF(i,i) = 1./RFF(i,i)
	do j=i+1,GT
	  INVRFF(i,j) = zero
	  do k= 1,j-1
	   INVRFF(i,j) = INVRFF(i,j) - RFF(k,j)*INVRFF(i,k)
	  enddo
 	  INVRFF(i,j) =INVRFF(i,j) /RFF(j,j)
	enddo
  enddo

 End subroutine INVERT
 


 
 
END MODULE PRESTORE
