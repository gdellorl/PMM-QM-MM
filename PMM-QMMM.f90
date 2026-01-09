program PMMQMMM
!gfortran PMM-QMMM.f90 liblapack.a librefblas.a libxdrf.a -o PMM-QMMM.x
        implicit none 

        ! Integer scalars
        integer :: ninte, ninte1, nninte
        integer :: n, i, j, k, g, m, ind
        integer :: nat, tot
        integer :: id1, ret, frame, nframe, step
        integer :: nsolv, nsolu, nrad
        integer :: nroot, fframe, lframe
        integer :: counter
        integer :: jarr        ! Final electronic state
        integer :: ipart       ! Initial electronic state
        integer :: ctype
        integer :: nch
        
        ! Integer arrays
        integer, allocatable :: a(:), asolu(:), num(:)
        integer, allocatable :: delta(:,:)
        
        ! Double precision scalars
        double precision :: ddelta
        double precision :: time, prec
        double precision :: box(9)
        double precision :: xcdm, ycdm, zcdm
        double precision :: xcdmqm,  ycdmqm,  zcdmqm
        double precision :: xcdmqm1, ycdmqm1, zcdmqm1
        double precision :: Enerx, Enery, Enerz
        double precision :: distx, disty, distz, dist
        double precision :: mu0x, mu0y, mu0z
        double precision :: Potential
        double precision :: QC_charge, QC_mass
        double precision :: the, phi, psi
        double precision :: mtot, mtot1
        double precision :: ExQMMM, EyQMMM, EzQMMM, PotQMMM
        double precision :: ExQMMM_MD, EyQMMM_MD, EzQMMM_MD
        double precision :: xcdmQC,  ycdmQC,  zcdmQC
        double precision :: xcdmQC1, ycdmQC1, zcdmQC1
        
        ! Double precision vectors and matrices
        double precision, allocatable :: xint(:), yint(:), zint(:)
        double precision, allocatable :: xintcm(:), yintcm(:), zintcm(:)
        double precision, allocatable :: ch(:), x(:)
        double precision, allocatable :: mass(:), massqm(:)
        double precision, allocatable :: mx(:,:), my(:,:), mz(:,:)
        double precision, allocatable :: mxp(:,:), myp(:,:), mzp(:,:)
        double precision, allocatable :: E0(:)
        double precision, allocatable :: ExMux(:,:), EyMuy(:,:), EzMuz(:,:)
        double precision, allocatable :: Hdiag(:,:), Hprim(:,:)
        double precision, allocatable :: Eivec(:,:), Eprim(:)
        double precision, allocatable :: Tvec(:,:)
        double precision, allocatable :: d0x(:), d0y(:), d0z(:)
        double precision, allocatable :: r0_1(:), r0_2(:), rfin(:)
        double precision, allocatable :: dist1(:)
        double precision, allocatable :: xQMMM(:), yQMMM(:), zQMMM(:), chQMMM(:)
        double precision, allocatable :: xQC(:), yQC(:), zQC(:), mQC(:)
        double precision, allocatable :: xQC1(:), yQC1(:), zQC1(:), mQC1(:)
        
        ! Rotation matrices and vectors
        double precision :: AA(3), BB(3), CC(3)
        
        ! Character variables
        character(len=5)  :: sdum
        character(len=60) :: rif_QC, solvent, solute, topol
        character(len=60) :: traj, dipoles, output, proj
        character(len=60) :: file_ch

        ! For the diagonalization subroutine
        character :: jobz, uplo
        integer   :: info, lda, lwork, lwmax
        double precision  :: work(1000)
        double precision, allocatable :: w(:)

        ! Read input parameters
        ! Open the PMM input file
        open(unit=2, file="input-pmm-QMMM")
        read(2,*)
        read(2,*) file_ch                   
        read(2,*) ! Point-charges file name 
        close(2)
        ! Open the point-charges file and read the charges
        open(unit=2, file=file_ch)         
        read(2,*) nch
        if (nch==0) then
                nch=1
                allocate(chQMMM(nch))
                allocate(xQMMM(nch))
                allocate(yQMMM(nch))
                allocate(zQMMM(nch))
                chQMMM(nch)=0.
                xQMMM(nch)=0.
                yQMMM(nch)=0.
                zQMMM(nch)=0.
        else
                allocate(chQMMM(nch))
                allocate(xQMMM(nch))
                allocate(yQMMM(nch))
                allocate(zQMMM(nch))
                do i=1,nch
                        read(2,*) chQMMM(i), xQMMM(i), yQMMM(i), zQMMM(i)
                enddo
        endif
        close(2)
        
        ! Open the QC file
        open(unit=2, file="input-pmm-QMMM")
        read(2,*)
        read(2,*)
        read(2,*) ! QC file name from the QM/MM calculation
        read(2,*) rif_QC
        close(2)
        open(unit=2, file=rif_QC)
        ! This refers to the whole QC system
        read(2,*) ninte
        allocate(mQC(ninte))
        allocate(xQC(ninte))
        allocate(yQC(ninte))
        allocate(zQC(ninte))
        xcdmQC=0
        ycdmQC=0
        zcdmQC=0
        mtot=0
        do i=1,ninte
                read(2,*) mQC(i), xQC(i), yQC(i), zQC(i)
                mtot=mtot+mQC(i)
                xcdmQC=xcdmQC+mQC(i)*xQC(i)
                ycdmQC=ycdmQC+mQC(i)*yQC(i)
                zcdmQC=zcdmQC+mQC(i)*zQC(i)
        enddo
        ! Center of mass of the whole QC system
        xcdmQC=xcdmQC/mtot
        ycdmQC=ycdmQC/mtot
        zcdmQC=zcdmQC/mtot
        ! This refers to the QC subgroup
        read(2,*) ninte1
        allocate(mQC1(ninte1))
        allocate(xQC1(ninte1))
        allocate(yQC1(ninte1))
        allocate(zQC1(ninte1))
        xcdmQC1=0
        ycdmQC1=0
        zcdmQC1=0
        mtot1=0
        do i=1,ninte1
                read(2,*) mQC1(i), xQC1(i), yQC1(i), zQC1(i)
                mtot1=mtot1+mQC1(i)
                xcdmQC1=xcdmQC1+mQC1(i)*xQC1(i)
                ycdmQC1=ycdmQC1+mQC1(i)*yQC1(i)
                zcdmQC1=zcdmQC1+mQC1(i)*zQC1(i)
        enddo
        ! Center of mass of the QC subgroup
        xcdmQC1=xcdmQC1/mtot1
        ycdmQC1=ycdmQC1/mtot1
        zcdmQC1=zcdmQC1/mtot1
        ! Redefine the atoms of the whole QC system with respect
        ! to the center of mass of the QC subgroup
        do i=1,ninte
                xQC(i)=xQC(i)-xcdmQC1
                yQC(i)=yQC(i)-ycdmQC1
                zQC(i)=zQC(i)-zcdmQC1
        enddo
        
        ! Redefine the charges with respect to the center of mass
        do i=1,nch
                xQMMM(i)=xQMMM(i)-xcdmQC1
                yQMMM(i)=yQMMM(i)-ycdmQC1
                zQMMM(i)=zQMMM(i)-zcdmQC1
        enddo
        
        ! Initialize the reference field
        Enerx=0.0
        Enery=0.0
        Enerz=0.0
        Potential=0.0
        
        ! Convert the center-of-mass coordinates of the QC subgroup to a.u.
        xcdmQC1=xcdmQC1 
        ycdmQC1=ycdmQC1 
        zcdmQC1=zcdmQC1 
        do m=1,nch
                ! Compute the distances of the charges from the center of mass
                ! of the QC subgroup
                distx=((xQMMM(m)/0.5291772109))**2 ! Convert the distance from Å to Bohr
                disty=((yQMMM(m)/0.5291772109))**2 ! Convert the distance from Å to Bohr
                distz=((zQMMM(m)/0.5291772109))**2 ! Convert the distance from Å to Bohr
                dist=sqrt(distx+disty+distz)
                ! Compute the electric field
                Enerx=Enerx+(chQMMM(m)*(0.-(xQMMM(m)/0.5291772109))/(dist)**3) ! Convert the distance from Å to Bohr
                Enery=Enery+(chQMMM(m)*(0.-(yQMMM(m)/0.5291772109))/(dist)**3) ! Convert the distance from Å to Bohr
                Enerz=Enerz+(chQMMM(m)*(0.-(zQMMM(m)/0.5291772109))/(dist)**3) ! Convert the distance from Å to Bohr
                ! Compute the electric potential
                Potential=Potential+chQMMM(m)/dist
        enddo
        
        ExQMMM=Enerx
        EyQMMM=Enery
        EzQMMM=Enerz
        PotQMMM=Potential
        close(2)
        open(unit=2, file="input-pmm-QMMM")
        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*) ! MD trajectory file name in XTC format
        read(2,*) traj
        read(2,*) ! First frame
        read(2,*) fframe
        read(2,*) ! Last frame
        read(2,*) lframe
        read(2,*) ! Number of atoms
        read(2,*) nat
        read(2,*) ! MD simulation charges file name
        read(2,*) topol
        read(2,*) ! QC index file name with masses
        read(2,*) solute
        read(2,*) ! Environment index file name
        read(2,*) solvent
        read(2,*) ! Charge of the QC
        read(2,*) QC_charge
        read(2,*) ! Number of states (including the ground state)
        read(2,*) nroot
        nrad=nroot
        
        ! Perform some operations
        tot=nat*3
        nframe=lframe-fframe+1 
        read(2,*) ! Number of QC atoms used for dipole calculations
        read(2,*) nninte
        
        allocate(xint(ninte))
        allocate(yint(ninte))
        allocate(zint(ninte))    
        allocate(xintcm(ninte))
        allocate(yintcm(ninte))
        allocate(zintcm(ninte))
        allocate(massqm(ninte))
        allocate(r0_1(3*ninte1))
        
        do n=1,ninte
                xint(n)=xQC(n)
                yint(n)=yQC(n)
                zint(n)=zQC(n)
                massqm(n)=mQC(n)
        end do
        ! Compute the QC center of mass
        QC_mass=0.0
        do i=1, ninte
                QC_mass=QC_mass+massqm(i)
        end do
        ! Center of mass of the entire QC
        xcdmqm=0.0                            
        ycdmqm=0.0
        zcdmqm=0.0
        do i=1,ninte
                xcdmqm=xcdmqm+xint(i)*massqm(i)
                ycdmqm=ycdmqm+yint(i)*massqm(i)
                zcdmqm=zcdmqm+zint(i)*massqm(i)
        end do
        xcdmqm=xcdmqm/QC_mass
        ycdmqm=ycdmqm/QC_mass
        zcdmqm=zcdmqm/QC_mass
        allocate(mx(nroot,nroot))
        allocate(my(nroot,nroot))
        allocate(mz(nroot,nroot))

        ! Now consider the QC subgroup
        ! Compute its center of mass
        QC_mass=0.0
        do i=1, ninte1
                QC_mass=QC_mass+mQC1(i)       
        end do
        xcdmqm1=0.0
        ycdmqm1=0.0
        zcdmqm1=0.0
        do i=1,ninte1
                xcdmqm1=xcdmqm1+xQC1(i)*mQC1(i)     
                ycdmqm1=ycdmqm1+yQC1(i)*mQC1(i)    
                zcdmqm1=zcdmqm1+zQC1(i)*mQC1(i)    
        end do
        xcdmqm1=xcdmqm1/QC_mass
        ycdmqm1=ycdmqm1/QC_mass
        zcdmqm1=zcdmqm1/QC_mass
        ! Reference structure used to fit all MD configurations
        j=1
        do i=1,ninte1
                r0_1(j)=xQC1(i)-xcdmQC1
                r0_1(j+1)=yQC1(i)-ycdmQC1
                r0_1(j+2)=zQC1(i)-zcdmQC1
                j=j+3
        enddo
        
        read(2,*) ! Dipole matrix
        do i=1,nroot
                do j=1,nroot
                        read(2,*) sdum, sdum, mx(i,j), my(i,j), mz(i,j)
                        if (j==i) then
                                mx(i,j)=mx(i,j)
                                my(i,j)=my(i,j)
                                mz(i,j)=mz(i,j)
                                ! Translate the dipole to the center of mass
                                ! of the QC subgroup
                                mx(i,j)=mx(i,j)+QC_charge*(xcdmQC1-xcdmQC)/0.5291772109 ! Convert the distance from Å to Bohr
                                my(i,j)=my(i,j)+QC_charge*(ycdmQC1-ycdmQC)/0.5291772109 ! Convert the distance from Å to Bohr
                                mz(i,j)=mz(i,j)+QC_charge*(zcdmQC1-zcdmQC)/0.5291772109 ! Convert the distance from Å to Bohr
                        end if
                end do
        end do
        allocate(E0(nroot))
        allocate(w(nroot))
        read(2,*) ! Ground-state energy
        read(2,*) E0(1) 
        read(2,*) ! Excitation energies (eV)
        counter=1
        do i=1,nroot-1
                counter=counter+1
                read(2,*) ddelta
                E0(counter)=E0(1)+ddelta*0.0367493 ! Convert excitation energies (eV) to absolute state energies (Hartree)
        end do
        read(2,*) ! Calculation type
        !   0 -> perturbed electronic state properties
        !   otherwise -> spectral absorption band
        read(2,*) ctype
        close(2)
        
        ! Read the indices of the system atoms in the environment
        open(unit=20, file=solvent, action='READ')
        read(20,*) nsolv
        allocate(a(nsolv))
        allocate(dist1(nsolv))
        
        ! Check that nsolv=nat-ninte
        if (nsolv.ne.(nat-ninte)) print*,'ERROR: WRONG NUMBER OF ENVIRONMENT ATOMS'
        do i=1,nsolv
                read(20,*) a(i)
        end do 
        close(20)
        open(unit=25, file=solute, action='READ') 
        read(25,*) nsolu
        allocate(asolu(nsolu))
        allocate(mass(nsolu))
        if (nsolu.ne.ninte1) print*,'ERROR: WRONG NUMBER OF QC ATOMS'
        QC_mass=0.0
        do i=1, nsolu
                read(25,*) asolu(i), mass(i) 
                QC_mass=QC_mass+mass(i)
        end do
        close(25)
        
        allocate(r0_2(3*nsolu))
        allocate(rfin(3*nsolu))
        
        ! Read the MD simulation charges file
        allocate(num(nat))
        allocate(ch(nat))
        
        open(unit=30, file=topol, action='READ')
        do j=1,nat
                read(30,*) sdum, ch(j)
        end do 
        close(30)

        ! Build a Dirac delta function to be used later
        allocate(delta(nroot,nroot))
        do i=1,nroot
                do j=1,nroot
                        if (i==j) then
                                delta(i,j)=1
                        else
                                delta(i,j)=0
                        end if
                end do
        end do
        
        allocate(x(tot))
        allocate(mxp(nroot,nroot))
        allocate(myp(nroot,nroot))
        allocate(mzp(nroot,nroot))
        allocate(ExMux(nroot,nroot))
        allocate(EyMuy(nroot,nroot))
        allocate(EzMuz(nroot,nroot))
        allocate(Hprim(nroot,nroot))       
        allocate(Hdiag(nrad,nrad))
        allocate(Eprim(nrad))
        allocate(Eivec(nrad,nrad))        
        allocate(Tvec(nrad,nrad))   
        allocate(d0x(nrad))
        allocate(d0y(nrad))
        allocate(d0z(nrad))   

        ! Loop over the electronic states
        do ind=1,nroot
                ! Calculation type:
                !   0 -> perturbed electronic state properties
                !   otherwise -> spectral absorption band
                if (ctype==0) then
                        ipart=ind
                        jarr=ind
                else
                        ipart=1
                        jarr=ind
                endif

                write(output, '(A,I0,A)') 'output_S', ind-1, '.xvg'
                open(unit=10, file=output)
                write(10,*) "#     frame   Energy (a.u.)             mux (a.u.)" // &
                "                muy (a.u.)                muz (a.u.)                |mu| (a.u.)"
                write(dipoles, '(A,I0,A)') 'dipole_fluctuations_S', ind-1, '.xvg'
                open(unit=9, file=dipoles)
                write(9,*) "#     frame   mux/mux0                 muy/muy0                   muz/muz0"
                write(proj, '(A,I0,A)') 'unpert_state_coefficients_to_S', ind-1, '_pert.xvg'
                open(unit=251, file=proj)

                call xdrfopen (id1,traj,'r',ret)
                if (ret.eq.0) then
                        print*, 'ERROR: unable to open XTC trajectory file'
                        print*, 'File:', traj
                        stop
                end if 
                
                print*, 'info: XTC trajectory file successfully opened'
                
                do frame=1,nframe
                        call readxtc(id1,nat,step,time,box,x,prec,ret) 
                        do j=1,tot
                                x(j)=x(j)*10 ! Convert from nm to Å
                        end do
                        ! Compute the center of mass of the quantum center used for fitting in the MD reference frame
                        xcdm=0.0
                        ycdm=0.0
                        zcdm=0.0
                        
                        do m=1,nsolu
                                g=asolu(m)
                                xcdm=xcdm+x(3*g-2)*mass(m)
                                ycdm=ycdm+x(3*g-1)*mass(m)
                                zcdm=zcdm+x(3*g)*mass(m)
                        end do
                        xcdm=xcdm/QC_mass
                        ycdm=ycdm/QC_mass
                        zcdm=zcdm/QC_mass
                        do m=1,nsolu
                                g=asolu(m)
                        end do
                        j=1
                        do m=1,nsolu
                                g=asolu(m)
                                r0_2(j)=x(3*g-2)-xcdm
                                r0_2(j+1)=x(3*g-1)-ycdm
                                r0_2(j+2)=x(3*g)-zcdm
                                j=j+3
                        enddo
                        
                        call mwsfit(nsolu,r0_2,r0_1,rfin,AA,BB,CC,mass,QC_mass,the,phi,psi)
                        
                        ! The AABBCC matrix performing the rotation from MD to QM
                        ! has been obtained; now apply its transpose to obtain dipoles from QM to MD
                        ! Compute the projections of dipoles onto the three vectors
                        
                        do i=1,nroot
                                do j=1,nroot
                                        mxp(i,j)=mx(i,j)*AA(1)+my(i,j)*AA(2)+mz(i,j)*AA(3)
                                        myp(i,j)=mx(i,j)*BB(1)+my(i,j)*BB(2)+mz(i,j)*BB(3) 
                                        mzp(i,j)=mx(i,j)*CC(1)+my(i,j)*CC(2)+mz(i,j)*CC(3)
                                end do
                        end do
                        ExQMMM_MD=ExQMMM*AA(1)+EyQMMM*AA(2)+EzQMMM*AA(3)
                        EyQMMM_MD=ExQMMM*BB(1)+EyQMMM*BB(2)+EzQMMM*BB(3)
                        EzQMMM_MD=ExQMMM*CC(1)+EyQMMM*CC(2)+EzQMMM*CC(3)
                        
                        j=1
                        do m=1,nsolu
                                g=asolu(m)
                                j=j+3
                        enddo
                        
                        Enerx=0.0
                        Enery=0.0
                        Enerz=0.0
                        Potential=0.0
                        
                        ! Loop over the solvent to compute the electric field
                        xcdm=xcdm/0.5291772109 ! Convert the distance from Å to Bohr
                        ycdm=ycdm/0.5291772109 ! Convert the distance from Å to Bohr
                        zcdm=zcdm/0.5291772109 ! Convert the distance from Å to Bohr
                        do m=1,nsolv
                                k=a(m)
                                distx=((x(3*k-2)/0.5291772109)-xcdm)**2 ! Convert the distance from Å to Bohr
                                disty=((x(3*k-1)/0.5291772109)-ycdm)**2 ! Convert the distance from Å to Bohr
                                distz=((x(3*k)/0.5291772109)-zcdm)**2   ! Convert the distance from Å to Bohr
                                dist=sqrt(distx+disty+distz)
                                dist1(m)=dist
                                Enerx=Enerx+(ch(k)*(xcdm-(x(3*k-2)/0.5291772109))/(dist)**3) ! Convert the distance from Å to Bohr
                                Enery=Enery+(ch(k)*(ycdm-(x(3*k-1)/0.5291772109))/(dist)**3) ! Convert the distance from Å to Bohr
                                Enerz=Enerz+(ch(k)*(zcdm-(x(3*k)/0.5291772109))/(dist)**3)   ! Convert the distance from Å to Bohr
                                ! Compute the electric potential at the center of mass
                                Potential=Potential+ch(k)/dist
                        end do
                        
                        Enerx=Enerx-ExQMMM_MD
                        Enery=Enery-EyQMMM_MD
                        Enerz=Enerz-EzQMMM_MD
                        Potential=Potential-PotQMMM
                               
                        do i=1,nroot
                                do j=1,nroot
                                        ExMux(i,j)=-Enerx*mxp(i,j)
                                        EyMuy(i,j)=-Enery*myp(i,j)
                                        EzMuz(i,j)=-Enerz*mzp(i,j)
                                end do
                        end do
                        
                        ! Build the perturbed Hamiltonian matrix to be diagonalized
                        do i=1,nroot
                                do j=1,nroot
                                        Hprim(i,j)=E0(i)*delta(i,j)+ExMux(i,j)+EyMuy(i,j)+EzMuz(i,j)
                                end do
                        end do
                        
                        do i=1,nrad
                                do j=1,nrad
                                        Hdiag(i,j)=Hprim(i,j)
                                end do
                        end do
                        
                        ! Call the diagonalization subroutine
                        lwork=-1
                        jobz='V'
                        uplo='U'
                        lda=nrad
                        lwmax=1000
                        
                        call DSYEV('V','U',nrad,Hprim,lda,w,work,lwork,info)
                        lwork=min(lwmax,int(work(1)))
                        ! Solve the eigenproblem
                        call DSYEV('V','U',nrad,Hprim,lda,w,work,lwork,info)
                        
                        ! Eprim are the eigenvalues and Eivec the eigenvectors
                        do i=1,nrad
                                Eprim(i)=w(i)+QC_charge*Potential
                        end do
                        
                        do i=1,nrad
                                do j=1,nrad
                                        Eivec(i,j)=Hprim(i,j)
                                end do
                        end do
                        
                        ! Compute the transpose of the Hdiag matrix
                        do i=1,nrad
                                do j=1,nrad  
                                        Tvec(i,j)=Eivec(j,i)
                                end do
                        end do
        
                        ! Write the coefficients of the unperturbed states contributing to the perturbed state
                        write(251,*) frame-1, (Eivec(i,jarr),i=1,nroot)
                        
                        ! Compute the components of the d0 vector
                        do i=1,nrad
                                d0x(i)=0.0 
                                d0y(i)=0.0 
                                d0z(i)=0.0 
                                do j=1,nrad
                                        d0x(i)=d0x(i)+mxp(i,j)*Eivec(j,ipart)           
                                        d0y(i)=d0y(i)+myp(i,j)*Eivec(j,ipart) 
                                        d0z(i)=d0z(i)+mzp(i,j)*Eivec(j,ipart)      
                                end do
                        end do
                        
                        mu0x=0.0
                        mu0y=0.0
                        mu0z=0.0
                        do i=1,nrad 
                                mu0x=mu0x+d0x(i)*Tvec(jarr,i)
                                mu0y=mu0y+d0y(i)*Tvec(jarr,i)  
                                mu0z=mu0z+d0z(i)*Tvec(jarr,i)   
                        end do
                        
                        write(10,*) frame-1, Eprim(jarr), mu0x, mu0y, mu0z, sqrt(mu0x**2+mu0y**2+mu0z**2)
                        write(9,*) frame-1, mu0x/mx(ipart,jarr), mu0y/my(ipart,jarr), mu0z/mz(ipart,jarr)
                end do
                close(9)
                close(10)
                close(251)
        enddo

99      format(a5,f13.3,2x,f13.3,2x,f13.3,f15.6,a8)
100     format(i6,f12.3,f12.3)
202     format(f12.3,f12.3,f12.3)
320     format(i6)
300     format(i10)
350     format(3i6)
400     format(i3,3f6.3)
420     format(i6,2x,e12.6,2x,e12.6,2x,e12.6)
421     format(i6,i3,i3,2x,e12.6,2x,e12.6,2x,e12.6)
450     format(i6,f8.3)
1010    format(15i5)
end program PMMQMMM

subroutine angles(the,phi,psi,A,B,C,Athe,Aphi,Apsi,Bthe,Bphi,Bpsi,Cthe,Cphi,Cpsi)
        implicit none

        !-------------------- Euler angles --------------------
        double precision :: the, phi, psi
        
        !-------------------- Trigonometric functions --------------------
        double precision :: cosphi, sinphi, costhe, sinthe, cospsi, sinpsi
        
        !-------------------- Rotation vectors --------------------
        double precision :: A(3), B(3), C(3)
        
        !-------------------- Derivatives of rotation vectors --------------------
        double precision :: Athe(3), Aphi(3), Apsi(3)
        double precision :: Bthe(3), Bphi(3), Bpsi(3)
        double precision :: Cthe(3), Cphi(3), Cpsi(3)

        cosphi=cos(phi)
        sinphi=sin(phi)
        costhe=cos(the)
        sinthe=sin(the)
        cospsi=cos(psi)
        sinpsi=sin(psi)
        
        A(1)=cosphi*cospsi-sinphi*costhe*sinpsi
        A(2)=-(cosphi*sinpsi+sinphi*costhe*cospsi)
        A(3)=sinphi*sinthe
        
        B(1)=sinphi*cospsi+cosphi*costhe*sinpsi
        B(2)=cosphi*costhe*cospsi-sinphi*sinpsi
        B(3)=-cosphi*sinthe
        
        C(1)=sinthe*sinpsi
        C(2)=sinthe*cospsi
        C(3)=costhe
        
        Athe(1)=sinphi*sinthe*sinpsi
        Aphi(1)=-sinphi*cospsi-cosphi*costhe*sinpsi
        Apsi(1)=-cosphi*sinpsi-sinphi*costhe*cospsi
        Bthe(1)=-cosphi*sinthe*sinpsi
        Bphi(1)=cosphi*cospsi-sinphi*costhe*sinpsi
        Bpsi(1)=-sinphi*sinpsi+cosphi*costhe*cospsi
        Cthe(1)=costhe*sinpsi
        Cphi(1)=0.
        Cpsi(1)=sinthe*cospsi
        
        Athe(2)=sinphi*sinthe*cospsi
        Aphi(2)=sinphi*sinpsi-cosphi*costhe*cospsi
        Apsi(2)=-cosphi*cospsi+sinphi*costhe*sinpsi
        Bthe(2)=-cosphi*sinthe*cospsi
        Bphi(2)=-sinphi*costhe*cospsi-cosphi*sinpsi
        Bpsi(2)=-cosphi*costhe*sinpsi-sinphi*cospsi
        Cthe(2)=costhe*cospsi
        Cphi(2)=0.
        Cpsi(2)=-sinthe*sinpsi
        
        Athe(3)=sinphi*costhe
        Aphi(3)=cosphi*sinthe
        Apsi(3)=0.
        Bthe(3)=-cosphi*costhe
        Bphi(3)=sinphi*sinthe
        Bpsi(3)=0.
        Cthe(3)=-sinthe
        Cphi(3)=0.
        Cpsi(3)=0.
end

subroutine mwsfit(natoms,r,rref,rfin,A,B,C,ATmass,ATmasstot,the,phi,psi)
        implicit none

        integer :: i, j, lk, natoms, kll
        
        !-------------------- Angles and updates --------------------
        double precision :: the, phi, psi
        double precision :: the1, phi1, psi1
        double precision :: dthe, dphi, dpsi
        
        !-------------------- First and second derivatives --------------------
        double precision :: DtheL,  DphiL,  DpsiL
        double precision :: DtheL1, DphiL1, DpsiL1
        double precision :: DDtheL, DDphiL, DDpsiL
        
        !-------------------- Steepest descent parameters --------------------
        double precision :: ulambda, pot, pot1, deltapot
        double precision :: sum, DGDr, DGDG
        
        !-------------------- Rotation vectors and derivatives --------------------
        double precision :: A(3), B(3), C(3)
        double precision :: Athe(3), Aphi(3), Apsi(3)
        double precision :: Bthe(3), Bphi(3), Bpsi(3)
        double precision :: Cthe(3), Cphi(3), Cpsi(3)
        
        !-------------------- Coordinates and gradients --------------------
        double precision :: r(3*natoms), rref(3*natoms), rfin(3*natoms)
        double precision :: V(natoms,3)
        double precision :: Wthe(natoms,3), Wphi(natoms,3), Wpsi(natoms,3)
        
        !-------------------- Masses --------------------
        double precision :: ATmass(natoms), ATmasstot

        ! The potential to be minimized using a steepest-descent procedure is L = sum_i m_i * (|r_i - r_i^ref|^2) / m_total
        ! V(i,1) = ∂L/∂x_i   V(i,2) = ∂L/∂y_i   V(i,3) = ∂L/∂z_i (i=atom index)
        ! The rotation matrix with respect to the Euler angles is ROT = (A, B, C), where A, B, and C are three-vectors
        
        ! Initial angles
        the=0.5
        phi=0.5
        psi=0.5

        ! Compute the rotation matrix (A, B, C).
        ! Athe is the derivative of A with respect to the angle theta,
        ! Aphi is the derivative of A with respect to the angle phi, etc.
        call angles(the,phi,psi,A,B,C,Athe,Aphi,Apsi,Bthe,Bphi,Bpsi,Cthe,Cphi,Cpsi)
        do i=1,natoms
                do j=1,3
                        if(j.eq.1)lk=2
                        if(j.eq.2)lk=1
                        if(j.eq.3)lk=0
                        V(i,j)=2.d0*ATmass(i)*(r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)-rref(3*i-lk))/ATmasstot
                        Wthe(i,j)=r(3*i-2)*Athe(j)+r(3*i-1)*Bthe(j)+r(3*i)*Cthe(j)
                        Wphi(i,j)=r(3*i-2)*Aphi(j)+r(3*i-1)*Bphi(j)+r(3*i)*Cphi(j)
                        Wpsi(i,j)=r(3*i-2)*Apsi(j)+r(3*i-1)*Bpsi(j)+r(3*i)*Cpsi(j)
                enddo
        enddo
        
        
        DtheL=0.0
        do i=1,natoms
                sum=0.
                do j=1,3
                        sum=sum+V(i,j)*Wthe(i,j)
                enddo
                DtheL=DtheL+sum
        enddo
        DphiL=0.0
        do i=1,natoms
                sum=0.
                do j=1,3
                        sum=sum+V(i,j)*Wphi(i,j)
                enddo
                DphiL=DphiL+sum
        enddo
        DpsiL=0.0
        do i=1,natoms
                sum=0.
                do j=1,3
                        sum=sum+V(i,j)*Wpsi(i,j)
                enddo
                DpsiL=DpsiL+sum
        enddo
        
        ! Set the initial step size
        ulambda=1.d-10
        
        ! Compute the potential L
        pot=0.d0
        do i=1,natoms
                do j=1,3
                        if(j.eq.1)lk=2
                        if(j.eq.2)lk=1
                        if(j.eq.3)lk=0
                        pot=ATmass(i)*(r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)-rref(3*i-lk))**2/ATmasstot+pot
                enddo
        enddo

        deltapot=1.
        kll=0
        do while (dabs(deltapot).gt.1.d-15)
                kll=kll+1
                
                if (kll.ne.1) pot=pot1

                ! New angles
                the1=the-ulambda*DtheL
                phi1=phi-ulambda*DphiL
                psi1=psi-ulambda*DpsiL
                
                call angles(the1,phi1,psi1,A,B,C,Athe,Aphi,Apsi,Bthe,Bphi,Bpsi,Cthe,Cphi,Cpsi)
                
                pot1=0.
                do i=1,natoms
                        do j=1,3
                                if(j.eq.1)lk=2
                                if(j.eq.2)lk=1
                                if(j.eq.3)lk=0
                                V(i,j)=2.d0*ATmass(i)*(r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)-rref(3*i-lk))/ATmasstot
                                Wthe(i,j)=r(3*i-2)*Athe(j)+r(3*i-1)*Bthe(j)+r(3*i)*Cthe(j)
                                Wphi(i,j)=r(3*i-2)*Aphi(j)+r(3*i-1)*Bphi(j)+r(3*i)*Cphi(j)
                                Wpsi(i,j)=r(3*i-2)*Apsi(j)+r(3*i-1)*Bpsi(j)+r(3*i)*Cpsi(j)
                                pot1=ATmass(i)*(r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)-rref(3*i-lk))**2/ATmasstot+pot1
                        enddo
                enddo
                
                deltapot=pot1-pot
                DtheL1=0.0
                do i=1,natoms
                        sum=0.
                        do j=1,3
                                sum=sum+V(i,j)*Wthe(i,j)
                        enddo
                        DtheL1=DtheL1+sum
                enddo
                
                DphiL1=0.0
                do i=1,natoms
                        sum=0.
                        do j=1,3
                                sum=sum+V(i,j)*Wphi(i,j)
                        enddo
                        DphiL1=DphiL1+sum
                enddo
                DpsiL1=0.0
                do i=1,natoms
                        sum=0.
                        do j=1,3
                                sum=sum+V(i,j)*Wpsi(i,j)
                        enddo
                        DpsiL1=DpsiL1+sum
                enddo
                do i=1,3
                        DDtheL=DtheL1-DtheL
                        DDphiL=DphiL1-DphiL
                        DDpsiL=DpsiL1-DpsiL
                        dthe=the1-the
                        dphi=phi1-phi
                        dpsi=psi1-psi
                enddo
                
                DGDr=0.d0
                DGDG=0.d0
                
                DGDr=DDtheL*dthe+DDphiL*dphi+DDpsiL*dpsi
                DGDG=DDtheL*DDtheL+DDphiL*DDphiL+DDpsiL*DDpsiL
                ! New step size
                ulambda=dabs(DGDr/DGDG)
                
                the=the1
                phi=phi1
                psi=psi1
                
                DtheL=DtheL1
                DphiL=DphiL1
                DpsiL=DpsiL1
        enddo

        do i=1,natoms
                do j=1,3
                        if(j.eq.1)lk=2
                        if(j.eq.2)lk=1
                        if(j.eq.3)lk=0
                        rfin(3*i-lk)=r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)
                enddo
        enddo
end
