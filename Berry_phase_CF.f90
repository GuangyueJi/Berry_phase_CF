program main
    implicit none
    
    interface
        real(8) function genrand_real1(random_number)
            integer(4):: random_number
        end function
    end interface
    
    real*8::t1,t2
    integer*4,parameter::Ne=108,q=2,Nphi=Ne*q,Ne2=Ne**2,Np=24,Nt=10000*Ne,Ns=1000
    complex*16,parameter::cj=(0.,1.)
    
    integer*4:: d(Np+1,Ne,2),p(Ne+1,2),ccN=5
    real*8::ratio,pi,e
    complex*16::L1,L2
    
    integer*4:: i,j,k,m,n,x(Np+1,Ne,2)
    complex*16:: x0,sigma0
    complex*16:: wei_sigma,f
    complex*16,ALLOCATABLE::f0(:,:),Cphi(:),Matrix11(:,:,:),Matrix(:,:,:,:)
    integer*4::ERR_MESSAGE
    complex*16::s(Ns,Np+1,8)=0,result(Ns,Np+1,8)=0,resultf(Ns,Np+1,4)=0
    
    
    call CPU_TIME(t1)
    ALLOCATE(f0(0:Ne2-1,0:Ne2-1),Cphi(0:2*Ne2-1),Matrix11(Np+1,Ne,Ne),Matrix(Np+1,5,Ne,Ne),stat=ERR_MESSAGE)    
    if(ERR_MESSAGE.NE.0) print*,'ALLOCATION ERROR'


    pi=asin(1.d0)*2
    L1=sqrt(Nphi*pi)
    L2=cj*L1

    e=10.d0**(-5.d0)
    !$OMP PARALLEL PRIVATE(x0)
    !$OMP DO       
    do i = 0, Ne2-1 ,1 
        do j = 0, Ne2-1 ,1
            x0=( i*L1 +j*L2 )/Ne2
            f0(i,j)=wei_sigma(x0,L1,L2,e)*exp(-1.d0/(2*Nphi)*x0* conjg(x0))
        end do
    end do
    !$OMP END DO
	!$OMP END PARALLEL
    
    !$OMP PARALLEL
    !$OMP DO
    do i =0,2*Ne2-1,1
        Cphi(i)=exp(i*cj*pi/Ne2)
    end do
    !$OMP END DO
	!$OMP END PARALLEL
    

    p(1,1:2)=(/-1,-1/)
    p(2,1:2)=(/-1,0/)
    p(3,1:2)=(/-1,1/)

    p(4,1:2)=(/1,-1/)
    p(5,1:2)=(/1,0/)
    p(6,1:2)=(/1,1/)
    p(7,1:2)=(/1,2/)
    p(8,1:2)=(/2,0/)
    p(9,1:2)=(/2,1/)    
    p(10,1:2)=(/2,2/)
    p(11,1:2)=(/2,3/)
    p(12,1:2)=(/3,0/)
    p(13,1:2)=(/3,1/)
    p(14,1:2)=(/3,2/)
    p(15,1:2)=(/3,3/)
    p(16,1:2)=(/3,4/)
    p(17,1:2)=(/4,0/)
    p(18,1:2)=(/4,1/)
    p(19,1:2)=(/4,2/) 
    p(20,1:2)=(/4,3/)
    p(21,1:2)=(/4,4/)
    p(22,1:2)=(/5,0/)
    p(23,1:2)=(/5,1/)
    p(24,1:2)=(/5,2/)
    p(25,1:2)=(/5,3/)
    
    
    p(26,1:2)=(/0,-1/)
    p(27,1:2)=(/0,0/)
    p(28,1:2)=(/0,1/)
    
    p(29,1:2)=(/-5,-3/)
    p(30,1:2)=(/-5,-2/)
    p(31,1:2)=(/-5,-1/)
    p(32,1:2)=(/-5,0/)
    p(33,1:2)=(/-5,1/)
    p(34,1:2)=(/-5,2/)
    p(35,1:2)=(/-5,3/) 
    p(36,1:2)=(/-4,-4/)
    p(37,1:2)=(/-4,-3/)
    p(38,1:2)=(/-4,-2/)
    p(39,1:2)=(/-4,-1/)
    p(40,1:2)=(/-4,0/)
    p(41,1:2)=(/-4,1/)
    p(42,1:2)=(/-4,2/)    
    p(43,1:2)=(/-4,3/)
    p(44,1:2)=(/-4,4/)
    p(45,1:2)=(/-3,-5/)
    p(46,1:2)=(/-3,-4/)
    p(47,1:2)=(/-3,-3/)
    p(48,1:2)=(/-3,-2/)
    p(49,1:2)=(/-3,-1/)
    p(50,1:2)=(/-3,0/)
    p(51,1:2)=(/-3,1/)
    p(52,1:2)=(/-3,2/) 
    p(53,1:2)=(/-3,3/)
    p(54,1:2)=(/-3,4/)
    p(55,1:2)=(/-3,5/)
    p(56,1:2)=(/-2,-5/) 
    p(57,1:2)=(/-2,-4/) 
    p(58,1:2)=(/-2,-3/) 
    p(59,1:2)=(/-2,-2/) 
    p(60,1:2)=(/-2,-1/) 
    p(61,1:2)=(/-2,0/) 
    p(62,1:2)=(/-2,1/) 
    p(63,1:2)=(/-2,2/) 
    p(64,1:2)=(/-2,3/) 
    p(65,1:2)=(/-2,4/) 
    p(66,1:2)=(/-2,5/) 
    p(67,1:2)=(/-1,-5/) 
    p(68,1:2)=(/-1,-4/) 
    p(69,1:2)=(/-1,-3/)
    p(70,1:2)=(/-1,-2/)
    p(71,1:2)=(/-1,2/)
    p(72,1:2)=(/-1,3/)    
    p(73,1:2)=(/-1,4/)
    p(74,1:2)=(/-1,5/)
    p(75,1:2)=(/0,-5/)
    p(76,1:2)=(/0,-4/)
    p(77,1:2)=(/0,-3/)
    p(78,1:2)=(/0,-2/)
    p(79,1:2)=(/0,2/)
    p(80,1:2)=(/0,3/)
    p(81,1:2)=(/0,4/)
    p(82,1:2)=(/0,5/) 
    p(83,1:2)=(/1,-5/)
    p(84,1:2)=(/1,-4/)
    p(85,1:2)=(/1,-3/)
    p(86,1:2)=(/1,-2/)
    p(87,1:2)=(/1,3/)
    p(88,1:2)=(/1,4/)
    p(89,1:2)=(/1,5/)
    p(90,1:2)=(/2,-5/)
    p(91,1:2)=(/2,-4/)
    p(92,1:2)=(/2,-3/)
    p(93,1:2)=(/2,-2/)
    p(94,1:2)=(/2,-1/)
    p(95,1:2)=(/2,4/) 
    p(96,1:2)=(/2,5/)
    p(97,1:2)=(/3,-5/)
    p(98,1:2)=(/3,-4/)
    p(99,1:2)=(/3,-3/)
    p(100,1:2)=(/3,-2/)
    p(101,1:2)=(/3,-1/)
    p(102,1:2)=(/3,5/)    
    p(103,1:2)=(/4,-4/)
    p(104,1:2)=(/4,-3/)
    p(105,1:2)=(/4,-2/)
    p(106,1:2)=(/4,-1/)
    p(107,1:2)=(/5,-3/)
    p(108,1:2)=(/5,-2/)
    p(109,1:2)=(/5,-1/)


    
    d(1,:,:)=p(2:109,:)
    do i=2,Np+1,1
        d(i,:,:)=d(i-1,:,:)
        d(i,i-1,:)=p(i-1,:)
    end do
    d=d*Ne
    
    open(11, file='data11.dat', status='old')  
    open(12, file='data12.dat', status='old')
    open(13, file='data13.dat', status='old')
    open(14, file='data14.dat', status='old')
    open(15, file='data15.dat', status='old')  
    open(16, file='data16.dat', status='old')
    open(17, file='data17.dat', status='old')  
    open(18, file='data18.dat', status='old')
    open(19, file='data19.dat', status='old')
    open(20, file='data20.dat', status='old')
    open(21, file='data21.dat', status='old')  
    open(22, file='data22.dat', status='old')
    
    
    !this is for the Markov chain
    !$OMP PARALLEL
    !$OMP DO   
    do i =1,Np+1,1
        call mainFunction1(i,d,Ne,q,Nphi,Ne2,Np,f0,Cphi,Nt,Matrix11(i,:,:),x(i,:,:))
    end do 
    !$OMP END DO
	!$OMP END PARALLEL
    
    
    !this is for the number of slices
    s=0
    do j=1,Ns,1
        
        !$OMP PARALLEL
        !$OMP DO   
        do i =1,Np+1,1
            call mainFunction2(i,j,d,s,Ne,q,Nphi,Ne2,Np,f0,Cphi,Nt,Ns,Matrix11(i,:,:),x(i,:,:))
        end do 
        !$OMP END DO
	    !$OMP END PARALLEL
    
    
        write(11,*)j
        write(12,*)j
        write(13,*)j 

        write(14,*)j
        write(15,*)j
        write(16,*)j
    
        write(17,*)j
        write(18,*)j
        write(19,*)j
    
        write(20,*)j
        write(21,*)j
        write(22,*)j


        do i =1,Np+1,1     

    
        if (i.ne.(Np+1)) then
            result(j,i,1) =s(j,i,1) / (Nt/10)
            result(j,i,2) = sqrt(  s(j,i,2)/(Nt/10)    )
            resultf(j,i,1)=result(j,i,1)/result(j,i,2)
            write(11,*)resultf(j,i,1) 
            write(12,*)abs(resultf(j,i,1))
            write(13,*)atan2(aimag(resultf(j,i,1)),real(resultf(j,i,1))) 

            result(j,i,5) =s(j,i,5) / (Nt/10)
            result(j,i,6) = sqrt(  s(j,i,6)/(Nt/10)    )
            resultf(j,i,3)=result(j,i,5)/result(j,i,6)
            write(17,*)resultf(j,i,3) 
            write(18,*)abs(resultf(j,i,3))
            write(19,*)atan2(aimag(resultf(j,i,3)),real(resultf(j,i,3)))
        end if
    
    
        if(i.ne.1) then
            result(j,i,3) =s(j,i,3) / (Nt/10) 
            result(j,i,4) = sqrt(  s(j,i,4)/(Nt/10)    )
            resultf(j,i,2)=result(j,i,3)/result(j,i,4)
            write(14,*)resultf(j,i,2) 
            write(15,*)abs(resultf(j,i,2))
            write(16,*)atan2(aimag(resultf(j,i,2)),real(resultf(j,i,2))) 
    
            result(j,i,7) =s(j,i,7) / (Nt/10) 
            result(j,i,8) = sqrt(  s(j,i,8)/(Nt/10)    )
            resultf(j,i,4)=result(j,i,7)/result(j,i,8)
            write(20,*)resultf(j,i,4) 
            write(21,*)abs(resultf(j,i,4))
            write(22,*)atan2(aimag(resultf(j,i,4)),real(resultf(j,i,4))) 
        end if

        
        end do

    end do 
    
    
    

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
    
    call cpu_time ( t2  )
    write ( *, *  ) 'Elapsed CPU time = ', t2 - t1   
    end program main  

    
    
    
    
subroutine mainFunction1(i,d,Ne,q,Nphi,Ne2,Np,f0,Cphi,Nt,Matrix11,x)
    implicit none
    integer*4::i,d(Np+1,Ne,2),Ne,q,Nphi,Ne2,Np,Nt
    complex*16::f0(0:Ne2-1,0:Ne2-1),Cphi(0:2*Ne2-1)


    integer*4:: j,k,m,n,ii,kk,jj
    integer*4::x(Ne,2),x2(Ne,2),alpha1(q, 2),alpha2(q, 2),alpha3(q, 2),num1,num2,num3,dbar1(2),dbar2(2),dbar3(2),ccN,ccr,Zt(2),ccd(2)
    real*8:: xt(Ne,2),genrand_real1
    real*8::m11,m12,xm(5),T
    complex*16::psi11,psi12,xpsi(5),Meigen11(Ne),Meigen12(Ne),xMeigen(5,Ne),Matrix11(Ne,Ne),Matrix12(Ne,Ne),Matrix(5,Ne,Ne)
    complex*16::xpsir(5),Txpsi(5),Txm(5),r1

    
    alpha1(:, :)=0
    dbar1=0
    do j=1,Ne,1
        alpha1(1, :)=alpha1(1, :)+d(i,j,:)/2
        dbar1(:)=dbar1(:)+d(i,j,:)/Ne
    end do
    
    if(i.ne.(Np+1)) then 
        alpha2(:, :)=0
        dbar2=0
        do j=1,Ne,1
            alpha2(1, :)=alpha2(1, :)+d(i+1,j,:)/2
            dbar2(:)=dbar2(:)+d(i+1,j,:)/Ne
        end do
    end if
    
    if(i.ne.1) then
        alpha3(:, :)=0
        dbar3=0
        do j=1,Ne,1
            alpha3(1, :)=alpha3(1, :)+d(i-1,j,:)/2
            dbar3(:)=dbar3(:)+d(i-1,j,:)/Ne
        end do 
    end if
    !generate the random number of x(i,j)
    do ii=1,Ne,1
        do j=1,2,1
            x(ii,j) = floor(Nphi*genrand_real1(0))*(Ne/2) !the range of x is from 0 to Nphi
        end do
    end do
    call PsiCFL0(psi11,Meigen11,m11,x,Matrix11,alpha1,dbar1,d(i,:,:),d(i,:,:),Ne,q,Nphi,Ne2,f0,Cphi)

    
    ccN=5
    num1=0
    num2=0
    num3=0
    do  m=1,Nt/ccN-1,1  
        
        ii=modulo(m,Ne)+1
        x2=x
        do j=1,2,1
            x2(ii,j)= modulo((x(ii,j)/(Ne/2)+floor(51*genrand_real1(0))-25),Nphi)*(Ne/2)  
        end do
        Matrix12=Matrix11
        call PsiCFL(ii,psi12,Meigen12,m12,x2,Matrix12,x,alpha1,dbar1,d(i,:,:),d(i,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
        
        T=0
        do kk=1,Ne,1
            T=T+2*LOG(abs(Meigen12(kk))/abs(Meigen11(kk)))
        end do
        T =T+ 2*Ne*LOG(m12/m11) + 2*LOG(abs(psi12)/abs(psi11)) 
        if (T>=0) then
            x=x2
            psi11=psi12
            m11=m12
            Meigen11=Meigen12
            Matrix11=Matrix12
            num1=num1+1
        else
            if (LOG(genrand_real1(0))<=T) then
                x=x2
                psi11=psi12
                m11=m12
                Meigen11=Meigen12
                Matrix11=Matrix12
                num2=num2+1
            else
                num3=num3+1
            end if
        end if
    end do  
    print*,num1,num2,num3

    end subroutine mainFunction1
        

    
    
    
    
subroutine mainFunction2(i,j,d,s,Ne,q,Nphi,Ne2,Np,f0,Cphi,Nt,Ns,Matrix11,x)
    implicit none
    integer*4::i,j,d(Np+1,Ne,2),Ne,q,Nphi,Ne2,Np,Nt,Ns
    complex*16::f0(0:Ne2-1,0:Ne2-1),Cphi(0:2*Ne2-1)

    integer*4:: k,m,n,ii,kk,jj
    integer*4::x(Ne,2),x2(Ne,2),alpha1(q, 2),alpha2(q, 2),alpha3(q, 2),num1,num2,num3,dbar1(2),dbar2(2),dbar3(2),ccN,ccr,Zt(2),ccd(2)
    real*8:: xt(Ne,2),genrand_real1
    real*8::m11,m12,xm(5),T
    complex*16::psi11,psi12,xpsi(5),Meigen11(Ne),Meigen12(Ne),xMeigen(5,Ne),Matrix11(Ne,Ne),Matrix12(Ne,Ne),Matrix(5,Ne,Ne)
    complex*16::xpsir(5),Txpsi(5),Txm(5),r1
    complex*16::s(Ns,Np+1,8)
       
       
        alpha1(:, :)=0
        dbar1=0
        do jj=1,Ne,1
            alpha1(1, :)=alpha1(1, :)+d(i,jj,:)/2
            dbar1(:)=dbar1(:)+d(i,jj,:)/Ne
        end do
    
        if(i.ne.(Np+1)) then 
            alpha2(:, :)=0
            dbar2=0
            do jj=1,Ne,1
                alpha2(1, :)=alpha2(1, :)+d(i+1,jj,:)/2
                dbar2(:)=dbar2(:)+d(i+1,jj,:)/Ne
            end do
        end if
    
        if(i.ne.1) then
            alpha3(:, :)=0
            dbar3=0
            do jj=1,Ne,1
                alpha3(1, :)=alpha3(1, :)+d(i-1,jj,:)/2
                dbar3(:)=dbar3(:)+d(i-1,jj,:)/Ne
            end do 
        end if

        call PsiCFL0(psi11,Meigen11,m11,x,Matrix11,alpha1,dbar1,d(i,:,:),d(i,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
        
        ccN=5
        
        m=Nt/ccN + Nt/10*(j-1)
        ii=modulo(m,Ne)+1
        x2=x
        do jj=1,2,1
            x2(ii,jj)= modulo((x(ii,jj)/(Ne/2)+floor(51*genrand_real1(0))-25),Nphi)*(Ne/2)  !x2 belongs to (0,Nphi-1)
        end do
        
        Matrix12=Matrix11
        call PsiCFL(ii,psi12,Meigen12,m12,x2,Matrix12,x,alpha1,dbar1,d(i,:,:),d(i,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
        T=0
        do kk=1,Ne,1
            T=T+2*LOG(abs(Meigen12(kk))/abs(Meigen11(kk)))
        end do
        T =T+ 2*Ne*LOG(m12/m11) + 2*LOG(abs(psi12)/abs(psi11))
        if (T>=0) then
            x=x2
            xpsi(1)=psi12
            xMeigen(1,:)=Meigen12
            xm(1)=m12
            Matrix(1,:,:)=Matrix12
            if(i.ne.(Np+1)) then
                call PsiCFL0(xpsi(2),xMeigen(2,:),xm(2),x2,Matrix(2,:,:),alpha2,dbar2,d(i,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                !call PsiCFL0(xpsi(3),xMeigen(3,:),xm(3),x2,Matrix(3,:,:),alpha2,dbar2,d(i+1,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
            end if
            
            if(i.ne.1) then
                call PsiCFL0(xpsi(4),xMeigen(4,:),xm(4),x2,Matrix(4,:,:),alpha3,dbar3,d(i,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                !call PsiCFL0(xpsi(5),xMeigen(5,:),xm(5),x2,Matrix(5,:,:),alpha3,dbar3,d(i-1,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
            end if
            
            !num1=num1+1
        else
            if (LOG(genrand_real1(0))<=T) then
                x=x2
                xpsi(1)=psi12
                xMeigen(1,:)=Meigen12
                xm(1)=m12
                Matrix(1,:,:)=Matrix12
                if(i.ne.(Np+1)) then
                    call PsiCFL0(xpsi(2),xMeigen(2,:),xm(2),x2,Matrix(2,:,:),alpha2,dbar2,d(i,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                    !call PsiCFL0(xpsi(3),xMeigen(3,:),xm(3),x2,Matrix(3,:,:),alpha2,dbar2,d(i+1,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                end if
            
                if(i.ne.1) then
                    call PsiCFL0(xpsi(4),xMeigen(4,:),xm(4),x2,Matrix(4,:,:),alpha3,dbar3,d(i,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                    !call PsiCFL0(xpsi(5),xMeigen(5,:),xm(5),x2,Matrix(5,:,:),alpha3,dbar3,d(i-1,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                end if
                !num2=num2+1
            else
                xpsi(1)=psi11
                xMeigen(1,:)=Meigen11
                xm(1)=m11
                Matrix(1,:,:)=Matrix11
                if(i.ne.(Np+1)) then
                    call PsiCFL0(xpsi(2),xMeigen(2,:),xm(2),x,Matrix(2,:,:),alpha2,dbar2,d(i,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                   ! call PsiCFL0(xpsi(3),xMeigen(3,:),xm(3),x,Matrix(3,:,:),alpha2,dbar2,d(i+1,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                end if
            
                if(i.ne.1) then
                    call PsiCFL0(xpsi(4),xMeigen(4,:),xm(4),x,Matrix(4,:,:),alpha3,dbar3,d(i,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                    !call PsiCFL0(xpsi(5),xMeigen(5,:),xm(5),x,Matrix(5,:,:),alpha3,dbar3,d(i-1,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                end if
                !num3=num3+1
            end if
        end if
        psi11=xpsi(1)
        Meigen11=xMeigen(1,:)
        m11=xm(1)
        Matrix11=Matrix(1,:,:)
        
        
        
        
            num1=0
            num2=0
            num3=0
            
           
            do  m=Nt/ccN + Nt/10*(j-1)+1,Nt/ccN + Nt/10*j, 1 
            
                ii=modulo(m,Ne)+1
                x2=x
                do jj=1,2,1
                    x2(ii,jj)= modulo((x(ii,jj)/(Ne/2)+floor(51*genrand_real1(0))-25),Nphi)*(Ne/2)
                end do
                
                Matrix12=Matrix11
                call PsiCFL(ii,psi12,Meigen12,m12,x2,Matrix12,x,alpha1,dbar1,d(i,:,:),d(i,:,:),Ne,q,Nphi,Ne2,f0,Cphi) 
                T=0
                do kk=1,Ne,1
                    T=T+2*LOG(abs(Meigen12(kk))/abs(Meigen11(kk)))
                end do
                T =T+ 2*Ne*LOG(m12/m11) + 2*LOG(abs(psi12)/abs(psi11))
                if (T>=0) then

                    if(i.ne.(Np+1)) then
                        call PsiCFL(ii,xpsi(2),xMeigen(2,:),xm(2),x2,Matrix(2,:,:),x,alpha2,dbar2,d(i,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                        !call PsiCFL(ii,xpsi(3),xMeigen(3,:),xm(3),x2,Matrix(3,:,:),x,alpha2,dbar2,d(i+1,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                    end if
            
                    if(i.ne.1) then
                        call PsiCFL(ii,xpsi(4),xMeigen(4,:),xm(4),x2,Matrix(4,:,:),x,alpha3,dbar3,d(i,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                        !call PsiCFL(ii,xpsi(5),xMeigen(5,:),xm(5),x2,Matrix(5,:,:),x,alpha3,dbar3,d(i-1,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                    end if
                    x=x2
                    xpsi(1)=psi12
                    xMeigen(1,:)=Meigen12
                    xm(1)=m12
                    Matrix(1,:,:)=Matrix12
                    
                    num1=num1+1
                else
                    if (LOG(genrand_real1(0))<=T) then

                        if(i.ne.(Np+1)) then
                            call PsiCFL(ii,xpsi(2),xMeigen(2,:),xm(2),x2,Matrix(2,:,:),x,alpha2,dbar2,d(i,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                            !call PsiCFL(ii,xpsi(3),xMeigen(3,:),xm(3),x2,Matrix(3,:,:),x,alpha2,dbar2,d(i+1,:,:),d(i+1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                        end if
            
                        if(i.ne.1) then
                            call PsiCFL(ii,xpsi(4),xMeigen(4,:),xm(4),x2,Matrix(4,:,:),x,alpha3,dbar3,d(i,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                            !call PsiCFL(ii,xpsi(5),xMeigen(5,:),xm(5),x2,Matrix(5,:,:),x,alpha3,dbar3,d(i-1,:,:),d(i-1,:,:),Ne,q,Nphi,Ne2,f0,Cphi)
                        end if
                        x=x2
                        xpsi(1)=psi12
                        xMeigen(1,:)=Meigen12
                        xm(1)=m12
                        Matrix(1,:,:)=Matrix12
                        
                        num2=num2+1
                    else
                        num3=num3+1
                    end if
                end if
            psi11=xpsi(1)
            Meigen11=xMeigen(1,:)
            m11=xm(1)
            Matrix11=Matrix(1,:,:)
            

            Txm=0
            if(i.ne.(Np+1)) then
                do kk=1,Ne,1
                    Txm(2)=Txm(2)+LOG(xMeigen(2,kk)/xMeigen(1,kk))
                end do
                Txm(2)=Txm(2)+Ne*LOG(xm(2)/xm(1))
                xpsir(2)=xpsi(2)/xpsi(1)*EXP(Txm(2))
                s(j,i,1)=s(j,i,1) + xpsir(2)
                s(j,i,2)=s(j,i,2) + abs(xpsir(2))**2
            

            end if
        
            if(i.ne.1) then
                do kk=1,Ne,1
                    Txm(4)=Txm(4)+LOG(xMeigen(4,kk)/xMeigen(1,kk))
                end do
                Txm(4)=Txm(4)+Ne*LOG(xm(4)/xm(1))
                xpsir(4)=xpsi(4)/xpsi(1)*EXP(Txm(4))
                s(j,i,3)=s(j,i,3) + xpsir(4)
                s(j,i,4)=s(j,i,4) + abs(xpsir(4))**2
            
            end if

            end do 
            print*,num1,num2,num3,i,j

end subroutine

    
    

    
        
function wei_sigma(z,L1,L2,e)
    implicit none
    complex*16:: z,L1,L2
    real*8::e
    integer*4:: n, k ,m,i,j
    complex*16:: xw,product1 ,product2,  L(2,2),  w1(2,2)
    complex*16::wei_sigma !the value of return
    n=0
    k=1 !can't initial it as when declare
    product1 = z !don,t change the initial value
    do while(.true.)
        do n = 0,k,1
            m = k-n
            product2 = product1
            do i =1,2,1
                do j =1,2,1
                    L(i,j)=(-1.d0)**i *m* L1 + (-1.d0)**j *n*L2
                    !xw stands for z/L
                    xw=z/L(i,j)
                    w1(i,j)=(1 - xw) * exp( xw + 1.d0/2.d0 * xw ** 2 )
                    product1 = product1 * w1(i,j)   
                end do
            end do
            
            if (n==0) then
                product1 =( product1 / w1(1,1) )/w1(2,1)
            end if
            
            if (m==0) then
                product1 =( product1 / w1(1,2) )/w1(1,1) 
            end if
        end do

        if (abs(product2-product1) < e) then
            wei_sigma=product1 !wei_sigma is the return value
            exit
        else
            k = k + 1
            product2 = product1
        end if 
        
    end do
    end function

function f(i,j,f0,Ne,Nphi,Ne2,Cphi)
    implicit none
    integer*4::i,j,Ne,Nphi,Ne2
    complex*16::f0(0:Ne2-1,0:Ne2-1), Cphi(0:2*Ne2-1)
    integer*4::m1,m2,n1,n2,ccf0
    complex*16::f

    m1=modulo(i,Ne2)
    n1=(i-m1)/(Ne2)
    m2=modulo(j,Ne2)
    n2=(j-m2)/(Ne2) 
    ccf0=modulo((n1*m2-n2*m1),(2*Ne2))
    if (modulo(n1*n2+n1+n2,2)==0) then
        f=f0(m1,m2)*Cphi(ccf0)
    else
        f=-f0(m1,m2)*Cphi(ccf0)
    end if            
    end function

    
    
subroutine PsiCFL0(psi,Meigen,m11,z,Matrix1,alpha,dbar,d1,d,Ne,q,Nphi,Ne2,f0,Cphi)
    implicit none
    integer*4::Ne,q,Nphi,Ne2,z(Ne,2),alpha(q,2),d(Ne,2),d1(Ne,2),dbar(2)
    complex*16::f0(0:Ne2-1,0:Ne2-1), Cphi(0:2*Ne2-1)
    integer*4::i,j,k,Zt(2),ccM
    real*8::m11
    complex*16::Matrix1(Ne,Ne),f,Meigen(Ne),psi

    Zt=0 !this is important
    do k =1,Ne,1
        Zt=Zt+z(k,:)
    end do
    
    do i =1,Ne,1
        do j =1,Ne,1              
            ccM= modulo( (d1(j,1)*z(i,2) -d1(j,2)*z(i,1) )*2/Ne , (2*Ne2) )
            Matrix1(i,j)= Cphi(ccM)
            ccM= modulo( -(d(j,1)*z(i,2) -d(j,2)*z(i,1) )/Ne , (2*Ne2) )
            Matrix1(i,j)=Matrix1(i,j)*Cphi(ccM)
            do k =1,Ne,1
                if (k/=i) then
                    Matrix1(i,j)=Matrix1(i,j)*f(z(i,1)-z(k,1)-d(j,1)+dbar(1),z(i,2)-z(k,2)-d(j,2)+dbar(2),f0,Ne,Nphi,Ne2,Cphi)
                end if
            end do
        end do
    end do
    m11=MAXVAL(abs(Matrix1))
    call determinant(Ne,Matrix1/m11,Meigen)
    psi = ( f(Zt(1) -alpha(1,1)  , Zt(2) -alpha(1,2) ,f0,Ne,Nphi,Ne2,Cphi)   )**q 
    end subroutine     
    
subroutine PsiCFL(ii,psi,Meigen,m11,z2,Matrix2,z1,alpha,dbar,d1,d,Ne,q,Nphi,Ne2,f0,Cphi)
    implicit none
    integer*4::ii,Ne,q,Nphi,Ne2,z2(Ne,2),z1(Ne,2),alpha(q,2),d(Ne,2),d1(Ne,2),dbar(2)
    complex*16::f0(0:Ne2-1,0:Ne2-1), Cphi(0:2*Ne2-1)
    integer*4::i,j,k,Zt(2),ccM
    real*8::m11
    complex*16::Matrix2(Ne,Ne),f,Meigen(Ne),psi

    Zt=0 !this is important
    do k =1,Ne,1
        Zt=Zt+z2(k,:)
    end do
    
    do i =1,Ne,1
        if (i==ii) then
            do j =1,Ne,1              
                ccM= modulo( (d1(j,1)*z2(i,2) -d1(j,2)*z2(i,1) )*2/Ne , (2*Ne2) )
                Matrix2(i,j)= Cphi(ccM)
                ccM= modulo( -(d(j,1)*z2(i,2) -d(j,2)*z2(i,1) )/Ne , (2*Ne2) )
                Matrix2(i,j)=Matrix2(i,j)*Cphi(ccM)
                do k =1,Ne,1
                    if (k/=i) then
                        Matrix2(i,j)=Matrix2(i,j)*f(z2(i,1)-z2(k,1)-d(j,1)+dbar(1),z2(i,2)-z2(k,2)-d(j,2)+dbar(2),f0,Ne,Nphi,Ne2,Cphi)
                    end if
                end do
            end do
        else
            do j =1,Ne,1 
                Matrix2(i,j)=Matrix2(i,j)*f(z2(i,1)-z2(ii,1)-d(j,1)+dbar(1),z2(i,2)-z2(ii,2)-d(j,2)+dbar(2),f0,Ne,Nphi,Ne2,Cphi)&
                    &/f(z2(i,1)-z1(ii,1)-d(j,1)+dbar(1),z2(i,2)-z1(ii,2)-d(j,2)+dbar(2),f0,Ne,Nphi,Ne2,Cphi)
            end do    
        end if
    end do
    m11=MAXVAL(abs(Matrix2))
    call determinant(Ne,Matrix2/m11,Meigen)
    psi = ( f(Zt(1) -alpha(1,1)  , Zt(2) -alpha(1,2) ,f0,Ne,Nphi,Ne2,Cphi)   )**q 
    end subroutine 
    
subroutine determinant(n,Mat,Meigen)
    implicit none
    integer*4::n
    complex*16::Mat(n,n),Meigen(n),a(n,n)
    integer*4 ::i,j,ipiv(n),info
    
    a=Mat !the mat is not changed in this way.
    call zgetrf( n, n, a, n,ipiv,info)
    !if(info/=0) write(0,*) 'Error occured in zgetrf!'
  
    do i = 1,n,1
         if(ipiv(i) .ne. i) then
             Meigen(i)=-a(i,i)
         else 
             Meigen(i)=a(i,i)
         end if
    end do

    end subroutine