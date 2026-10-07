!  EuLagPATS.f90 
!
!  FUNCTIONS:
!  EuLagPATS - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: EuLagPATS
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program EuLagPATS
    
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none
    
    real(kind=4):: Sp0
    integer:: i, j, k, n, m, nt, ni, ii, iii, jj, kk, nl, nr, np, istart, iRe
    integer:: n_data1, n_data2, n_data3, it
    integer:: mode, nsmall, modecaco3, datamode, mode_eps, mode_init
    integer:: flagtot, flag4, flagRe, im, lr, tflag
    integer:: dminint, dmaxint, nRelim, nf, ntime
    integer, allocatable:: vartime(:), flag(:)
    integer:: countzj
    integer, allocatable::  countij(:) 
    
    character(len = 1000):: fileepc100, fileepcalc100, fileepsi100, filethetao, fileo2, fileco3, fileco3satcalc
    character(len = 1000):: path, root
    character(len = 5):: nmonth, nfrac
    character(len = 100), allocatable:: frac(:)
    
    real(kind=4):: fluxtot, romid, k1i, k2i, k3i, k4i, gammai, timez
    real(kind=4), allocatable:: k1(:), k2(:), k3(:), k4(:)
    real(kind=4), allocatable:: pM(:), pV(:), masspi(:), mai(:)
    real(kind=4), allocatable:: mlimi(:), roi(:), Vi(:)
    real(kind=4), allocatable:: ai(:), bi(:), Cri(:),  Trefi(:), Q10i(:), K_O2i(:)
    real(kind=4), allocatable:: Spi_curij(:), Spi_curij1(:), Spi_curij2(:)
    
    real(kind=4):: aj, bj, Re, nu, g, Cd, Recur, Sp_measur, Fd_measur, Z_measur 
    real(kind=4)::  wp2, dwp2,  H, dz, temp, dd0,  Flux0
    real(kind=4):: wp11, wpi, wpi1, wpi2, wpn, wpn1, wpn2
    real(kind=4):: Aij, Aij1, Aij2, dt2, dz2, dz22, Spn
    real(kind=4):: Aijwpi2, Aijwpi, Aij1wpi, Aij1wpi1
    
    real(kind=4):: M0xCm, dM0xCm, masscur, massnext, dcur, dnext
    real(kind=4):: m1, m2, m3, m4
    real(kind=4):: lll, ma, massp, teta
    real(kind=4):: x, y, dp, dlim, dmin, dmax, Cm, ratio, dstart, dend
    real(kind=4):: droa, da, row, rop, V, pi
    real(kind=4):: CO3z1, CO3z2, CO3z3, CO3z4
    real (kind=4):: time, dt,  ttime, dtsave, dtime, Time_end, null
    real(kind=4)::  time1, time2, time3 
    real (kind=4):: t_max, delta_t1, delta_t2
    real (kind=4):: T1, O1,  T2, O2,  T3, O3, T4, O4, Tz
    real (kind=4):: Cw, sigma, eta, M0, gamma0, psi  
    real (kind=4):: epsilon, eps_max, eps_min
    real (kind=4):: CO3sat1, CO3sat2, CO3sat3, CO3sat4,  SIG
    real(kind=4)::  totflux, totconc
    real(kind=4):: gamma1, gamma2, gamma3, gamma4
    
    real (kind=4), allocatable:: O2_conc(:), T(:), d2Tdt(:), d2O2dt(:), d2wpdt(:)
    real (kind=4), allocatable::  T_data(:), O2_data(:)
    real (kind=4), allocatable:: Sp_norm(:), Fd_norm(:), wp(:),  Aai(:, :)
    real (kind=4), allocatable:: wp_interp(:, :), Total_Sp(:), Total_Fd(:)
    real (kind=4), allocatable:: Total_Spi(:, :), Total_Fdi(:, :)
    real (kind = 4), allocatable::  gamma_zi(:), massicur(:), massinext(:), temp_arr(:, :)
    real (kind=4), allocatable:: Spi_norm(:, :), Fdi_norm(:, :), flux(:), conc(:)
    
    real (kind = 4), allocatable::  Spi_cur(:, :, :), Spi_next(:, :, :), time_z(:), d2timedt(:)
    real (kind = 4), allocatable::  d2CO3dt(:), d2CO3satdt(:), time_interp(:, :)
    real (kind = 4), allocatable::  SIGz(:), CO3z(:), CO3sat(:), CO3_data(:), CO3sat_data(:)
    real (kind = 4), allocatable:: d0(:), z(:), Zgrid(:),Zgrid1(:), Z_data1(:), Z_data2(:), Z_data3(:)
    real (kind = 4), allocatable::  wp0(:)
    integer, allocatable:: counttot(:), counti(:, :)
    common /cm1/ Cw, sigma, eta, n_data1, n_data2
 
    Cw=0.
    Cm=0.
    

    namelist /Inp1D/ nf, n_data1, n_data2, n_data3, filethetao, fileo2, fileco3satcalc, fileco3, root
    namelist /Inp2D/ H, n, ni, frac, vartime, roi, ai, bi,  &
                    Cri, Q10i, Cw, K_O2i, sigma, eta, teta, epsilon, eps_max, eps_min, &
                    time1, time2, delta_t1, delta_t2, Time_end, dtime, dtsave, nu, dstart, dend, dp, & 
                    Trefi, mode, modecaco3, mode_eps, mode_init, datamode, pM, &
                    Z_measur, Sp_measur, Fd_measur, flux, conc
    
    open(1111,file='Input1.nml',form='formatted')
    read(1111,Inp1D)

    allocate(k1(nf), k2(nf), k3(nf), k4(nf))
    allocate(pM(nf), pV(nf), masspi(nf), mai(nf))
    allocate(mlimi(nf), roi(nf), Vi(nf))
    allocate(ai(nf), bi(nf), Cri(nf),  Trefi(nf), Q10i(nf), K_O2i(nf))
    allocate(gamma_zi(nf))
    allocate(vartime(nf), flag(nf), flux(nf), conc(nf))
    allocate(countij(nf), frac(nf))
    allocate(Spi_curij(nf), Spi_curij1(nf), Spi_curij2(nf))
    
    
    
    open(1111,file='Input2.nml',form='formatted')
    read(1111,Inp2D)
    
    
    time = 0.
    
    allocate(massicur(nf), massinext(nf))
    allocate(SIGz(n+1))
    allocate(CO3z(n+1), CO3sat(n+1), CO3_data(n_data3), CO3sat_data(n_data3))
    
    allocate(Spi_next(nf, n+1, ni+1), Spi_cur(nf, n+1, ni+1))
    allocate(Zgrid(n+1), Zgrid1(n+1), wp_interp(n+1, ni+1), Aai(n+1, ni+1))
    allocate(Total_Sp(n+1), Total_Fd(n+1))
    allocate(Total_Spi(nf, n+1), Total_Fdi(nf, n+1))
    
     allocate(O2_data(n_data2), T_data(n_data1)) 
    allocate(O2_conc(n+1), T(n+1), d2Tdt(n_data1), d2O2dt(n_data2))
    allocate(d2CO3dt(n_data3), d2CO3satdt(n_data3))
    
    allocate(d2wpdt(int(40*n)), z(int(40*n)),wp(int(40*n)))  
    allocate(Sp_norm(n+1), Fd_norm(n+1), Spi_norm(nf, n+1), Fdi_norm(nf, n+1))
    allocate(Z_data1(n_data1), Z_data2(n_data2), Z_data3(n_data3)) 
    allocate(counttot(ni+1), counti(nf, ni+1), wp0(ni+1))
     
    allocate(time_z(int(40*n)), d2timedt(int(40*n)), time_interp(n+1, ni+1))
    allocate(d0(ni+1))

    dz = (H-Z_measur)/n
    
    dd0 = (dend-dstart)/ni
    ntime = floor(Time_end/dtime) + 1
    time3= 0
    
    d0 = 0.     !particle size spectrum from d_min to d_max
    
    do i = 1, ni+1
        d0(i) = dstart + dd0*(i-1)
    enddo    
    
    pi = 3.14159265358979323846

    do i = 1, nf, 1
    flag(i) = 1.
    enddo
    flagtot = 1.
    Cm=0.
time = 0.

path = trim(root) !output path

    
    dz = (H-Z_measur)/n
    dd0 = (dend-dstart)/ni
    
    d0 = 0.     !particle size spectrum from d_min to d_max
    
    do i = 1, ni+1
        d0(i) = dstart + dd0*(i-1)
    enddo
        
    do i = 1, nf, 1
        Vi(i) = 0.
    enddo
    if (datamode.eq.3) then
            totflux = 0.
            do i = 1, nf, 1
                totflux = totflux + flux(i)
            enddo
            
            do i = 1, nf, 1
                pM(i) = flux(i)/totflux !component mass fraction
            enddo
     elseif (datamode.eq.2) then
     
            totconc = 0.
            do i = 1, nf, 1
                totconc = totconc + conc(i)
            enddo
            
            do i = 1, nf, 1
                pM(i) = conc(i)/totconc
            enddo
     endif
     
    romid = 0.
    do jj = 1, nf, 1
        romid = romid + pM(jj)/roi(jj)
    enddo
    romid = 1./romid
    
    do jj = 1, nf, 1
        pV(jj) = pM(jj)*romid/roi(jj) !component volume fraction
    enddo
    rop = 0.
    do i = 1, nf, 1
       masspi(i) = pi*roi(i)*dp**3/6. !mass of the primary particle
       rop = rop + pV(i)*roi(i)
    enddo
                
    dlim = (0.01*dstart) !diameter limit
    if (dlim.lt.dp) then
        dlim = dp  !diameter limit cannot be less then diameter of primary particle
    endif
    
    do i = 1, nf, 1
        mlimi(i) = pV(i)*masspi(i)*(dlim/dp)**sigma !component mass limit 
    enddo

    null = 0
    row = 1027.
    g = 9.8
    Zgrid = 0. !regular Z coordinate, m
    dcur = 0.     !diameter d_i of some size class i, m
    dnext = 0.     !diameter d_i of some size class i, m
     
    do i = 1, nf, 1
        massicur(i) = 0.  !mass of the aggregate component
        massinext(i) = 0.       
    enddo
    masscur = 0.
    massnext = 0.
    z = 0.     !irregular z-coordinate z*, m
    wp = 0.    !particle settling velocity, wp, m/s
    
    
    Spi_cur = 0.    !Particulate organic matter concentration of some size class i, Sp_i, kg/m^3
    Spi_next = 0.
    
    do i = 1, nf, 1
        gamma_zi(i) = 0. !degradation rate gamma, 1/s
    enddo
    
    !Normalized POM concentration Sp and flux Fd (kg/(m^2*s))
    Sp_norm = 0.
    Fd_norm = 0
    
    do i = 1, nf, 1
        Spi_norm(i, :) = 0. !Total POM concentration Sp, kg/m^3
        Fdi_norm(i, :) = 0. !Total POM flux Fd, kg/(m^2*s)
    enddo
    Zgrid = 0. !regular Z coordinate, m
    
    O2_conc = 0. !Oxygen concentration, kg/m^3 
    T = 0.       !Temperature, C
    !Measurements data
    O2_data = 0.
    T_data = 0.
    !Derivatives
    d2Tdt = 0. 
    d2O2dt = 0.
    d2wpdt = 0.
    time_z = 0. !age of the aggregate
    d2timedt = 0.
    
    !interpolated over regular Z-grid values
    wp_interp = 0. 
    time_interp = 0.
    
    !Total POM concentration and flux (calculated using numerical integration)
    Total_Sp = 0. !Total POM concentration Sp, kg/m^3
    Total_Fd = 0. !Total POM flux Fd, kg/(m^2*s)
    
    do i = 1, nf, 1
        Total_Spi(i, :) = 0. !Total POM concentration Sp, kg/m^3
        Total_Fdi(i, :) = 0. !Total POM flux Fd, kg/(m^2*s)
    enddo
    
    counti = n+1
    counttot = n+1
    
    open(1, file = filethetao, status='old', form='formatted', action='read')

    do i = 1, n_data1
        read(1, *) Z_data1(i), T_data(i)
    enddo
    close(1)

    open(1, file = fileo2, status='old', form='formatted', action='read')
    do i = 1, n_data2
        read(1, *) Z_data2(i), O2_data(i)
    enddo
    close(1)
   
    if (modecaco3.eq.1) then   
        open(1, file = fileco3satcalc, status='old', form='formatted', action='read')
        
        do i = 1, n_data3
            read(1, *) Z_data3(i), CO3sat_data(i)
        enddo
        close(1)  
        
        open(1, file = fileco3, status='old', form='formatted', action='read')
        do i = 1, n_data3
            read(1, *) Z_data3(i), CO3_data(i)
        enddo
        close(1)
    else
        CO3_data  = 2.
        CO3sat_data = 1.
        
    endif
    
    call spline(Z_data1, T_data, n_data1, null, null,d2Tdt)
    call spline(Z_data2, O2_data, n_data2, null, null,d2O2dt)
    call spline(Z_data3, CO3_data, n_data3, null, null,d2CO3dt)
    call spline(Z_data3, CO3sat_data, n_data3, null, null,d2CO3satdt)
    
   
    do i = 1, n+1
        Zgrid(i) = Z_measur+(i-1)*dz
        call splint(Z_data1, T_data, d2Tdt, n_data1, Zgrid(i), T(i))
        call splint(Z_data2, O2_data, d2O2dt, n_data2, Zgrid(i), O2_conc(i))
        
        if (modecaco3.eq.1) then                 
            call splint(Z_data3, CO3_data, d2CO3dt, n_data3, Zgrid(i), CO3z(i))
            call splint(Z_data3, CO3sat_data, d2CO3satdt, n_data3, Zgrid(i), CO3sat(i))
        else
            CO3z = 2.
            CO3sat = 1.
            SIG = 2.
            
        endif
    enddo
    
    do i = 1, n+1
        Zgrid1(i) = (i-1)*dz
    enddo
    
    rop = 0.
    do i = 1, nf, 1
       masspi(i) = pi*roi(i)*dp**3/6. !mass of the primary particle
       rop = rop + pV(i)*roi(i) !average density of the primary particle
    enddo
    massp = pi*rop*dp**3/6.
    
    Aai = 0.
    z = 0.  

    
 do j = 1, ni+1, 1
 
    do ii = 1,nf, 1
        mai(ii) = 0.
    enddo
    ma = 0.
    
    z = 0.
    time_z = 0.
    wp(:) = 0.
    do ii = 1,nf, 1
        gamma_zi(ii) = 0.
        massicur(ii) = 0. 
        massinext(ii) = 0. 
    enddo
    
    dcur = 0. 
    z(1) = Zgrid(1)
    dcur = d0(j)
    dnext = 0. 
    masscur = 0.
    massnext = 0.
    rop = 0.
    time = 0.
    do ii = 1,nf, 1
        massicur(ii) = pV(ii)*masspi(ii)*(dcur/dp)**sigma
        masscur = masscur + massicur(ii)
        rop =  rop + pV(ii)*roi(ii)
    enddo
    droa = (rop - row)*(dp/dcur)**(3.-sigma)
    
    if (Q10i(1).gt.1) then
        nu = 5.*10**(-10.)*exp(2250./(T(1)+273.15)) !viscosity
    endif
    
    if (mode.eq.0) then
        if (j.eq.1) then
            wp(1) = (rop - row)*g*(dp**(3.-sigma))*(dcur**(sigma-1.))/18./row/nu
            wp0(1) = wp(1)
        Recur = 0.
        else 
            Re = wp0(j-1)*d0(j-1)/nu
            Recur = Re
            Cd = 24./Re + 6./(1+Re**0.5) + 0.4
            wp(1) = ((4*droa*g*dcur)/(3.*row*Cd))**0.5
            wp0(j) = wp(1)
        endif    
        Re =  wp(1)*dcur/nu 
        iRe = 0
        
        do while (abs(Recur - Re)/Recur.gt.0.0)   
            Recur = Re
            Cd = 24./Re + 6./(1+Re**0.5) + 0.4
            wp(1) = ((4*droa*g*dcur)/(3.*row*Cd))**0.5
            Re = wp(1)*dcur/nu
            iRe = iRe+1
        enddo
  
        Cd = 24./Re + 6./(1+Re**0.5) + 0.4
        
    elseif (mode.eq.1) then
        Cw = (rop - row)*g*(dp**(3.-sigma))/18./row/nu
        eta = (sigma-1.)
        wp(1) =Cw*dcur**eta 
    else
        wp(1) =Cw*dcur**eta 
    endif
    
        wp0(j) = wp(1)

    
    i = 0
    call splint(Z_data1, T_data, d2Tdt, n_data1, z(1), T1)
    call splint(Z_data2, O2_data, d2O2dt, n_data2, z(1), O1)
    
    if (modecaco3.eq.1) then       
        call splint(Z_data3, CO3_data, d2CO3dt, n_data3, z(1), CO3z1)
        call splint(Z_data3, CO3sat_data, d2CO3satdt, n_data3, z(1), CO3sat1)
        SIG = CO3z1/CO3sat1
    else
        SIG = 2.
    endif
    
    do ii = 1, nf, 1
        call gamma_func(frac(ii), vartime(ii), ai(ii), bi(ii), Cri(ii),  Trefi(ii), Q10i(ii), K_O2i(ii), &
         SIG, teta, time, T1, O1, gamma_zi(ii))
        
        if ((massicur(ii)-dz*gamma_zi(ii)*massicur(ii)/wp(1)).gt.0) then
            dt = dz/wp(1)
        else
            dt = int(massicur(ii)/gamma_zi(ii)/massicur(ii)/10.)
        endif
    enddo
        
    do ii = 1, nf, 1 
        counti(ii, j) = n+1
    enddo
    counttot(j) = n+1

do ii = 1, nf, 1 
    flag(ii) = 1.
enddo
flagtot = 1.
flag4 = 0
tflag = 1
ii = 1

time_z(1) = time
i = 0
do while ((z(i+1).le.H))
    i = i+1   
    call splint(Z_data1, T_data, d2Tdt, n_data1,z(i), T1)
    call splint(Z_data2, O2_data, d2O2dt, n_data2, z(i), O1)
    if (modecaco3.eq.1) then       
        call splint(Z_data3, CO3_data, d2CO3dt, n_data3, z(i), CO3z1)
        call splint(Z_data3, CO3sat_data, d2CO3satdt, n_data3, z(i), CO3sat1)

        SIG = CO3z1/CO3sat1
    else
        SIG = 2.
    endif
    
    V = 0.
    rop = 0.
    do jj = 1, nf, 1
        call gamma_func(frac(jj), vartime(jj), ai(jj), bi(jj), Cri(jj),  Trefi(jj), Q10i(jj), K_O2i(jj), & 
         SIG, teta, time, T1, O1, gamma_zi(jj))
       
        k1(jj) = -gamma_zi(jj)*massicur(jj)
        Vi(jj) = massicur(jj)/roi(jj)
        V = V + Vi(jj)
        rop = rop + roi(jj)*Vi(jj)
    enddo
    rop = rop/V
    da = dcur
    m1 = wp(i) 
    
    call splint(Z_data1, T_data, d2Tdt, n_data1, z(i)+dt*m1*0.5, T2)
    call splint(Z_data2, O2_data, d2O2dt, n_data2, z(i)+dt*m1*0.5, O2)
    
    
    if (modecaco3.eq.1) then               
        call splint(Z_data3, CO3_data, d2CO3dt, n_data3, z(i)+dt*m1*0.5, CO3z2)
        call splint(Z_data3, CO3sat_data, d2CO3satdt, n_data3, z(i)+dt*m1*0.5, CO3sat2)
        SIG = CO3z2/CO3sat2
        
    else
        SIG = 2.
    endif
    
    V = 0.
    rop = 0.
    ma = 0.
    do jj = 1, nf, 1
        call gamma_func(frac(jj), vartime(jj), ai(jj), bi(jj), Cri(jj), Trefi(jj), Q10i(jj), K_O2i(jj), SIG, & 
        teta, time+dt*0.5, T2, O2, gamma2)
        
        k2(jj) = -gamma2*(massicur(jj)+dt*k1(jj)*0.5)
        
        if (flag(jj).gt.0) then            
            mai(jj) = (massicur(jj)+dt*k1(jj)*0.5)
            Vi(jj) = (massicur(jj)+dt*k1(jj)*0.5)/roi(jj)
        else
            Vi(jj) = 0.
        endif
       
        V = V + Vi(jj)
        rop = rop + roi(jj)*Vi(jj)
        ma = ma + mai(jj)
    enddo
    
    rop = rop/V
    da =  (ma/massp)**(1/sigma)*dp 
    droa = (rop - row)*(dp/da)**(3.-sigma)
    
    if (da.eq.0) then
        m2 = 0.
    else
        if (mode.eq.0) then
            m2 = ((4*droa*g*da)/(3.*row*Cd))**0.5  
        else
            m2 = Cw*da**eta
        endif
    endif
    call splint(Z_data1, T_data, d2Tdt, n_data1, z(i)+dt*m2*0.5, T3)
    call splint(Z_data2, O2_data, d2O2dt, n_data2, z(i)+dt*m2*0.5, O3)
      
    if (modecaco3.eq.1) then   
        call splint(Z_data3, CO3_data, d2CO3dt, n_data3, z(i)+dt*m2*0.5, CO3z3)
        call splint(Z_data3, CO3sat_data, d2CO3satdt, n_data3, z(i)+dt*m2*0.5, CO3sat3)

        SIG = CO3z3/CO3sat3
        
    else
        SIG = 2.
    endif
    V = 0.
    rop = 0.
    ma = 0.
    do jj = 1, nf, 1
        call gamma_func(frac(jj), vartime(jj), ai(jj), bi(jj), Cri(jj),  Trefi(jj), Q10i(jj), K_O2i(jj), SIG, &
        teta, time+dt*0.5, T3, O3, gamma3)      
        
        k3(jj) = -gamma3*(massicur(jj)+dt*k2(jj)*0.5)
    
        if (flag(jj).gt.0) then
            mai(jj) = (massicur(jj)+dt*k2(jj)*0.5)
            Vi(jj) = (massicur(jj)+dt*k2(jj)*0.5)/roi(jj)
        else
            Vi(jj) = 0.
        endif
        
        V = V + Vi(jj)
        rop = rop + roi(jj)*Vi(jj)
        ma = ma + mai(jj)
    enddo
    
    rop = rop/V
    da =  (ma/massp)**(1/sigma)*dp 
    droa = (rop - row)*(dp/da)**(3.-sigma)
    
    if (da.eq.0) then
        m3 = 0.
    else
        if (mode.eq.0) then
            m3 = ((4*droa*g*da)/(3.*row*Cd))**0.5  
        else
            m3 = Cw*da**eta
        endif 
    endif
    
    call splint(Z_data1, T_data, d2Tdt, n_data1, z(i)+dt*m3, T4)
    call splint(Z_data2, O2_data, d2O2dt, n_data2, z(i)+dt*m3, O4)
    
     
    if (modecaco3.eq.1) then 
        call splint(Z_data3, CO3_data, d2CO3dt, n_data3, z(i)+dt*m3, CO3z4)
        call splint(Z_data3, CO3sat_data, d2CO3satdt, n_data3, z(i)+dt*m3, CO3sat4)
        SIG = CO3z4/CO3sat4
    else
        SIG = 2.
    endif
  
   
    V = 0.
    rop = 0.
    ma = 0.
    do jj = 1, nf, 1
        call gamma_func(frac(jj), vartime(jj), ai(jj), bi(jj), Cri(jj),  Trefi(jj), Q10i(jj), K_O2i(jj), SIG, & 
        teta, time+dt, T4, O4, gamma4)      
     
        k4(jj) = -gamma4*(massicur(jj)+dt*k3(jj))
        if (flag(jj).gt.0) then
            mai(jj) = (massicur(jj)+dt*k3(jj))
            Vi(jj) = (massicur(jj)+dt*k3(jj))/roi(jj)
        else
            Vi(jj) = 0.
        endif
        
        V = V + Vi(jj)
        rop = rop + roi(jj)*Vi(jj)
        ma = ma + mai(jj)
    enddo
    
    rop = rop/V
    da =  (ma/massp)**(1/sigma)*dp 
    droa = (rop - row)*(dp/da)**(3.-sigma)

    if (da.eq.0) then
        m4 = 0.
    else
        if (mode.eq.0) then
            m4 = ((4*droa*g*da)/(3.*row*Cd))**0.5  
        else
            m4 = Cw*da**eta
        endif 
    endif
    V = 0. 
    rop = 0.
    massnext = 0. 
    do jj = 1, nf, 1
        if (flag(jj).gt.0) then
            massinext(jj) = massicur(jj) + dt*(k1(jj) + 2.*k2(jj) + 2.*k3(jj) + k4(jj))/6.
        else
            massinext(jj) = 0.
        endif
    
        if (massinext(jj).lt.mlimi(jj).and.flag(jj).ne.0.) then
            flag(jj) = 0.
            massinext(jj) = 0.
            counti(jj, j) = i
        endif
        Vi(jj) = massinext(jj)/roi(jj)
        V = V + Vi(jj)
        rop = rop + roi(jj)*Vi(jj)
        massnext = massnext + massinext(jj)
    enddo
  
    rop = rop/V
    dnext = (massnext/massp)**(1/sigma)*dp
    droa = (rop - row)*(dp/dnext)**(3.-sigma)
    time = time + dt
   
    time_z(i+1) = time
    
    if (dnext.lt.dlim) then
        flagtot = 0.
        dnext = 0.
        counttot(j) = i+1
        exit
    endif


    z(i+1) = z(i) + dt*(m1 + 2.*m2 + 2.*m3 + m4)/6.
    
    
    call splint(Z_data1, T_data, d2Tdt, n_data1,z(i+1), T1)
    call splint(Z_data2, O2_data, d2O2dt, n_data2, z(i+1), O1)
    
     
    if (modecaco3.eq.1) then
        call splint(Z_data3, CO3_data, d2CO3dt, n_data3, z(i+1), CO3z1)
        call splint(Z_data3, CO3sat_data, d2CO3satdt, n_data3, z(i+1), CO3sat1)
        SIG = CO3z1/CO3sat1
    else
        SIG = 2.
    endif
    

    do jj = 1, nf, 1
        call gamma_func(frac(jj), vartime(jj), ai(jj), bi(jj), Cri(jj),  Trefi(jj), Q10i(jj), K_O2i(jj), SIG, & 
        teta, time, T4, O4, gamma_zi(jj))   
         
    enddo
       
    if (mode.eq.0) then
        Recur = Re
          wp(i+1) =((4*droa*g*dnext)/(3.*row*Cd))**0.5 

        call splint(Z_data1, T_data, d2Tdt, n_data1, z(i+1), Tz)
        
        if (Q10i(1).gt.1) then
            nu = 5.*10**(-10.)*exp(2250./(Tz+273.15))
        endif
        Re = wp(i+1)*dnext/nu
        iRe = 0
        
        do while (abs(Recur - Re)/Recur.gt.0.0)   
            Recur = Re
            Cd = 24./Re + 6./(1+Re**0.5) + 0.4
            wp(i+1) = ((4*droa*g*dnext)/(3.*row*Cd))**0.5
            Re = wp(i+1)*dnext/nu
            iRe = iRe+1
        enddo
            
        Cd = 24./Re + 6./(1+Re**0.5) + 0.4  
    else
        wp(i+1) = Cw*dnext**eta
    endif
      
    if (ABS(z(i+1)-z(i)).lt.0.01*dz) then
        flagtot = 0.
        z(i+1) = 0. 
        wp(i+1) = 0. 
        dnext = 0.
        counttot(j) = i
        exit
    endif

    do jj = 1, nf, 1
        if ((massinext(jj)-dz*gamma_zi(jj)*massinext(jj)/wp(i+1)).le.0) then
            tflag = 0
        endif
    enddo
    
    if (tflag.ne.0) then
        dt = dz/wp(i+1)
    endif
    
    do jj = 1, nf, 1
        massicur(jj) = massinext(jj)
    enddo
    
    masscur = massnext
    dcur = dnext
enddo

if (ieee_is_nan(z(i+1)).or.z(i+1).eq.0) then
    m = i
else
    m = i+1
endif

time = 0 

if (z(m).gt.Zgrid(n+1)) then
    counttot(j) = n+1
else
    do ii = 1, n+1, 1
        if (Zgrid(ii).gt.z(m)) then
            counttot(j)  = ii-1
            exit
        endif
    enddo 
endif

do jj = 1, nf, 1
    if (flag(jj).eq.1) then
        counti(jj, j) = counttot(j)
    else
        if (counti(jj, j).gt.m) then
            counti(jj, j) = counttot(j)
        else
            do ii = 1, n+1, 1
                if (Zgrid(ii).gt.z(counti(jj, j))) then
                    counti(jj, j)  = ii-1
                    exit
                endif
            enddo 
        endif 
    endif
enddo


if (m.gt.1) then 
    call spline(z(:), wp(:), m, null, null,d2wpdt(:))
    call spline(z(:), time_z(:), m, null, null, d2timedt(:))
endif


do i = 1, counttot(j), 1

    if (m.eq.1) then
        wp_interp(1, j) = wp(i)
        Aai(1, j) = wp(i)
        time_interp(i, j) = time_z(i)
    else
        call splint(z(:), wp(:), d2wpdt(:), m, Zgrid(i), wp_interp(i, j))
        call splint(z(:), wp(:), d2wpdt(:), m, Zgrid(i) + dz/2, Aai(i, j))
        call splint(z(:), time_z(:), d2timedt(:), m, Zgrid(i), time_interp(i, j))
    endif
    
    if (Zgrid(i).gt.z(m)) then
        exit
    endif
enddo

enddo
 

deallocate(wp, d2wpdt, z, time_z, d2timedt)
do i = 1, n+1
    Zgrid(i) = Z_measur+(i-1)*dz
enddo    
do i = 1, n+1
    Zgrid1(i) = (i-1)*dz
enddo   


ttime = 0.
Flux0 = 0
!M0xCm = M0xCm0


if (mode_eps.eq.1) then
    if ((ttime/86400).le.time1) then
        epsilon = eps_min
    else if ((ttime/86400).le.time2) then
        epsilon = (eps_max-eps_min)*ttime/(time2-time1)/86400 + eps_min - time1*(eps_max-eps_min)/(time2-time1)
    else 
        epsilon = eps_max
    endif
endif

if (mode_eps.eq.2) then
    time2 = time1+delta_t1
    time3 = time1+delta_t2
    if ((ttime/86400).le.time1) then
        epsilon = eps_max
        
    else if ((ttime/86400).le.time2) then
        epsilon = (eps_min-eps_max)*ttime/(time2-time1)/86400 + eps_max- time1*(eps_min-eps_max)/(time2-time1)
    else if  (((ttime/86400).gt.time2).and.((ttime/86400).le.time3))  then
        epsilon = (eps_max-eps_min)*ttime/(time3-time2)/86400 + eps_min - time2*(eps_max-eps_min)/(time3-time2)
    else
        epsilon = eps_max
    endif
endif
Total_Sp = 0.
Total_Fd = 0.
Total_Spi = 0.
Total_Fdi = 0.

   do j = 1, ni+1, 1
        do jj = 1, nf, 1
            Spi_cur(jj, 1, j) =  dd0*(pi/6.)*roi(jj)*pV(jj)*dp**(3.-sigma)*d0(j)**(sigma-epsilon)
        enddo
    enddo

      
     do j = 1, ni+1, 1
        do jj = 1, nf, 1
            if (j.eq.1.or.j.eq.ni+1) then
                    Total_Spi(jj, 1) = Total_Spi(jj, 1) + Spi_cur(jj, 1, j)
                    Total_Fdi(jj, 1) = Total_Fdi(jj, 1) + Spi_cur(jj, 1, j)*wp_interp(1, j) 
            else 
                if (mod(j,2).eq.0) then
                    Total_Spi(jj, 1) = Total_Spi(jj, 1) + 4*Spi_cur(jj, 1, j)
                    Total_Fdi(jj, 1) = Total_Fdi(jj, 1) + 4*Spi_cur(jj, 1, j)*wp_interp(1, j) 
                else
                    Total_Spi(jj, 1) = Total_Spi(jj, 1) + 2*Spi_cur(jj, 1, j)
                    Total_Fdi(jj, 1) = Total_Fdi(jj, 1) + 2*Spi_cur(jj, 1, j)*wp_interp(1, j) 
                endif
            endif
            
        enddo
     enddo
     
        do jj = 1, nf, 1
            Total_Spi(jj, 1) = Total_Spi(jj, 1)/3.
            Total_Fdi(jj, 1) = Total_Fdi(jj, 1)/3.
        enddo
     
      
        do jj = 1, nf, 1
            Total_Sp(1) = Total_Sp(1) +  Total_Spi(jj, 1)
            Total_Fd(1) = Total_Fd(1) +  Total_Fdi(jj, 1)
        enddo

   

         if (datamode.eq.3) then
            M0xCm = totflux/Total_Fd(1) 
         elseif (datamode.eq.2) then
            M0xCm = totconc/Total_Sp(1)
         elseif (datamode.eq.0.) then
             M0xCm = Sp_measur/Total_Sp(1)
         elseif (datamode.eq.1) then
             M0xCm = Fd_measur/Total_Fd(1)
         endif
                        
        do i = 1, n+1, 1 
            Total_Sp(i) = Total_Sp(i)*M0xCm
            Total_Fd(i) = Total_Fd(i)*M0xCm
        enddo
            
        do jj = 1, nf, 1   
            do i = 1, n+1, 1  
                Total_Spi(jj, i) = Total_Spi(jj, i)*M0xCm
                Total_Fdi(jj, i) = Total_Fdi(jj, i)*M0xCm
            enddo
        enddo

    do j= 1, ni+1, 1 
        do jj = 1, nf, 1
            Spi_cur(jj, 1, j) = Spi_cur(jj, 1, j)*M0xCm
        enddo
    enddo  
if ((nf.eq.1).and.(mode_init.eq.1)) then   
    allocate(temp_arr(n+1, ni+1))
    temp_arr = 0. 
  
  !analytical solution can be used as a starting profile for 1-component model  
  call solution(vartime(1), epsilon, roi(1), pV(1), dp, ai(1), bi(1), Cri(1), n, ni, M0xCm, dd0, Zgrid1, d0, temp_arr, Total_Sp, Total_Fd)
   
   Spi_cur(1, :, :) = temp_arr
   Total_Spi(1, :) =  Total_Sp
   Total_Fdi(1, :) = Total_Fd
   deallocate(temp_arr)
    
endif

    do i = 1, n+1, 1
        Sp_norm(i) = Total_Sp(i)/Total_Sp(1)
        Fd_norm(i) = Total_Fd(i)/Total_Fd(1)
        
        do jj = 1, nf, 1  
            Spi_norm(jj, i) = Total_Spi(jj, i)/Total_Spi(jj, 1)
            Fdi_norm(jj, i) = Total_Fdi(jj, i)/Total_Fdi(jj, 1)
        enddo
    enddo   
  
    Flux0 = Total_Fd(1)

print *,"time =", int(ttime/86400), "days"

        
do jj = 1, nf, 1
    write(nfrac, '(I0)') jj
    call write_arr(Zgrid1(:),Total_Spi(jj, :),n+1,path,'Sp'//trim(nfrac)//'(z)_',7,int(ttime/86400))
    call write_arr(Zgrid1(:),Total_Fdi(jj, :),n+1,path,'Fd'//trim(nfrac)//'(z)_',7,int(ttime/86400))
    call write_arr(Zgrid1(:),Spi_norm(jj,:),n+1,path,'Sp'//trim(nfrac)//'(z)_norm_',12,int(ttime/86400))
    call write_arr(Zgrid1(:),Fdi_norm(jj, :),n+1,path,'Fd'//trim(nfrac)//'(z)_norm_',12,int(ttime/86400))
enddo    
    call write_arr(Zgrid1(:),Total_Sp(:),n+1,path,'Sp(z)_',6,int(ttime/86400))
    call write_arr(Zgrid1(:),Total_Fd(:),n+1,path,'Fd(z)_',6,int(ttime/86400))
    call write_arr(Zgrid1(:),Sp_norm(:),n+1,path,'Sp(z)_norm_',11,int(ttime/86400))
    call write_arr(Zgrid1(:),Fd_norm(:),n+1,path,'Fd(z)_norm_',11,int(ttime/86400))
            

dt2 = dtime**2
dz2 = 2.*dz
dz22 = 2.*dz**2
do  it = 2, ntime, 1
    ttime = (it-1)*dtime
        
    Total_Sp = 0.
    Total_Fd = 0.
        
    Total_Spi = 0.
    Total_Fdi = 0.

    if (mode_eps.eq.1) then
        if ((ttime/86400).le.time1) then
            epsilon = eps_min
        else if ((ttime/86400).le.time2) then
            epsilon = (eps_max-eps_min)*ttime/(time2-time1)/86400 + eps_min - time1*(eps_max-eps_min)/(time2-time1)
        else 
            epsilon = eps_max
        endif
    endif

    if (mode_eps.eq.2) then
        time2 = time1+delta_t1
        time3 = time1+delta_t2
        if ((ttime/86400).le.time1) then
            epsilon = eps_max
            
        else if ((ttime/86400).le.time2) then
            epsilon = (eps_min-eps_max)*ttime/(time2-time1)/86400 + eps_max - time1*(eps_min-eps_max)/(time2-time1)
        else if  (((ttime/86400).gt.time2).and.((ttime/86400).le.time3))  then
            epsilon = (eps_max-eps_min)*ttime/(time3-time2)/86400 + eps_min - time2*(eps_max-eps_min)/(time3-time2)
        else
            epsilon = eps_max
        endif
    endif
    
    do j = 1, ni+1, 1
        do jj = 1, nf, 1
            Spi_next(jj, 1, j) = M0xCm*dd0*(pi/6.)*roi(jj)*pV(jj)*dp**(3.-sigma)*d0(j)**(sigma-epsilon)
        enddo
        
    wp11 = wp_interp(1, j)
        do jj = 1, nf, 1
            if (j.eq.1.or.j.eq.ni+1) then
                Total_Spi(jj, 1) = Total_Spi(jj, 1) + Spi_next(jj, 1, j)
                Total_Fdi(jj, 1) = Total_Fdi(jj, 1) + Spi_next(jj, 1, j)*wp11
            else 
                if (mod(j,2).eq.0) then
                    Total_Spi(jj, 1) = Total_Spi(jj, 1) + 4*Spi_next(jj, 1, j)
                    Total_Fdi(jj, 1) = Total_Fdi(jj, 1) + 4*Spi_next(jj, 1, j)*wp11
                else
                    Total_Spi(jj, 1) = Total_Spi(jj, 1) + 2*Spi_next(jj, 1, j)
                    Total_Fdi(jj, 1) = Total_Fdi(jj, 1) + 2*Spi_next(jj, 1, j)*wp11
                endif
            endif
        enddo
        countzj = counttot(j)
        
        
        do jj = 1, nf, 1
            countij(jj) = counti(jj, j)
        enddo
        
        do i = 2, countzj-1, 1
            wpi = wp_interp(i, j)
            wpi1 = wp_interp(i-1, j)
            wpi2 = wp_interp(i+1, j)
            Aij = Aai(i, j)
            Aij1 = Aai(i-1, j)
            
            do jj = 1, nf, 1
                Spi_curij(jj) = Spi_cur(jj, i, j)
                Spi_curij1(jj) = Spi_cur(jj, i-1, j)
                Spi_curij2(jj) = Spi_cur(jj, i+1, j)
            enddo      
            
            Aijwpi2 = Aij*wpi2
            Aijwpi = Aij*wpi
            Aij1wpi = Aij1*wpi
            Aij1wpi1=Aij1*wpi1
        
            do jj = 1, nf, 1
                if (i.le.countij(jj)) then
                        
                    if (modecaco3.eq.1) then
                        SIG = CO3z(i)/CO3sat(i)
                    else
                        SIG = 2.
                    endif
                    
                    call gamma_func(frac(jj), vartime(jj), ai(jj), bi(jj), Cri(jj),  Trefi(jj), Q10i(jj), K_O2i(jj), SIG, &
                    teta, time_interp(i, j), T(i), O2_conc(i), gammai)
                    
                    Spn = Spi_curij(jj) - dtime*(wpi2*Spi_curij2(jj)-wpi1*Spi_curij1(jj))/dz2
                    Spn = Spn + dt2*((Aijwpi2*Spi_curij2(jj)-Aijwpi*Spi_curij(jj))-(Aij1wpi*Spi_curij(jj)-Aij1wpi1*Spi_curij1(jj)))/dz22
                    Spn = Spn - dtime*gammai*Spi_curij(jj)
                    Spi_next(jj, i, j) = Spn 
                endif
                            
        
                if (mod(ttime, dtsave).eq.0) then 
                    if (j.eq.1.or.j.eq.ni+1) then
                            Total_Spi(jj, i) = Total_Spi(jj, i) + Spi_next(jj, i, j)
                            Total_Fdi(jj, i) = Total_Fdi(jj, i) + Spi_next(jj, i, j)*wpi
                    else 
                        if (mod(j,2).eq.0) then
                            Total_Spi(jj, i) = Total_Spi(jj, i) + 4*Spi_next(jj, i, j)
                            Total_Fdi(jj, i) = Total_Fdi(jj, i) + 4*Spi_next(jj, i, j)*wpi
                        else
                            Total_Spi(jj, i) = Total_Spi(jj, i) + 2*Spi_next(jj, i, j)
                            Total_Fdi(jj, i) = Total_Fdi(jj, i) + 2*Spi_next(jj, i, j)*wpi
                        endif 
                    endif
                endif
          enddo
        enddo
        wpn = wp_interp(countzj, j)
        wpn1 = wp_interp(countzj-1, j)
        
        
        do jj = 1, nf, 1
            if (countzj.le.countij(jj)) then
                        
                if (modecaco3.eq.1) then
                    SIG = CO3z(countzj)/CO3sat(countzj)
                else
                    SIG = 2.
                endif
                
                call gamma_func(frac(jj), vartime(jj), ai(jj), bi(jj), Cri(jj),  Trefi(jj), Q10i(jj), K_O2i(jj), SIG, &
                teta, time_interp(countzj, j), T(countzj), O2_conc(countzj), gammai)
                    
                Spn = Spi_cur(jj, countzj, j) - dtime*(wpn*Spi_cur(jj, countzj, j)-wpn1*Spi_cur(jj, countzj-1, j))/dz 
                Spn = Spn - dtime*gammai*Spi_cur(jj, countzj, j)
                Spi_next(jj, countzj, j) = Spn
      
                if (mod(ttime, dtsave).eq.0) then 
                    if (j.eq.1.or.j.eq.ni+1) then
                            Total_Spi(jj, countzj) = Total_Spi(jj, countzj) + Spi_next(jj, countzj, j)
                            Total_Fdi(jj, countzj) = Total_Fdi(jj, countzj) + Spi_next(jj, countzj, j)*wpn
                    else 
                        if (mod(j,2).eq.0) then
                            Total_Spi(jj, countzj) = Total_Spi(jj, countzj) + 4*Spi_next(jj, countzj, j)
                            Total_Fdi(jj, countzj) = Total_Fdi(jj, countzj) + 4*Spi_next(jj, countzj, j)*wpn
                        else
                            Total_Spi(jj, countzj) = Total_Spi(jj, countzj) + 2*Spi_next(jj, countzj, j)
                            Total_Fdi(jj, countzj) = Total_Fdi(jj, countzj) + 2*Spi_next(jj, countzj, j)*wpn
                        endif 
                    endif     
                endif 
            endif
        enddo
        
            
    enddo  
    
    Spi_cur = Spi_next
    
    Total_Sp = 0. 
    Total_Fd = 0.
    
    
    if (mod(ttime, dtsave).eq.0) then  
        do jj = 1, nf, 1     
            do i = 1, n+1, 1  
                Total_Spi(jj, i) = Total_Spi(jj, i)/3.
                Total_Fdi(jj, i) = Total_Fdi(jj, i)/3.
                Total_Sp(i) = Total_Sp(i) + Total_Spi(jj, i)
                Total_Fd(i) = Total_Fd(i) + Total_Fdi(jj, i)
            enddo
        enddo
    else
        do jj = 1, nf, 1     
            Total_Spi(jj, 1) = Total_Spi(jj, 1)/3.
            Total_Fdi(jj, 1) = Total_Fdi(jj, 1)/3.
            
            Total_Sp(1) = Total_Sp(1) + Total_Spi(jj, 1)
            Total_Fd(1) = Total_Fd(1) + Total_Fdi(jj, 1)
        enddo
    endif 
    
  

        dM0xCm = Flux0/Total_Fd(1)
        M0xCm = M0xCm*dM0xCm
        
        if (mod(ttime, dtsave).eq.0) then 
        
            do i = 1, n+1, 1  
                Total_Sp(i) = Total_Sp(i)*dM0xCm
                Total_Fd(i) = Total_Fd(i)*dM0xCm
            enddo  
            
            do jj = 1, nf, 1     
                do i = 1, n+1, 1  
                    Total_Spi(jj, i) = Total_Spi(jj, i)*dM0xCm
                    Total_Fdi(jj, i) = Total_Fdi(jj, i)*dM0xCm
                enddo  
            enddo 
            
            do i = 1, n+1, 1
                Sp_norm(i) = Total_Sp(i)/Total_Sp(1)
                Fd_norm(i) = Total_Fd(i)/Total_Fd(1) 
                
                do jj = 1, nf, 1  
                    Spi_norm(jj, i) = Total_Spi(jj, i)/Total_Spi(jj, 1)
                    Fdi_norm(jj, i) = Total_Fdi(jj, i)/Total_Fdi(jj, 1)
                enddo
            enddo   
            
       print *,"time =", int(ttime/86400), "days"
               
            do jj = 1, nf, 1
                write(nfrac, '(I0)') jj
                call write_arr(Zgrid1(:),Total_Spi(jj, :),n+1,path,'Sp'//trim(nfrac)//'(z)_',7,int(ttime/86400))
                call write_arr(Zgrid1(:),Total_Fdi(jj, :),n+1,path,'Fd'//trim(nfrac)//'(z)_',7,int(ttime/86400))
                call write_arr(Zgrid1(:),Spi_norm(jj, :),n+1,path,'Sp'//trim(nfrac)//'(z)_norm_',12,int(ttime/86400))
                call write_arr(Zgrid1(:),Fdi_norm(jj, :),n+1,path,'Fd'//trim(nfrac)//'(z)_norm_',12,int(ttime/86400))
            enddo
                
            call write_arr(Zgrid1(:),Total_Sp(:),n+1,path,'Sp(z)_',6,int(ttime/86400))
            call write_arr(Zgrid1(:),Total_Fd(:),n+1,path,'Fd(z)_',6,int(ttime/86400))
            call write_arr(Zgrid1(:),Sp_norm(:),n+1,path,'Sp(z)_norm_',11,int(ttime/86400))
            call write_arr(Zgrid1(:),Fd_norm(:),n+1,path,'Fd(z)_norm_',11,int(ttime/86400))
        endif

    enddo
     
    
    deallocate(k1, k2, k3, k4)
    deallocate(pM, pV, masspi, mai)
    deallocate(mlimi, roi, Vi)
    deallocate(ai, bi, Cri,  Trefi, Q10i, K_O2i)
    deallocate(frac, gamma_zi)
    deallocate(vartime, flag, flux, conc)
    deallocate(countij, frac)
    deallocate(Spi_curij, Spi_curij1, Spi_curij2)
    deallocate(massicur, massinext)
    deallocate(SIGz)
    deallocate(CO3z, CO3sat, CO3_data, CO3sat_data)
    
    deallocate(Spi_next, Spi_cur)
    deallocate(Zgrid, Zgrid1, wp_interp, Aai)
    deallocate(Total_Sp, Total_Fd)
    deallocate(Total_Spi, Total_Fdi)
    
     deallocate(O2_data, T_data) 
    deallocate(O2_conc, T, d2Tdt, d2O2dt)
    deallocate(d2CO3dt, d2CO3satdt)
      
    deallocate(Sp_norm, Fd_norm, Spi_norm, Fdi_norm)
    deallocate(Z_data1, Z_data2, Z_data3) 
    deallocate(counttot, counti, wp0)
     
    deallocate(time_interp)
    deallocate(d0)


    end program EuLagPATS


    subroutine gamma_func(frac, vartime, a, b, Cr, Tref, Q10, K_O2, SIG, teta, time, T, O2_conc, gamma_z)
        real (kind = 4):: time, gamma_z, T, O2_conc 
        real (kind = 4):: a, b, Cr, Tref, Q10, Cw, K_O2, sigma, eta, teta, SIG 
        character(len = 100):: frac
        
        integer n_data1,n_data2, vartime
        common /cm1/ Cw, sigma, eta, n_data1, n_data2
        
        if (vartime.eq.0) then
            gamma_z = Cr
        else
            gamma_z = b/(a + time)
        endif
               
        gamma_z = gamma_z * (Q10**((T-Tref)/10.))
        if ((O2_conc.eq.0.).or.(K_O2.eq.0.)) then
            gamma_z = gamma_z
        else
            gamma_z = gamma_z * O2_conc/(K_O2 + O2_conc)
        endif
        
        if (trim(frac).eq."caco3") then
        if (SIG.lt.1.) then
            gamma_z = gamma_z*(1.-SIG)**teta
        else
            gamma_z = 0.
        endif
        endif
   
        
    end subroutine
    
    subroutine spline(x,y,n,yp1,ypn,y2)
        integer n,NMAX
        real(kind = 4)::x(n), yp1,ypn,y(n),y2(n)
        integer i,k
        real(kind = 4):: p,qn,sig,un,u(n)
        if (yp1.gt..99e30) then 
            y2(1)=0.
            u(1)=0.
        else 
            y2(1)=-0.5
            u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif
        
        do i=2,n-1 
            sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
            p=sig*y2(i-1)+2.
            y2(i)=(sig-1.)/p
            u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        enddo 
    
        if (ypn.gt..99e30) then
            qn=0.
            un=0.
        else 
            qn=0.5
            un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
        
        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
        do k=n-1,1,-1 
            y2(k)=y2(k)*y2(k+1)+u(k)
        enddo 
    end

   SUBROUTINE splint(xa,ya,y2a,n,x,y)
    INTEGER n
    REAL(kind = 4):: x,y,xa(n),y2a(n),ya(n)
    INTEGER k,khi,klo
    REAL(kind = 4)::  a,b,h
    klo=1 
    khi=n
    1 if (khi-klo.gt.1)then
    k=(khi+klo)/2
    if(xa(k).gt.x)then
    khi=k
    else
    klo=k
    endif
    goto 1
    endif 
    h=xa(khi)-xa(klo)
    if (h.eq.0.) pause 
    a=(xa(khi)-x)/h 
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

    return
    END
   
       
    subroutine write_arr(arr1,arr2,N,path,fname1,len,tlev)
    implicit none
    ! writes arr1 and arr2 as columns in file
    integer N,tlev,len
    real(kind = 4) arr1(N)
    real(kind = 4) arr2(N)
    character(len) fname1
    integer len1,curlen,i
    character*1000 fname, fname2, path
    character*10 str1,tstr
    
    str1= ''
    tstr = ''
    fname = ''
    len1=len_trim(fname1)
    fname(1:len1)=fname1(1:len1)
    curlen=len1
    write(str1,'(I10)')tlev
    tstr=ADJUSTL(str1)
    len1=len_trim(tstr)
    fname(curlen+1:curlen+len1)=tstr(1:len1)
    curlen=curlen+len1
    fname(curlen+1:curlen+4)='.dat'
    fname2 = trim(path)//trim(fname)    
    OPEN(2024,FILE=fname2)
    do i=1,N
    write(2024,*)arr1(i),arr2(i)
    end do
    !print *, "file is written ", fname2
    CLOSE(2024)
    end subroutine


!subroutine to calculate analytical solution for 1-fraction model; can be used as starticg profile
    subroutine solution(vartime, epsilon, ro, pV, dp, a, b, Cr, n, ni, M0xCm, dd0, Zgrid, d0, Sp_cur, Total_Sp, Total_Fd)
    implicit none
    
    integer:: n, ni, i, j, vartime,  n_data1, n_data2
    real (kind = 4)::  epsilon, a, b, Cr, Cw, sigma, eta, ro, pV, dp
    real (kind = 4):: Zgrid(n+1), d0(ni+1),  Sp_cur(n+1, ni+1), Total_Sp(n+1), Total_Fd(n+1)
    real (kind = 4):: a123,  M0xCm, dd0, fi, pi
    real (kind = 4), allocatable:: Fd_cur(:, :)
    common /cm1/ Cw, sigma, eta, n_data1, n_data2
    
    allocate(Fd_cur(n+1, ni+1))
    
    pi = 3.14159265358979323846
    Fd_cur = 0.    
    if (vartime.eq.0) then
        do i=1, n+1, 1
            Sp_cur(i, 1) = 0
            Total_Sp(i) = 0
                a123 = (sigma*Cw*(d0(1)**eta)/eta/Cr)
                if (Zgrid(i).lt.a123) then
                    Sp_cur(i, 1) =  M0xCm*dd0*(pi/6.)*ro*pV*dp**(3.-sigma)
                    Sp_cur(i, 1) = Sp_cur(i, 1)*(d0(1)**(sigma-epsilon))*(1-(eta*Cr*(Zgrid(i)))/sigma/Cw/(d0(1)**eta))**((sigma-eta)/eta) 
                    Total_Sp(i) = Sp_cur(i, 1)
                endif
                do j=2, ni, 1
                    a123 = (sigma*Cw*(d0(j)**eta)/eta/Cr)
                    if (Zgrid(i).lt.a123) then
                    Sp_cur(i, j) =  M0xCm*dd0*(pi/6.)*ro*pV*dp**(3.-sigma)
                    Sp_cur(i, j) =  Sp_cur(i, j)*(d0(j)**(sigma-epsilon))*(1-(eta*Cr*Zgrid(i))/sigma/Cw/(d0(j)**eta))**((sigma-eta)/eta)
                    if (mod(j,2).eq.0) then
                        Total_Sp(i) = Total_Sp(i) + 4*Sp_cur(i, j)
                    else
                        Total_Sp(i) = Total_Sp(i) + 2*Sp_cur(i, j)
                    endif    
                    endif
                enddo
                a123 = (sigma*Cw*(d0(ni+1)**eta)/eta/Cr)
                if (Zgrid(i).lt.a123) then
                    Sp_cur(i, ni+1) =  M0xCm*dd0*(pi/6.)*ro*pV*dp**(3.-sigma)
                    Sp_cur(i, ni+1) =  Sp_cur(i, ni+1)*(d0(ni+1)**(sigma-epsilon))*(1-(eta*Cr*Zgrid(i))/sigma/Cw/(d0(ni+1)**eta))**((sigma-eta)/eta)
                    Total_Sp(i) = Total_Sp(i) + Sp_cur(i, ni+1)
                endif
           Total_Sp(i) =  Total_Sp(i)/3
        enddo


        do i=1, n+1, 1
            Fd_cur(i,1) = 0
            Total_Fd(i) = 0
                a123 = (sigma*Cw*(d0(1)**eta)/eta/Cr)
                if (Zgrid(i).lt.a123) then
                    Fd_cur(i,1) = M0xCm * Cw * dd0*(pi/6.)*ro*pV*dp**(3.-sigma)
                    Fd_cur(i,1) = Fd_cur(i,1)*((d0(1)**(eta+sigma-epsilon))*(1-(eta*Cr*Zgrid(i))/sigma/Cw/(d0(1)**eta))**(sigma/eta))
                    Total_Fd(i) =Fd_cur(i,1)
                endif
                do j=2, ni, 1
                    a123 = (sigma*Cw*(d0(j)**eta)/eta/Cr)
                    if (Zgrid(i).lt.a123) then
                    Fd_cur(i,j) = M0xCm * Cw * dd0*(pi/6.)*ro*pV*dp**(3.-sigma)
                    Fd_cur(i,j) = Fd_cur(i,j)*(d0(j)**(eta+sigma-epsilon))*(1-(eta*Cr*Zgrid(i))/sigma/Cw/(d0(j)**eta))**(sigma/eta)
                    if (mod(j,2).eq.0) then
                        Total_Fd(i) = Total_Fd(i) + 4*Fd_cur(i,j)
                    else
                        Total_Fd(i) = Total_Fd(i) + 2*Fd_cur(i,j)
                    endif    
                    endif
                enddo
                a123 = (sigma*Cw*(d0(ni+1)**eta)/eta/Cr)
                if (Zgrid(i).lt.a123) then
                    Fd_cur(i,ni+1) =  M0xCm * Cw * dd0*(pi/6.)*ro*pV*dp**(3.-sigma)
                    Fd_cur(i,ni+1) = Fd_cur(i,ni+1)*(d0(ni+1)**(eta+sigma-epsilon))*(1-(eta*Cr*Zgrid(i))/sigma/Cw/(d0(ni+1)**eta))**(sigma/eta)
                    Total_Fd(i) =Total_Fd(i) + Fd_cur(i,ni+1)
                endif
             Total_Fd(i) =  Total_Fd(i)/3
                
        enddo
    
    
    else
    
    do i=1, n+1, 1
            Sp_cur(i, 1) = 0
            Total_Sp(i) = 0
            
            fi = (1-eta*b/sigma)/a/Cw/(d0(1)**eta)
            Sp_cur(i, 1) = M0xCm * dd0 * (pi/6.)*ro*pV*dp**(3.-sigma)
            Sp_cur(i, 1) = Sp_cur(i, 1)*(d0(1)**(sigma-epsilon))*(1+fi*Zgrid(i))**((eta-sigma)*b/(sigma-eta*b))
            Total_Sp(i) = Sp_cur(i, 1)
            do j=2, ni, 1
                fi = (1-eta*b/sigma)/a/Cw/(d0(j)**eta)
                Sp_cur(i, j) = M0xCm * dd0 * (pi/6.)*ro*pV*dp**(3.-sigma)
                Sp_cur(i, j) = Sp_cur(i, j)*(d0(j)**(sigma-epsilon))*(1+fi*Zgrid(i))**((eta-sigma)*b/(sigma-eta*b))
                if (mod(j,2).eq.0) then
                    Total_Sp(i) = Total_Sp(i) + 4*Sp_cur(i, j)
                else
                    Total_Sp(i) = Total_Sp(i) + 2*Sp_cur(i, j)
                endif 
            enddo
            fi = (1-eta*b/sigma)/a/Cw/(d0(ni+1)**eta)
            Sp_cur(i, ni+1) = M0xCm * dd0 *(pi/6.)*ro*pV*dp**(3.-sigma)
            Sp_cur(i, ni+1) = Sp_cur(i, ni+1)*(d0(ni+1)**(sigma-epsilon))*(1+fi*Zgrid(i))**((eta-sigma)*b/(sigma-eta*b))
            Total_Sp(i) = Total_Sp(i) + Sp_cur(i, ni+1) 
            Total_Sp(i) = Total_Sp(i)/3
    enddo


    
    do i=1, n+1, 1
            Fd_cur(i,1) = 0
            Total_Fd(i) = 0
            fi = (1-eta*b/sigma)/a/Cw/(d0(1)**eta)
            Fd_cur(i,1) = M0xCm * Cw * dd0*(pi/6.)*ro*pV*dp**(3.-sigma)
            Fd_cur(i,1) = Fd_cur(i,1)*(d0(1)**(eta+sigma-epsilon))*(1+fi*Zgrid(i))**(-sigma*b/(sigma-eta*b))
            Total_Fd(i) = Fd_cur(i,1)
            do j=2, ni, 1
                fi = (1-eta*b/sigma)/a/Cw/(d0(j)**eta)
                Fd_cur(i,j) = M0xCm * Cw * dd0*(pi/6.)*ro*pV*dp**(3.-sigma)
                Fd_cur(i,j) = Fd_cur(i,j) *(d0(j)**(eta+sigma-epsilon))*(1+fi*Zgrid(i))**(-sigma*b/(sigma-eta*b))
                if (mod(j,2).eq.0) then
                    Total_Fd(i) = Total_Fd(i) + 4*Fd_cur(i,j)
                else
                    Total_Fd(i) = Total_Fd(i) + 2*Fd_cur(i,j)
                endif 
            enddo
            fi = (1-eta*b/sigma)/a/Cw/(d0(ni+1)**eta)
            Fd_cur(i,ni+1) = M0xCm * Cw * dd0*(pi/6.)*ro*pV*dp**(3.-sigma)
            Fd_cur(i,ni+1) = Fd_cur(i,ni+1)*(d0(ni+1)**(eta+sigma-epsilon))*(1+fi*Zgrid(i))**(-sigma*b/(sigma-eta*b)) 
            Total_Fd(i) = Total_Fd(i) + Fd_cur(i,ni+1)   
       
        Total_Fd(i) =  Total_Fd(i)/3
    enddo

    endif
  deallocate(Fd_cur)
   end subroutine


