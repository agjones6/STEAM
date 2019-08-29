! Equations of State Module
! Courtesy of Dr. Doster, NCSU Department of Nuclear Engineering
MODULE state

    IMPLICIT NONE

CONTAINS

    !ALL VARIABLES IN IMPERIAL UNITS

    !Function to compute the density of a saturated liquid
    REAL(8) FUNCTION rhof(P)
        IMPLICIT NONE

        REAL(8):: P

        rhof = 1./vf(P)

        RETURN
    END


    !Function to compute the derivative of saturated liquid density
    !with respect to pressure
    REAL(8) FUNCTION drhofdP(P)
        IMPLICIT NONE

        REAL(8):: P

        drhofdP = -(1./vf(P)**2)*dvfdP(P)

        RETURN
    END


    !Function to compute the density of a saturated vapor
    REAL(8) FUNCTION rhog(P)
        IMPLICIT NONE

        REAL(8) P

        rhog = 1./vg(P)

        RETURN
    END


    !Function to compute the derivative of saturated vapor density
    !with respect to pressure
    REAL(8) FUNCTION drhogdP(P)
        IMPLICIT NONE

        REAL(8):: P

        drhogdP = -(1./vg(P)**2)*dvgdP(P)

        RETURN
    END


    !Function to compute the difference between saturated vapor and
    !saturated liquid specific volume
    REAL(8) FUNCTION vfg(P)
        IMPLICIT NONE

        REAL(8):: P

        vfg = vg(P)-vf(P)

        RETURN
    END


    !Function to compute the derivative of the difference between
    !saturated vapor and saturated liquid specIFic volume
    REAL(8) FUNCTION dvfgdP(P)
        IMPLICIT NONE

        REAL(8):: P

        dvfgdP = dvgdP(P)-dvfdP(P)

        RETURN
    END


	!Function to compute saturation temperature
    REAL(8) FUNCTION Tsat(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: a

        DIMENSION a(8)
        DATA a/-127.45099,.0736883,-5.127918e-5,2.941807e-8, &
               -8.968781e-12,1.066619e-15,228.4795,.1463839/
        Tsat = a(1)+a(2)*P+a(3)*P**2+a(4)*P**3+a(5)*P**4+a(6)*P**5 &
              +a(7)*P**a(8)

        RETURN
    END


	!Function to compute saturated liquid specific volume
    REAL(8) FUNCTION vf(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: a

        DIMENSION a(8)
        DATA a/.0158605,1.436698e-6,-6.546245e-10,1.2621567e-12, &
               -6.106028e-16,1.17416e-19,3.004294e-4,.3809203/
        vf = a(1)+a(2)*P+a(3)*P**2+a(4)*P**3+a(5)*P**4+a(6)*P**5 &
            +a(7)*P**a(8)

        RETURN
    END


	!Function to compute derivative of saturated liquid specific
	!volume with respect to pressure
    REAL(8) FUNCTION dvfdP(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: a

        DIMENSION a(8)
        DATA a/.0158605,1.436698e-6,-6.546245e-10,1.2621567e-12, &
               -6.106028e-16,1.17416e-19,3.004294e-4,.3809203/
        dvfdP = a(2)+2.*a(3)*P+3.*a(4)*P**2+4.*a(5)*P**3+5.*a(6)*P**4 &
               +a(8)*a(7)*P**(a(8)-1.)

        RETURN
    END


	!Function to compute saturated vapor specific volume
    REAL(8) FUNCTION vg(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: sum
        REAL(8):: a, b
        INTEGER:: i

        DIMENSION a(7),b(7)
        DATA a/5931.557,1142.2341,171.5671,41.76546, &
               11.64542,3.264609,.8898603/
        DATA b/11.60044,1.990131,.3299698,.0806798, &
               .0200894,4.596498e-3,7.761257e-4/
        sum = 0.

        DO i = 1,7
            sum = sum+a(i)*exp(-b(i)*P)
        ENDDO

        vg = sum

        RETURN
    END


	!Function to compute derivative of saturated vapor specific
	!volume with respect to pressure
    REAL(8) FUNCTION dvgdP(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: sum
        REAL(8):: a, b
        INTEGER:: i

        DIMENSION a(7),b(7)
        DATA a/5931.557,1142.2341,171.5671,41.76546, &
               11.64542,3.264609,.8898603/
        DATA b/11.60044,1.990131,.3299698,.0806798, &
               .0200894,4.596498e-3,7.761257e-4/
        sum = 0.

        DO i = 1,7
            sum = sum-a(i)*b(i)*exp(-b(i)*P)
        ENDDO

        dvgdP = sum

        RETURN
    END


	!Function to compute saturated liquid internal energy
    REAL(8) FUNCTION uf(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: a

        DIMENSION a(8)
        DATA a/-158.61531184,0.10213234,-5.88930771e-5,3.77381711e-8, &
               -1.22429068e-11,1.65122388e-15,227.55166589,0.14666680/
        uf = a(1)+a(2)*P+a(3)*P**2+a(4)*P**3+a(5)*P**4+a(6)*P**5 &
            +a(7)*P**a(8)

        RETURN
    END


	!Function to compute the derivative of saturated liquid
	!internal energy with respect to pressure
    REAL(8) FUNCTION dufdP(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: a

        DIMENSION a(8)
        DATA a/-158.61531184,0.10213234,-5.88930771e-5,3.77381711e-8, &
               -1.22429068e-11,1.65122388e-15,227.55166589,0.14666680/
        dufdP = a(2)+2.*a(3)*P+3.*a(4)*P**2+4.*a(5)*P**3+5.*a(6)*P**4 &
               +a(8)*a(7)*P**(a(8)-1.)

        RETURN
    END


	!Function to compute saturated vapor internal energy
    REAL(8) FUNCTION ug(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: a

        DIMENSION a(8)
        DATA a/961.57632930,-0.06346313,2.69645643e-5,-2.46758641e-8, &
               9.45803668e-12,-1.53574346e-15,82.19290877,0.13028136/
        ug = a(1)+a(2)*P+a(3)*P**2+a(4)*P**3+a(5)*P**4+a(6)*P**5 &
             +a(7)*P**a(8)

        RETURN
    END


	!Function to compute the derivative of saturated vapor
	!internal energy with respect to pressure
    REAL(8) FUNCTION dugdP(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: a

        DIMENSION a(8)
        DATA a/961.57632930,-0.06346313,2.69645643e-5,-2.46758641e-8, &
               9.45803668e-12,-1.53574346e-15,82.19290877,0.13028136/
        dugdP = a(2)+2.*a(3)*P+3.*a(4)*P**2+4.*a(5)*P**3+5.*a(6)*P**4 &
               +a(8)*a(7)*P**(a(8)-1.)

        RETURN
    END


    !Function to compute liquid density as a funtion of Internal Energy
	!and pressure
    REAL(8) FUNCTION rhol(u,P)
        IMPLICIT NONE

        REAL(8):: u, P
        REAL(8):: au_P, bu_P, cu_P, delu
        REAL(8):: a, b, c

        DIMENSION a(5),b(5),c(5)
        DATA a/-16.6292550e-3,-88.8920769e-6,80.7871571e-9,-33.8804030e-12, &
               4.68248456e-15/
        DATA b/42.3254971e-6,-269.512210e-9,289.924067e-12,-126.247514e-15, &
               18.2130057e-18/
        DATA c/100.115638e-9,-297.848894e-12,302.240803e-15, &
               -127.204073e-18,18.1566424e-21/
        au_P  = a(1)+a(2)*P+a(3)*P**2+a(4)*P**3+a(5)*P**4
        bu_P  = b(1)+b(2)*P+b(3)*P**2+b(4)*P**3+b(5)*P**4
        cu_P  = c(1)+c(2)*P+c(3)*P**2+c(4)*P**3+c(5)*P**4
        delu  = u-uf(P)
        rhol  = rhof(P)+au_P*delu+bu_P*delu**2+cu_P*delu**3

        RETURN
    END


	!Function to compute the derivative of liquid density with respect to
	!Internal energy
    REAL(8) FUNCTION drholdu(u,P)
        IMPLICIT NONE

        REAL(8):: u, P
        REAL(8):: au_P, bu_P, cu_P, delu
        REAL(8):: a, b, c

        DIMENSION a(5),b(5),c(5)
        DATA a/-16.6292550e-3,-88.8920769e-6,80.7871571e-9,-33.8804030e-12, &
               4.68248456e-15/
        DATA b/42.3254971e-6,-269.512210e-9,289.924067e-12,-126.247514e-15, &
               18.2130057e-18/
        DATA c/100.115638e-9,-297.848894e-12,302.240803e-15, &
               -127.204073e-18,18.1566424e-21/
        au_P    = a(1)+a(2)*P+a(3)*P**2+a(4)*P**3+a(5)*P**4
        bu_P    = b(1)+b(2)*P+b(3)*P**2+b(4)*P**3+b(5)*P**4
        cu_P    = c(1)+c(2)*P+c(3)*P**2+c(4)*P**3+c(5)*P**4
        delu    = u-uf(P)
        drholdu = au_P+2.*bu_P*delu+3.*cu_P*delu**2

        RETURN
    END


	!Function to compute the derivative of liquid density with respect to
	!pressure
    REAL(8) FUNCTION drholdP(u,P)
        IMPLICIT NONE

        REAL(8):: u, P
        REAL(8):: au_P, bu_P, cu_P, delu, daudP, dbudP, dcudP
        REAL(8):: a, b, c

        DIMENSION a(5),b(5),c(5)
        DATA a/-16.6292550e-3,-88.8920769e-6,80.7871571e-9,-33.8804030e-12, &
               4.68248456e-15/
        DATA b/42.3254971e-6,-269.512210e-9,289.924067e-12,-126.247514e-15, &
               18.2130057e-18/
        DATA c/100.115638e-9,-297.848894e-12,302.240803e-15, &
               -127.204073e-18,18.1566424e-21/
        au_P    = a(1)+a(2)*P+a(3)*P**2+a(4)*P**3+a(5)*P**4
        bu_P    = b(1)+b(2)*P+b(3)*P**2+b(4)*P**3+b(5)*P**4
        cu_P    = c(1)+c(2)*P+c(3)*P**2+c(4)*P**3+c(5)*P**4
        delu    = u-uf(P)
       !drholdu = au_P+2.*bu_P*delu+3.*cu_P*delu**2
        daudP   = a(2)+2.*a(3)*P+3.*a(4)*P**2+4.*a(5)*P**3
        dbudP   = b(2)+2.*b(3)*P+3.*b(4)*P**2+4.*b(5)*P**3
        dcudP   = c(2)+2.*c(3)*P+3.*c(4)*P**2+4.*c(5)*P**3
        drholdP = drhofdP(P)-dufdP(P)*drholdu(u,P)+ &
                  delu*daudP+delu**2*dbudP+delu**3*dcudP

        RETURN
    END


    !Function to compute surface tension
    REAL(8) FUNCTION sigma(P) ![lbf/ft]
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: T
        REAL(8):: a
        INTEGER:: i

        DIMENSION a(6)
        DATA a/5.244774e-3,-1.409364e-6,-3.573894e-8,1.087337e-10, &
               -1.700704e-13,1.000511e-16/
        T = Tsat(P)
        sigma = 0.

        DO i = 1,6
            sigma = sigma+a(i)*T**(i-1)
        ENDDO

        RETURN
    END


    !Function to compute Liquid Dynamic Viscosity
    REAL(8) FUNCTION Viscosity(T,P)
        IMPLICIT NONE

        REAL(8):: T, P
        REAL(8):: DeltaT, Vsat, R
        REAL(8):: a, b, aa, bb
        INTEGER:: i

        DIMENSION a(3),b(3),aa(4),bb(4)
        DATA a/3.69971196,4.27115194,0.75003508/
        DATA b/-0.01342834,-0.03890983,-2.19455284e-3/
        DATA aa/8.52917235e-4,-4.17979848e-5,2.6043459e-7,-2.20531928e-11/
        DATA bb/-1.13658775,0.01495184,-2.86548888e-5,2.17440064e-9/
        DeltaT = Tsat(P)-T

        IF(DeltaT .LT. 0.)THEN
            DeltaT = 0.
        END IF

        Vsat = 0.

        DO i = 1,3
            Vsat = Vsat+a(i)*exp(b(i)*T)
        ENDDO

        IF(DeltaT .GT. 0.)THEN
            R = (aa(1)+aa(2)*DeltaT+aa(3)*DeltaT**2+aa(4)*DeltaT**3)* &
                (bb(1)+bb(2)*T+bb(3)*T**2+bb(4)*T**3)
        ELSE
            R = 0.
        END IF

        Viscosity = Vsat+R

        RETURN
    END


    !Function Subprogram to Calculate Temperature from
    !Liquid Internal Energy
    REAL(8) FUNCTION Temp(u)
        IMPLICIT NONE

        REAL(8):: u
        REAL(8):: a
        INTEGER:: k

        DIMENSION a(6)
        DATA a /32.180814,.9858671,1.8576575e-4,-8.0930376e-7, &
                1.0831764e-9,-8.1894562e-13/
        Temp = 0.

        DO k = 1,6
            Temp = Temp+a(k)*u**(k-1)
        ENDDO

        RETURN
    END


    !Function to compute saturated vapor viscosity
    REAL(8) FUNCTION Mug(P)
        IMPLICIT NONE

        REAL(8):: P
        REAL(8):: T
        REAL(8):: a

        DIMENSION a(4)
        DATA a/0.01790175,5.22632124e-5,4.52957731e-19,0.05532552/
        T   = Tsat(P)
        Mug = a(1)+a(2)*T+a(3)*exp(a(4)*T)

        RETURN
    END


	!Function to compute density of a superheated vapor
    REAL(8) FUNCTION rhov(u,P)
        IMPLICIT NONE

        REAL(8):: u, P
        REAL(8):: delu, a0, b0, c0, v
        REAL(8):: a, b, c

        DIMENSION a(4),b(4),c(3)
        DATA a/1.63874718,2.35848674,-15.1292824,12.9097851/
        DATA b/-145.604704e-6,-5.97575487e-3,36.6261462e-3,-30.9658047e-3/
        DATA c/-6.06404047e-9,14.4121205e-12,-7.54565037e-15/
        delu = u-ug(P)
        a0   = a(1)/P+a(2)/P**2+a(3)/P**3+a(4)/P**4
        b0   = b(1)/P+b(2)/P**2+b(3)/P**3+b(4)/P**4
        c0   = c(1)+c(2)*P+c(3)*P**2
        v    = vg(P)+a0*delu+b0*delu**2+c0*delu**3
        rhov = 1/v

        RETURN
    END


    !Function to compute the derivative of superheated vapor density
    !with respect to internal energy
    REAL(8) FUNCTION drhovdu(u,P)
        IMPLICIT NONE

        REAL(8):: u, P
        REAL(8):: delu, a0, b0, c0, v, dvdu
        REAL(8):: a, b, c

        DIMENSION a(4),b(4),c(3)
        DATA a/1.63874718,2.35848674,-15.1292824,12.9097851/
        DATA b/-145.604704e-6,-5.97575487e-3,36.6261462e-3,-30.9658047e-3/
        DATA c/-6.06404047e-9,14.4121205e-12,-7.54565037e-15/
        delu = u-ug(P)
        a0   = a(1)/P+a(2)/P**2+a(3)/P**3+a(4)/P**4
        b0   = b(1)/P+b(2)/P**2+b(3)/P**3+b(4)/P**4
        c0   = c(1)+c(2)*P+c(3)*P**2
        v    = vg(P)+a0*delu+b0*delu**2+c0*delu**3

        dvdu    = a0+2.*b0*delu+3.*c0*delu**2
        drhovdu = (-1./v**2)*dvdu

        RETURN
    END


    !Function to compute the derivative of superheated vapor density
    !with respect to pressure
    REAL(8) FUNCTION drhovdP(u,P)
        IMPLICIT NONE

        REAL(8):: u, P
        REAL(8):: delu, a0, b0, c0, v, dadP, dbdP, dcdP, dvdP
        REAL(8):: a, b, c

        DIMENSION a(4),b(4),c(3)
        DATA a/1.63874718,2.35848674,-15.1292824,12.9097851/
        DATA b/-145.604704e-6,-5.97575487e-3,36.6261462e-3,-30.9658047e-3/
        DATA c/-6.06404047e-9,14.4121205e-12,-7.54565037e-15/
        delu = u-ug(P)
        a0   = a(1)/P+a(2)/P**2+a(3)/P**3+a(4)/P**4
        b0   = b(1)/P+b(2)/P**2+b(3)/P**3+b(4)/P**4
        c0   = c(1)+c(2)*P+c(3)*P**2
        v    = vg(P)+a0*delu+b0*delu**2+c0*delu**3

        dadP = -a(1)/P**2-2.*a(2)/P**3-3.*a(3)/P**4-4.*a(4)/P**5
        dbdP = -b(1)/P**2-2.*b(2)/P**3-3.*b(3)/P**4-4.*b(4)/P**5
        dcdP = c(2)+2.*c(3)*P

        dvdP = dvgdP(P)-dugdP(P)*(a0+2.*delu*b0+3.*delu**2*c0)+ &
            delu*dadP+delu**2*dbdP+delu**3*dcdP
        drhovdP = (-1./v**2)*dvdP

        RETURN
    END


    !Function to compute saturation pressure from temperature
    REAL(8) FUNCTION Psat(T)
        IMPLICIT NONE

        REAL(8):: T
        REAL(8):: F, a0, b0, c0, d0, d1, d2, d3, T1, T2, T3, a1, b1, c1, &
                  a2, b2, c2, a3, b3, c3, P0, epsilon, DeltaP, dFdP, DelP, P

        !F(P) = T-Tsat(P)

        !Spline coefficients

        a0 = -4.2783969
        b0 = 0.1700619
        c0 = -1.8445571e-3
        d0 = 6.8826857e-6
        d1 = 2.1727190e-5
        d2 = 3.5465308e-5
        d3 = 8.4908055e-5
        T1 = 262.9033023
        T2 = 441.5206305
        T3 = 602.5002447

        a1 = a0+b0*T1+c0*T1**2+d0*T1**3
        b1 = b0+2.*c0*T1+3.*d0*T1**2
        c1 = c0+3.*d0*T1

        a2 = a1+b1*(T2-T1)+c1*(T2-T1)**2+d1*(T2-T1)**3
        b2 = b1+2.*c1*(T2-T1)+3.*d1*(T2-T1)**2
        c2 = c1+3.*d1*(T2-T1)

        a3 = a2+b2*(T3-T2)+c2*(T3-T2)**2+d2*(T3-T2)**3
        b3 = b2+2.*c2*(T3-T2)+3.*d2*(T3-T2)**2
        c3 = c2+3.*d2*(T3-T2)

        IF(T .LT. T1)THEN
            Psat = a0+b0*T+c0*T**2+d0*T**3
        ELSEIF(T .LT. T2)THEN
            Psat = a1+b1*(T-T1)+c1*(T-T1)**2+d1*(T-T1)**3
        ELSEIF(T .LT. T3)THEN
            Psat = a2+b2*(T-T2)+c2*(T-T2)**2+d2*(T-T2)**3
        ELSE
            Psat = a3+b3*(T-T3)+c3*(T-T3)**2+d3*(T-T3)**3
        END IF

        P0 = Psat

        epsilon = 1.

        DO WHILE(epsilon.gt.1.e-4)
            DeltaP  = abs(0.001*P0)
            dFdP    = ((T-Tsat(P0+DeltaP))-(T-Tsat(P0-DeltaP)))/(2.*DeltaP)
            DelP    = (T-Tsat(P0))/dFdP
            P       = P0-DelP
            epsilon = abs(DelP/P0)
            P0      = P
        ENDDO

        Psat = P

        RETURN
    END


    !Function Subprogram to Calculate Liquid Internal Energy
    REAL(8) FUNCTION uliq(T)
        IMPLICIT NONE

        REAL(8):: T
        REAL(8):: F, u0, epsilon, Deltau, dFdu, Delu, u
        REAL(8):: a
        INTEGER:: k

        DIMENSION a(6)
        DATA a /-42.1658015,1.3305375,-3.0673856e-3,1.1675009e-5, &
                -1.9395597e-8,1.2214095e-11/
        !F(u) = T-Temp(u)

        u0 = 0.

        DO k = 1,6
            u0 = u0+a(k)*T**(k-1)
        ENDDO

        epsilon = 1.

        DO WHILE(epsilon .GT. 1.e-4)
            Deltau  = abs(0.001*u0)
            dFdu    = ((T-Temp(u0+Deltau))-(T-Temp(u0-Deltau)))/(2.*Deltau)
            Delu    = (T-Temp(u0))/dFdu
            u       = u0-Delu
            epsilon = abs(Delu/u0)
            u0      = u
        ENDDO

        uliq = u

        RETURN
    END

    REAL(8) FUNCTION hg(P)
        IMPLICIT NONE

        REAL(8):: P

        hg = ug(P)+P*vg(P)*144./778.

        RETURN
    END

    REAL(8) FUNCTION hf(P)
        IMPLICIT NONE

        REAL(8):: P

        hf = uf(P)+P*vf(P)*144./778.

        RETURN
    END

    REAL(8) FUNCTION hliq(T,P)
        IMPLICIT NONE

        REAL(8):: T, P

        hliq = uliq(T)+(144./778)*P/rhol(uliq(T),P)

        RETURN
    END

    !Function to calculate u from h (liquid) **NOT FROM DR. DOSTER
    REAL(8) FUNCTION u_given_h(h_liquid,T,P)
        IMPLICIT NONE

        REAL(8):: h_liquid, T, P

        u_given_h = h_liquid-(144./778)*P/rhol(uliq(T),P)

    END FUNCTION

    REAL(8) FUNCTION hfg(P)
        IMPLICIT NONE

        REAL(8):: P

        hfg = hg(P)-hf(P)

        RETURN
    END

    REAL(8) FUNCTION F1(T)
        IMPLICIT NONE

        REAL(8):: T
        REAL(8):: hl, P

        COMMON/F1Parameters/hl,P
        F1 = hl-hliq(T,P)

        RETURN
    END


    !Function to compute Specific Heat
    REAL(8) FUNCTION Cpl(T,P)
        IMPLICIT NONE

        REAL(8):: T, P
        REAL(8):: DeltaT, Cpsat, R
        REAL(8):: a, b, aa, bb
        INTEGER:: i

        DIMENSION a(3),b(3),aa(3),bb(4)
        DATA a/0.98850267,3.11434479e-4,9.79793383e-27/
        DATA b/1.00787645e-4,0.01223534,0.08836906/
        DATA aa/0.03500856,7.40539427e-4,-1.30297916e-6/
        DATA bb/0.41844219,-7.71336906e-3,3.23610762e-5,-3.94022105e-8/
        DeltaT = Tsat(P)-T

        IF(DeltaT .LT. 0.)THEN
            DeltaT = 0.
        END IF

        Cpsat = 0.

        DO i = 1,3
            Cpsat = Cpsat+a(i)*exp(b(i)*T)
        ENDDO

        IF(DeltaT .GT. 0.)THEN
            R = (aa(1)+aa(2)*DeltaT+aa(3)*DeltaT**2)* &
                (bb(1)+bb(2)*T+bb(3)*T**2+bb(4)*T**3)
        ELSE
            R = 0.
        END IF

        Cpl   = Cpsat+R

        RETURN
    END


	!Function to compute Thermal Conductivity of Water
    REAL(8) FUNCTION kl(T,P)
        IMPLICIT NONE

        REAL(8):: T, P
        REAL(8):: DeltaT, Ksat, R
        REAL(8):: a, aa, bb
        INTEGER:: i

        DIMENSION a(6),aa(6),bb(6)
        DATA a/0.28956598,9.98373531e-4,-2.76514034e-6,1.31610616e-9, &
               3.99581573e-12,-5.18550975e-15/
        DATA aa/-3.51256646e-3,6.04273366e-5,2.48976537e-7,3.85754267e-11, &
                -1.59857317e-13,2.20172921e-16/
        DATA bb/-0.01305876,9.88477177e-4,-5.52334508e-6,6.66724984e-9, &
                3.03459927e-11,-3.78351489e-14/
        DeltaT = Tsat(P)-T

        IF(DeltaT .LT. 0.)THEN
            DeltaT = 0.
        END IF

        Ksat = 0.

        DO i = 1,6
            Ksat = Ksat+a(i)*T**(i-1)
        ENDDO

        IF(DeltaT .GT. 0.)THEN
            R = (aa(1)+aa(2)*DeltaT+aa(3)*DeltaT**2+aa(4)*DeltaT**3+ &
                 aa(5)*DeltaT**4+aa(6)*DeltaT**5)* &
                (bb(1)+bb(2)*T+bb(3)*T**2+bb(4)*T**3+bb(5)*T**4+bb(6)*T**5)
        ELSE
            R = 0.
        END IF

        kl = Ksat+R

        RETURN
    END


END MODULE
