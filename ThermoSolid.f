       Program ThemoSolid
       
C     Defini‡Æo das Vari veis

      implicit none
      
      character*30 fin, fout, fres, ext1, ext2
      character*80 texto(4), titulo
      character*9 eos(3)

      real*8 :: h, s, sol, ex
      real*8 :: a1, a2, a3, pr, Vid
      integer i, j, k, m, maxit, t

      integer N, np, nm
      real*8 :: Temp(3,200), V0(3,200), B0(3,200), B1(3,200), Voc(200)
      real*8 :: P(200), E0(3,200)
      real*8 :: V(200,200), E(200,200), G(200,200)
      real*8 :: PV(200,200)
      real*8 :: fat
       
      fat=1312.7496997450642
      ext1='.out'
      ext2='.res'

      eos(1)='Murnaghan'
      eos(2)='BirchMurn'
      eos(3)='Vinet_exp'
      
C     Defini‡Æo dos Formatos
12    format (' ',125('*'))
14    format (A175)
15    format (A50)
16    format (A17, F5.2, A18, ES8.2)
38    format (A65)

      write(*,*) 'Name of the Input File ?'
      read(*,'(a)') fin

      open(5,file=fin,status='old')

C     manual: input
C     Primeira Linha: n£mero de valores de T, n£mero de valores de pressÆo, n£mero de mol‚culas por c‚lula unit ria
C     Segunda Linha em diante: T(K), Vo(ua^3), B0(kbar), B1, E0(ryd/cell)
C     éltima Linha: Valores de PressÆo (kbar)

      read(5,14) titulo
      read(5,*) N, np, nm
      read(5,*) (P(i), i=1, np)
      do t=1,3
      read(5,*)
      do i=1,N
      read(5,*) Temp(t,i), V0(t,i), B0(t,i), B1(t,i), E0(t,i)
      enddo
      enddo

      do t=1,3
      
      fout=eos(t)//ext1
      fres=eos(t)//ext2
      
      open(6,file=fout)
      open(7,file=fres)

      do i=1,N
      write(6,*) Temp(t,i), V0(t,i), B0(t,i), B1(t,i), E0(t,i)
      enddo

      write(6,*) 'This program calculates the molar volume of a solid',
     *' at given pressure and temperature value with the ', eos(t),
     *' EOS'
      write(6,*)

      write(6,14) titulo
      write(6,*) 'Newton-Raphson Calculation of the Volume'
      write(6,*)
C     bloco da equa‡Æo de estado
      do i=1,N
      Voc(i)=V0(t,i)*(0.5292**3)*((0.00000001)**3)*6.02D+23/nm
      a1=B0(t,i)
      a2=B1(t,i)
      a3=Voc(i)
      do j=1,np
      pr=P(j)
      if(j.eq.1) then
      Vid=Voc(i)*1.05
      else
      Vid=V(i,j-1)
      endif
      
      write(6,16) 'Temperature (K): ', Temp(t,i), ' Pressure (kbar): ',
     *P(j)
      write(6,*)
      
      CALL NR(Vid,a1,a2,a3,pr,t,h,s,sol)
      V(i,j)=sol
      enddo
      enddo

74    format (A8,200(ES14.4))
75    format (F8.2,200(F14.7))
76    format (F8.2,200(F14.4))

      ex=0.3333333333333333
      
      do i=1,N
      do j=1,np
      if (t.eq.1) then
      E(i,j)=E0(t,i)*fat/nm+((B0(t,i)*V(i,j)*0.1/B1(t,i))*(((Voc(i)/
     *V(i,j))**B1(t,i))*(1/(B1(t,i)-1))+1))-(B0(t,i)*Voc(i)*0.1/
     *(B1(t,i)-1))
      elseif (t.eq.2) then
      E(i,j)=E0(t,i)*fat/nm+(9*B0(t,i)*Voc(i)*0.1/16)*
     *((((((Voc(i)/V(i,j))**(2*ex))-1)**3)*B1(t,i))+(((((Voc(i)/
     *V(i,j))**(2*ex))-1)**2)*(6-4*((Voc(i)/V(i,j))**(2*ex)))))
      else
      E(i,j)=E0(t,i)*fat/nm+((2*B0(t,i)*Voc(i)*0.1/((B1(t,i)-1)**2))*
     *(2-(5+3*((V(i,j)/Voc(i))**(ex))*(B1(t,i)-1)-3*B1(t,i))*
     *EXP(-1.5*(B1(t,i)-1)*(((V(i,j)/Voc(i))**(ex))-1))))
      endif
      PV(i,j)=V(i,j)*P(j)/10
      G(i,j)=E(i,j)+PV(i,j)
      enddo
      enddo

      write(7,*)
      write(7,*) '                   ThermoSolid Program - 2021'
      write(7,*) '            Laboratorio de Cinetica Quimica - UFRRJ'
      write(7,*) '                      Version: May, 2021'
      write(7,*)
      write(7,*) 'Please cite: '
      write(7,38) ' ',
     *'*************************************************************',
     *'  Xavier, Jr, N.F., da Silva, Jr., A.M.; Bauerfeldt, G.F.    ',
     *'  What Rules the Relative Stability of a-, b-, and           ',
     *'  g-Glycine Polymorphs? Crystal Growth & Design, v. 20,      ',
     *'  p. 4695-4706, 2021. doi.org/10.1021/acs.cgd.0c00489        ',
     *'*************************************************************',
     *' '
      write(7,14) titulo
      write(7,*) 'Results from the ', eos(t), ' EOS'
      write(7,12)
      write(7,*) 'Volume (cm**3/mol)'
      write(7,15) 'P(kbar)'
      write(7,74) 'T(K)', (P(i),i=1,np)
      do i=1,N
      write(7,75) Temp(t,i), (V(i,j), j=1,np)
      enddo

      write(7,12)
      write(7,*) 'Energy (kJ/mol)'
      write(7,15) 'P(kbar)'
      write(7,74) 'T(K)', (P(i),i=1,np)
      do i=1,N
      write(7,76) Temp(t,i), (E(i,j), j=1,np)
      enddo

      write(7,12)
      write(7,*) 'Gibbs Free Energy (kJ/mol)'
      write(7,15) 'P(kbar)'
      write(7,74) 'T(K)', (P(i),i=1,np)
      do i=1,N
      write(7,76) Temp(t,i), (G(i,j), j=1,np)
      enddo
      enddo
      end program

C     Inicia a subrotina NR

      SUBROUTINE NR(Vid,a1,a2,a3,pr,t,h,s,sol)

      implicit none
      
      real*8 :: a(500), termo(500), tmais(500), tmenos(500)
      real*8 :: terro(500), termoxtres(500)
      real*8 :: f, fmais, fmenos, dp, ds
      real*8 :: h, y, yr, ex
      real*8 :: s, smais, smenos, sol
      real*8 :: somaterro, erro
      real*8 :: xum, xdois, xtres
      real*8 :: Vid, a1, a2, a3, pr
      integer i, j, k, m, maxit, t

      maxit=50
      h=pr/1000
      y=1e-6

43    format (A2,A5,A15,A13,A12,A9,2(A13),3(A12))
      write(6,43) 'k', 's', 'f(s)', 'f+(s)', 'f-(s)', 'df', 'd2f',
     *'x(1)', 'x(2)', 'x(3)', 'erro'
      s=Vid
      smais=s+s*h
      smenos=s-s*h

      erro=1000
      k=0
      do while (erro.gt.y)
      i=1
      f=0
      fmais=0
      fmenos=0
      somaterro=0
      ex=0.3333333333333333

C     bloco de constru‡Æo da fun‡Æo
      if (t.eq.1) then
      f=pr-((a1/a2)*(((a3/s)**(a2))-1))
      fmais=pr-((a1/a2)*(((a3/smais)**(a2))-1))
      fmenos=pr-((a1/a2)*(((a3/smenos)**(a2))-1))
      elseif (t.eq.2) then
      f=pr-(1.5*a1*(((a3/s)**(7*ex))-((a3/s)**(5*ex)))*(1+0.75*(a2-4)*
     *(((a3/s)**(2*ex))-1)))
      fmais=pr-(1.5*a1*(((a3/smais)**(7*ex))-((a3/smais)**(5*ex)))*
     *(1+0.75*(a2-4)*(((a3/smais)**(2*ex))-1)))
      fmenos=pr-(1.5*a1*(((a3/smenos)**(7*ex))-((a3/smenos)**(5*ex)))*
     *(1+0.75*(a2-4)*(((a3/smenos)**(2*ex))-1)))
      else
      f=pr-(3*a1*((a3/s)**(2*ex))*(1-((s/a3)**(ex)))*
     *exp(-1.5*(a2-1)*(((s/a3)**(ex))-1)))
      fmais=pr-(3*a1*((a3/smais)**(2*ex))*
     *(1-((smais/a3)**(ex)))*
     *exp(-1.5*(a2-1)*(((smais/a3)**(ex))-1)))
      fmenos=pr-(3*a1*((a3/smenos)**(2*ex))*
     *(1-((smenos/a3)**(ex)))*
     *exp(-1.5*(a2-1)*(((smenos/a3)**(ex))-1)))
      endif
C     fim do bloco

      yr=s*h
      dp=(fmais-fmenos)/(2*yr)
      ds=(fmais-2*f+fmenos)/(yr**2)
      xum=s-f/dp
      xdois=s-f/(dp+(ds*(xum-s)/2))
      xtres=s-f/(dp+(ds*(xdois-s)/2))

      if (t.eq.1) then
      somaterro=pr-((a1/a2)*((xtres/a3)**(-a2)-1))
      elseif (t.eq.2) then
      somaterro=pr-(1.5*a1*(((a3/xtres)**(7*ex))-((a3/s)**(5*ex)))*
     *(1+0.75*(a2-4)*(((a3/xtres)**(2*ex))-1)))
      else
      somaterro=pr-(3*a1*((a3/xtres)**(2*ex))*
     *(1-((xtres/a3)**(ex)))*
     *exp(-1.5*(a2-1)*(((xtres/a3)**(ex))-1)))
      endif
      
      k=k+1
      
      if(k.gt.maxit) then
80    format (A30, I2, A11)
      write(6,80) 'Convergence was not achieved after ', maxit,
     *' iteractions'
      sol=0
      go to 122
      
      else
      
      erro=abs(somaterro)

85    format (I2, 10(ES12.2))

      write(6,85) k, s, f, fmais, fmenos, dp, ds, xum, xdois, xtres,
     *erro


      if ((f/somaterro).gt.0) then
      
      s=xtres
      
      else
      
      s=(s+xtres)/2

      endif

      smais=s+s*h
      smenos=s-s*h
      endif
      sol=s
90    end do

      write(6,*)
95    format (' ',54('*'))
96    format (A21, ES12.4, A10, ES12.4)
      write(6,95)
      write(6,96) 'Final Solution: V = ', sol, ' ; error = ',erro
      write(6,95)
      write(6,*)
122   END SUBROUTINE NR

