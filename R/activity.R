
TEOS<-function(V0,alpha,K0,dK0,d2K0,Pbar,TK,S0,n){
  #' THERMOCALC 2011 (Tait) EOS FOR SOLIDS
  #' Adapted from matlab code by Jesse Walters
  #' @export

  aT   <- (1+dK0)/(1+dK0+K0*d2K0)
  bT   <- (dK0/K0)-(d2K0/(1+dK0))
  cnum <- 1+dK0+K0*d2K0
  cden <- ((dK0*dK0)+dK0)-(K0*d2K0)
  cT   <- cnum/cden

  # Einstein Temperature approximation
  theta <- 10636/(S0/n+6.44)
  x <- theta/TK
  x0 <- theta/298.15
  ex <- exp(x)
  ex0 <- exp(x0)

  # Einstein thermal energy
  Eth  <- 3*n*8.3144598*theta*(0.5+(1/(ex-1)))
  Eth0 <- 3*n*8.3144598*theta*(0.5+(1/(ex0-1)))

  # Einstein heat capacity
  Cv0 <- 3*n*8.3144598*((x0*x0*ex0)/((ex0-1)*(ex0-1)))

  # thermal pressures
  Pt <- (alpha*K0*Eth)/Cv0
  Pt0 <- ((alpha*K0*Eth0)/Cv0)

  # thermal pressure relative to standard state
  Pth <- Pt-Pt0

  # intVdP
  P0 <- 100000;
  psubpth <- Pbar-P0-Pth

  intVdP <- (Pbar-P0)*V0*(1-aT+(aT*(((1-bT*Pth)^(1-cT))-((1+bT*(psubpth))^(1-cT)))/(bT*(cT-1)*(Pbar-P0))))
  # cat("TEOS:",intVdP,"\n",
  #     "n:",n,"\n",
  #     "Cv0:",Cv0,"\n",
  #     "Pth:",Pth,"\n"
  #     )
  return(intVdP)
}

Landau <- function(Tc0,Vmax,Smax,Pbar,TK){
  #' Landau Model - Holland and Powell 2011
  #' @details the HP2011 database uses a corrected Landau model for some
  #' minerals. This correction was not incorporated into the HP2011 paper;
  #' however, this is code uses the same formulism as Thermocalc. For the
  #' original HP1998 approach, please see the other code. In addition, thermal
  #' expansivity and bulk modulus are not included. In practice, HP set the Vf
  #' term of the volume integral to 1. (JW)
  #' Adapted from matlab code by Jesse Walters
  #' @export

  P0 <- 1e5 # pressure of experimental data in pascals
  T0 <- 298.15 # temperature of experimental data

  # Critical Temperature

  Tc <- Tc0+(Vmax*(Pbar-P0))/Smax # critical T of landau transition at P of interest

  # Q - Order Parameter
  # When T>Tc then Q=0
  if(TK>Tc){ # for T's greater than the critical T
    Q <- 0
  } else {
    Q <- ((Tc-TK)/Tc0)^(0.25)
  }
  if(T0>Tc0){# for T's greater than the critical T
    # CHANGED Tc -> Tc0 ???
    Q0 <- 0
  } else {
    Q0 <- ((Tc0-T0)/Tc0)^(1/4) # Q at reference T (298K)
  }

  # Excess Gibbs free energies
  Term1 <- Tc0*Smax*((Q0^2)-(1/3)*(Q0^6))
  Term2 <- Smax*(Tc*Q*Q-Tc0*(Q^6)*(1/3))
  Term3 <- TK*Smax*(Q0*Q0-Q*Q)
  Term4 <- (Pbar-P0)*Vmax*Q0*Q0

  # To deal with minerals that do not have a Landau correction, the Gibbs free
  # energy contribution of the Landau correction is set to zero for every row
  # with a critical temperature of zero. This requires the input data to have
  # zeros even if a critical temperature is not reported.

  if(Tc0==0){
    Gl <- 0
  }else{
    Gl <- Term1-Term2-Term3+Term4
  }


  dGdT <- Smax*(Q*Q-Q0*Q0)  # molar entropy contribution
  dGdP <- -Vmax*(Q*Q-Q0*Q0) # molar volume contribution
  d2GdP2 <- -(Vmax^2)/(2*Smax*Tc0*Q*Q) # bulk modulus contribution
  d2GdT2 <- -Smax/(2*Tc0*Q*Q) # heat capacity contribution
  d2GdPdT <- Vmax/(2*Tc0*Q*Q) # thermal expansivity contribution


  return(list(Gl=Gl,
              dGdT=dGdT,
              dGdP=dGdP,
              d2GdP2=d2GdP2,
              d2GdT2=d2GdT2,
              d2GdPdT=d2GdPdT))
}

GibbsSolid_landau<-function(Hf,S0,V0,a,b,c,d,alpha,K0,dK0,d2K0,PGPa,TK,n,Vmax,Smax,Tc0,doLandau){
  #' This function calculates the molar gibbs free energy for solid phases in
  #' two parts: G=G0+RTlnK
  #' where
  #' G0=Hf-TS0-integral(Cp)dT-T*integral(Cp/T)dT+integral(Vsolid)dP
  #'
  #' Adapted from matlab code by Jesse Walters
  #' @export

  ## Integral(Cp)dT
  T0 <- 298.15
  Pbar <- PGPa*1e9 # converts GPa to pascals
  intCpdT <- (a*TK+0.5*b*TK*TK-c/TK+2*d*sqrt(TK))-(a*T0+0.5*b*T0*T0-c/T0+2.0*d*sqrt(T0))

  ## Integral(Cp/T)dT
  intCpoverTdT <- a*log(TK/298.15)+b*(TK-298.15)-(c/2)*(1/(TK*TK)-1/(298.15*298.15))-2*d*(1/sqrt(TK)-1/(sqrt(298.15)))

  # calls the TEOS to calculate the
  # volume contribution to the Gibbs free energy
  intVdP <- TEOS(V0,alpha,K0,dK0,d2K0,Pbar,TK,S0,n)

  # excess gibbs free energy from Landau
  ### ONLY if the mineral has a Landau transition !
  if(doLandau){
    Gl <- Landau(Tc0,Vmax,Smax,Pbar,TK)$Gl
  }else{Gl <-0}

  Gs <- Hf+intCpdT-TK*(S0+intCpoverTdT)+intVdP+Gl;

  return(Gs)
}

Gibbs<-function(mineral,PGPa,TK){
  #' Wrapper
  #' @export
  G <- GibbsSolid_landau(mineral["Hf"],
                         mineral["S0"],
                         mineral["V0"],
                         mineral["aCp"],
                         mineral["bCp"],
                         mineral["cCp"],
                         mineral["dCp"],
                         mineral["alpha"],
                         mineral["K0"],
                         mineral["dK0"],
                         mineral["d2K0"],
                         PGPa,
                         TK,
                         mineral["nAt"]+mineral["nOx"],
                         mineral["Vmax"],
                         mineral["Smax"],
                         mineral["Tc0"],
                         doLandau=mineral["Tc0"]>0)

  names(G)<-"G"
  return(G)
}

vGibbs <- Vectorize(Gibbs,c("PGPa","TK"))
  #' vectorized version of Gibbs
  # Export tag is somewhere else !

###### fO2 #####

CORK<-function(Ca0,Ca1,Cb0,Cc0,Cc1,Cd0,Cd1,Tcf,Pcf,Pbar,TK){
  #' Compensated-Redlich-Kwong (CORK) Equation
  #' Following Holland and Powell 1991
  #' Jesse Walters
  #' Rtlnf
  #' Equations 9
  #' @export

  Ca <- ((Ca0*(Tcf^(2.5)))/Pcf)+((Ca1*(TK*Tcf^(1.5)))/(Pcf))
  Cb <- (Cb0*Tcf)/Pcf
  Cc <- ((Cc0*Tcf)/(Pcf^(1.5)))+((TK*Cc1)/((Pcf^(1.5))))
  Cd <- ((Cd0*Tcf)/(Pcf*Pcf))+((TK*Cd1)/(Pcf*Pcf))

  ## Fugacity equation
  Pr <- Pbar-10000 # relative pressure (Pa)

  if(Tcf==0){RTlnf <- 0}else{
    RTlnf <- R*TK*log(1e-5*Pr)+(Cb*Pr)+Ca/(Cb*sqrt(TK))*(log((R*TK)+(Cb*Pr))-log((R*TK)+(2*Cb*Pr)))+((2/3)*Cc*Pr*sqrt(Pr))+((Cd/2)*Pr*Pr)
  }

  return(RTlnf)
}


GibbsPure<-function(TK,PGPa,Hff,S0f,af,bf,cf,df,Ca0,Ca1,Cb0,Cc0,Cc1,Cd0,Cd1,Tcf,Pcf){
  #' Molar Gibbs Free Energy Calculator for Pure Solid Phases
  #' Jesse Walters
  #' Gf
  #' This function calculates the molar gibbs free energy for solid phases in
  #' two parts: G=G0+RTlnK
  #' where
  #' G0=Hf-TS0-integral(Cp)dT-T*integral(Cp/T)dT+integral(Vsolid)dP
  #' @export

  Pbar <- PGPa*1e9 # converts GPa to pascals

  ## Integral(Cp)dT
  T0 <- 298.15
  intCpdTf <- (af*TK+0.5*bf*TK*TK-cf/TK+2*df*sqrt(TK))-
    (af*T0+0.5*bf*T0*T0-cf/T0+2*df*sqrt(T0))

  ## Integral(Cp/T)dT
  intCpoverTdTf <- af*log(TK/298.15)+bf*(TK-298.15)-(cf/2)*(1/(TK*TK)-
                                                              1/(298.15*298.15))-2*df*(1/sqrt(TK)-1/(sqrt(298.15)))

  ## calls the CORK EoS to calculate the volume contribution to the Gibbs free energy
  RTlnf <- CORK(Ca0,Ca1,Cb0,Cc0,Cc1,Cd0,Cd1,Tcf,Pcf,Pbar,TK);

  Gf<- Hff + intCpdTf - TK*(S0f+intCpoverTdTf)+RTlnf

  return(Gf)

}

# Wrapper
GibbsO2<-function(PGPa,TK){
  #' Wrapper
  #' @export
  G <- GibbsPure(TK,PGPa,
                 Hff=O2thermo["Hf"],
                 S0f=O2thermo["S0"],
                 af=O2thermo["a"],
                 bf=O2thermo["b"],
                 cf=O2thermo["c"],
                 df=O2thermo["d"],
                 Ca0=O2thermo["Ca0"],
                 Ca1=O2thermo["Ca1"],
                 Cb0=O2thermo["Cb0"],
                 Cc0=O2thermo["Cc0"],
                 Cc1=O2thermo["Cc1"],
                 Cd0=O2thermo["Cd0"],
                 Cd1=O2thermo["Cd1"],
                 Tcf=O2thermo["Tc"],
                 Pcf=O2thermo["Pc"]
  )


  names(G)<-"G"
  return(G)
}

fO2<-function(mu,PGpa,TK){
  #' Oxygen fugacity
  #' @export
  G0<-GibbsO2(0.0001,TK)
  return( exp((mu - G0) / R / TK )  )
}

FMQ<-function(Pbar,TK){
  #' FMQ
  #' @export
  A <- 5.5976
  B <- 24505
  C <- 0.8099
  D <- 0.0937

  logfO2<-A - B / TK + C * log10(TK) + D * (Pbar - 1) / TK
  return(logfO2)
}

aTiO2<-function(Pbar,TK,muTiO2){
  #' Final user-facing function !
  #' @param Pbar Pressure, bars
  #' @param TK Temperature, Kelvin
  #' @param muTiO2 Chemical potential of TiO2
  #' @return The activity of TiO2
  #' @export

  aTiO2 <- exp( - ( vGibbs(Rt,Pbar*1e-4,TK ) - muTiO2) / R / TK )

  #The above code is rather brittle, better check...
  if( any(aTiO2 > 1)  ){
    cat("Something is wrong here, activities shouldn't be greater than 1 ! \n
         Maximum value reached for aTiO2:",
        round(max(aTiO2,na.rm=T),3),
        "\n
          Probable causes:
          1. (if you are very far from 1): this has to do with HSC_conversion tag (should be OFF)\n
          2. (if you are just above 1): you are not using dataset 6.2\n")
  }

  return(aTiO2)
}

aSiO2<-function(Pbar,TK,muSiO2){
  #' Final user-facing function !
  #' @param Pbar Pressure, bars
  #' @param TK Temperature, Kelvin
  #' @param muSiO2 Chemical potential of SiO2
  #' @return The activity of SiO2
  #' @export

  aSiO2 <- exp( - ( vGibbs(Qtz,Pbar*1e-4,TK ) - muSiO2) / R / TK )

  #The above code is rather brittle, better check...
  if( any(aSiO2 > 1)  ){
    cat("Something is wrong here, activities shouldn't be greater than 1 ! \n
         Maximum value reached for aSiO2:",
         round(max(aSiO2,na.rm=T),3),
         "\n
          Probable causes:
          1. (if you are very far from 1): this has to do with HSC_conversion tag (should be OFF)\n
          2. (if you are just above 1): you are not using dataset 6.2\n")

  }

  return(aSiO2)
}


fO2<-function(Pbar,TK,muO2){
  #' Final user-facing function !
  #' @param Pbar Pressure, bars
  #' @param TK Temperature, Kelvin
  #' @param muSiO2 Chemical potential of SiO2
  #' @return a data frame with fO2 and delta FMQ
  #' @export

  fO2 <- exp((muO2 - GibbsO2(0.0001,TK) ) / R / TK )
  Delta_FMQ = log10(fO2) - FMQ(Pbar,TK)

  #The above code is rather brittle, better check...
  if(any(fO2 > 1)  ){
    cat("Something is wrong here, activities shouldn't be greater than 1 ! \n
          Probable causes:
          1. (if you are very far from 1): this has to do with HSC_conversion tag (should be OFF)\n
          2. (if you are just above 1): you are not using dataset 6.2\n")

    cat("Maximum value reached for fO2:",
        round(max(fO2,na.rm=T),3),
        "\n")
  }

  return(data.frame(fO2=fO2,
                    delta_FMQ = Delta_FMQ))
}


