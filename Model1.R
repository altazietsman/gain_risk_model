#input parameters
carbmax = function(dfe) {
  
  enddf = length(dfe[[1]])
  GW = list()
  Evaporation = list()
  Assimilation = list()
  Canoppypresure = list()
  GC = list()
  VPDLEAF = list()
  
  
  for (i in 1:enddf) {
    
    
    airtemp = dfe[[4]][i] #airtemperature in degree celcius
    RH = dfe[[6]][i] #humidity in %
    patm = dfe[[8]][i] #atmospheric presure in kpa
    wind = dfe[[5]][i] #wind speed in m s-2
    solar = dfe[[3]][i] #total incoming solar radiation in W m-2
    date = dfe[[1]][i] #read date from csv or excel as number and not date time
    time = dfe[[2]][i] #read time as hours (12h)
    P = dfe[[7]][i] #soil presure 
    Tleaf = dfe[[9]][i] #leaf temperature in degree celcius
    kmax = dfe[[13]][i] #kmax
    b = dfe[[14]][i] #vulnerability curve parameter
    c = dfe[[15]][i] #vulnerability curve parameter
    vmax25 = dfe[[16]][i] 
    jmax25 = dfe[[17]][i] 
    r = dfe[[18]][i] #CO2 compensation point rate at 25 degree 
    recovery = dfe[[19]][i] #xylem recovery index
    LABA = dfe[[20]][i] #leaf area: basel area ratio
    
    
    
    
    ##Hydralic cost model
    
    #weibull curve to get p-crit
    Pint = 0 #soil water presure
    wbr = kmax*(exp(-1*(Pint/b)^c)) # weibull function
    PLC = 1 - wbr/kmax #PLC function
    
    ##get p critical
    while (PLC <= 0.99) {
      wbr = kmax*(exp(-1*(Pint/b)^c))
      PLC = 1 - wbr/kmax
      Pint = Pint + 0.1
    }
    
    #get critical P value
    Pcrit = Pint
    
    
    #get new PLC based on soil P value
    #weibull curve from P to Pcrit
    if (i<2){
      #get new PLC based on soil P value
      #weibull curve from P to Pcrit
      P = dfe[[7]][i]
      wbr = kmax*(exp(-1*(P/b)^c)) # weibull function
      kcrit = kmax*(exp(-1*(Pcrit/b)^c)) #kcrit at P99
      kmax = max(wbr) #instantanous kmax
      PLC = (kmax - wbr)/(kmax - kcrit)
      l = list() #initialize empty list for PLC
      pl = list() #initialize empty list for soil presures
      wbrl = list() #initialize empty list for k
      inc = 0.1
      
      
      while (P <=Pcrit){
        wbr = kmax*(exp(-1*(P/b)^c))
        PLC = 1 - wbr/kmax
        l = c(l,PLC)
        wbrl = c(wbrl, wbr)
        pl = c(pl,P)
        P = P + inc
      }
      df = do.call(rbind, Map(data.frame, A=pl, B=l,C=wbrl))
      colnames(df)[1] = 'canopy_presure'
      colnames(df)[2] = 'PLC'
      colnames(df)[3] = 'K'
      plcprev = df$PLC[1]
      
      
    } else  if (i > 1 && (dfe[[7]][i] >= dfe[[7]][i-1])) {
      
      #get new PLC based on soil P value
      #weibull curve from P to Pcrit
      P = dfe[[7]][i]
      wbr = kmax*(exp(-1*(P/b)^c)) # weibull function
      kcrit = kmax*(exp(-1*(Pcrit/b)^c)) #kcrit at P99
      kmax = max(wbr) #instantanous kmax
      PLC = (kmax - wbr)/(kmax - kcrit)
      l = list() #initialize empty list for PLC
      pl = list() #initialize empty list for soil presures
      wbrl = list() #initialize empty list for k
      inc = 0.1
      
      
      while (P <=Pcrit){
        wbr = kmax*(exp(-1*(P/b)^c))
        PLC = 1 - wbr/kmax
        l = c(l,PLC)
        wbrl = c(wbrl, wbr)
        pl = c(pl,P)
        P = P + inc
      }
      df = do.call(rbind, Map(data.frame, A=pl, B=l, C=wbrl))
      colnames(df)[1] = 'canopy_presure'
      colnames(df)[2] = 'PLC'
      colnames(df)[3] = 'K'
      plcnew = plcprev-(plcprev-df$PLC[1])*recovery
      df$PLCnew = ((1-plcnew)/(1-df$PLC[1]))*(df$PLC-1)+1
      df$knew = (-1*(df$PLCnew-1))*kmax
      
      #counter loop to keep original PLC if recoverd PLC is smaller than PLC
      
      if (df$PLCnew[1] > df$PLC[1]){
        df = subset(df, select = -c(PLC,K))
        df$PLC = df$PLCnew
        df$K =df$knew
      }
      plcprev = df$PLC[1]
      
      
    } else if (i > 1 && (dfe[[7]][i] <= dfe[[7]][i-1])){
      
      #get new PLC based on soil P value
      #weibull curve from P to Pcrit
      P = dfe[[7]][i]
      wbr = kmax*(exp(-1*(P/b)^c)) # weibull function
      kcrit = kmax*(exp(-1*(Pcrit/b)^c)) #kcrit at P99
      kmax = max(wbr) #instantanous kmax
      PLC = (kmax - wbr)/(kmax - kcrit)
      l = list() #initialize empty list for PLC
      pl = list() #initialize empty list for soil presures
      wbrl = list() #initialize empty list for k
      inc = 0.1
      
      
      while (P <=Pcrit){
        wbr = kmax*(exp(-1*(P/b)^c))
        PLC = (1 - wbr/kmax)
        l = c(l,PLC)
        wbrl = c(wbrl, wbr)
        pl = c(pl,P)
        P = P + inc
      }
      df = do.call(rbind, Map(data.frame, A=pl, B=l, C=wbrl))
      colnames(df)[1] = 'canopy_presure'
      colnames(df)[2] = 'PLC'
      colnames(df)[3] = 'K'
      plcprev = df$PLC[1]
      
      
    }
    
    df$integral_area = df$K/inc
    df$integral_sum = cumsum(df$integral_area)
    #E needs to be converted from kg h-1 m-2 Mpa -1 per basel area to mmol s-1 m-2 p leaf area
    df$E = (df$integral_sum - df$integral_sum[1])*55.4*(1/3600)*(1/LABA)*(1000)*(1/100)
    
    
    
    #Assimilation gain model
    ca = 40 #atmospheric CO2 partial presure in kpa
    oa = 21000 #atmospheric O2 partial presure in kpa
    
    satvp = 0.61365*exp((17.502*Tleaf)/(240.97+Tleaf)) # Buck 1981
    
    
    VPchm = ((RH/100))*satvp
    
    VPleaf = (100/100)*satvp
    
    df$VPDleaf = VPleaf - VPchm
    
    df$gw = (df$E/df$VPDleaf)*100 #stomatal conductance in mmol s-1 m-2
    df$gc = df$gw/1.6 #CO2 stomatal conductance in mmol s-1 m-2
    
    
    ko25 = 28202 #michales-menton constant for oxylation (pa) at 25 degrees celcious
    kc25 = 41 #michales-menton constant for carboxylation (pa) at 25 degrees celcious
    qe = 0.3 #quatum yield of electron transport rate (mol photon mol-1 e)
    cc = 0.9 #curvature of light response curve
    cj = 0.98 #curvature of factor je v. jc limited photosynthesis (collatz et al 1991 and in Campbell and Norman)
    
    
    df$vnumerator = vmax25*(1+exp((-4424)/(298*8.314)))*exp((73637/(8.314*298))*(1-(298/(273+Tleaf))))
    df$vdenominator = 1+exp((486*(Tleaf+273)-149252)/(8.314*(273+Tleaf)))
    df$vmax = df$vnumerator/df$vdenominator #Vmax adjusted according to temperature (Leuning 2002 temperature adjusted)
    
    df$kc = kc25*exp((79430*((Tleaf+273)-298))/(298*8.314*(Tleaf+273))) #Kc temperature corrected (Bernacchi et al. 2001 temperature adjusted)
    df$ko = ko25*exp((36380*((Tleaf+273)-298))/(298*8.314*(Tleaf+273))) #Ko temperature corrected (Bernacchi et al. 2001 temperature adjusted)
    df$r = r*exp((37830*((Tleaf+273)-298))/(298*8.314*(Tleaf+273))) #temperature adjusted (Bernacchi et al. 2001)
    
    par = solar*0.45*4.57 #PAR (calculated from Wm-2) (umols-1m-2) Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf)
    
    df$jnumerator = jmax25*(1+exp(-4534/(298*8.314)))*exp((50300/(8.314*298))*(1-(298/(Tleaf+273))))
    df$jdenominator = 1+exp((495*(Tleaf+273)-152044)/(8.314*(273+Tleaf)))
    df$jmax = df$jnumerator/df$jdenominator #Jmax adjusted according to temperature (Leuning 2002 temperature adjusted)
    
    df$resp = (0.01*vmax25*2^(((Tleaf+273)-298)/10))* (1 + exp(1.3 * (Tleaf - 55))) #respiration rate (umol s-1m-2 (collatz et a; 1991/ also in Campbell and Norman) and medlyn 2002
    
    #get assimilation rate and ci
    ci = 0
    A = list(0)
    Ag = 0.1
    An = 0
    end = length((df$gc))
    
    for (i in 2:end){
      An = 0
      Ag = 0.1
      while (Ag > An) {
        jc =(df$vmax[i]*(ci-r))/(ci+df$kc[i]*(oa/df$ko[i]))-df$resp[i]
        j =(qe*par+df$jmax[i]-((((qe*par+df$jmax[i])^2)-
                                  (par*df$jmax[i]*cc*qe*4))^0.5))/(2*cc)
        je = (j/4)*((ci-df$r[i])/(ci+2*df$r[i]))
        An = (je +jc-((((je+jc)^2)-4*cj*je*jc)^0.5))/(2*cj)
        Ag = (df$gc[i]*(ca-ci))/101.3
        ci = ci + inc
      }
      A = c(A, An)
      i = i + 1
    }
    
    df[["An"]] = as.numeric(as.character(unlist(A))) #add assimilation to dataframe (umol s-1 m-2)
    Amax = df$An[[end]] #get maximum An
    df$gain = df$An/Amax #normalize An
    df$profit = df$gain - df$PLC #get profit
    
    index = which.max(df$profit) #get index of maximum profit
    
    g = df$gw[index]
    E = df$E[index]
    As = df$An[index]
    Tleaf = df$Tleaf[index]
    Cpresure = df$canopy_presure[index]
    gc = df$gc[index]
    VPDL = df$VPDleaf[index]
    
    #create list of all the results
    GW = c(GW,g)
    Evaporation =  c(Evaporation,E)
    Assimilation =  c(Assimilation,As)
    Canoppypresure =  c(Canoppypresure,Cpresure)
    GC = c(GC,gc)
    VPDLEAF = c(VPDLEAF, VPDL)
    
    
  }
  
  #create dataframe of results
  #create dataframe of results
  results = do.call(rbind, Map(data.frame, A=GW, B=Evaporation, C=Assimilation,
                               E=Canoppypresure))
  colnames(results)[1] = 'Gw'
  colnames(results)[2] = 'E'
  colnames(results)[3] = 'An'
  colnames(results)[4] = 'Cpresure'
  
  
  results <<- results #move results df to global variable
  
  
}
