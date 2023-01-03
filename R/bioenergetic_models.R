#' Gavrilchuk bioenergetic model
#'
#' @param input_file dataframe
#' @param coeff.capt depending on species. see lehoux et al. 2020
#' @param mean.ed mean energy density of the taxa consumed
#' @param sd.ed standard deviation of the energy density of the  taxa consumed
#' @param change_150m logicial does the behaviour change below 150m
#' @param ID unique ID of sample
#' @importFrom dplyr mutate arrange
#' @importFrom magrittr %>%
#' @return a data frame
#' @export

#'
#' @examples
#' #todo
NARW_bionerg<-function(input_file, coeff.capt=1, mean.ed,sd.ed, change_150m,ID){
  ## Import  prey data:
  df<-input_file

  #using output en energy, convertir en biomass avec la valeur moyenne,
  df <-  df %>%  dplyr::mutate(DW_Z =(Energy_kJ_g_Z/mean.ed) *1000,
                        DW_Z_gm3=DW_Z*0.0001) %>%  dplyr::arrange(Zsampled)# Add a column to convert Calanus_DW_Z (mg/10m3) into g/m3


  if(change_150m==T){
  # Split df according to depth:
  df1<- df %>%  dplyr::filter(Zsampled <= 150)
  df2<- df %>%  dplyr::filter(Zsampled > 150)
  }

  # NARW morphometric parameters:
  length.cm<-1400 # average adult length
  length.m<-length.cm*0.01 # convert length to m
  mass<-35000 # Moore et al. 2004.Didn't change mass of pregnant female; calf mass at birth only represents 3% of mother mass (1101 kg/36258 kg)
  width.m<-3.15  # Nousek McGregor 2010: Body width measurements from necropsy reports of stranded right whales animals (Moore et al. 2004) & aerial photogrammetric measurements of live whales at sea (Miller Angell, 2006; Perryman & Lynn, 2002)

  # NARW swim speed (m/s):
  U.s.b<-1.0 # surface, bottom
  U.d.a<-1.45   # descent; Baumgartner et al. 2017
  # Body pitch during descent/ascent, used to calculate distance covered:
  desc_ang<-sin(74*pi/180) # angle of descent=74; Nousek McGregor et al. 2014
  asc_ang<-sin(62*pi/180) # angle of ascent=62; Nousek McGregor et al. 2014

  if(change_150m==T){
  # 0-150m dive parameters:
  df1 <- df1 %>%  dplyr::mutate( bottom_t = (0.0704 * Zsampled)*60) # DBA_t = submerged time


  # 160-500m dive parameters:
  # Set max bottom time (BT) constant & equal to max predicted by Baumgartner's (2017) eq. at 150m
  if(nrow(df2) > 0){
    df2 <- df2 %>%  dplyr::mutate( bottom_t = (0.0704 * 150)*60) # DBA_t = submerged time
    dfbt<-rbind(df1, df2)
  }

  if(nrow(df2)==0) dfbt<-df1
  }

  if(change_150m ==F){
    dfbt<- df %>% dplyr::mutate(bottom_t = (0.0704 * 150)*60)
  }

  dfbt <-  dfbt %>%  dplyr::mutate(bottom_d = U.s.b * bottom_t , # distance (m) = speed * time
                              desc_d = Zsampled/desc_ang,
                              asc_d = Zsampled/asc_ang,
                              desc_t = desc_d/U.d.a,
                              asc_t = asc_d/U.d.a,
                              DBA_t = desc_t + asc_t + bottom_t)

  # Surface time:
  # A fraction of total time spent submerged; Baumgartner & Mate 2003: resting females (n=11): 21.1%, pregnant (n=1, only 4 dives recorded from tag data): 51.1%, lactating females (n=4): 34.2%
  # set % surface time to 34.2 for both pregnant and lactating
  dfbt <-  dfbt %>%  dplyr::mutate(surface_t.rest = 0.211*DBA_t, #surface time
                              surface_t.prla = 0.342*DBA_t,
                              surface_d.rest = surface_t.rest*U.s.b, #surface distance
                              surface_d.prla = surface_t.prla*U.s.b,
                              dive.cycle.t_rest = DBA_t + surface_t.rest,# Total dive cycle time (all 4 phases):
                              dive.cycle.t_prla = DBA_t + surface_t.prla,
                              total.dives.h_rest = 3600/dive.cycle.t_rest, # Total no. of dives per hour:
                              total.dives.h_prla = 3600/dive.cycle.t_prla
                              )


  # Total time spent foraging per day (converted to sec)
  forage.s.min<-15.1*3600
  forage.s.max<-17.2*3600
  forage.s.mean<-16.15*3600

  # Total no. of dives per day
  dfbt <-  dfbt %>%  dplyr::mutate(total.dives.d_rest.forage.min = forage.s.min/dive.cycle.t_rest, # Resting, min & max foraging time per day
                              total.dives.d_rest.forage.max = forage.s.max/dive.cycle.t_rest,
                              total.dives.d_rest.forage.mean = forage.s.mean/dive.cycle.t_rest,
                              total.dives.d_prla.forage.min = forage.s.min/dive.cycle.t_prla,# Pregnant & Lactating, min & max foraging time per day
                              total.dives.d_prla.forage.max = forage.s.max/dive.cycle.t_prla,
                              total.dives.d_prla.forage.mean = forage.s.mean/dive.cycle.t_prla,
                              pctBT_rest = bottom_t/dive.cycle.t_rest, # Percent bottom time per dive:
                              pctBT_prla = bottom_t/dive.cycle.t_prla,
                              forage_t_rest.min = pctBT_rest * forage.s.min,  # Total bottom time in sec (ingestion time) per day:
                              forage_t_rest.max = pctBT_rest * forage.s.max,
                              forage_t_rest.mean = pctBT_rest * forage.s.mean,
                              forage_t_prla.min = pctBT_prla * forage.s.min,
                              forage_t_prla.max = pctBT_prla * forage.s.max,
                              forage_t_prla.mean = pctBT_prla * forage.s.mean
  )



  # Shreer & Kovacs 1997: Maximum dive duration is significantly correlated with mass for several species of air-breathing vertebrates:
  # Max duration = 0.04*mass^0.61 (= 23.65 min)

  # Halsey et al. 2006: For a 35000kg whale,
  # Mean maximum dive duration= 35.5*mass^0.326 = 17.92 min
  # Mean dive depth (mammals)= 3.8*mass^0.389 = 223 m
  # Mean maximum dive depth= 9.4*mass^0.327 = 287 m


  ## Maintenance costs (J/s):
  RMR<-(581*mass^0.68)/1000 # 95% of BMR; Marine mammal resting metabolism (Castellini & Mellish 2016 p.51) in kJ/day; where 1 kJ/day = 0.0115741 J/s; divide by 1000 for MJ/d

  ## Activity costs (foraging, traveling):
  # Expressed as Propulsive power (Pp, J/s)= Fdrag * U / nm * np
  # Fdrag = 0.5*rho*Cd*A*U^2*gamma*g*lambda
  # Cd (Drag coef) = Cf*(1 + k) (skin friction coef*form factor)
  # Cf = 0.072*Re^-1/5
  # Re = L*U/v

  # Propulsive power parameters:

  U.s.b<-1.0 # surface & bottom speed (m/s)
  U.d.a<-1.45 # descent & ascent speed (m/s)
  U.trav<-2.0 # traveling/migration speed (m/s)
  v<-1.83*10^-6  # seawater kinematic viscosity (m2/s)
  rho<-1029 # density of seawater; range 1026-1032 (kg/m3)
  gamma<-1 # Surface/wave drag (increase if animal spends extended periods of time at surface)
  g<-1.3 # Appendage drag (flippers & flukes); van der Hoop et al. 2014
  lambda<-1 # Body oscillation drag (range 0.5-2.0)
  A<-0.08*mass^0.65 # body surface area (m2)
  A.preg<-A*1.05 # Increase surface area of pregnant whale by 5% (Nousek McGregor 2010, p.107)
  nm<-0.25 # metabolic efficiency (nm)
  np<-0.51 # propulsive efficiency (np)

  # Propulsive power formulae:

  # Skin friction coefficient (Cf):
  Cf.s.b<- 0.072*((length.m*U.s.b)/v)^(-1/5) # surface, bottom
  Cf.d.a<- 0.072*((length.m*U.d.a)/v)^(-1/5) # descent, ascent
  Cf.trav<-0.072*((length.m*U.trav)/v)^(-1/5) # traveling/migration

  # Drag coefficient (Cd)
  Cd.s.b<- Cf.s.b* (1 + 1.5*(length.m/width.m)^(-3/2) + 7*(length.m/width.m)^-3)
  Cd.d.a<- Cf.d.a* (1 + 1.5*(length.m/width.m)^(-3/2) + 7*(length.m/width.m)^-3)
  Cd.trav<-Cf.trav*(1 + 1.5*(length.m/width.m)^(-3/2) + 7*(length.m/width.m)^-3)

  # Drag force (Fdrag, N)
  Drag.s<-   0.5 * rho * Cd.s.b * A * U.s.b^2 * gamma * g * lambda # surface
  Drag.b.2<- Drag.s * 2 # bottom, low drag (x2)
  Drag.b.2.5<- Drag.s * 2.5 # bottom, high drag (x3)
  Drag.b.3<- Drag.s * 3 # bottom, high drag (x3)
  Drag.d.a<- 0.5 * rho * Cd.d.a * A * U.d.a^2 * gamma * g * lambda # descent, ascent
  Drag.trav<-0.5 * rho * Cd.trav * A * U.trav^2 * gamma * g * lambda # trav

  # Drag force for pregnant females
  Drag.s.preg<-   0.5 * rho * Cd.s.b * A.preg * U.s.b^2 * gamma * g * lambda # surface
  Drag.b.2.preg<- Drag.s.preg * 2 # bottom, low drag (x2)
  Drag.b.2.5.preg<- Drag.s.preg * 2.5 # bottom, high drag (x3)
  Drag.b.3.preg<- Drag.s.preg * 3 # bottom, high drag (x3)
  Drag.d.a.preg<- 0.5 * rho * Cd.d.a * A.preg * U.d.a^2 * gamma * g * lambda # descent, ascent
  Drag.trav.preg<-0.5 * rho * Cd.trav * A.preg * U.trav^2 * gamma * g * lambda # trav

  # Propulsive power (J/s)
  Pp.s<-   (Drag.s    * U.s.b) / (nm*np) # surface
  Pp.b.2<- (Drag.b.2  * U.s.b) / (nm*np) # bottom, low drag
  Pp.b.2.5<- (Drag.b.2.5  * U.s.b) / (nm*np) # bottom, high drag
  Pp.b.3<- (Drag.b.3  * U.s.b) / (nm*np) # bottom, high drag
  Pp.d.a<- (Drag.d.a  * U.d.a) / (nm*np) # descent, ascent
  Pp.trav<-(Drag.trav * U.trav)/ (nm*np) # traveling/migration

  # Propulsive power for pregnant females
  Pp.s.preg<-   (Drag.s.preg    * U.s.b) / (nm*np) # surface
  Pp.b.2.preg<- (Drag.b.2.preg  * U.s.b) / (nm*np) # bottom, low drag
  Pp.b.2.5.preg<- (Drag.b.2.5.preg  * U.s.b) / (nm*np) # bottom, high drag
  Pp.b.3.preg<- (Drag.b.3.preg  * U.s.b) / (nm*np) # bottom, high drag
  Pp.d.a.preg<- (Drag.d.a.preg  * U.d.a) / (nm*np) # descent, ascent
  Pp.trav.preg<-(Drag.trav.preg * U.trav)/ (nm*np) # traveling/migration

  # Percentage of time spent gliding per dive phase
  s.gl<- 0.09
  d.gl<- 0.36
  b.gl<- 0.09
  a.gl<- 0.30


  # Propulsive power (J/s) required per dive
                            # RESITNG
                            # Low forage drag
  dfbt <-  dfbt %>% dplyr::mutate(Pp.forage.2_BT.rest=(Pp.d.a * (desc_t    * (1-d.gl))) +
                             (Pp.b.2 * (bottom_t  * (1-b.gl))) +
                             (Pp.d.a * (asc_t     * (1-a.gl))) +
                             (Pp.s   * (surface_t.rest * (1-s.gl))),
                           # High forage drag
                           Pp.forage.3_BT.rest=(Pp.d.a * (desc_t    * (1-d.gl))) +
                             (Pp.b.3 * (bottom_t  * (1-b.gl))) +
                             (Pp.d.a * (asc_t     * (1-a.gl))) +
                             (Pp.s   * (surface_t.rest * (1-s.gl))),
                          # Median forage drag
                           Pp.forage.2.5_BT.rest=(Pp.d.a * (desc_t    * (1-d.gl))) +
                            (Pp.b.2.5 * (bottom_t  * (1-b.gl))) +
                            (Pp.d.a * (asc_t     * (1-a.gl))) +
                            (Pp.s   * (surface_t.rest * (1-s.gl))),
                          # PREGNANT
                          # Low forage drag
                           Pp.forage.2_BT.preg=(Pp.d.a.preg * (desc_t * (1-d.gl))) +
                            (Pp.b.2.preg * (bottom_t  * (1-b.gl))) +
                            (Pp.d.a.preg * (asc_t     * (1-a.gl))) +
                            (Pp.s.preg   * (surface_t.prla * (1-s.gl))),
                          # High forage drag
                          Pp.forage.3_BT.preg=(Pp.d.a.preg * (desc_t * (1-d.gl))) +
                            (Pp.b.3.preg * (bottom_t  * (1-b.gl))) +
                            (Pp.d.a.preg * (asc_t     * (1-a.gl))) +
                            (Pp.s.preg   * (surface_t.prla * (1-s.gl))),

                          # Median forage drag
                          Pp.forage.2.5_BT.preg=(Pp.d.a.preg * (desc_t * (1-d.gl))) +
                            (Pp.b.2.5.preg * (bottom_t  * (1-b.gl))) +
                            (Pp.d.a.preg * (asc_t     * (1-a.gl))) +
                            (Pp.s.preg   * (surface_t.prla * (1-s.gl))),
                          # LACTATING
                          # Low forage drag
                          Pp.forage.2_BT.lact=(Pp.d.a * (desc_t    * (1-d.gl))) +
                            (Pp.b.2 * (bottom_t  * (1-b.gl))) +
                            (Pp.d.a * (asc_t     * (1-a.gl))) +
                            (Pp.s   * (surface_t.prla * (1-s.gl))),

                          # High forage drag
                          Pp.forage.3_BT.lact=(Pp.d.a * (desc_t    * (1-d.gl))) +
                            (Pp.b.3 * (bottom_t  * (1-b.gl))) +
                            (Pp.d.a * (asc_t     * (1-a.gl))) +
                            (Pp.s   * (surface_t.prla * (1-s.gl))),

                          # Median forage drag
                          Pp.forage.2.5_BT.lact=(Pp.d.a * (desc_t    * (1-d.gl))) +
                            (Pp.b.2.5 * (bottom_t  * (1-b.gl))) +
                            (Pp.d.a * (asc_t     * (1-a.gl))) +
                            (Pp.s   * (surface_t.prla * (1-s.gl)))

  )


  # Summer foraging cost (MJ/d):
  dfbt <-  dfbt %>% dplyr::mutate(# RESTING
                            E.forage.2_BT.rest = (Pp.forage.2_BT.rest*1e-6) * total.dives.d_rest.forage.min,
                            E.forage.3_BT.rest = (Pp.forage.3_BT.rest*1e-6) * total.dives.d_rest.forage.max,
                            E.forage.2.5_BT.rest = (Pp.forage.2.5_BT.rest*1e-6) * total.dives.d_rest.forage.mean,
                            # PREGNANT
                            E.forage.2_BT.preg = (Pp.forage.2_BT.preg*1e-6) * total.dives.d_prla.forage.min,
                            E.forage.3_BT.preg = (Pp.forage.3_BT.preg*1e-6) * total.dives.d_prla.forage.max,
                            E.forage.2.5_BT.preg = (Pp.forage.2.5_BT.preg*1e-6) *  total.dives.d_prla.forage.mean,
                            # LACTATING
                            E.forage.2_BT.lact = (Pp.forage.2_BT.lact*1e-6) * total.dives.d_prla.forage.min,
                            E.forage.3_BT.lact = (Pp.forage.3_BT.lact*1e-6) * total.dives.d_prla.forage.max,
                            E.forage.2.5_BT.lact = (Pp.forage.2.5_BT.lact*1e-6) * total.dives.d_prla.forage.mean)


  #plot(E.forage.2_BT.rest ~dfbt$Zsampled, pch=2, cex=1.2, ylab='Energy cost of foraging (MJ/day)',xlab='Depth (m)')

  # Traveling cost (MJ/d):
  # 6.2 to 8.3 h per day (summer)
  # 19.9 to 21.4 h per day (winter)

  travel.summ.min<-6.2 * 3600 # convert to sec
  travel.summ.max<-8.3 * 3600
  travel.summ.mean<-7.25 * 3600

  travel.mig<- 20 * 3600
  travel.wint.min<-19.9 * 3600
  travel.wint.max<-21.4 * 3600
  travel.wint.mean<-20.65 * 3600

  # RESTING & LACTATING
  E.trav.summ.min<-((Pp.trav * (1-s.gl)) * travel.summ.min) * 1e-6
  E.trav.summ.max<-((Pp.trav * (1-s.gl)) * travel.summ.max) * 1e-6
  E.trav.summ.mean<-((Pp.trav * (1-s.gl)) * travel.summ.mean) * 1e-6

  E.trav.migr<-((Pp.trav * (1-s.gl)) * travel.mig) * 1e-6

  E.trav.wint.min<-((Pp.trav * (1-s.gl)) * travel.wint.min) * 1e-6
  E.trav.wint.max<-((Pp.trav * (1-s.gl)) * travel.wint.max) * 1e-6
  E.trav.wint.mean<-((Pp.trav * (1-s.gl)) * travel.wint.mean) * 1e-6

  # PREGNANT
  E.trav.summ.preg.min<-((Pp.trav.preg * (1-s.gl)) * travel.summ.min)  * 1e-6
  E.trav.summ.preg.max<-((Pp.trav.preg * (1-s.gl)) * travel.summ.max)  * 1e-6
  E.trav.summ.preg.mean<-((Pp.trav.preg * (1-s.gl)) * travel.summ.mean)  * 1e-6

  E.trav.migr.preg<-((Pp.trav.preg * (1-s.gl)) * travel.mig) * 1e-6

  E.trav.wint.preg.min<-((Pp.trav.preg * (1-s.gl)) * travel.wint.min) * 1e-6
  E.trav.wint.preg.max<-((Pp.trav.preg * (1-s.gl)) * travel.wint.max) * 1e-6
  E.trav.wint.preg.mean<-((Pp.trav.preg * (1-s.gl)) * travel.wint.mean) * 1e-6


  ## TOTAL daily summer activity cost (MJ/d)
  dfbt <-  dfbt %>%  dplyr::mutate(# RESTING
                            E.summ.2_BT.rest = E.forage.2_BT.rest + E.trav.summ.min + RMR,
                            E.summ.3_BT.rest = E.forage.3_BT.rest + E.trav.summ.max + RMR ,
                            E.summ.2.5_BT.rest = E.forage.2.5_BT.rest + E.trav.summ.mean + RMR,

                            # PREGNANT
                            E.summ.2_BT.preg = E.forage.2_BT.preg + E.trav.summ.preg.min + RMR,
                            E.summ.3_BT.preg = E.forage.3_BT.preg + E.trav.summ.preg.max + RMR,
                            E.summ.2.5_BT.preg = E.forage.2.5_BT.preg + E.trav.summ.preg.mean + RMR,

                            # LACTATING
                            E.summ.2_BT.lact = E.forage.2_BT.lact + E.trav.summ.min + RMR,
                            E.summ.3_BT.lact = E.forage.3_BT.lact + E.trav.summ.max + RMR,
                            E.summ.2.5_BT.lact = E.forage.2.5_BT.lact + E.trav.summ.mean + RMR
  )

  ## TOTAL daily spring/fall activity cost (MJ/d)

  # RESTING & LACTATING
  E.migr<- E.trav.migr + RMR

  # PREGNANT
  E.migr.preg<- E.trav.migr.preg + RMR


  ## TOTAL daily winter activity cost (MJ/d)

  # RESTING & LACTATING
  E.wint.min<- E.trav.wint.min + RMR
  E.wint.max<- E.trav.wint.max + RMR
  E.wint.mean<- E.trav.wint.mean + RMR

  # PREGNANT
  E.wint.preg.min<- E.trav.wint.preg.min + RMR
  E.wint.preg.max<- E.trav.wint.preg.max + RMR
  E.wint.preg.mean<- E.trav.wint.preg.mean + RMR



  ## Reproduction costs:

  # Gestation:
  mass.calf.min<-1101-311 # kg (mean 1101 +/- SD 311;  Fortune et al. 2012)
  mass.calf.max<-1101+311
  mass.calf.mean<-1101

  #age<-0
  #mass.calf<-3169.39+1773.666*age # kg; weight vs age (Moore et al 2005)


  # Heat of gestation over period ~ 12 months (MJ)
  E.gest.min <- 0.001*(18421.9*mass.calf.min^1.2) # MJ
  E.gest.max <- 0.001*(18421.9*mass.calf.max^1.2) # MJ
  E.gest.mean <- 0.001*(18421.9*mass.calf.mean^1.2) # MJ

  #(E.gest.min <- 0.001*(18421.9*mass.calf^1.2)) # MJ

  # Lactation:

  # Milk energy transfer efficiency (90%, so mother provides 110% of calf's energy req.)
  lact.eff<-1.1

  # Calf daily energy requirements (1767 +/-SD 261 MJ/d; Fortune et al. 2013)
  E.calf.min<-1767-261
  E.calf.max<-1767+261
  E.calf.mean<-1767

  t.calf<- seq(from = 1, to = 182.5, by = 1) # Days in 2nd half of lactation
  p.mom<- 1 - (t.calf/182.5) # Proportion of calf's daily energy requirements provided by mother during the 2nd phase of lactation (assuming a linear decline in contribution)

  l1.min<-(lact.eff*E.calf.min) * 182.5 # lactation cost (MJ) during 1st half of lactation (0-6 months)
  l1.max<-(lact.eff*E.calf.max) * 182.5
  l1.mean<-(lact.eff*E.calf.mean) * 182.5

  l2.min<-sum(p.mom * (lact.eff*E.calf.min)) # lactation cost (MJ) during weaning phase (6-12 months)
  l2.max<-sum(p.mom * (lact.eff*E.calf.max))
  l2.mean<-sum(p.mom * (lact.eff*E.calf.mean))

  E.lact.min.d<-(l1.min + l2.min) / 365 # lactation cost per day
  E.lact.max.d<-(l1.max + l2.max) / 365
  E.lact.mean.d<-(l1.mean + l2.mean) / 365

  E.lact.min<-l1.min + l2.min # lactation cost per year
  E.lact.max<-l1.max + l2.max
  E.lact.mean<-l1.mean + l2.mean



  ## Seasonal time-activity budget (days per season/year):

  #  Migration time (d, one way); Firestone et al. 2008
  t.migr.min<-21*2
  t.migr.max<-24*2
  t.migr.mean<-22.5*2


  # Winter time (d)

  # Fortune et al. 2013: Minimum residency times on winter grounds
  # Non-lactating adult females (mean ? SD)= 23.75 ? 18.60 d
  # Lactating females = 46.32 ? 14.60 d

  # Krzystan et al. 2018: Modeled residency times on winter grounds
  # Non-calving adult females range of means from 2004-2011 = 26.1-55.5 d
  # Lactating females range of means from 2004-2011 = 78.9-99.6 d

  # Resting/Pregnant females
  t.wint.repr.min.F<-23.75-18.6  # Fortune et al. 2013
  t.wint.repr.max.F<-23.75+18.6
  t.wint.repr.mean.F<-23.75

  t.wint.repr.min.K<-26.1  # Krzystan et al. 2018
  t.wint.repr.max.K<-55.5
  t.wint.repr.mean.K<-41.1

  # Lactating
  t.wint.la.min.F<-46.32-14.6
  t.wint.la.max.F<-46.32+14.6
  t.wint.la.mean.F<-46.32

  t.wint.la.min.K<-78.9
  t.wint.la.max.K<-99.6
  t.wint.la.mean.K<-87.5

  # Summer time (d):

  # Resting/Pregnant females
  t.sum.repr.min.F<-365 - (t.migr.max + t.wint.repr.max.F)
  t.sum.repr.max.F<-365 - (t.migr.min + t.wint.repr.min.F)
  t.sum.repr.mean.F<-365 - (t.migr.mean + t.wint.repr.mean.F)
  t.sum.repr.min.K<-365 - (t.migr.max + t.wint.repr.max.K)
  t.sum.repr.max.K<-365 - (t.migr.min + t.wint.repr.min.K)
  t.sum.repr.mean.K<-365 - (t.migr.mean + t.wint.repr.mean.K)

  # Lactating
  t.sum.la.min.F<-365 - (t.migr.max + t.wint.la.max.F)
  t.sum.la.max.F<-365 - (t.migr.min + t.wint.la.min.F)
  t.sum.la.mean.F<-365 - (t.migr.mean + t.wint.la.mean.F)
  t.sum.la.min.K<-365 - (t.migr.max + t.wint.la.max.K)
  t.sum.la.max.K<-365 - (t.migr.min + t.wint.la.min.K)
  t.sum.la.mean.K<-365 - (t.migr.mean + t.wint.la.mean.K)


  ###########################
  # Annual daily energy output: daily seasonal costs * total days / no. feeding days
  ###########################

  ### Resting: ###
  dfbt <- dfbt %>% dplyr::mutate(
                          # Min, with Fortune's (F) estimates of winter residency times
                          Eo.rest.min_BT.F = ((E.summ.2_BT.rest * t.sum.repr.max.F) +                    # min. foraging drag (2), max time spent on feeding ground
                                                     (E.migr           * t.migr.min) +                          # min time spent migrating
                                                     (E.wint.min       * t.wint.repr.min.F)) / t.sum.repr.max.F, # min winter enery costs, min time spent on breeding ground, all divided by the max time spent on feeding ground
                          # Min, with Krzystan's (K) estimates of winter residency times
                          Eo.rest.min_BT.K = ((E.summ.2_BT.rest * t.sum.repr.max.K) +
                                                     (E.migr           * t.migr.min) +
                                                     (E.wint.min       * t.wint.repr.min.K)) / t.sum.repr.max.K,
                          # Max, with Fortune's (F) estimates of winter residency times
                          Eo.rest.max_BT.F = ((E.summ.3_BT.rest * t.sum.repr.min.F) +                    # max foraging drag, min time spent on feeding ground
                                                     (E.migr           * t.migr.max) +                          # max time spent migrating
                                                     (E.wint.max       * t.wint.repr.max.F)) / t.sum.repr.min.F, # max winter enery costs, max time spent on breeding ground, all divided by the min time spent on feeding ground
                          # Max, with Krzystan's (K) estimates of winter residency times
                          Eo.rest.max_BT.K = ((E.summ.3_BT.rest * t.sum.repr.min.K) +
                                                     (E.migr           * t.migr.max) +
                                                     (E.wint.max       * t.wint.repr.max.K)) / t.sum.repr.min.K,
                          # Mean, with Fortune's (F) estimates of winter residency times
                          Eo.rest.mean_BT.F = ((E.summ.2.5_BT.rest * t.sum.repr.mean.F) +                    # max foraging drag, min time spent on feeding ground
                                                      (E.migr           * t.migr.mean) +                          # max time spent migrating
                                                      (E.wint.mean       * t.wint.repr.mean.F)) / t.sum.repr.mean.F, # max winter enery costs, max time spent on breeding ground, all divided by the min time spent on feeding ground
                          # Mean, with Krzystan's (K) estimates of winter residency times
                          Eo.rest.mean_BT.K = ((E.summ.2.5_BT.rest * t.sum.repr.mean.K) +
                                                      (E.migr           * t.migr.mean) +
                                                      (E.wint.mean       * t.wint.repr.mean.K)) / t.sum.repr.mean.K,
                          ### Pregnant: ###

                          # Min, with Fortune's (F) estimates of winter residency times
                          Eo.preg.min_BT.F = ((E.summ.2_BT.preg * t.sum.repr.max.F) +                    # min. foraging drag (2), max time spent on feeding ground
                                                     (E.migr.preg      * t.migr.min) +                          # min time spent migrating
                                                     (E.wint.min       * t.wint.repr.min.F) + E.gest.min) / t.sum.repr.max.F, # min winter energy costs, min time spent on breeding ground, all divided by the max time spent on feeding ground
                          # Min, with Krzystan's (K) estimates of winter residency times
                          Eo.preg.min_BT.K = ((E.summ.2_BT.preg * t.sum.repr.max.K) +
                                                     (E.migr.preg      * t.migr.min) +
                                                     (E.wint.min       * t.wint.repr.min.K) + E.gest.min) / t.sum.repr.max.K,
                          # Max, with Fortune's (F) estimates of winter residency times
                          Eo.preg.max_BT.F = ((E.summ.3_BT.preg * t.sum.repr.min.F) +                    # max foraging drag (3), min time spent on feeding ground
                                                     (E.migr.preg      * t.migr.max) +                          # max time spent migrating
                                                     (E.wint.max       * t.wint.repr.max.F) + E.gest.max) / t.sum.repr.min.F, # max winter energy costs, max time spent on breeding ground, all divided by the min time spent on feeding ground
                          # Max, with Krzystan's (K) estimates of winter residency times
                          Eo.preg.max_BT.K = ((E.summ.3_BT.preg * t.sum.repr.min.K) +
                                                     (E.migr.preg      * t.migr.max) +
                                                     (E.wint.max       * t.wint.repr.max.K) + E.gest.max) / t.sum.repr.min.K,
                          # Mean, with Fortune's (F) estimates of winter residency times
                          Eo.preg.mean_BT.F = ((E.summ.2.5_BT.preg * t.sum.repr.mean.F) +                    # max foraging drag (3), min time spent on feeding ground
                                                      (E.migr.preg      * t.migr.mean) +                          # max time spent migrating
                                                      (E.wint.mean       * t.wint.repr.mean.F) + E.gest.mean) / t.sum.repr.mean.F, # max winter energy costs, max time spent on breeding ground, all divided by the min time spent on feeding ground
                          # Mean, with Krzystan's (K) estimates of winter residency times
                          Eo.preg.mean_BT.K = ((E.summ.2.5_BT.preg * t.sum.repr.mean.K) +
                                                      (E.migr.preg      * t.migr.mean) +
                                                      (E.wint.mean       * t.wint.repr.mean.K) + E.gest.mean) / t.sum.repr.mean.K,
                          ### Lactating: ###

                          # Min, with Fortune's (F) estimates of winter residency times
                          Eo.lact.min_BT.F = ((E.summ.2_BT.lact * t.sum.la.max.F) +
                                                     (E.migr           * t.migr.min) +
                                                     (E.wint.min       * t.wint.la.min.F) + E.lact.min) / t.sum.la.max.F,
                          # Min, with Krzystan's (K) estimates of winter residency times
                          Eo.lact.min_BT.K = ((E.summ.2_BT.lact * t.sum.la.max.K) +
                                                     (E.migr           * t.migr.min) +
                                                     (E.wint.min       * t.wint.la.min.K) + E.lact.min) / t.sum.la.max.K,
                          # Max, with Fortune's (F) estimates of winter residency times
                          Eo.lact.max_BT.F = ((E.summ.3_BT.lact * t.sum.la.min.F) +
                                                     (E.migr           * t.migr.max) +
                                                     (E.wint.max       * t.wint.la.max.F) + E.lact.max) / t.sum.la.min.F,
                          # Max, with Krzystan's (K) estimates of winter residency times
                          Eo.lact.max_BT.K = ((E.summ.3_BT.lact * t.sum.la.min.K) +
                                                     (E.migr           * t.migr.max) +
                                                     (E.wint.max       * t.wint.la.max.K) + E.lact.max) / t.sum.la.min.K,
                          # Mean, with Fortune's (F) estimates of winter residency times
                          Eo.lact.mean_BT.F = ((E.summ.2.5_BT.lact * t.sum.la.mean.F) +
                                                      (E.migr           * t.migr.mean) +
                                                      (E.wint.mean       * t.wint.la.mean.F) + E.lact.mean) / t.sum.la.mean.F,
                          # Mean, with Krzystan's (K) estimates of winter residency times
                          Eo.lact.mean_BT.K = ((E.summ.2.5_BT.lact * t.sum.la.mean.K) +
                                                      (E.migr           * t.migr.mean) +
                                                      (E.wint.mean       * t.wint.la.mean.K) + E.lact.mean) / t.sum.la.mean.K
  )



  #####################
  ## Annual energy input
  ## or Ingestion rate:
  #####################
  # Parameters:
  # Mouth area (m2)
  Am.min<-1.7 # For NARW measuring between 13.4 and 14.1m in length; van der Hoop et al. 2018 (WP3; Table 1)
  Am.max<-1.9
  Am.mean<-1.8

  # Bottom speed in m/s
  Ub<-1.0

  dfbt <- dfbt %>% dplyr::mutate(#Time spent ingesting (sec per day); summary presented as hours/day
                          Tb_BT.rest.min = forage_t_rest.min,
                          Tb_BT.rest.max = forage_t_rest.max,
                          Tb_BT.rest.mean = forage_t_rest.mean,

                          Tb_BT.preg.min = forage_t_prla.min,
                          Tb_BT.preg.max = forage_t_prla.max,
                          Tb_BT.preg.mean = forage_t_prla.mean,

                          Tb_BT.lact.min = forage_t_prla.min,# same time spent foraging as pregnant
                          Tb_BT.lact.max = forage_t_prla.max,
                          Tb_BT.lact.mean = forage_t_prla.mean
  )

  # Prey energy density (27.9 ? 5.0 kJ/g (Davies et al. 2012); converted to MJ/g)
  min.ed <- mean.ed - sd.ed
  max.ed <- mean.ed + sd.ed

  Ep.min<- min.ed * 0.001
  Ep.max<-max.ed * 0.001
  Ep.mean<-mean.ed * 0.001

  # Assimilation efficiency (after fecal & urinary loss; Swaim et al. 2009)
  eps.a.min<-0.80
  eps.a.max<-0.92
  eps.a.mean<-0.86


  #efficacit? de capture, environ 50% sur les petite esp?ces
  if(unique(input_file$Taxa) %in% c("psca","temo")) dfbt$eps.c <-ifelse(dfbt$MeshSize==333, 1, 0.5)

  if(!unique(input_file$Taxa) %in% c("psca","temo", "Tras","Mnor")) dfbt$eps.c <- 1

  dfbt <- dfbt %>%  dplyr::mutate( # Prey density (g/m3)
                            Dp_BT = DW_Z_gm3,
                            # Energy intake (MJ/day):
                            # RESTING
                            Ei.min_BT.rest = (Am.min * Ub * Tb_BT.rest.min * Ep.min * Dp_BT) * eps.a.min * eps.c,
                            Ei.max_BT.rest = (Am.max * Ub * Tb_BT.rest.max * Ep.max * Dp_BT) * eps.a.max * eps.c,
                            Ei.mean_BT.rest = (Am.mean * Ub * Tb_BT.rest.mean * Ep.mean * Dp_BT) * eps.a.mean * eps.c,
                            # PREGNANT
                            Ei.min_BT.preg = (Am.min * Ub * Tb_BT.preg.min * Ep.min * Dp_BT) * eps.a.min * eps.c,
                            Ei.max_BT.preg = (Am.max * Ub * Tb_BT.preg.max * Ep.max * Dp_BT) * eps.a.max * eps.c,
                            Ei.mean_BT.preg = (Am.mean * Ub * Tb_BT.preg.mean * Ep.mean * Dp_BT) * eps.a.mean * eps.c,
                            # LACTATING
                            Ei.min_BT.lact = (Am.min * Ub * Tb_BT.lact.min * Ep.min * Dp_BT) * eps.a.min * eps.c,
                            Ei.max_BT.lact = (Am.max * Ub * Tb_BT.lact.max * Ep.max * Dp_BT) * eps.a.max * eps.c,
                            Ei.mean_BT.lact = (Am.mean * Ub * Tb_BT.lact.mean * Ep.mean * Dp_BT) * eps.a.mean * eps.c
                            )
  ###################
  # NET ENERGY RATIO: (Ei-Eo)/Eo
  # expressed as a proportion of Eo
  ###################

  dfbt <- dfbt %>%  dplyr::mutate(# RESTING:
                            Enet.rest.min.F = (Ei.min_BT.rest - Eo.rest.max_BT.F) / Eo.rest.max_BT.F,
                            Enet.rest.max.F = (Ei.max_BT.rest - Eo.rest.min_BT.F) / Eo.rest.min_BT.F,
                            Enet.rest.mean.F = (Ei.mean_BT.rest - Eo.rest.mean_BT.F) / Eo.rest.mean_BT.F,

                            Enet.rest.min.K = (Ei.min_BT.rest - Eo.rest.max_BT.K) / Eo.rest.max_BT.K,
                            Enet.rest.max.K = (Ei.max_BT.rest - Eo.rest.min_BT.K) / Eo.rest.min_BT.K,
                            Enet.rest.mean.K = (Ei.mean_BT.rest - Eo.rest.mean_BT.K) / Eo.rest.mean_BT.K,
                           # PREGNANT:
                            Enet.preg.min.F = (Ei.min_BT.preg - Eo.preg.max_BT.F) /Eo.preg.max_BT.F,
                            Enet.preg.max.F = (Ei.max_BT.preg - Eo.preg.min_BT.F) / Eo.preg.min_BT.F,
                            Enet.preg.mean.F = (Ei.mean_BT.preg - Eo.preg.mean_BT.F) / Eo.preg.mean_BT.F,

                            Enet.preg.min.K = (Ei.min_BT.preg - Eo.preg.max_BT.K) / Eo.preg.max_BT.K,
                            Enet.preg.max.K = (Ei.max_BT.preg - Eo.preg.min_BT.K) / Eo.preg.min_BT.K,
                            Enet.preg.mean.K = (Ei.mean_BT.preg - Eo.preg.mean_BT.K) / Eo.preg.mean_BT.K,
                            # LACTATING:
                            Enet.lact.min.F = (Ei.min_BT.lact - Eo.lact.max_BT.F) / Eo.lact.max_BT.F,
                            Enet.lact.max.F = (Ei.max_BT.lact - Eo.lact.min_BT.F) / Eo.lact.min_BT.F,
                            Enet.lact.mean.F = (Ei.mean_BT.lact - Eo.lact.mean_BT.F) / Eo.lact.mean_BT.F,

                            Enet.lact.min.K = (Ei.min_BT.lact - Eo.lact.max_BT.K) / Eo.lact.max_BT.K,
                            Enet.lact.max.K = (Ei.max_BT.lact - Eo.lact.min_BT.K) / Eo.lact.min_BT.K,
                            Enet.lact.mean.K = (Ei.mean_BT.lact - Eo.lact.mean_BT.K) / Eo.lact.mean_BT.K
                            )

  ####################################
  # Minimum prey density requirement:
  # replace Ei with Eo, and solve for Dp
  ####################################
  dfbt <- dfbt %>%  dplyr::mutate(# RESTING: Dp max refers to the density that maximize Enet it is lower thatn Dp min.
                            Dp.rest.max.F = (Eo.rest.min_BT.F) / (Am.max*Ub*Tb_BT.rest.max*Ep.max*eps.a.max),
                            Dp.rest.min.F = (Eo.rest.max_BT.F) / (Am.min*Ub*Tb_BT.rest.min*Ep.min*eps.a.min),
                            Dp.rest.mean.F = (Eo.rest.mean_BT.F) / (Am.mean*Ub*Tb_BT.rest.mean*Ep.mean*eps.a.mean),

                            Dp.rest.max.K = (Eo.rest.min_BT.K) / (Am.max*Ub*Tb_BT.rest.max*Ep.max*eps.a.max),
                            Dp.rest.min.K = (Eo.rest.max_BT.K) / (Am.min*Ub*Tb_BT.rest.min*Ep.min*eps.a.min),
                            Dp.rest.mean.K = (Eo.rest.mean_BT.K) / (Am.mean*Ub*Tb_BT.rest.mean*Ep.mean*eps.a.mean),
                            # PREGNANT:
                            Dp.preg.max.F = (Eo.preg.min_BT.F) / (Am.max*Ub*Tb_BT.preg.max*Ep.max*eps.a.max),
                            Dp.preg.min.F = (Eo.preg.max_BT.F) / (Am.min*Ub*Tb_BT.preg.min*Ep.min*eps.a.min),
                            Dp.preg.mean.F = (Eo.preg.mean_BT.F) / (Am.mean*Ub*Tb_BT.preg.mean*Ep.mean*eps.a.mean),

                            Dp.preg.max.K = (Eo.preg.min_BT.K) / (Am.max*Ub*Tb_BT.preg.max*Ep.max*eps.a.max),
                            Dp.preg.min.K = (Eo.preg.max_BT.K) / (Am.min*Ub*Tb_BT.preg.min*Ep.min*eps.a.min),
                            Dp.preg.mean.K = (Eo.preg.mean_BT.K) / (Am.mean*Ub*Tb_BT.preg.mean*Ep.mean*eps.a.mean),
                            # LACTATING:
                            Dp.lact.max.F = (Eo.lact.min_BT.F) / (Am.max*Ub*Tb_BT.lact.max*Ep.max*eps.a.max),
                            Dp.lact.min.F = (Eo.lact.max_BT.F) / (Am.min*Ub*Tb_BT.lact.min*Ep.min*eps.a.min),
                            Dp.lact.mean.F = (Eo.lact.mean_BT.F) / (Am.mean*Ub*Tb_BT.lact.mean*Ep.mean*eps.a.mean),

                            Dp.lact.min.K = (Eo.lact.min_BT.K) / (Am.max*Ub*Tb_BT.lact.max*Ep.max*eps.a.max),
                            Dp.lact.max.K = (Eo.lact.max_BT.K) / (Am.min*Ub*Tb_BT.lact.min*Ep.min*eps.a.min),
                            Dp.lact.mean.K = (Eo.lact.mean_BT.K) / (Am.mean*Ub*Tb_BT.lact.mean*Ep.mean*eps.a.mean))


   dfbt<-dfbt[, c(ID, "Zsampled",
                   "Energy_kJ_g_Z",     "Taxa",
                   "DW_Z",              "DW_Z_gm3",
                   "Eo.rest.min_BT.K",  "Eo.rest.max_BT.K",  "Eo.rest.mean_BT.K",
                   "Eo.preg.min_BT.K",  "Eo.preg.max_BT.K",  "Eo.preg.mean_BT.K",
                   "Eo.lact.min_BT.K",  "Eo.lact.max_BT.K",  "Eo.lact.mean_BT.K",
                   "Ei.min_BT.rest",    "Ei.max_BT.rest",    "Ei.mean_BT.rest",
                   "Ei.min_BT.preg",    "Ei.max_BT.preg",    "Ei.mean_BT.preg",
                   "Ei.min_BT.lact",    "Ei.max_BT.lact",    "Ei.mean_BT.lact",
                   "Enet.rest.min.K",   "Enet.rest.max.K",   "Enet.rest.mean.K" ,
                   "Enet.preg.min.K",   "Enet.preg.max.K",   "Enet.preg.mean.K" ,
                   "Enet.lact.min.K",   "Enet.lact.max.K",   "Enet.lact.mean.K",
                  "Dp.rest.min.K",     "Dp.rest.max.K",    "Dp.rest.mean.K",
                  "Dp.preg.min.K",     "Dp.preg.max.K",    "Dp.preg.mean.K",
                  "Dp.lact.min.K",     "Dp.lact.max.K",    "Dp.lact.mean.K")]



  # Clean up Inf values, convert to NA
  do.call(data.frame,lapply(dfbt,  function(x) replace(x, is.infinite(x),NA)))

  return(dfbt)

}
