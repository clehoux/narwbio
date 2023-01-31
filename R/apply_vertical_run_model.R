
#' Transform Calanus biomass into Enet for NARW
#'
#' @param dat dataframe to import
#' @param cfin name of the column (quoted) with Cfin biomass
#' @param chyp name of the column (quoted) with Chyp biomass
#' @param psca name of the column (quoted) with Pseudocalanus spp. biomass
#' @param temo  name of the column (quoted) with temora biomass
#' @param cfinZ default to NULL, optionnaly if the vertical distribution doesn't not need to be explored, name of the column for depth of the maximum of Cfin biomass.
#' @param chypZ default to NULL, optionnaly if the vertical distribution doesn't not need to be explored, name of the column for depth of the maximum of Chyp biomass.
#' @param pscaZ default to NULL, optionnaly if the vertical distribution doesn't not need to be explored, name of the column for depth of the maximum of pseudocalanus biomass.
#' @param temoZ default to NULL, optionnaly if the vertical distribution doesn't not need to be explored, name of the column for depth of the maximum of temora biomass.
#' @param units units of the biomass, one of "ug", "mg", "g", "kJ", default to mg
#' @param ID name of the column (quoted) for the unique ID of samples
#' @param change_150m does the NARW behavior change when it dives below 150m. Default to T  see details
#' @param state decide to extract results for "resting", "pregnant" or "lactating" females. default to "pregnant". At the moment does not accommodate more that one state.
#' @param param  decide to extract results for "min", "mean" or "max" enet. default to "mean". At the moment does not accommodate more that one state.
#' @param split_80 choose whether the model is applied on the whole water column or on the layer above and below 80 m separatly
#'
#'
#'
#' @details Z columns applies for Nicolas Lecorre outputs. I don't see other possible usage.
#' The parameter change_150m is exploratory. In Gavrilchuk et al. and Lehoux et al. the whale change its behaviour below 150m , the dive has a maximum duration and thus time spent feedding decrease with feeding depth. If FALSE, then feeding duration decrease with depth for all depth. I strongly suggest using TRUE except for exploration.
#' @return the same data frame as the input, with added columns for NARW bioenergy.
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select mutate rename left_join full_join filter group_by summarise slice ungroup
#' @importFrom magrittr %>%
#' @importFrom utils menu
#'
#' @export
#'
#' @examples
#' #TO DO
 execute_bioenergy <- function(dat, cfin = "Z_Cf.Cg4.6_mgm2",
                              chyp = "Z_Ch4.6_mgm2",
                              psca = "Z_Pcal.Micro.Scol_mgm2",
                              temo = "Z_Tem.Eury_mgm2",
                              cfinZ = NULL,
                              chypZ = NULL,
                              pscaZ = NULL,
                              temoZ = NULL,
                              units = "mg",
                              ID = "Label",
                              change_150m = T,
                              state="pregnant",
                              param="mean",
                              split_80=F) {
  dat <- as.data.frame(dat) # en tibble fait des trucs bizarre

  if(!state %in% c("resting", "pregnant", "lactating")) stop("state must be one of resting, pregnant or lactating")
  if(!param %in% c("min", "mean", "max")) stop("param must be one of min, mean, max")


  if (length(unique(dat[, ID])) != nrow(dat)) {
    warning("IDs are not unique, a new ID column was created for adequate join")
    dat$IDcol <- seq(1, nrow(dat), 1)

    ID <- "IDcol"
  }

  colnames_out <- colnames(dat)


  if (!units %in% c("ug", "mg", "g", "kJ")) stop("units should be in ug, mg, g, kJ")

  if (is.null(cfinZ) & is.null(chypZ) & is.null(pscaZ) & is.null(temoZ)) {
    get_vert <- T
  }
  if (!is.null(cfinZ) | !is.null(chypZ) | !is.null(pscaZ) | !is.null(temoZ)) {
    get_vert <- F
  }

  if (is.null(cfin)) {
    dat$cfin <- NA
    dat$cfinz <- NA
  }
  if (!is.null(cfin) & !is.null(cfinZ)) {
    dat$cfin <- dat[, cfin]
    dat$cfinz <- dat[, cfinZ]
  }
  if (!is.null(cfin) & is.null(cfinZ)) {
    dat$cfin <- dat[, cfin]
    dat$cfinz <- NA
  }
  if (is.null(chyp)) {
    dat$chyp <- NA
    dat$chypz <- NA
  }
  if (!is.null(chyp) & !is.null(chypZ)) {
    dat$chyp <- dat[, chyp]
    dat$chypz <- dat[, chypZ]
  }
  if (!is.null(chyp) & is.null(chypZ)) {
    dat$chyp <- dat[, chyp]
    dat$chypz <- NA
  }
  if (is.null(psca)) {
    dat$psca <- NA
    dat$pscaz <- NA
  }
  if (!is.null(psca) & !is.null(pscaZ)) {
    dat$psca <- dat[, psca]
    dat$pscaz <- dat[, pscaZ]
  }
  if (!is.null(psca) & is.null(pscaZ)) {
    dat$psca <- dat[, psca]
    dat$pscaz <- NA
  }
  if (is.null(temo)) {
    dat$temo <- NA
    dat$temoz <- NA
  }
  if (!is.null(temo) & !is.null(temoZ)) {
    dat$temo <- dat[, temo]
    dat$temoz <- dat[, temoZ]
  }
  if (!is.null(temo) & is.null(temoZ)) {
    dat$temo <- dat[, temo]
    dat$temoz <- NA
  }


  if (units == "ug") {
    dat <- dat %>% dplyr::mutate(
      cfin_kj = ifelse(!is.na(cfin), cfin / 1000 / 1000 * 27.9, NA),
      chyp_kj = ifelse(!is.na(chyp), chyp / 1000 / 1000 * 27.9, NA),
      psca_kj = ifelse(!is.na(psca), psca / 1000 / 1000 * 22.73, NA),
      temo_kj = ifelse(!is.na(temo), temo / 1000 / 1000 * 27.9, NA)
    )
  }


  if (units == "mg") {
    dat <- dat %>% dplyr::mutate(
      cfin_kj = ifelse(!is.na(cfin), (cfin / 1000) * 27.9, NA),
      chyp_kj = ifelse(!is.na(chyp), (chyp / 1000) * 27.9, NA),
      psca_kj = ifelse(!is.na(psca), (psca / 1000) * 22.73, NA),
      temo_kj = ifelse(!is.na(temo), (temo / 1000) * 19.72, NA)
    )
  }

  if (units == "g") {
    dat <- dat %>% dplyr::mutate(
      cfin_kj = ifelse(!is.na(cfin), cfin * 27.9, NA),
      chyp_kj = ifelse(!is.na(chyp), chyp * 27.9, NA),
      psca_kj = ifelse(!is.na(psca), psca * 22.73, NA),
      temo_kj = ifelse(!is.na(temo), temo * 19.72, NA)
    )
  }

  if (units == "kJ") {
    dat <- dat %>% dplyr::mutate(
      cfin_kj = ifelse(!is.na(cfin), cfin, NA),
      chyp_kj = ifelse(!is.na(chyp), chyp, NA),
      psca_kj = ifelse(!is.na(psca), psca, NA),
      temo_kj = ifelse(!is.na(temo), temo, NA)
    )
  }
  if (!"MeshSize" %in% colnames(dat)) test <- utils::menu(c("200", "333"), title = "There is no mesh size in the data. Which mesh size should the function use? This affects the filtering efficiency for small calanus.")

  if (!"MeshSize" %in% colnames(dat)) dat <- dat %>% mutate(MeshSize = ifelse(test == 1, 200, 333))





  dat2 <- dat
  # do you want to estimate the vertical distribution
  if (get_vert == T) {
    dat3 <- dat2
    if (!"season" %in% colnames(dat)) {
      warning("To estimate the vertical distribution, you must include column for season or Month. Season is missing and set to NA and can lead to missing results")
      dat3 <- dat3 %>% mutate(season = NA)

    }

 if (!"Month" %in% colnames(dat)) {
   warning("To estimate the vertical distribution, you must include column for season or Month. Month is missing and set to NA and can lead to missing results")
   dat3 <- dat3 %>% mutate(Month = NA)

 }
    if (!"Region" %in% colnames(dat)) {
      warning("To estimate the vertical distribution outside of the GSL, you must include column for Region. Region is set to Gulf of St Lawrence (gsl) and can lead to misleading results")
      dat3 <- dat3 %>% mutate(Region = "gsl")
    }

    if (!"Day_Night" %in% colnames(dat)) {
      warning("To estimate the vertical distribution according to Lehoux et al. 2020, you must include column for Day_Night. Day_Night is set to NA and can lead to missing results")
      dat3 <- dat3 %>% mutate(Day_Night = NA)
    }
    if (!"Zstation" %in% colnames(dat)) {
      warning("To estimate the vertical distribution, you must include column for Zstation. Zstation is set to NA and can lead to missing results")
      dat3 <- dat3 %>% dplyr::mutate(Zstation = NA)
    }

    if(!"Zdev_Ch4.6_mgm2" %in% colnames(dat3)) dat3$Zdev_Ch4.6_mgm2 <- NA

    dat3 <- dat3 %>%
      dplyr::mutate(
        dev_index = ifelse(!"dev_index" %in% colnames(dat) &
          !"civ.hyp_dev.index_station" %in% colnames(dat) &
          !"civ.hyp_dev.index_transect" %in% colnames(dat), NA,
        ifelse("dev_index" %in% colnames(dat), dev_index,
          ifelse("civ.hyp_dev.index_station" %in% colnames(dat) &
            !"civ.hyp_dev.index_transect" %in% colnames(dat), civ.hyp_dev.index_station,
          ifelse(!"civ.hyp_dev.index_station" %in% colnames(dat) &
            "civ.hyp_dev.index_transect" %in% colnames(dat), civ.hyp_dev.index_transect,
          ifelse("civ.hyp_dev.index_station" %in% colnames(dat) &
            "civ.hyp_dev.index_transect" %in% colnames(dat) &
            is.na(civ.hyp_dev.index_station), civ.hyp_dev.index_transect, civ.hyp_dev.index_station)
          )
          )
        )
        ),
        dev_index=ifelse(dev_index < 0.65, "Active",
                         ifelse(dev_index >= 0.65 & dev_index < 0.85, yes ="Transition",
                                ifelse(dev_index >= 0.85, "Diapause", NA)))# Zdev_Ch4.6_mgm2 = ifelse(!"Zdev_Ch4.6_mgm2" %in% colnames(dat), NA, Zdev_Ch4.6_mgm2)
      ) %>%
      dplyr::select(all_of(ID), season, Day_Night, Region,Month,Zstation, dev_index, cfin_kj, chyp_kj, psca_kj, temo_kj, Zdev_Ch4.6_mgm2)

    dat3 <- dplyr::rename(dat3, D_N = Day_Night)

    dat_out <- kj_df(datkj = dat3, ID=ID)



    dat2 <- suppressMessages(dplyr::full_join(dat2 %>% dplyr::select(-cfinz, -chypz, -temoz, -pscaz), dplyr::mutate(dat_out,
      Zsampled = as.numeric(Zsampled),
      cfin_kj_str = as.numeric(cfin_kj_str),
      chyp_kj_str = as.numeric(chyp_kj_str),
      temo_kj_str = as.numeric(temo_kj_str),
      psca_kj_str = as.numeric(psca_kj_str)
    )))

    dat2[dat2 == "Inf"] <- NA
  } # end of vertical distribution

  if (get_vert == F) {
    datz <- suppressMessages(dplyr::full_join(
      dat2 %>% dplyr::select(ID, MeshSize, cfinz, chypz, pscaz, temoz) %>% tidyr::pivot_longer(cols = c(cfinz, chypz, pscaz, temoz), names_to = "Taxa", values_to = "Zsampled") %>% dplyr::mutate(Taxa = substring(Taxa, 1, 4)),
      dat2 %>% dplyr::select(ID, cfin_kj, chyp_kj, psca_kj, temo_kj) %>% pivot_longer(cols = c(cfin_kj, chyp_kj, psca_kj, temo_kj), names_to = "Taxa", values_to = "Energy_kJ_g_Z") %>% dplyr::mutate(Taxa = substring(Taxa, 1, 4))
    ))
  }

  if (get_vert == T) {
    datz <- suppressMessages(dplyr::full_join(
      dat2 %>% dplyr::select(all_of(ID), MeshSize, Zsampled),
      dat2 %>% dplyr::select(all_of(ID), Zsampled, cfin_kj_str, chyp_kj_str, psca_kj_str, temo_kj_str) %>% tidyr::pivot_longer(cols = c(cfin_kj_str, chyp_kj_str, psca_kj_str, temo_kj_str), names_to = "Taxa", values_to = "Energy_kJ_g_Z") %>% dplyr::mutate(Taxa = substring(Taxa, 1, 4))
    ))
  }

  # dat_long<- dplyr::left_join(dat_long, dat[,c("Label", "MeshSize")])

  cfin_b <- NARW_bionerg(input_file = subset(datz, Taxa == "cfin"), mean.ed = 27.9, sd.ed = 5.0, change_150m = change_150m, ID = ID)
  chyp_b <- NARW_bionerg(input_file = subset(datz, Taxa == "chyp"), mean.ed = 27.9, sd.ed = 5.0, change_150m = change_150m, ID = ID)
  psca_b <- NARW_bionerg(input_file = subset(datz, Taxa == "psca"), mean.ed = 22.73, sd.ed = 0.65, change_150m = change_150m, ID = ID)
  temo_b <- NARW_bionerg(input_file = subset(datz, Taxa == "temo"), mean.ed = 19.72, sd.ed = 1.05, change_150m = change_150m, ID = ID)



  # join data
  sp_b <- rbind(cfin_b, chyp_b, temo_b, psca_b)
if(any(duplicated(sp_b[,c(ID,"Taxa", "Zsampled")]))) stop("duplicates found rework vertical distribution code")
  #sp_b<- sp_b[!duplicated(sp_b[,c(ID,"Taxa", "Zsampled")]),]# in some cases when the vertical distribution creates 2 times the same layer, should be resolve with +9

  # function for summarize
  sp.out <- sp_b %>% ungroup() %>%
    dplyr::filter(!is.na(Zsampled)) %>%
    dplyr::group_by(get(ID), Zsampled) %>%
    dplyr::summarise(
      across(starts_with("Eo"), ~ mean(.x, na.rm = TRUE)),
      across(starts_with("Ei"), ~ sum(.x, na.rm = TRUE)),
      Energy_kJ_g_Z = sum(Energy_kJ_g_Z, na.rm = T),
      DW_Z_gm3 = sum(DW_Z_gm3, na.rm = T),
      across(starts_with("Dp"), ~ mean(.x, na.rm = TRUE))
    ) %>%
    dplyr::mutate(Enet.preg.mean.K = (Ei.mean_BT.preg - Eo.preg.mean_BT.K) / Eo.preg.mean_BT.K,
                  Enet.rest.mean.K = (Ei.mean_BT.rest - Eo.rest.mean_BT.K) / Eo.rest.mean_BT.K,
                  Enet.lact.mean.K = (Ei.mean_BT.lact - Eo.lact.mean_BT.K) / Eo.lact.mean_BT.K,
                  Enet.preg.max.K = (Ei.max_BT.preg - Eo.preg.min_BT.K) / Eo.preg.min_BT.K,
                  Enet.rest.max.K = (Ei.max_BT.rest - Eo.rest.min_BT.K) / Eo.rest.min_BT.K,
                  Enet.lact.max.K = (Ei.max_BT.lact - Eo.lact.min_BT.K) / Eo.lact.min_BT.K,
                  Enet.preg.min.K = (Ei.min_BT.preg - Eo.preg.max_BT.K) / Eo.preg.max_BT.K,
                  Enet.rest.min.K = (Ei.min_BT.rest - Eo.rest.max_BT.K) / Eo.rest.max_BT.K,
                  Enet.lact.min.K = (Ei.min_BT.lact - Eo.lact.max_BT.K) / Eo.lact.max_BT.K,


                  )
  sp.out[, ID] <- sp.out$`get(ID)`


  #this is how each state and param shopuld be
  if (get_vert == T) {
    if(state=="pregnant" & param=="mean"){
      pat <-  "preg.mean" # for Enet and Eo and Dp
      pat2 <-  "mean_BT.preg" # for Ei
    }
    if(state=="pregnant" & param=="min"){
      pat <-  "preg.min" # for Enet and Eo and Dp
      pat2 <-  "min_BT.preg" # for Ei
     }
    if(state=="pregnant" & param=="max"){
      pat <-  "preg.max" # for Enet and Eo and Dp
      pat2 <-  "max_BT.preg" # for Ei
    }
    if(state=="resting" & param=="mean"){
      pat <-  "rest.mean" # for Enet and Eo and Dp
      pat2 <-  "mean_BT.rest" # for Ei
     }
    if(state=="resting" & param=="min"){
      pat <-  "rest.min" # for Enet and Eo and Dp
      pat2 <-  "min_BT.rest" # for Ei
      }
    if(state=="resting" & param=="max"){
      pat <-  "rest.max" # for Enet and Eo and Dp
      pat2 <-  "rest_BT.max" # for Ei
      }

    if(state=="lactating" & param=="mean"){
      pat <-  "lact.mean" # for Enet and Eo and Dp
      pat2 <-  "mean_BT.lact" # for Ei
     }
    if(state=="lactating" & param=="min"){
      pat <-  "lact.min" # for Enet and Eo and Dp
      pat2 <-  "min_BT.lact" # for Ei
     }
    if(state=="lactating" & param=="max"){
      pat <-  "lact.max" # for Enet and Eo and Dp
      pat2 <-  "max_BT.lact" # for Ei
    }


#general way of subsetting
      stEnet<-  stringr::str_subset(colnames(sp.out),pattern=paste("Enet", pat, sep="."))
      sp.outa <- sp.out %>%
      dplyr::group_by(get(ID)) %>%
      dplyr::slice(., which.max(get(stEnet)))
      cols_to_keepa=c(paste0("Ei.", pat2),paste0("Eo.",pat,"_BT.K"),paste0("Enet.",pat,".K"),paste0("Dp.",pat,".K"))

      if(split_80){
      sp.outm80 <- sp.out %>% dplyr::filter(Zsampled <= 80) %>%
      dplyr::group_by(get(ID)) %>%
      dplyr::slice(., which.max(get(stEnet))) %>%  dplyr::rename_with(starts_with(c("E", "Dp")), .fn=~ paste0(.x, "_m80")) %>%  dplyr::rename(ED.max_Zm80=Zsampled)

      sp.outp80 <- sp.out %>% dplyr::filter(Zsampled > 80) %>%
        dplyr::group_by(get(ID)) %>%
        dplyr::slice(., which.max(get(stEnet))) %>%  dplyr::rename_with(starts_with(c("E", "Dp")), .fn=~ paste0(.x, "_p80")) %>%  dplyr::rename(ED.max_Zp80=Zsampled)

      sp.outb <- dplyr::full_join(sp.outm80, sp.outp80, by="get(ID)")
      cols_to_keepb=c("ED.max_Zm80","ED.max_Zp80", c(paste0("Ei.", pat2,"_m80"),paste0("Eo.",pat,"_BT.K_m80"),paste0("Enet.",pat,".K_m80"),paste0("Dp.",pat,".K_m80")),
                     c(paste0("Ei.", pat2,"_p80"),paste0("Eo.",pat,"_BT.K_p80"),paste0("Enet.",pat,".K_p80"),paste0("Dp.",pat,".K_p80")))

      sp.outb[, ID] <- sp.outb$`get(ID)`
      }

    sp.outa[, ID] <- sp.outa$`get(ID)`
  }

# 80m not done from here, start here...................

if(state=="pregnant"){
  bioena <- sp.outa %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(ID), Energy_kJ_g_Z, Zsampled, all_of(cols_to_keepa)) %>%
    dplyr::rename(ED.max = "Energy_kJ_g_Z", ED.max_Z = "Zsampled", NARW_Eout.preg = grep(cols_to_keepa, pattern="Eo", value=T), NARW_Ein =grep(cols_to_keepa, pattern="Ei", value=T), NARW_Enet.preg = grep(cols_to_keepa, pattern="Enet", value=T), Minimum_prey_density_gm3=grep(cols_to_keepa, pattern="Dp", value=T))
if(split_80){
  bioenb <- sp.outb %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(ID), Energy_kJ_g_Z_m80,  Energy_kJ_g_Z_p80,all_of(cols_to_keepb)) %>%
    dplyr::rename(ED.max_m80 = "Energy_kJ_g_Z_m80",ED.max_p80 = "Energy_kJ_g_Z_p80", NARW_Eout.preg_m80 = grep(cols_to_keepb, pattern=paste0("Eo.",pat,"_BT.K_m80"), value=T, fixed=T) ,NARW_Eout.preg_p80 = grep(cols_to_keepb, pattern=paste0("Eo.",pat,"_BT.K_p80"), value=T, fixed=T),
                                                                                     NARW_Ein_m80 =grep(cols_to_keepb, pattern=paste0("Ei.", pat2,"_m80"), value=T),
                                                                                     NARW_Ein_p80 =grep(cols_to_keepb, pattern=paste0("Ei.", pat2,"_p80"), value=T),
                               NARW_Enet.preg_m80 = grep(cols_to_keepb, pattern=paste0("Enet.",pat,".K_m80"), value=T),
                               NARW_Enet.preg_p80 = grep(cols_to_keepb, pattern=paste0("Enet.",pat,".K_p80"), value=T),
                               Minimum_prey_density_gm3_m80=grep(cols_to_keepb, pattern=paste0("Dp.",pat,".K_m80"), value=T),
                                Minimum_prey_density_gm3_p80=grep(cols_to_keepb, pattern=paste0("Dp.",pat,".K_p80"), value=T))
 }
}

  if(state!="pregnant"){
    bioena <- sp.outa %>%
      dplyr::ungroup() %>%
      dplyr::select(all_of(ID), Energy_kJ_g_Z, Zsampled, all_of(cols_to_keepa)) %>%
      dplyr::rename(ED.max = "Energy_kJ_g_Z", ED.max_Z = "Zsampled", NARW_Eout = grep(cols_to_keepa, pattern="Eo", value=T), NARW_Ein =grep(cols_to_keepa, pattern="Ei", value=T), NARW_Enet = grep(cols_to_keepa, pattern="Enet", value=T), Minimum_prey_density_gm3=grep(cols_to_keepa, pattern="Dp", value=T))

    if(split_80){
      bioenb <- sp.outb %>%
        dplyr::ungroup() %>%
        dplyr::select(all_of(ID), Energy_kJ_g_Z_m80,  Energy_kJ_g_Z_p80,Zsampled, all_of(cols_to_keepb)) %>%
        dplyr::rename(ED.max_m80 = "Energy_kJ_g_Z_m80",ED.max_p80 = "Energy_kJ_g_Z_p80", NARW_Eout_m80 = grep(cols_to_keepb, pattern=paste0("Eo.",pat,"_BT.K_m80"), value=T, fixed=T) ,NARW_Eout_p80 = grep(cols_to_keepb, pattern=paste0("Eo.",pat,"_BT.K_p80"), value=T, fixed=T),
                      NARW_Ein_m80 =grep(cols_to_keepb, pattern=paste0("Ei.", pat2,"_m80"), value=T),
                      NARW_Ein_p80 =grep(cols_to_keepb, pattern=paste0("Ei.", pat2,"_p80"), value=T),
                      NARW_Enet_m80 = grep(cols_to_keepb, pattern=paste0("Enet.",pat,".K_m80"), value=T),
                      NARW_Enet_p80 = grep(cols_to_keepb, pattern=paste0("Enet.",pat,".K_p80"), value=T),
                      Minimum_prey_density_gm3_m80=grep(cols_to_keepb, pattern=paste0("Dp.",pat,".K_m80"), value=T),
                      Minimum_prey_density_gm3_p80=grep(cols_to_keepb, pattern=paste0("Dp.",pat,".K_p80"), value=T))
    }

    }


  if(split_80) bioen <-  full_join(bioena, bioenb)
  if(!split_80) bioen <-  bioena

  # si contient d?j? les colonnes
  if ("NARW_Enet.preg" %in% colnames_out) dat_final <- suppressMessages(dplyr::left_join(dat[, -which(colnames(dat) %in% colnames(bioen)[2:ncol(boen)])], bioen))
  # si ne contient pas les colonnes
  if (!"NARW_Enet.preg" %in% colnames_out) dat_final <- suppressMessages(dplyr::left_join(dat, bioen))
  # check nrow NA

  # si ne contient pas tous les colonnes
  if (!"NARW_Enet.preg" %in% colnames_out) dat_final <- dat_final[, c(colnames_out, colnames(bioen)[2:ncol(bioen)])]

  # si contient d?j? les colonnes
  if ("NARW_Enet.preg" %in% colnames_out) dat_final <- dat_final[, colnames_out]


  if (nrow(dat_final[duplicated(dat_final[, ID]), ]) != 0) message("duplicated labels")
  if (nrow(dat_final) != nrow(dat)) stop("join not done correctly")

  if(state=="pregnant") {
    dat_final$NARW_Enet.preg[dat_final$NARW_Enet.preg == -1] <- NA
    if(split_80){
      dat_final$NARW_Enet.preg_m80[dat_final$NARW_Enet.preg_m80 == -1] <- NA
      dat_final$NARW_Enet.preg_p80[dat_final$NARW_Enet.preg_p80 == -1] <- NA
    }

    }
  if(state!="pregnant"){
    dat_final$NARW_Enet[dat_final$NARW_Enet == -1] <- NA
    if(split_80){
      dat_final$NARW_Enet.preg_m80[dat_final$NARW_Enet.preg_m80 == -1] <- NA
      dat_final$NARW_Enet.preg_p80[dat_final$NARW_Enet.preg_p80 == -1] <- NA
    }

  }

  return(dat_final)
}






#' Apply vertical distribution to integrated energy
#'
#' @param datkj  dataframe
#' @param ID name of the column for unique identifier.
#' @return a dataframe with Zsampled for 3D distribution
#' @importFrom mgcv predict.gam
#' @export
#'
#' @examples
#'#to do
 kj_df <- function(datkj, ID) {
   output <- data.frame(matrix(ncol = 6))

   pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                        max = nrow(datkj), # Maximum value of the progress bar
                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                        width = 50,   # Progress bar width. Defaults to getOption("width")
                        char = "=")

   colnames(output) <- c(ID, "Zsampled", "cfin_kj_str", "chyp_kj_str", "temo_kj_str", "psca_kj_str")
   for (i in 1:nrow(datkj)) {
    datkji <- datkj[i,]

    setTxtProgressBar(pb, i)
     if ((is.na(datkji$season) | is.na(datkji$D_N)) &  is.na(datkji$Region)) {
       naout <- cbind(datkj[i, ID], NA, NA, NA, NA, NA)
       colnames(naout) <- colnames(output)
       output <- rbind(output, naout)
     }
     if ((is.na(datkji$season) | is.na(datkji$D_N)) &  tolower(datkji$Region) %in% c("gom", "wss", "ess", "gsl")) {
       # mod
       if (tolower(datkji$Region)== "gom") {
         datkj$fMonth<- as.factor(datkj$Month)
         mod.cfin <- Cfin_CV_CVI_GOM
         mod.chyp <- Cfin_CV_CVI_GOM #it is the same, we don't have one for GoM and using wss causes missing values.

       }
       if (tolower(datkji$Region) == "wss") {
         mod.cfin <- mcfin_glac_c_wss
         mod.chyp <- mchyp_c_wss
        }
       if (tolower(datkji$Region)=="ess") {
         mod.cfin <- mcfin_glac_c_ess
         mod.chyp <- mchyp_c_ess
        }
       if (tolower(datkji$Region)=="gsl") {
         mod.cfin <- mcfin_glac_c_gsl
         mod.chyp <- mchyp_c_gsl
        }

       ## cr?ation des pr?dictions th?oriques#
       Zsampled <- seq(0, datkji$Zstation + 9, 10) # Z sampled en intervalle de 10 #
       tjn <- as.data.frame(Zsampled)
       tjn[, ID] <- rep(datkj[i, ID], nrow(tjn))

       datx <- suppressMessages(dplyr::left_join(datkji, tjn))
       datx[nrow(datx), "Zsampled"] <- datx[nrow(datx), "Zstation"]


       datx$percZ_stn <- datx$Zsampled / datx$Zstation * 100

       datx$predcfin <- tryCatch( # if months are missing
         {mgcv::predict.gam(mod.cfin, datx, type = "response", exclude = "s(ID)", newdata.guaranteed = T) }, error = function(e) {return(NA)})
       maxcfin<- max(datx$predcfin)
       if(!is.na(maxcfin)){if(maxcfin < 0.85) {warning(paste0(ID,"maxcfin < 0.85"))}}
       datx$predcfin <-  datx$predcfin/maxcfin
       datx$propcfin <- c(NA, diff(datx$predcfin, lag = 1)) # lag de 1 parce que j'ai entr? en bloc de 10 donc valeur -1 est un bloc de 10
       datx$propcfin[datx$propcfin < 0] <- 0
       datx$cfin_kj_str <- datx$cfin_kj * datx$propcfin

       datx$predchyp <- tryCatch( # if months are missing
         {mgcv::predict.gam(mod.chyp, datx, type = "response", exclude = "s(ID)", newdata.guaranteed = T) }, error = function(e) {return(NA)})
       maxchyp<- max(datx$predchyp)
       if(!is.na(maxchyp)){if(maxchyp < 0.85) {warning(paste0(ID,"maxchyp < 0.85"))}}
       datx$predchyp <-  datx$predchyp/maxchyp
       datx$propchyp <- c(NA, diff(datx$predchyp, lag = 1)) # lag de 1 parce que j'ai entr? en bloc de 10 donc valeur -1 est un bloc de 10
       datx$propchyp[datx$propchyp < 0] <- 0
       datx$chyp_kj_str <- datx$chyp_kj * datx$propchyp

       datx$predtemo <- NA
       maxtemo<- max(datx$predtemo)
       datx$predtemo <-  datx$predtemo/maxtemo
       datx$proptemo <- c(NA, diff(datx$predtemo, lag = 1)) # lag de 1 parce que j'ai entr? en bloc de 10 donc valeur -1 est un bloc de 10
       datx$proptemo[datx$proptemo < 0] <- 0
       datx$temo_kj_str <- datx$temo_kj * datx$proptemo

       datx$predpsca <- NA
       maxpsca<- max(datx$predpsca)
       datx$predpsca <-  datx$predpsca/maxpsca
       datx$proppsca <- c(NA, diff(datx$predpsca, lag = 1)) # lag de 1 parce que j'ai entr? en bloc de 10 donc valeur -1 est un bloc de 10
       datx$proppsca[datx$proppsca < 0] <- 0
       datx$psca_kj_str <- datx$psca_kj * datx$proppsca
       output <- rbind(output, datx[, c(ID, "Zsampled", "cfin_kj_str", "chyp_kj_str", "temo_kj_str", "psca_kj_str")])
       }


     if (!is.na(datkji$season) & !is.na(datkji$D_N) & tolower(datkji$Region) %in% c("gsl", "sgsl", "ngsl")) {
       # mod
       if (datkji$season == "esum" & datkji$D_N == "D") {
         mod.cfin <- mcfin_glac_c_gsl_esum_D
         if (is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_c_gsl_esum_D
         if (!is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_dev_catD450
       }
       if (datkji$season == "esum" & datkji$D_N == "N") {
         mod.cfin <- mcfin_glac_c_gsl_esum_N
         if (is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_c_gsl_esum_N
         if (!is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_dev_catN450
       }
       if (datkji$season == "lsum" & datkji$D_N == "D") {
         mod.cfin <- mcfin_glac_c_gsl_lsum_D
         if (is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_c_gsl_lsum
         if (!is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_dev_catD450
       }
       if (datkji$season == "lsum" & datkji$D_N == "N") {
         mod.cfin <- mcfin_glac_c_gsl_lsum_N
         if (is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_c_gsl_lsum
         if (!is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_dev_catN450
       }
       if (datkji$D_N == "D") {
         mod.temo <- mtem_c_gsl_D
         mod.psca <- mpseudo_c_gsl_D
       }
       if (datkji$D_N == "N") {
         mod.temo <- mtem_c_gsl_N
         mod.psca <- mpseudo_c_gsl_N
       }
       if (datkji$season == "fall") {
         mod.cfin <- mcfin_glac_c_gsl_fall
         mod.chyp <- mchyp_c_gsl_fall
       }
       if (datkji$season == "fall" & datkji$D_N == "D") {
         if (!is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_dev_catD450
       }

       if (datkji$season == "fall" & datkji$D_N == "N") {
         if (!is.na(datkji$Zdev_Ch4.6_mgm2)) mod.chyp <- mchyp_dev_catN450
       }
       ## cr?ation des pr?dictions th?oriques#
       Zsampled <- seq(0, datkji$Zstation + 9, 10) # Z sampled en intervalle de 10 #
       tjn <- as.data.frame(Zsampled)
       tjn[, ID] <- rep(datkj[i, ID], nrow(tjn))

       datx <- suppressMessages(dplyr::left_join(datkji, tjn))
       datx[nrow(datx), "Zsampled"] <- datx[nrow(datx), "Zstation"]


       datx$percZ_stn <- datx$Zsampled / datx$Zstation * 100

       datx$predcfin <- mgcv::predict.gam(mod.cfin, datx, type = "response", exclude = "s(ID)", newdata.guaranteed = T)
       maxcfin<- max(datx$predcfin)
       if(maxcfin < 0.85) warning(paste0(ID,"maxcfin < 0.85"))
       datx$predcfin <-  datx$predcfin/maxcfin
       datx$propcfin <- c(NA, diff(datx$predcfin, lag = 1)) # lag de 1 parce que j'ai entr? en bloc de 10 donc valeur -1 est un bloc de 10
       datx$propcfin[datx$propcfin < 0] <- 0
       datx$cfin_kj_str <- datx$cfin_kj * datx$propcfin

       datx$predchyp <- mgcv::predict.gam(mod.chyp, datx, type = "response", exclude = "s(ID)", newdata.guaranteed = T)
       maxchyp<- max(datx$predchyp)
       if(maxchyp < 0.85) warning(paste0(ID,"maxchyp < 0.85"))
       datx$predchyp <-  datx$predchyp/maxchyp
       datx$propchyp <- c(NA, diff(datx$predchyp, lag = 1)) # lag de 1 parce que j'ai entr? en bloc de 10 donc valeur -1 est un bloc de 10
       datx$propchyp[datx$propchyp < 0] <- 0
       datx$chyp_kj_str <- datx$chyp_kj * datx$propchyp

       datx$predtemo <- mgcv::predict.gam(mod.temo, datx, type = "response")
       maxtemo<- max(datx$predtemo)
       if(maxtemo < 0.85) warning(paste0(ID,"maxtemo < 0.85"))
       datx$predtemo <-  datx$predtemo/maxtemo
       datx$proptemo <- c(NA, diff(datx$predtemo, lag = 1)) # lag de 1 parce que j'ai entr? en bloc de 10 donc valeur -1 est un bloc de 10
       datx$proptemo[datx$proptemo < 0] <- 0
       datx$temo_kj_str <- datx$temo_kj * datx$proptemo

       datx$predpsca <- mgcv::predict.gam(mod.psca, datx, type = "response")
       maxpsca<- max(datx$predpsca)
       if(maxpsca < 0.85) warning(paste0(ID,"maxpsca < 0.85"))
       datx$predpsca <-  datx$predpsca/maxpsca
       datx$proppsca <- c(NA, diff(datx$predpsca, lag = 1)) # lag de 1 parce que j'ai entr? en bloc de 10 donc valeur -1 est un bloc de 10
       datx$proppsca[datx$proppsca < 0] <- 0
       datx$psca_kj_str <- datx$psca_kj * datx$proppsca
       output <- rbind(output, datx[, c(ID, "Zsampled", "cfin_kj_str", "chyp_kj_str", "temo_kj_str", "psca_kj_str")])
     }

   }
    return(output)
 }
