library("survival")
#datos.citom <- read.delim2("~/Documentos/articulos/LPL/figuras/originales/graficas_lpl/surv_lpl.csv", sep = ",") #Load data table
datos.citom <- read.delim2("~/Documentos/articulos/LPL/figuras/originales/graficas_lpl/ConHFresco_corr12-5-17multi_0dead_1alive.csv", sep = ",") #Load data table

datos.citom$X <- NULL
View(datos.citom) #Shows data frame
surv.lpl <- Surv(datos.citom$OS.months., datos.citom$status)
#head(surv.lpl)
lpl.km <- survfit(Surv(datos.citom$OS.months., datos.citom$status) ~ lpl.fc, data = datos.citom)
plot(lpl.km, lty = 1:2, col = c("darkblue", "darkred"))
legend("topright", legend = c("LPL+", "LPL-"), lty = 1:2, col = c("darkblue", "darkred"))
##
#Compute confidence intervals and plot them
lpl.km2 <- survfit(Surv(datos.citom$TFT.months., datos.citom$treat.status) ~ lpl.fc, data = datos.citom, conf.type = "log-log")
#lpl.km2 <- survfit(Surv(datos.citom$OS.months., datos.citom$status) ~ lpl.fc, data = datos.citom, conf.type = "log-log", type="fh")
summary(lpl.km2)
plot(lpl.km2, mark.time = FALSE, conf.int = TRUE, lty = 1, col = c("darkgreen", "darkred"))
legend("topright", legend = c("LPL-", "LPL+"), lty = 1, col = c("darkgreen", "darkred"))
##
#Test for difference (log-rank test)
survdiff(Surv(datos.citom$OS.months., datos.citom$status) ~ lpl.fc, data = datos.citom)
##
#Now a little more pro
library(survminer)
ggsurvplot(
  lpl.km2,                     # survfit object with calculated statistics.
  data = datos.citom,  # data used to fit survival curves. 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,100),        # present narrower X axis, but not affect
  xlab = "time (months)",
  # survival estimates.
  break.time.by = 5,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE, # show bars instead of names in text annotations
  # in legend of risk table
  surv.median.line = "hv"
)
##
#With the median surv
ggsurv <- ggsurvplot(
  lpl.km,                     # survfit object with calculated statistics.
  data = datos.citom,             # data used to fit survival curves.
  #risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,100),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time (months)",   # customize X axis label.
  break.time.by = 5,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  #risk.table.y.text.col = T,# colour risk table text annotations.
  #risk.table.height = 0.25, # the height of the risk table
  #risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  #ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv"  # add the median survival pointer.
  #legend.labs = 
   # c("Male", "Female")    # change legend labels.
)

# Apply custom color palettes and print
ggpar(ggsurv, palette = c("#E7B800", "#2E9FDF"))
##
## Using specific fonts for risk table and ncensor plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Font for Risk Table
ggsurv$table <- ggpar(ggsurv$table,
                      font.title = c(13, "bold.italic", "green"),
                      font.subtitle = c(15, "bold", "pink"),
                      font.caption = c(11, "plain", "darkgreen"),
                      font.x = c(8, "bold.italic", "orange"),
                      font.y = c(11, "bold.italic", "darkgreen"),
                      font.tickslab = c(9, "bold", "red")
)


# Font for ncensor plot
ggsurv$ncensor.plot <- ggpar(ggsurv$ncensor.plot,
                             font.title = c(13, "bold.italic", "green"),
                             font.subtitle = c(15, "bold", "pink"),
                             font.caption = c(11, "plain", "darkgreen"),
                             font.x = c(8, "bold.italic", "orange"),
                             font.y = c(11, "bold.italic", "darkgreen"),
                             font.tickslab = c(9, "bold", "red")
)

print(ggsurv)
##
## Changing Labels
# %%%%%%%%%%%%%%%%%%%%%%%%%%
# Labels for Survival Curves (plot)
ggsurv$plot <- ggsurv$plot + labs(
  title = "Survival curves",                     
  subtitle = "Based on Kaplan-Meier estimates"  
  #caption = "created with survminer"             
)

# Labels for Risk Table 
ggsurv$table <- ggsurv$table + labs(
  title = "Note the risk set sizes",          
  subtitle = "and remember about censoring.", 
  caption = "source code: website.com"        
)

# Labels for ncensor plot 
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + labs( 
  title = "Number of censorings", 
  subtitle = "over the time.",
  caption = "source code: website.com"
)

# Changing the font size, style and color
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Applying the same font style to all the components of ggsurv:
# survival curves, risk table and censor part

ggsurv <- ggpar(ggsurv,
                font.title = c(16, "bold", "darkblue"),         
                font.subtitle = c(15, "bold.italic", "gray"), 
                font.caption = c(14, "plain", "orange"),        
                font.x = c(14, "bold", "black"),          
                font.y = c(14, "bold", "black"),      
                font.tickslab = c(12, "plain", "darkgreen"),
                legend = "top"
)

# Apply custom color palettes and print
ggpar(ggsurv, palette = "lancet")#jco,jem,lancet
##
#
#Distribution of events’ times
#The function ggsurvevents() [in survminer] calculates and plots the distribution for events (both status = 0 and status = 1). It helps to notice when censoring is more common (@pbiecek, #116). This is an alternative to cumulative events and censored tables, described in the previous section.
#
ggsurvevents(surv.lpl)
#
##
# Data preparation and computing cox model
library(survival)
datos.citom <- read.delim2("~/Documentos/articulos/LPL/figuras/originales/graficas_lpl/ConHFresco_corr12-5-17multi_0dead_1alive.csv", sep = ",") #Load data table
datos.citom$lpl.fc <- factor(datos.citom$lpl.fc, levels = c(1,2), labels = c("LPL-", "LPL+"))
#survival_cens<- (as.numeric(as.character survival_cens))
#res.cox <- coxph(Surv(datos.citom$OS.months., datos.citom$status) ~ datos.citom$lpl.fc, data = datos.citom)
res.cox <- coxph(Surv(datos.citom$OS.months., datos.citom$status.residuals) ~ datos.citom$lpl.fc, data = datos.citom)
# Plot the baseline survival function
# with showing all individual predicted surv. curves
ggcoxadjustedcurves(res.cox, data = datos.citom, individual.curves = TRUE)
##
#Forest plot
ggforest(res.cox)
##
#Cut-off determination
library("survivalMPL")
#hist(datos.citom$FC, breaks = 20, col = "blue")
densp <- density(datos.citom$FC)
plot(densp)
#fit_mpl <- coxph_mpl(Surv(Multi_corr$OS.months., Multi_corr$dead1.alive0) ~ FC, data = Multi_corr)
fit_mpl <- coxph_mpl(Surv(datos.citom$TFT.months., datos.citom$status.residuals) ~ datos.citom$FC, data = datos.citom)
#plot(fit_mpl, ask=TRUE, which=1:3, upper.quantile=.95)
plot(residuals(fit_mpl), which=1, ask=FALSE)
#plot(fitted(fit_mpl), which=1, ask=T)
