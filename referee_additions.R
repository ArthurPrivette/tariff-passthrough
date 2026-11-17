# referee_additions.R
# Additional analyses reported in Appendix F of "Who Pays for a Tariff?"
#   1. Rebuild the monthly panel from the FGKK file (changes, fixed-effect ids)
#   2. Variation table (Table A6)
#   3. Event study and pre-trend tests (Exhibit 3)
#   4. First stage with a one-month lag (Table A8)
#   5. China-only pass-through (Table A7)
#
# Data: m_flow_hs10_fm_new.dta from the FGKK replication package,
# Harvard Dataverse, https://doi.org/10.7910/DVN/KSOVSE
# Set FILE to its location before running.

library(haven); library(data.table); library(fixest); library(ggplot2)

FILE  <- "extracted/rtp/data/analysis/m_flow_hs10_fm_new.dta"   # change to your path
CHINA <- 5700
M0 <- (2017 - 1960) * 12          # 2017m1 in Stata months
M1 <- (2019 - 1960) * 12 + 3      # 2019m4

# ---------------------------------------------------------------------------
# 1. Rebuild the panel
# ---------------------------------------------------------------------------
dt <- as.data.table(read_dta(FILE, col_select = c(
  "hs10", "hs8", "cty_code", "naics_str", "year", "mdate", "m_val",
  "lm_p", "lm_pduty", "lm_val", "lm_q1",
  "m_stattariff2", "m_apptariff", "m_status2", "m_effective_mdate2")))
dt <- dt[cty_code > 0 & mdate >= (M0 - 1) & mdate <= M1]   # one extra month so Jan 2017 has a change
gc()

dt[, id := .GRP, by = .(hs10, cty_code)]
dt[, naics4 := substr(naics_str, 1, 4)]
dt[, ct := .GRP, by = .(cty_code, mdate)]       # country x month
dt[, ht := .GRP, by = .(hs10, mdate)]           # product x month
dt[, cs := .GRP, by = .(cty_code, naics4)]      # country x industry

# m_stattariff2 is day-weighted within the month a tariff takes effect
dt[is.na(m_stattariff2), m_stattariff2 := 0]
dt[, lstattf := log1p(m_stattariff2)]
dt[, lapptf  := log1p(m_apptariff)]

setorder(dt, id, mdate)
d1 <- function(x, m) fifelse(m - shift(m) == 1, x - shift(x), NA_real_)
dt[, `:=`(d_lm_p     = d1(lm_p, mdate),
          d_lm_pduty = d1(lm_pduty, mdate),
          d_lm_val   = d1(lm_val, mdate),
          d_lm_q1    = d1(lm_q1, mdate),
          d_lstattf  = d1(lstattf, mdate),
          d_lapptf   = d1(lapptf, mdate)), by = id]
dt <- dt[mdate >= M0]
gc()

# Check: reduced form on the duty-paid price should be about 0.57-0.58
feols(d_lm_pduty ~ d_lstattf | ct + ht + cs, data = dt, cluster = ~ hs8 + cty_code)

# ---------------------------------------------------------------------------
# 2. Variation table
# ---------------------------------------------------------------------------
hit <- dt[d_lstattf != 0 & !is.na(d_lstattf)]
var_tab <- data.table(
  stat = c("Observations, 2017m1-2019m4",
           "HS10 products in sample",
           "Origin countries in sample",
           "Observations with a statutory tariff change",
           "HS10 products with at least one tariff change",
           "Countries with at least one tariff change",
           "Share of tariff-change observations from China"),
  value = c(nrow(dt), uniqueN(dt$hs10), uniqueN(dt$cty_code),
            nrow(hit), uniqueN(hit$hs10), uniqueN(hit$cty_code),
            round(mean(hit$cty_code == CHINA), 3)))
print(var_tab)
fwrite(var_tab, "table_variation.csv")

# ---------------------------------------------------------------------------
# 3. Event study: targeted Chinese varieties vs. other sellers of the same product
# ---------------------------------------------------------------------------
dt[, ever := max(m_status2 == 2, na.rm = TRUE), by = id]
dt[, eff  := suppressWarnings(max(m_effective_mdate2, na.rm = TRUE)), by = id]
dt[!is.finite(eff), eff := NA]
dt[, targeted := as.integer(cty_code == CHINA & ever == 1 & !is.na(eff))]
dt[, etime := fifelse(targeted == 1, mdate - eff, NA_real_)]
dt[, ebin  := fifelse(targeted == 1, pmin(pmax(etime, -12), 9), -1)]   # bin the tails
dt[targeted == 0, ebin := -1]                                            # controls at reference

es_run <- function(y) feols(as.formula(paste0(y, " ~ i(ebin, targeted, ref = -1) | id + ht + ct")),
                            data = dt, cluster = ~ hs8 + cty_code, lean = TRUE, mem.clean = TRUE)
es_p <- es_run("lm_p");   gc()     # pre-duty price
es_v <- es_run("lm_val"); gc()     # import value

grab <- function(m, lab) {
  b <- coef(m); s <- se(m)
  data.table(outcome = lab,
             k   = as.integer(sub(".*ebin::(-?[0-9]+):targeted", "\\1", names(b))),
             est = unname(b), se = unname(s))
}
es <- rbind(grab(es_p, "Pre-duty price"), grab(es_v, "Import value"),
            data.table(outcome = c("Pre-duty price", "Import value"), k = -1L, est = 0, se = 0))
setorder(es, outcome, k)
fwrite(es, "event_study_coefs.csv")

# Pre-trend: weighted linear trend through pre-period coefficients
trend_test <- function(lab) {
  pre <- es[outcome == lab & k < -1]
  fit <- lm(est ~ k, data = pre, weights = 1 / se^2)
  s <- summary(fit)$coefficients["k", ]
  c(slope = unname(s["Estimate"]), se = unname(s["Std. Error"]), p = unname(s["Pr(>|t|)"]))
}
print(rbind(price = trend_test("Pre-duty price"), value = trend_test("Import value")))
print(es[k != -1, .(mean_est = round(mean(est), 4)),
         by = .(outcome, period = fifelse(k < -1, "pre", "post"))])

es[, panel := factor(outcome, levels = c("Import value", "Pre-duty price"))]
g <- ggplot(es, aes(k, est)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = -0.5, colour = "grey40", linewidth = 0.3, linetype = "dashed") +
  geom_ribbon(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), fill = "#2a78d6", alpha = 0.15) +
  geom_line(colour = "#2a78d6", linewidth = 0.6) + geom_point(colour = "#2a78d6", size = 1.6) +
  facet_wrap(~ panel) +
  scale_x_continuous(breaks = seq(-12, 9, 3)) +
  labs(x = "Months relative to the tariff taking effect",
       y = "Log points vs. other sellers") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0))
ggsave("fig_event_study.png", g, width = 8, height = 3.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# 4. First stage with a one-month lag
# ---------------------------------------------------------------------------
setorder(dt, id, mdate)
dt[, d_lstattf_l1 := fifelse(mdate - shift(mdate) == 1, shift(d_lstattf), NA_real_), by = id]

fs0 <- feols(d_lapptf ~ d_lstattf                | ct + ht + cs, data = dt, cluster = ~ hs8 + cty_code)
fs1 <- feols(d_lapptf ~ d_lstattf + d_lstattf_l1 | ct + ht + cs, data = dt, cluster = ~ hs8 + cty_code)
etable(fs0, fs1, headers = c("Same month", "With one lag"))

b <- coef(fs1)[c("d_lstattf", "d_lstattf_l1")]
V <- vcov(fs1)[c("d_lstattf", "d_lstattf_l1"), c("d_lstattf", "d_lstattf_l1")]
cat("Cumulative first stage (b0 + b1):", round(sum(b), 4), "\n")
z <- (sum(b) - 1) / sqrt(sum(V))
cat("Test b0 + b1 = 1: z =", round(z, 3), " p =", round(2 * pnorm(-abs(z)), 4), "\n")

# ---------------------------------------------------------------------------
# 5. China-only pass-through
#    Product-by-month effects are not identified with one country, so the
#    comparison is across HS10 products within the same HS6 group and month.
# ---------------------------------------------------------------------------
cn <- dt[cty_code == CHINA]
cn[, hs6m := .GRP, by = .(floor(hs10 / 1e4), mdate)]
iv_cn_pre  <- feols(d_lm_p     ~ 1 | ct + hs6m + cs | d_lapptf ~ d_lstattf, data = cn, cluster = ~ hs8)
iv_cn_duty <- feols(d_lm_pduty ~ 1 | ct + hs6m + cs | d_lapptf ~ d_lstattf, data = cn, cluster = ~ hs8)
iv_cn_val  <- feols(d_lm_val   ~ 1 | ct + hs6m + cs | d_lapptf ~ d_lstattf, data = cn, cluster = ~ hs8)
etable(iv_cn_pre, iv_cn_duty, iv_cn_val,
       headers = c("Pre-duty price", "Duty-paid price", "Import value"))
