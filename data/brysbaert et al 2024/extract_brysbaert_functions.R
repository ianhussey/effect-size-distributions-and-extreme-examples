# Load Brysbaert et al. (2024) corpus and compute participant-level sum scores
#
# Iterates over a catalog of (scale folder x data file) combinations under
#   data/brysbaert et al 2024/Data and analysis code/
# loads each, slices the item columns from the corresponding analysis script,
# and computes a per-participant sum score = mean(observed items) * n_items so
# missing data is handled by projecting the partial mean to the full-scale total.
#
# Outputs:
#   dat_brysbaert_sumscores  one row per (kept participant x dataset)
#     id            unique participant identifier (source_rownum)
#     scale         canonical instrument name (e.g. "BFI-2", "BISBAS", "BDI-II")
#     source        data file basename (lab/study identifier)
#     n_items       total items in the scale
#     n_observed    items the participant actually answered
#     sumscore      mean(items) * n_items
#     min_response  catalog response minimum
#     max_response  catalog response maximum
#     observed_min  empirical min for this dataset (kept participants only)
#     observed_max  empirical max for this dataset (kept participants only)
#
#   dat_brysbaert_long       one row per (kept participant x item); schema
#                            mirrors the Bainbridge dataset for compatibility
#     id, subscale, item, response, scale, min_response, max_response
#     (subscale carries the source/dataset identifier since Brysbaert has no
#     genuine subscales after the per-instrument collapse)
#
# At write time both are filtered to scales with prop_excluded < 0.10 (i.e.,
# fewer than 10% of participants had at least one out-of-range item response).
# This drops datasets where the catalog response_min/max is wrong; the in-
# memory variables stay unfiltered for diagnostics.
#
# Response ranges are populated from standard psychometric conventions where
# known. Where the catalog max disagrees with the observed data max, the
# script prints a warning so the user can refine. Instrument names come from
# INSTRUMENT_RULES (see below); folder->instrument mapping is one-to-one
# except for multi-instrument folders (Personality, Personality disorders,
# Impulsivity UPPS and BIS, Emotional intelligence, Interoceptive accuracy,
# OCD, Compulsivity Muela), which dispatch on filename patterns.
#
# To use:
#   source("extract_brysbaert_functions.R")
#   # produces `dat_brysbaert_sumscores` in the calling environment
# or run directly to also save a CSV.

suppressPackageStartupMessages({
  library(tidyverse)
  library(haven)
  library(readxl)
})

# --------------------------------------------------------------------- catalog

# Tab-separated. item_cols is an R expression to be eval'd in the data frame's
# scope; can be a range like "257:276" or a vector like "c(3:6,8:11)".
# Set valid = FALSE to skip an entry without deleting it.
brysbaert_catalog <- read.delim(
  text = "scale_folder\tdata_file\titem_cols\tn_items\tresponse_min\tresponse_max\tvalid
BISBAS\tBISBAS_Binter.xlsx\t10:29\t20\t1\t4\tTRUE
BISBAS\tBISBAS_Jonker.sav\t(23:46)[-c(1,6,11,17)]\t20\t1\t4\tTRUE
BISBAS\tBISBAS_Khaliq.sav\t(298:321)[-c(1,6,11,17)]\t20\t1\t4\tTRUE
BISBAS\tBISBAS_Molenda.sav\t(34:57)[-c(1,6,11,17)]\t20\t1\t4\tTRUE
BISBAS\tBISBAS_Rutten.csv\t257:276\t20\t1\t4\tTRUE
BISBAS\tBISBASdata_Dierickx.sav\tc(3:6,8:11,13:17,19:25)\t20\t1\t4\tTRUE
Beliefs\tBeliefs_Barnby.csv\tseq(10, ncol(.df), 3)\tNA\t1\t7\tFALSE
Compulsivity Muela\tcompulsivity_gambling_Muela.xlsx\t10:99\t90\t1\t5\tTRUE
Compulsivity Muela\tcompulsivity_videogames_Muela.xlsx\t10:99\t90\t1\t5\tTRUE
CovidStress\tFinal_COVIDiSTRESS_Vol2_cleaned.csv\t111:117\t7\t1\t6\tTRUE
DASS-21\tdass21_Thiyagarajan.xlsx\t9:29\t21\t0\t3\tTRUE
DOSPERT\tDOSPERT_Frey.csv\tc(3:8,21:26,39:44,57:62,75:80)\t30\t1\t7\tTRUE
Depression\tBDI-II_Garcia_Batista.sav\t3:23\t21\t0\t3\tTRUE
Depression\tBDI-II_Rabeya.xlsx\t17:37\t21\t0\t3\tTRUE
Depression\tBDI-II_Smith.sav\t49:68\t20\t0\t3\tTRUE
Depression\tdepression_Sandoval.xlsx\t27:47\t21\t0\t3\tTRUE
Emotion regulation DERS\tDERS_Valencia.sav\t8:43\t36\t1\t5\tTRUE
Emotional intelligence\tAES_Anglim.csv\t254:286\t33\t1\t5\tTRUE
Emotional intelligence\tAES_Demetrovics.sav\t38:70\t33\t1\t5\tTRUE
Emotional intelligence\tAES_Demetrovics2.sav\t44:76\t33\t1\t5\tTRUE
Emotional intelligence\tAES_Rua.xlsx\t5:37\t33\t1\t5\tTRUE
Emotional intelligence\tBainbridge_s2_IRI.csv\t224:239\t16\t1\t5\tTRUE
Emotional intelligence\tEI_Brienza.sav\t283:298\t16\t1\t5\tTRUE
Emotional intelligence\tEI_Robinson.xlsx\t2:17\t16\t1\t5\tTRUE
Emotional intelligence\tIRI_Stosic.xlsx\t2:29\t28\t0\t4\tTRUE
Emotional intelligence\tMartingano_IRI.sav\t187:214\t28\t0\t4\tTRUE
Emotional intelligence\tTEIQue_Perazzo.sav\t4:33\t30\t1\t7\tTRUE
Emotional intelligence\tTEIQue_Perez-Diaz.csv\t1:30\t30\t1\t7\tTRUE
FFMQ\tChile.sav\t22:60\t39\t1\t5\tTRUE
FFMQ\tGermany.sav\t16:54\t39\t1\t5\tTRUE
FFMQ\tHong Kong.sav\t17:55\t39\t1\t5\tTRUE
FFMQ\tIndia.sav\t8:46\t39\t1\t5\tTRUE
FFMQ\tNorway.sav\t5:43\t39\t1\t5\tTRUE
FFMQ\tPoland.xlsx\t1:39\t39\t1\t5\tTRUE
FFMQ\tPortugal.sav\t2:40\t39\t1\t5\tTRUE
FFMQ\tRomania.xlsx\t7:45\t39\t1\t5\tTRUE
FFMQ\tSpain.sav\t3:41\t39\t1\t5\tTRUE
FFMQ\taustralia.sav\t3:41\t39\t1\t5\tTRUE
FFMQ\taustria.sav\t5:43\t39\t1\t5\tFALSE
FFMQ\tcroatia.sav\t4:42\t39\t1\t5\tTRUE
FFMQ\tsweden.sav\t6:44\t39\t1\t5\tTRUE
FFMQ\tus1.xlsx\t3:41\t39\t1\t5\tTRUE
Face recognition\tDal_Lago_2023_PI20_raw_data.xlsx\t4:23\t20\t1\t5\tTRUE
Face recognition\tDanish PI20 validation study.xlsx\t15:34\t20\t1\t5\tTRUE
Face recognition\tEstudillo_2021_PI20_JCP.csv\t2:21\t20\t1\t5\tFALSE
Face recognition\tEstudillo_2021_PI20_PEERJ.csv\t2:21\t20\t1\t5\tFALSE
Face recognition\tKramer_2023_PI20.xlsx\t1:20\t20\t1\t5\tTRUE
Face recognition\tPI20_questionnaire_Jlowes.xlsx\t4:23\t20\t1\t5\tTRUE
Handedness\tDragan_handedness.dat\t5:14\t10\t1\t5\tTRUE
Handedness\tMcManus_handedness.sav\tc(2,14,3,8,12,17,4,5,7,27)\t10\t1\t5\tTRUE
Hope\tBoutilier_RRB_Hope_Optimism_Data.sav\t89:100\t12\t1\t8\tTRUE
Impulsivity UPPS and BIS\tBIS_Demetrovics1.sav\t3:32\t30\t1\t4\tTRUE
Impulsivity UPPS and BIS\tBIS_Demetrovics2.sav\t10:30\t21\t1\t4\tTRUE
Impulsivity UPPS and BIS\tBillieux.xlsx\t7:26\t20\t1\t4\tTRUE
Impulsivity UPPS and BIS\tBothe.sav\t1:20\t20\t1\t4\tTRUE
Impulsivity UPPS and BIS\tEben2020.xlsx\t2:21\t20\t1\t4\tTRUE
Impulsivity UPPS and BIS\tFinley.sav\t129:162\t30\t1\t4\tTRUE
Impulsivity UPPS and BIS\tFlayelle.xlsx\t27:46\t20\t1\t4\tTRUE
Impulsivity UPPS and BIS\tHaines.csv\t26:55\t30\t1\t4\tTRUE
Impulsivity UPPS and BIS\tLittrell.csv\t44:73\t30\t1\t4\tTRUE
Impulsivity UPPS and BIS\tMeule.sav\t5:19\t15\t1\t4\tTRUE
Impulsivity UPPS and BIS\tRacine_EMA_Data.sav\t74:132\t59\t1\t4\tTRUE
Impulsivity UPPS and BIS\tSteiner.xlsx\t2:46\t45\t1\t4\tFALSE
Impulsivity UPPS and BIS\tvalidation_sample_Wuellhorst.xlsx\t2:60\t59\t1\t4\tTRUE
Interoceptive accuracy\tIAS_Brand1.sav\t13:33\t21\t1\t5\tTRUE
Interoceptive accuracy\tIAS_Brand2.sav\t11:31\t21\t1\t5\tTRUE
Interoceptive accuracy\tIAS_Brand3.sav\t4:24\t21\t1\t5\tTRUE
Interoceptive accuracy\tIAS_Campos.sav\t76:96\t21\t1\t5\tTRUE
Interoceptive accuracy\tIAS_Murphy1.sav\t3:23\t21\t1\t5\tTRUE
Interoceptive accuracy\tIAS_Todd.sav\t80:100\t21\t1\t5\tTRUE
Interoceptive accuracy\tMAIA_Desmedt.csv\t27:63\t37\t0\t5\tTRUE
Interoceptive accuracy\tMAIA_Ferentzi.xlsx\t9:40\t32\t0\t5\tTRUE
Interoceptive accuracy\tMAIA_Randelovic.sav\t4:40\t37\t0\t5\tTRUE
Interoceptive accuracy\tMAIA_Reis.csv\t1:32\t32\t0\t5\tTRUE
Interoceptive accuracy\tMAIA_Rogowska.xlsx\t6:42\t37\t0\t5\tTRUE
Loneliness\tLardone.xlsx\t109:128\t20\t1\t4\tTRUE
Loneliness\tMcDanal.csv\t30:49\t20\t1\t4\tTRUE
Loneliness\tPanayiotou.xlsx\t2:21\t20\t1\t4\tTRUE
Loneliness\tWeinstein.csv\t50:69\t20\t1\t4\tTRUE
Moral foundations\tMoralF_Harper.csv\t25:65\t41\t0\t5\tTRUE
Moral foundations\tMoralF_Zakharin.csv\tc(3:7,9:23,25:34)\t30\t0\t5\tTRUE
OCD\tOCI-R dataset Ignatova.xlsx\t2:19\t18\t0\t4\tTRUE
OCD\tObsession_Ozcanli.sav\t2:51\t50\t0\t4\tTRUE
Optimism\tOptimism_Coelho_Bra.sav\t2:10\t9\t0\t4\tTRUE
Optimism\tOptimism_Coelho_UK.sav\t11:19\t9\t0\t4\tTRUE
PANAS\tPANAS_Diaz-Garcia.sav\t12:31\t20\t1\t5\tTRUE
Perfectionism\tMacKinnon_Perfectionism1.csv\t259:304\t46\t1\t7\tTRUE
Perfectionism\tPerfectionism_Kacar-Basaran.sav\t2:46\t45\t1\t7\tFALSE
Perfectionism\tPerfectionism_Linnett.sav\t7:36\t30\t1\t7\tTRUE
Perfectionism\tWorkye_Perfectionism.csv\t37:74\t38\t1\t7\tTRUE
Personality\tBFAS_Zajenkowski.sav\t49:148\t100\t1\t5\tTRUE
Personality\tBFI-2_Gallardo-Pujol.csv\t5:64\t60\t1\t5\tTRUE
Personality\tBFI2-Andrejevic.csv\t5:64\t60\t1\t5\tTRUE
Personality\tBFI2_Vermeiren.xlsx\t2:61\t60\t1\t5\tTRUE
Personality\tBig_Five_Mezquita.sav\t5:54\t50\t1\t5\tTRUE
Personality\tHEXACO_Gallagher.csv\t23:82\t60\t1\t5\tTRUE
Personality\tHEXACO_Henry_time1.xlsx\t1:96\t96\t1\t5\tTRUE
Personality\tIPIP_Silvia.csv\t107:156\t50\t1\t5\tTRUE
Personality\tTIPI_Coelho_Bra.sav\t11:20\t10\t1\t7\tTRUE
Personality\tTIPI_Coelho_UK.sav\t1:10\t10\t1\t7\tTRUE
Personality\tdata_RIASEC_openpsychometrics.csv\t1:48\t48\t1\t5\tFALSE
Personality disorders\tDark_triad_Nielsen.csv\t3:29\t27\t1\t5\tTRUE
Personality disorders\tICD11_Carnovale.csv\t6:65\t60\t0\t3\tTRUE
Personality disorders\tMachiavellism_Grabovac.csv\t1:52\t52\t1\t5\tTRUE
Personality disorders\tNatoli_Pers_Functioning.csv\t6:17\t12\t1\t4\tTRUE
Personality disorders\tPsychopathy_Knack.sav\t26:89\t64\t1\t5\tTRUE
Post Traumatic Stress Disorder\tAshbaugh_PCL-5.FR.sav\t4:23\t20\t0\t4\tTRUE
Post Traumatic Stress Disorder\tAshbaugh_PCL-5.sav\t4:23\t20\t0\t4\tTRUE
Post Traumatic Stress Disorder\tOrovou_PCL-5.sav\t26:45\t20\t0\t4\tTRUE
QEWB\tQEWB_Ishii.xlsx\t5:25\t21\t1\t5\tTRUE
Quality of Life WHO\tWHOQOL_McConachie.xlsx\t3:26\t24\t1\t5\tTRUE
Reinforcement sensitivity\tNuel_RST-PQ.txt\t2:55\t54\t1\t4\tTRUE
Reinforcement sensitivity\tNuel_RST-PQ2.txt\t2:55\t54\t1\t4\tTRUE
Reinforcement sensitivity\tRST_Sozer.sav\t3:67\t65\t1\t4\tFALSE
Reinforcement sensitivity\tSatchell_RST-PQ.csv\t5:77\t73\t1\t4\tTRUE
STAI-T\tSTAI-T_Sundelin.xlsx\t3:22\t20\t1\t4\tTRUE
STAI-T\tSTAI_Joshi.xlsx\t22:41\t20\t1\t4\tTRUE
Self control\tBSC1_Paap.sav\t63:75\t13\t1\t5\tTRUE
Self control\tBSC2_Paap.sav\t78:90\t13\t1\t5\tTRUE
Self control\tBSC_Sjastad1.xlsx\t14:26\t13\t1\t5\tTRUE
Self control\tBSC_Sjastad2.xlsx\t14:26\t13\t1\t5\tTRUE
Sense of coherence\tCoherence_Zajenkowska.sav\t86:114\t29\t1\t7\tTRUE
Sense of coherence\tSense_of_coherence_Lelek.xlsx\t2:30\t29\t1\t7\tTRUE
Sexual desire\tSDI-2_Jones.txt\t5:18\t14\t0\t8\tTRUE
Social dominance\tSDO_Fleeson1.sav\t37:52\t16\t1\t7\tTRUE
Social dominance\tSDO_Fleeson2.sav\t85:100\t16\t1\t7\tTRUE
Social dominance\tSDO_Perkins.csv\t8:23\t16\t1\t7\tTRUE
Social dominance\tSDO_Roets1.sav\t15:30\t16\t1\t7\tTRUE
Social dominance\tSDO_Roets2.sav\t16:31\t16\t1\t7\tTRUE
Social dominance\tSDO_Simon1.xlsx\t585:600\t16\t1\t7\tTRUE
Social dominance\tSDO_Simon3.xlsx\t579:594\t16\t1\t7\tTRUE
Symptom checklist 90R\tSCL90R_Chen.csv\t55:144\t90\t0\t4\tTRUE
Time perspective\tTime_perspective_Bodecka.sav\t41:60\t20\t1\t5\tTRUE
Time perspective\tZTPI_Akirmak.csv\t19:74\t56\t1\t5\tTRUE
Time perspective\tZTPI_Ceccato.xlsx\t8:25\t18\t1\t5\tTRUE",
  sep = "\t", stringsAsFactors = FALSE, na.strings = c("", "NA"),
  encoding = "UTF-8"
)

# ----------------------------------------------------------- instrument map --

# Map (scale_folder, data_file) -> canonical instrument name. Single-instrument
# folders use the catch-all `.*` rule. Multi-instrument folders try each rule
# in order; first regex match (case-insensitive) wins.
INSTRUMENT_RULES <- list(
  "BISBAS"                          = list(c(".*", "BISBAS")),
  "CovidStress"                     = list(c(".*", "COVIDiSTRESS Covid-19 stress")),
  "DASS-21"                         = list(c(".*", "DASS-21")),
  "DOSPERT"                         = list(c(".*", "DOSPERT risk-taking")),
  "Depression"                      = list(c(".*", "BDI-II depression")),
  "Emotion regulation DERS"         = list(c(".*", "DERS emotion regulation")),
  "FFMQ"                            = list(c(".*", "FFMQ mindfulness")),
  "Face recognition"                = list(c(".*", "PI20 prosopagnosia")),
  "Handedness"                      = list(c(".*", "EHI handedness")),
  "Hope"                            = list(c(".*", "Hope (Snyder)")),
  "Loneliness"                      = list(c(".*", "UCLA Loneliness")),
  "Moral foundations"               = list(c(".*", "MFQ moral foundations")),
  "Optimism"                        = list(c(".*", "LOT-R optimism")),
  "PANAS"                           = list(c(".*", "PANAS affect")),
  "Perfectionism"                   = list(c(".*", "Perfectionism")),
  "Post Traumatic Stress Disorder"  = list(c(".*", "PCL-5 PTSD")),
  "QEWB"                            = list(c(".*", "QEWB eudaimonic wellbeing")),
  "Quality of Life WHO"             = list(c(".*", "WHOQOL quality of life")),
  "Reinforcement sensitivity"       = list(c(".*", "RST-PQ reinforcement sensitivity")),
  "STAI-T"                          = list(c(".*", "STAI-T trait anxiety")),
  "Self control"                    = list(c(".*", "BSC self-control")),
  "Sense of coherence"              = list(c(".*", "SOC-29 sense of coherence")),
  "Sexual desire"                   = list(c(".*", "SDI-2 sexual desire")),
  "Social dominance"                = list(c(".*", "SDO social dominance")),
  "Symptom checklist 90R"           = list(c(".*", "SCL-90R symptom checklist")),
  "Time perspective"                = list(c(".*", "ZTPI time perspective")),

  "Compulsivity Muela" = list(
    c("gambling",   "Compulsivity (gambling)"),
    c("videogames", "Compulsivity (videogames)")
  ),
  "Emotional intelligence" = list(
    c("AES",     "AES emotional intelligence"),
    c("TEIQue",  "TEIQue emotional intelligence"),
    c("IRI",     "IRI empathy"),
    c("EI_",     "EI (short) emotional intelligence")
  ),
  "Impulsivity UPPS and BIS" = list(
    c("Meule",                                  "BIS-15 impulsivity"),
    c("BIS_Demetrovics|Finley|Haines|Littrell", "BIS-11 impulsivity"),
    c("Racine|Wuellhorst",                      "UPPS-P impulsivity"),
    c("Billieux|Bothe|Eben|Flayelle",           "UPPS impulsivity")
  ),
  "Interoceptive accuracy" = list(
    c("IAS",  "IAS interoception"),
    c("MAIA", "MAIA interoception")
  ),
  "OCD" = list(
    c("OCI-R",     "OCI-R OCD"),
    c("Obsession", "Obsession (50-item)")
  ),
  "Personality" = list(
    c("BFAS",     "BFAS-100"),
    c("BFI-?2",   "BFI-2"),
    c("Big_Five", "Big Five (50-item)"),
    c("HEXACO",   "HEXACO"),
    c("IPIP",     "IPIP-50"),
    c("TIPI",     "TIPI"),
    c("RIASEC",   "RIASEC")
  ),
  "Personality disorders" = list(
    c("Dark_triad",              "Dark Triad"),
    c("ICD11",                   "PiCD personality disorders"),
    c("Machiavellism",           "Mach IV Machiavellianism"),
    c("Natoli_Pers_Functioning", "LPFS-BF personality functioning"),
    c("Psychopathy",             "SRP-III psychopathy")
  )
)

# Subscale mappings: data_file -> named list of subscale_name -> 1-indexed
# positions WITHIN THE SLICED ITEM COLUMNS. Datasets present here get split
# into per-subscale sumscores; datasets absent stay as a single scale.
#
# Positions are derived from each lab's analysis script (the canonical-name
# rename combined with the lavaan model), or from canonical published item
# orderings where the lab uses the standard form.
SUBSCALE_MAPPINGS <- list(
  # BISBAS — standard rename order: B2,B3,B4,B5,B7,B8,B9,B10,B12,B13,B14,
  # B15,B16,B18,B19,B20,B21,B22,B23,B24. BIS items: B2,B8,B13,B16,B19,B22,B24.
  "BISBAS_Rutten.csv"       = list(
    "BIS behavioral inhibition" = c(1, 6, 10, 13, 15, 18, 20),
    "BAS behavioral activation" = c(2, 3, 4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 19)
  ),
  "BISBASdata_Dierickx.sav" = list(
    "BIS behavioral inhibition" = c(1, 6, 10, 13, 15, 18, 20),
    "BAS behavioral activation" = c(2, 3, 4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 19)
  ),
  "BISBAS_Jonker.sav" = list(
    "BIS behavioral inhibition" = c(1, 6, 10, 13, 15, 18, 20),
    "BAS behavioral activation" = c(2, 3, 4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 19)
  ),
  "BISBAS_Khaliq.sav" = list(
    "BIS behavioral inhibition" = c(1, 6, 10, 13, 15, 18, 20),
    "BAS behavioral activation" = c(2, 3, 4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 19)
  ),
  "BISBAS_Molenda.sav" = list(
    "BIS behavioral inhibition" = c(1, 6, 10, 13, 15, 18, 20),
    "BAS behavioral activation" = c(2, 3, 4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 19)
  ),
  # Binter uses a different rename order (BAS first, then BIS):
  # B3,B9,B12,B21,B5,B10,B15,B20,B4,B7,B14,B18,B23,B2,B8,B13,B16,B19,B22,B24
  "BISBAS_Binter.xlsx" = list(
    "BIS behavioral inhibition" = c(14, 15, 16, 17, 18, 19, 20),
    "BAS behavioral activation" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
  ),

  # SD3 (Short Dark Triad): items 1-9 Narcissism, 10-18 Mach, 19-27 Psychopathy.
  # Prefixed with "SD3" to disambiguate from Mach IV and Psychopathy_Knack.
  "Dark_triad_Nielsen.csv" = list(
    "SD3 narcissism"       = 1:9,
    "SD3 Machiavellianism" = 10:18,
    "SD3 psychopathy"      = 19:27
  ),

  # DASS-21: standard 21-item ordering per Lovibond & Lovibond (1995)
  "dass21_Thiyagarajan.xlsx" = list(
    "dass anxiety"    = c(2, 4, 7, 9, 15, 19, 20),
    "dass depression" = c(3, 5, 10, 13, 17, 18, 21),
    "dass stress"     = c(1, 6, 8, 11, 12, 14, 16)
  ),

  # BFI-2 (Soto & John 2017): canonical 5-item cycle E,A,C,N,O.
  # NB the analysis script's lavaan model used a different ordering that
  # contradicts the published BFI-2; using the canonical mapping here.
  "BFI-2_Gallardo-Pujol.csv" = list(
    "BFI-2 extraversion"      = c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56),
    "BFI-2 agreeableness"     = c(2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57),
    "BFI-2 conscientiousness" = c(3, 8, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58),
    "BFI-2 neuroticism"       = c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59),
    "BFI-2 openness"          = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
  ),
  "BFI2-Andrejevic.csv" = list(
    "BFI-2 extraversion"      = c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56),
    "BFI-2 agreeableness"     = c(2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57),
    "BFI-2 conscientiousness" = c(3, 8, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58),
    "BFI-2 neuroticism"       = c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59),
    "BFI-2 openness"          = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
  ),
  "BFI2_Vermeiren.xlsx" = list(
    "BFI-2 extraversion"      = c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56),
    "BFI-2 agreeableness"     = c(2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57),
    "BFI-2 conscientiousness" = c(3, 8, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58),
    "BFI-2 neuroticism"       = c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59),
    "BFI-2 openness"          = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
  ),

  # IPIP-50 (Goldberg's IPIP Big Five): blocks of 10 items per factor.
  "IPIP_Silvia.csv" = list(
    "IPIP-50 extraversion"      = 1:10,
    "IPIP-50 agreeableness"     = 11:20,
    "IPIP-50 conscientiousness" = 21:30,
    "IPIP-50 neuroticism"       = 31:40,
    "IPIP-50 openness"          = 41:50
  ),
  # Big_Five_Mezquita uses the same IPIP-50 layout per the analysis script.
  "Big_Five_Mezquita.sav" = list(
    "Big Five (50-item) extraversion"      = c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46),
    "Big Five (50-item) agreeableness"     = c(2, 7, 12, 17, 22, 27, 32, 37, 42, 47),
    "Big Five (50-item) conscientiousness" = c(3, 8, 13, 18, 23, 28, 33, 38, 43, 48),
    "Big Five (50-item) neuroticism"       = c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49),
    "Big Five (50-item) openness"          = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  ),

  # HEXACO-60 (Ashton & Lee 2009): canonical 6-item cycle O,C,A,X,E,H.
  "HEXACO_Gallagher.csv" = list(
    "HEXACO openness"          = c(1, 7, 13, 19, 25, 31, 37, 43, 49, 55),
    "HEXACO conscientiousness" = c(2, 8, 14, 20, 26, 32, 38, 44, 50, 56),
    "HEXACO agreeableness"     = c(3, 9, 15, 21, 27, 33, 39, 45, 51, 57),
    "HEXACO extraversion"      = c(4, 10, 16, 22, 28, 34, 40, 46, 52, 58),
    "HEXACO emotionality"      = c(5, 11, 17, 23, 29, 35, 41, 47, 53, 59),
    "HEXACO honesty-humility"  = c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60)
  ),
  # HEXACO-96 (Henry): blocked layout (16 items per factor, A,C,E,H,O,X).
  "HEXACO_Henry_time1.xlsx" = list(
    "HEXACO agreeableness"     = 1:16,
    "HEXACO conscientiousness" = 17:32,
    "HEXACO emotionality"      = 33:48,
    "HEXACO honesty-humility"  = 49:64,
    "HEXACO openness"          = 65:80,
    "HEXACO extraversion"      = 81:96
  ),

  # TIPI (Gosling et al. 2003): pairs.
  "TIPI_Coelho_Bra.sav" = list(
    "TIPI extraversion"      = c(1, 6),
    "TIPI agreeableness"     = c(2, 7),
    "TIPI conscientiousness" = c(3, 8),
    "TIPI neuroticism"       = c(4, 9),
    "TIPI openness"          = c(5, 10)
  ),
  "TIPI_Coelho_UK.sav" = list(
    "TIPI extraversion"      = c(1, 6),
    "TIPI agreeableness"     = c(2, 7),
    "TIPI conscientiousness" = c(3, 8),
    "TIPI neuroticism"       = c(4, 9),
    "TIPI openness"          = c(5, 10)
  ),

  # BFAS-100 / IPIP-100: 100 items, interleaved 5-item cycle.
  "BFAS_Zajenkowski.sav" = list(
    "BFAS-100 extraversion"      = seq(1, 100, by = 5),
    "BFAS-100 agreeableness"     = seq(2, 100, by = 5),
    "BFAS-100 conscientiousness" = seq(3, 100, by = 5),
    "BFAS-100 neuroticism"       = seq(4, 100, by = 5),
    "BFAS-100 openness"          = seq(5, 100, by = 5)
  )
)

resolve_instrument <- function(scale_folder, data_file) {
  base <- tools::file_path_sans_ext(data_file)
  rules <- INSTRUMENT_RULES[[scale_folder]]
  if (is.null(rules)) return(scale_folder)
  for (r in rules) {
    if (grepl(r[1], base, ignore.case = TRUE, perl = TRUE)) return(r[2])
  }
  scale_folder
}

# ---------------------------------------------------------------- file readers

# Try to read a data file based on its extension. Returns a data.frame on
# success, NULL on failure (with a warning).
read_dataset <- function(path) {
  ext <- tolower(tools::file_ext(path))
  # Try multiple readers; pick the one that yields the most columns
  # (delimited files collapse to 1 column when the wrong separator is used).
  try_readers <- function(readers) {
    best <- NULL
    for (rd in readers) {
      x <- tryCatch(rd(), error = function(e) NULL,
                    warning = function(w) tryCatch(rd(), error = function(e) NULL))
      if (!is.null(x) && (is.null(best) || ncol(x) > ncol(best))) best <- x
    }
    best
  }
  csv_readers <- list(
    function() utils::read.csv(path, stringsAsFactors = FALSE),
    function() utils::read.csv2(path, stringsAsFactors = FALSE),
    function() utils::read.delim(path, stringsAsFactors = FALSE)
  )
  txt_readers <- list(
    function() utils::read.delim(path, stringsAsFactors = FALSE),
    function() utils::read.table(path, header = TRUE, sep = "",
                                 stringsAsFactors = FALSE, fill = TRUE),
    function() utils::read.csv(path, stringsAsFactors = FALSE)
  )
  d <- tryCatch({
    switch(ext,
      csv  = try_readers(csv_readers),
      sav  = haven::read_sav(path),
      xlsx = readxl::read_excel(path),
      dat  = try_readers(txt_readers),
      txt  = try_readers(txt_readers),
      stop("unknown extension: ", ext)
    )
  }, error = function(e) {
    warning(sprintf("Failed to read %s: %s", basename(path), conditionMessage(e)))
    NULL
  })
  if (is.null(d)) return(NULL)
  as.data.frame(d, stringsAsFactors = FALSE)
}

# Coerce a slice of columns to a numeric matrix, stripping haven labels.
items_to_numeric <- function(items) {
  if (requireNamespace("haven", quietly = TRUE)) {
    items <- haven::zap_labels(items)
  }
  as.data.frame(lapply(items, function(x) suppressWarnings(as.numeric(x))))
}

# ---------------------------------------------------------------- ingestion --

# Resolve base_dir whether the script is run/sourced from:
#   - this folder (data/brysbaert et al 2024/), e.g. via extract_brysbaert.Rmd
#   - the repo root
#   - somewhere a few levels deep
find_base_dir <- function() {
  candidates <- c(
    "Data and analysis code",
    file.path("data", "brysbaert et al 2024", "Data and analysis code"),
    file.path("..", "data", "brysbaert et al 2024", "Data and analysis code"),
    file.path("..", "..", "data", "brysbaert et al 2024", "Data and analysis code"),
    file.path("..", "..", "..", "data", "brysbaert et al 2024", "Data and analysis code")
  )
  for (p in candidates) if (dir.exists(p)) return(p)
  stop("Could not locate 'Data and analysis code/' relative to working directory: ",
       getwd())
}
base_dir <- find_base_dir()

ingest_one <- function(row) {
  path <- file.path(base_dir, row$scale_folder, row$data_file)
  if (!file.exists(path)) {
    return(list(ok = FALSE, reason = "file not found", df = NULL))
  }
  d <- read_dataset(path)
  if (is.null(d)) return(list(ok = FALSE, reason = "read failed", df = NULL))

  # Eval the column expression in a sandbox where `.df` references the data
  cols <- tryCatch(
    eval(parse(text = row$item_cols), envir = list(.df = d)),
    error = function(e) NULL
  )
  if (is.null(cols)) {
    return(list(ok = FALSE, reason = "item_cols parse failed", df = NULL))
  }
  if (is.numeric(cols)) cols <- as.integer(cols)
  if (is.numeric(cols) && (any(cols < 1) || any(cols > ncol(d)))) {
    return(list(ok = FALSE, reason = "item_cols out of range", df = NULL))
  }
  items <- tryCatch(d[, cols, drop = FALSE], error = function(e) NULL)
  if (is.null(items)) {
    return(list(ok = FALSE, reason = "slicing failed", df = NULL))
  }

  items_num <- items_to_numeric(items)

  # NA out values that are clearly missing-data sentinels: anything more
  # than 10 below response_min or 10 above response_max can't be a Likert
  # response. This catches conventions like -99, -999, 999, 503506, etc.,
  # without disturbing data when the catalog response_max is slightly off.
  if (!is.na(row$response_min) && !is.na(row$response_max)) {
    floor_  <- as.numeric(row$response_min) - 10
    ceil_   <- as.numeric(row$response_max) + 10
    sentinel_mask <- !is.na(items_num) &
                     ((items_num < floor_) | (items_num > ceil_))
    if (any(as.matrix(sentinel_mask))) {
      items_num[as.matrix(sentinel_mask)] <- NA_real_
    }
  }

  n_items   <- ncol(items_num)
  n_obs     <- rowSums(!is.na(items_num))

  # Per-participant data quality gate: any non-NA item that falls outside
  # the catalog [response_min, response_max] flags the participant for
  # exclusion. Sentinels (>10 outside the range) have already been NA'd
  # above so they don't trigger this; only plausible-but-out-of-range
  # values (e.g. a 6 on a 1-5 catalog scale, often signalling a different
  # scale variant) cause exclusion.
  participant_n_total <- nrow(items_num)
  if (!is.na(row$response_min) && !is.na(row$response_max)) {
    oor <- !is.na(items_num) &
           ((items_num < as.numeric(row$response_min)) |
            (items_num > as.numeric(row$response_max)))
    n_oor_per_row <- rowSums(as.matrix(oor))
    keep <- n_oor_per_row == 0L
  } else {
    keep <- rep(TRUE, participant_n_total)
  }
  participant_n_kept     <- sum(keep)
  participant_n_excluded <- participant_n_total - participant_n_kept

  means    <- rowMeans(items_num, na.rm = TRUE)
  means[!is.finite(means)] <- NA_real_
  sumscore <- means * n_items

  # Raw observed min/max across ALL participants (after sentinel NA, before
  # per-participant range exclusion). Useful for the range_check diagnostic
  # since the kept-participants min/max trivially lies inside the catalog.
  raw_observed_min <- suppressWarnings(min(items_num, na.rm = TRUE))
  raw_observed_max <- suppressWarnings(max(items_num, na.rm = TRUE))
  if (!is.finite(raw_observed_min)) raw_observed_min <- NA_real_
  if (!is.finite(raw_observed_max)) raw_observed_max <- NA_real_

  # Observed min/max among kept participants (always within catalog range).
  observed_min <- suppressWarnings(min(items_num[keep, , drop = FALSE], na.rm = TRUE))
  observed_max <- suppressWarnings(max(items_num[keep, , drop = FALSE], na.rm = TRUE))
  if (!is.finite(observed_min)) observed_min <- NA_real_
  if (!is.finite(observed_max)) observed_max <- NA_real_

  source_id <- tools::file_path_sans_ext(row$data_file)
  instrument <- resolve_instrument(row$scale_folder, row$data_file)
  participant_ids <- paste(source_id, seq_len(participant_n_total), sep = "_")

  # Produce sumscore + long rows for a (subset of) item columns under a given
  # scale label. Reuses the dataset-level `keep` mask so the per-participant
  # exclusion is applied consistently across all subscales of a dataset.
  build_outputs <- function(item_subset, scale_label) {
    sub_n_items <- ncol(item_subset)
    sub_n_obs   <- rowSums(!is.na(item_subset))
    sub_means   <- rowMeans(item_subset, na.rm = TRUE)
    sub_means[!is.finite(sub_means)] <- NA_real_
    sub_sumscore <- sub_means * sub_n_items

    sub_obs_min <- suppressWarnings(min(item_subset[keep, , drop = FALSE], na.rm = TRUE))
    sub_obs_max <- suppressWarnings(max(item_subset[keep, , drop = FALSE], na.rm = TRUE))
    if (!is.finite(sub_obs_min)) sub_obs_min <- NA_real_
    if (!is.finite(sub_obs_max)) sub_obs_max <- NA_real_

    ss <- tibble::tibble(
      id           = participant_ids,
      scale        = scale_label,
      source       = source_id,
      n_items      = sub_n_items,
      n_observed   = sub_n_obs,
      sumscore     = sub_sumscore,
      min_response = row$response_min,
      max_response = row$response_max,
      observed_min = sub_obs_min,
      observed_max = sub_obs_max
    ) |>
      dplyr::filter(keep & sub_n_obs > 0L)

    keep_with_data <- keep & (sub_n_obs > 0L)
    if (sum(keep_with_data) > 0L) {
      items_kept <- item_subset[keep_with_data, , drop = FALSE]
      ids_kept   <- participant_ids[keep_with_data]
      n_kept_p   <- nrow(items_kept)
      lng <- tibble::tibble(
        id           = rep(ids_kept, times = sub_n_items),
        subscale     = source_id,
        item         = rep(seq_len(sub_n_items), each = n_kept_p),
        response     = as.numeric(unlist(items_kept)),
        scale        = scale_label,
        min_response = as.numeric(row$response_min),
        max_response = as.numeric(row$response_max)
      ) |>
        dplyr::filter(!is.na(response))
    } else {
      lng <- tibble::tibble(
        id = character(0), subscale = character(0), item = integer(0),
        response = numeric(0), scale = character(0),
        min_response = numeric(0), max_response = numeric(0)
      )
    }
    list(sumscores = ss, long = lng)
  }

  mapping <- SUBSCALE_MAPPINGS[[row$data_file]]
  if (is.null(mapping)) {
    parts <- list(build_outputs(items_num, instrument))
  } else {
    # Sanity-check that every position is in range
    max_pos <- max(unlist(mapping))
    if (max_pos > ncol(items_num)) {
      stop(sprintf(
        "Subscale mapping for %s references position %d but slice has only %d cols",
        row$data_file, max_pos, ncol(items_num)
      ))
    }
    parts <- lapply(names(mapping), function(sub_name) {
      build_outputs(items_num[, mapping[[sub_name]], drop = FALSE], sub_name)
    })
  }

  out_sumscores <- dplyr::bind_rows(lapply(parts, `[[`, "sumscores"))
  out_long      <- dplyr::bind_rows(lapply(parts, `[[`, "long"))

  list(ok = TRUE, reason = NA_character_,
       df = out_sumscores,
       df_long = out_long,
       n_total = participant_n_total,
       n_kept = participant_n_kept,
       n_excluded = participant_n_excluded,
       raw_observed_min = raw_observed_min,
       raw_observed_max = raw_observed_max,
       scale = instrument,
       source = source_id,
       min_response = row$response_min,
       max_response = row$response_max)
}

# ----------------------------------- accessor for raw item-level data ------

# Re-extract item-level data for a single dataset, applying the same sentinel
# filter and per-participant range exclusion as the main ingestion. Useful for
# downstream analyses that need raw items, not just sumscores (e.g. bootstrap
# validation of MOM vs MLE Beta fits at small samples).
get_dataset_items <- function(scale_folder_name, data_file_name,
                              catalog = brysbaert_catalog, base = base_dir) {
  row <- catalog[catalog$scale_folder == scale_folder_name &
                 catalog$data_file == data_file_name, , drop = FALSE]
  if (nrow(row) == 0) stop("Dataset not in catalog: ",
                           scale_folder_name, "/", data_file_name)
  if (nrow(row) > 1)  row <- row[1, ]
  if (!isTRUE(row$valid)) stop("Dataset is marked valid = FALSE")

  path <- file.path(base, row$scale_folder, row$data_file)
  if (!file.exists(path)) stop("File not found: ", path)
  d <- read_dataset(path)
  if (is.null(d)) stop("Failed to read: ", path)

  cols <- eval(parse(text = row$item_cols), envir = list(.df = d))
  items_num <- items_to_numeric(d[, cols, drop = FALSE])

  if (!is.na(row$response_min) && !is.na(row$response_max)) {
    floor_ <- as.numeric(row$response_min) - 10
    ceil_  <- as.numeric(row$response_max) + 10
    mask <- !is.na(items_num) & ((items_num < floor_) | (items_num > ceil_))
    items_num[as.matrix(mask)] <- NA_real_

    oor <- !is.na(items_num) &
           ((items_num < as.numeric(row$response_min)) |
            (items_num > as.numeric(row$response_max)))
    keep <- rowSums(as.matrix(oor)) == 0L
    items_num <- items_num[keep, , drop = FALSE]
  }

  list(items        = items_num,
       min_response = as.numeric(row$response_min),
       max_response = as.numeric(row$response_max),
       n_items      = ncol(items_num),
       scale        = resolve_instrument(row$scale_folder, row$data_file),
       source       = tools::file_path_sans_ext(row$data_file))
}

# Per-participant POMP score: (mean(items) - min) / (max - min). Returns a
# numeric vector on (approximately) [0, 1] with one entry per kept participant.
participant_pomp <- function(dataset) {
  rng <- dataset$max_response - dataset$min_response
  m   <- rowMeans(dataset$items, na.rm = TRUE)
  m   <- m[is.finite(m)]
  (m - dataset$min_response) / rng
}

# --------------------------------- bootstrap MOM vs MLE phi at small N -----

# For a single scale dataset, bootstrap subsamples at multiple sample sizes
# and compute MOM and MLE phi on each subsample. Returns a long data frame
# (scale, source, N, b, phi_mom, phi_mle, phi_star) where phi_star is the
# full-N MLE used as the ground-truth reference.
#
# Inputs:
#   pomp_full : numeric vector of participant POMP scores on (0, 1)
#   scale, source : labels (carried through to output)
#   N_grid : sample sizes to bootstrap (e.g. c(50, 100, 200, 500))
#   B : bootstrap reps per N
bootstrap_phi_finite_sample <- function(pomp_full, scale, source,
                                        N_grid = c(50, 100, 200, 500),
                                        B = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("betareg", quietly = TRUE)) {
    stop("Install betareg")
  }
  N_full <- length(pomp_full)
  if (N_full < max(N_grid)) {
    warning(sprintf("%s/%s: N_full = %d < max(N_grid) = %d; trimming",
                    scale, source, N_full, max(N_grid)))
    N_grid <- N_grid[N_grid <= N_full]
  }
  if (length(N_grid) == 0) return(NULL)

  # Ground-truth phi: MLE on full data (after SK transform)
  sk_full <- (pomp_full * (N_full - 1) + 0.5) / N_full
  fit_full <- tryCatch(betareg::betareg(sk_full ~ 1), error = function(e) NULL)
  phi_star <- if (is.null(fit_full)) NA_real_ else
    unname(fit_full$coefficients$precision)

  rows <- vector("list", length(N_grid) * B)
  k <- 1L
  for (N in N_grid) {
    for (b in seq_len(B)) {
      samp <- sample(pomp_full, N, replace = TRUE)
      mu   <- mean(samp)
      v    <- var(samp)
      phi_mom <- if (v > 0 && mu > 0 && mu < 1) {
                   mu * (1 - mu) / v - 1
                 } else NA_real_
      sk      <- (samp * (N - 1) + 0.5) / N
      fit     <- tryCatch(betareg::betareg(sk ~ 1), error = function(e) NULL)
      phi_mle <- if (is.null(fit)) NA_real_ else
        unname(fit$coefficients$precision)
      rows[[k]] <- data.frame(scale = scale, source = source,
                              N = N, b = b,
                              phi_mom = phi_mom, phi_mle = phi_mle,
                              phi_star = phi_star,
                              stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

# --------------------------------------------------- range-violation audit --

# For each dataset whose observed item range falls outside the catalog
# response_min/max, classify the offending values as:
#   - "rare" (< rare_pct of observations): likely sentinels or stray bad rows
#     that can be safely NA'd
#   - "common" (>= rare_pct): suggests the catalog response_min/max is wrong
#     for this dataset (different scale variant, off-by-one indexing, etc.)
audit_response_ranges <- function(catalog = brysbaert_catalog,
                                  base = base_dir,
                                  rare_pct = 1.0) {
  active <- subset(catalog, valid)
  rows <- vector("list", nrow(active))
  for (i in seq_len(nrow(active))) {
    row <- active[i, ]
    path <- file.path(base, row$scale_folder, row$data_file)
    if (!file.exists(path)) next
    d <- read_dataset(path)
    if (is.null(d)) next
    cols <- tryCatch(eval(parse(text = row$item_cols), envir = list(.df = d)),
                     error = function(e) NULL)
    if (is.null(cols)) next
    items <- tryCatch(d[, cols, drop = FALSE], error = function(e) NULL)
    if (is.null(items)) next
    items_num <- items_to_numeric(items)

    # Apply the same gross-sentinel filter the ingestion does (so we audit the
    # data the way it'd be scored, not raw with -99 codes)
    if (!is.na(row$response_min) && !is.na(row$response_max)) {
      floor_ <- as.numeric(row$response_min) - 10
      ceil_  <- as.numeric(row$response_max) + 10
      mask <- !is.na(items_num) &
              ((items_num < floor_) | (items_num > ceil_))
      items_num[as.matrix(mask)] <- NA_real_
    }

    vals <- as.numeric(unlist(items_num))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) next

    tab <- sort(table(vals), decreasing = TRUE)
    fr  <- data.frame(value = as.numeric(names(tab)),
                      n     = as.integer(tab),
                      pct   = 100 * as.integer(tab) / sum(as.integer(tab)))
    in_range <- fr$value >= row$response_min & fr$value <= row$response_max
    rare_oor   <- fr[!in_range & fr$pct <  rare_pct, ]
    common_oor <- fr[!in_range & fr$pct >= rare_pct, ]

    rows[[i]] <- data.frame(
      scale_folder = row$scale_folder,
      source = tools::file_path_sans_ext(row$data_file),
      catalog_min = row$response_min,
      catalog_max = row$response_max,
      observed_min = min(vals),
      observed_max = max(vals),
      n_obs = length(vals),
      pct_oor_total = round(sum(fr$pct[!in_range]), 3),
      rare_oor_pct  = round(sum(rare_oor$pct), 3),
      rare_oor      = paste0(rare_oor$value, "(", sprintf("%.2f%%", rare_oor$pct),
                             ")", collapse = "; "),
      common_oor    = paste0(common_oor$value,
                             "(", sprintf("%.1f%%", common_oor$pct), ")",
                             collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows[!sapply(rows, is.null)])
  out$verdict <- with(out, ifelse(
    pct_oor_total < 0.001, "in_range",
    ifelse(nchar(common_oor) == 0, "rare_only_drop_them", "catalog_likely_wrong")
  ))
  out[order(out$verdict, -out$pct_oor_total), ]
}

# Iterate
catalog_active     <- subset(brysbaert_catalog, valid)
results_sumscores  <- vector("list", nrow(catalog_active))
results_long       <- vector("list", nrow(catalog_active))
log_rows           <- vector("list", nrow(catalog_active))

for (i in seq_len(nrow(catalog_active))) {
  row <- catalog_active[i, ]
  res <- ingest_one(row)
  results_sumscores[[i]] <- res$df
  results_long[[i]]      <- res$df_long
  log_rows[[i]] <- tibble::tibble(
    scale_folder     = row$scale_folder,
    data_file        = row$data_file,
    scale            = if (is.null(res$scale))            NA_character_ else res$scale,
    source           = if (is.null(res$source))           NA_character_ else res$source,
    ok               = res$ok,
    reason           = res$reason,
    n_total          = if (is.null(res$n_total))          NA_integer_   else res$n_total,
    n_kept           = if (is.null(res$n_kept))           NA_integer_   else res$n_kept,
    n_excluded       = if (is.null(res$n_excluded))       NA_integer_   else res$n_excluded,
    n_rows           = if (is.null(res$df))               0L            else nrow(res$df),
    min_response     = if (is.null(res$min_response))     NA_real_      else res$min_response,
    max_response     = if (is.null(res$max_response))     NA_real_      else res$max_response,
    raw_observed_min = if (is.null(res$raw_observed_min)) NA_real_      else res$raw_observed_min,
    raw_observed_max = if (is.null(res$raw_observed_max)) NA_real_      else res$raw_observed_max
  )
}

dat_brysbaert_sumscores <- dplyr::bind_rows(results_sumscores)
dat_brysbaert_long      <- dplyr::bind_rows(results_long)
brysbaert_load_log      <- dplyr::bind_rows(log_rows)

# --------------------------------------------------------------- diagnostics --

n_total_part    <- sum(brysbaert_load_log$n_total, na.rm = TRUE)
n_excluded_part <- sum(brysbaert_load_log$n_excluded, na.rm = TRUE)
cat(sprintf(
  "Brysbaert ingestion: %d/%d datasets loaded (%d skipped/failed).\n",
  sum(brysbaert_load_log$ok), nrow(brysbaert_load_log),
  sum(!brysbaert_load_log$ok)
))
cat(sprintf(
  "Participants: %d kept / %d total (%d excluded for out-of-range item responses).\n",
  nrow(dat_brysbaert_sumscores), n_total_part, n_excluded_part
))

failures <- subset(brysbaert_load_log, !ok)
if (nrow(failures) > 0) {
  cat("\nFailed entries (refine catalog or fix file paths to recover):\n")
  print(failures, n = nrow(failures))
}

# Datasets where any participants were excluded (sorted by exclusion rate)
excl <- subset(brysbaert_load_log,
               ok & !is.na(n_excluded) & n_excluded > 0)
if (nrow(excl) > 0) {
  excl$pct_excluded <- round(100 * excl$n_excluded / excl$n_total, 2)
  excl <- excl[order(-excl$pct_excluded), ]
  cat("\nDatasets with participant exclusions (>10% suggests catalog response_min/max ",
      "is wrong rather than data quality):\n", sep = "")
  print(excl[, c("scale_folder", "data_file", "n_total", "n_excluded", "pct_excluded")],
        n = nrow(excl))
}

# Catalog vs observed range mismatch (suggests wrong response_min/max guess).
# Uses the RAW observed range (after sentinel NA, before per-participant
# exclusion); this is what's informative for diagnosing catalog errors,
# since the post-exclusion range is always inside the catalog by construction.
range_check <- brysbaert_load_log |>
  dplyr::filter(ok) |>
  dplyr::transmute(
    scale, source, min_response, max_response,
    observed_min = raw_observed_min,
    observed_max = raw_observed_max,
    min_off = !is.na(min_response) & is.finite(observed_min) &
              observed_min < min_response,
    max_off = !is.na(max_response) & is.finite(observed_max) &
              observed_max > max_response
  ) |>
  dplyr::filter(min_off | max_off)

if (nrow(range_check) > 0) {
  cat("\nDatasets where observed range exceeds catalog response_min/max ",
      "(catalog guess may be wrong):\n", sep = "")
  print(range_check, n = nrow(range_check))
}

# --------------------------------------------------------------- save output --

if (sys.nframe() == 0) {
  # Scale-level exclusion: drop datasets where >=10% of participants had
  # at least one out-of-range item response. These are mostly catalog errors
  # (different scale variant, off-by-one indexing) rather than data quality
  # issues, but until they're fixed the whole dataset is suspect.
  exclusions <- brysbaert_load_log |>
    dplyr::mutate(prop_excluded = dplyr::if_else(
      is.na(n_total) | n_total == 0L,
      NA_real_,
      n_excluded / n_total
    ))

  # Join on `source` only: it's unique per dataset, and each dataset may now
  # produce multiple split scales (BIS/BAS, dass depression/anxiety/stress,
  # etc.) whose `scale` values won't match the load_log's instrument-level
  # scale. The dataset-level prop_excluded applies uniformly to all subscales.
  clean_sources <- exclusions |>
    dplyr::filter(ok, !is.na(prop_excluded), prop_excluded < 0.10) |>
    dplyr::pull(source) |>
    unique()

  dat_sumscores_clean <- dat_brysbaert_sumscores |>
    dplyr::filter(source %in% clean_sources)

  # Long file uses subscale = source.
  dat_long_clean <- dat_brysbaert_long |>
    dplyr::filter(subscale %in% clean_sources)

  # Write as parquet (compact + fast). Read back with arrow::read_parquet().
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Install arrow: install.packages('arrow')")
  }
  out_dir <- file.path(dirname(base_dir), "processed")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  sumscores_file <- file.path(out_dir, "dat_sumscores_brysbaert.parquet")
  long_file      <- file.path(out_dir, "dat_long_brysbaert.parquet")
  arrow::write_parquet(dat_sumscores_clean, sumscores_file)
  arrow::write_parquet(dat_long_clean,      long_file)

  cat(sprintf(
    "\nApplied scale-level exclusion (prop_excluded < 0.10): %d / %d datasets kept.\n",
    length(clean_sources),
    sum(brysbaert_load_log$ok, na.rm = TRUE)
  ))
  cat(sprintf("Wrote %s (%d rows, %s).\n",
              sumscores_file, nrow(dat_sumscores_clean),
              format(structure(file.size(sumscores_file), class = "object_size"),
                     units = "auto")))
  cat(sprintf("Wrote %s (%d rows, %s).\n",
              long_file, nrow(dat_long_clean),
              format(structure(file.size(long_file), class = "object_size"),
                     units = "auto")))
}
