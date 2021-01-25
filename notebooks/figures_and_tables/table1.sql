WITH cohort AS (
SELECT DISTINCT
  c.anon_id, c.pat_enc_csn_id_coded, c.index_time, a.anon_id a_anon_id
FROM 
  `mining-clinical-decisions.abx.cohort_not_infected_rules` c
LEFT JOIN  
  `mining-clinical-decisions.abx.final_ast_labels` a
USING
  (pat_enc_csn_id_coded)
WHERE EXTRACT(YEAR FROM c.index_time) BETWEEN 2009 AND 2019
),

adt_dep as (
SELECT DISTINCT
  adt.pat_enc_csn_id_coded, 
  FIRST_VALUE(dm.department_name) OVER 
  (PARTITION BY adt.pat_enc_csn_id_coded ORDER BY adt.effective_time_jittered_utc) department_name,
FROM 
  `shc_core.adt` adt
INNER JOIN
  `som-nero-phi-jonc101.shc_core.dep_map` dm
USING
  (department_id)
)

SELECT DISTINCT 
  dep.department_name,
  DATE_DIFF(CAST(c.index_time as DATE), d.BIRTH_DATE_JITTERED, year) age,
  c.pat_enc_csn_id_coded,
  d.ANON_ID, d.GENDER, d.CANONICAL_RACE, d.CANONICAL_ETHNICITY,
  CASE WHEN d.LANGUAGE = "English" THEN "English"
  ELSE "Non-English" END LANGUAGE,
  CASE WHEN d.INSURANCE_PAYOR_NAME = "MEDICARE" THEN "Medicare"
  WHEN d.INSURANCE_PAYOR_NAME = "MEDI-CAL" THEN "Medi-Cal"
  ELSE "Other" END INSURANCE_PAYOR_NAME,
  CASE WHEN c.a_anon_id IS NULL THEN 0 ELSE 1 END culture_growth
FROM 
  `som-nero-phi-jonc101.shc_core.demographic` d
INNER JOIN
  cohort c
ON
  d.ANON_ID = c.anon_id
INNER JOIN
  adt_dep dep
ON
  c.pat_enc_csn_id_coded = dep.pat_enc_csn_id_coded

