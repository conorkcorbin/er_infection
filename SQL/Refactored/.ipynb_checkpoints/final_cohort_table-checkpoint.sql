CREATE OR REPLACE TABLE `mining-clinical-decisions.abx.final_cohort_table` AS
(
SELECT anon_id, pat_enc_csn_id_coded, index_time,
CASE WHEN abx_stopped_2_days = 1
AND not_discharged_with_orals = 1
AND no_pos_cult_2_weeks_later = 1
AND no_pos_cult_within_day = 1
AND no_inf_dx_codes = 1
THEN 1 ELSE 0 END not_infected
FROM `mining-clinical-decisions.abx.cohort_not_infected_rules`
);

CREATE OR REPLACE TABLE `mining-clinical-decisions.abx.final_cohort_table` AS (
WITH joined_table AS (
SELECT c.*, l.Cefazolin, l.Ceftriaxone, l.Cefepime, l.Zosyn, l.Vancomycin, l.Meropenem 
FROM `mining-clinical-decisions.abx.final_cohort_table` c
LEFT JOIN `mining-clinical-decisions.abx.final_ast_labels` l
USING (pat_enc_csn_id_coded)
),

labels as (SELECT anon_id, pat_enc_csn_id_coded, index_time,
CASE WHEN Cefazolin = 1 OR (Cefazolin IS NULL AND not_infected = 1) THEN 1 # if suscept or not infected
WHEN Cefazolin = 0 THEN 0 # if resistant
ELSE NULL END Cefazolin, # if no bug growth but not necessarily not infected
CASE WHEN Ceftriaxone = 1 OR (Ceftriaxone IS NULL AND not_infected = 1) THEN 1
WHEN Ceftriaxone = 0 THEN 0
ELSE NULL END Ceftriaxone,
CASE WHEN Cefepime = 1 OR (Cefepime IS NULL AND not_infected = 1) THEN 1
WHEN Cefepime = 0 THEN 0
ELSE NULL END Cefepime,
CASE WHEN Zosyn = 1 OR (Zosyn IS NULL AND not_infected = 1) THEN 1
WHEN Zosyn = 0 THEN 0
ELSE NULL END Zosyn,
CASE WHEN Vancomycin = 1 OR (Vancomycin IS NULL AND not_infected = 1) THEN 1
WHEN Vancomycin = 0 THEN 0
ELSE NULL END Vancomycin,
CASE WHEN Meropenem = 1 OR (Meropenem IS NULL AND not_infected = 1) THEN 1
WHEN Meropenem = 0 THEN 0
ELSE NULL END Meropenem
FROM joined_table)

SELECT l.*, CASE WHEN l.Cefazolin IS NULL THEN 1 ELSE 0 END label_unobserved
FROM labels l
)


-- SELECT 
-- AVG(Cefazolin) Cefazolin,
-- AVG(Ceftriaxone) Ceftriaxone,
-- AVG(Cefepime) Cefepime,
-- AVG(Zosyn) Zosyn,
-- AVG(Vancomycin) Vancomycin,
-- AVG(Meropenem) Meropenem,
-- FROM `mining-clinical-decisions.abx.final_ast_labels` 
