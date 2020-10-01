# Gets all antibiotic orders that happen within 24 hours
WITH er_admits AS (
SELECT pat_enc_csn_id_coded, min(effective_time_jittered_utc) as er_admit_time, max(effective_time_jittered_utc) as er_transfer_out_time
FROM `shc_core.adt`
WHERE pat_class_c = "112" AND pat_service = "Emergency"
GROUP BY pat_enc_csn_id_coded),

-- Prelim cohort before filtering
cohort_and_index_times AS (
SELECT abx.anon_id, abx.pat_enc_csn_id_coded, min(order_time) index_time
FROM `mining-clinical-decisions.abx.abx_orders_within_24_hrs` abx
INNER JOIN `mining-clinical-decisions.abx.culture_orders_within_24_hrs` cult 
USING (pat_enc_csn_id_coded)
GROUP BY anon_id, pat_enc_csn_id_coded
),

-- CSNS in cohort that have no abx ordered in two weeks prior to er admit time
has_no_prior_abx_2_weeks AS (
SELECT DISTINCT cit.pat_enc_csn_id_coded
FROM cohort_and_index_times cit
LEFT JOIN
  # Finds all csns' where abx were ordered between 1 hours and 14 weeks befor admit. 12 hours
  # gives a little bit of wiggle room because I'm not sure how accurate order_start_time is. 
  (SELECT DISTINCT cit.pat_enc_csn_id_coded
  FROM cohort_and_index_times cit
  LEFT JOIN er_admits ea
  USING (pat_enc_csn_id_coded)
  INNER JOIN (
     -- Gets columns we need from order med and finds order time from 2018 extract when missing in 2020 extract 
     SELECT om.anon_id, om.pat_enc_csn_id_coded, om.med_description, om.order_med_id_coded, om.discon_time_jittered_utc,
     CASE WHEN om.order_inst_utc IS NULL THEN omm.order_time_jittered_utc ELSE om.order_inst_utc END order_time
     FROM `shc_core.order_med` om
     LEFT JOIN `starr_datalake2018.order_med` omm
     USING (order_med_id_coded)
  ) om
  USING (anon_id)
  INNER JOIN `mining-clinical-decisions.abx.abx_types` abx_types
  USING (med_description)
  WHERE TIMESTAMP_DIFF(ea.er_admit_time, om.order_time, HOUR) BETWEEN 1 AND 24*14 
  AND abx_types.affects_not_infected_label = 1) has_prior_order
USING (pat_enc_csn_id_coded)
WHERE has_prior_order.pat_enc_csn_id_coded IS NULL
),

-- CSNs in cohort that have no cultures ordered within 2 weeks leading up until er admit
has_no_prior_cultures AS (
SELECT DISTINCT cit.pat_enc_csn_id_coded
FROM cohort_and_index_times cit
LEFT JOIN
  # Finds all csns where cultures were ordered between 1 hours and 14 weeks before admit
  (SELECT DISTINCT cit.pat_enc_csn_id_coded
   FROM cohort_and_index_times cit
   LEFT JOIN er_admits ea
   USING (pat_enc_csn_id_coded)
   INNER JOIN `shc_core.order_proc` op
   USING (anon_id)
   INNER JOIN `mining-clinical-decisions.abx.culture_types` culture_types
   USING (description)
   WHERE TIMESTAMP_DIFF(ea.er_admit_time, op.order_time_jittered_utc, HOUR) BETWEEN 1 AND 24*14
   AND culture_types.affects_not_infected_label = 1) has_prior_order
USING (pat_enc_csn_id_coded)
WHERE has_prior_order.pat_enc_csn_id_coded IS NULL
),

-- Removes examples from cohort where cultures or abx were ordered in 2 weeks prior to er admit time
filtered_examples AS (
SELECT cit.*
FROM cohort_and_index_times cit
INNER JOIN has_no_prior_abx_2_weeks 
USING (pat_enc_csn_id_coded)
INNER JOIN has_no_prior_cultures 
USING (pat_enc_csn_id_coded)
),

-- Gets table of max abx_stop_times for abx that were ordered within first two weeks and given
last_abx_stop_time_if_started as (
SELECT anon_id, pat_enc_csn_id_coded, max(abx_stop_time) abx_stop_time
FROM `mining-clinical-decisions.abx.abx_orders_given_and_stopped_2_weeks`
WHERE was_given = 1
GROUP BY anon_id, pat_enc_csn_id_coded
),

-- Gets CSN's where all abx ordered within 2 weeks were either not started or were stopped within a 2 days of er admit time. 
started_abx_stopped_within_day_after_cult_results_or_not_started AS (
SELECT DISTINCT fe.pat_enc_csn_id_coded
FROM filtered_examples fe
LEFT JOIN last_abx_stop_time_if_started las
USING (pat_enc_csn_id_coded)
LEFT JOIN er_admits ea
USING (pat_enc_csn_id_coded)
WHERE las.abx_stop_time IS NULL -- its null for csn's where abx were never started
OR las.abx_stop_time <= TIMESTAMP_ADD (ea.er_admit_time, INTERVAL 48 HOUR) -- abx that were started must be stopped 1 day after last culture result time
),

-- Gets CSNS where we say patients were discharge with oral abx defined as there exists an order for oral abx within the interval
-- of 12 hours before their discharge and 5 days after. 
discharged_with_orals AS (
SELECT DISTINCT fe.pat_enc_csn_id_coded
FROM filtered_examples fe
INNER JOIN 
  (SELECT fe.pat_enc_csn_id_coded, max(adt.effective_time_jittered_utc) discharge_time
  FROM filtered_examples fe
  INNER JOIN `shc_core.adt` adt
  USING (pat_enc_csn_id_coded)
  GROUP BY fe.pat_enc_csn_id_coded) discharge
USING (pat_enc_csn_id_coded)
LEFT JOIN (
     -- Gets columns we need from order med and finds order time from 2018 extract when missing in 2020 extract 
     SELECT om.anon_id, om.pat_enc_csn_id_coded, om.med_description, om.order_med_id_coded, om.discon_time_jittered_utc,
     CASE WHEN om.order_inst_utc IS NULL THEN omm.order_time_jittered_utc ELSE om.order_inst_utc END order_time
     FROM `shc_core.order_med` om
     LEFT JOIN `starr_datalake2018.order_med` omm
     USING (order_med_id_coded)
) om
USING (anon_id)
INNER JOIN `mining-clinical-decisions.abx.abx_types` abx_types
USING (med_description)
WHERE abx_types.is_oral = 1
AND om.order_time BETWEEN TIMESTAMP_SUB(discharge.discharge_time, INTERVAL 12 HOUR) AND TIMESTAMP_ADD(discharge.discharge_time, INTERVAL 24*5 HOUR)
),

-- Finds csns with a positive culture (that affects label) ordered within 2 weeks of er admit
has_pos_cult_ordered_2_weeks AS (
SELECT DISTINCT fe.pat_enc_csn_id_coded
FROM filtered_examples fe
LEFT JOIN er_admits ea
USING (pat_enc_csn_id_coded)
LEFT JOIN `shc_core.order_proc` op
USING (anon_id)
LEFT JOIN `mining-clinical-decisions.abx.culture_types` culture_types
USING (description)
INNER JOIN `shc_core.culture_sensitivity` cs
USING (order_proc_id_coded)
WHERE culture_types.affects_not_infected_label = 1 
AND TIMESTAMP_DIFF(op.order_time_jittered_utc, ea.er_admit_time, HOUR) BETWEEN 0 AND 14*24
AND cs.organism <> "COAG NEGATIVE STAPHYLOCOCCUS" 
AND cs.organism NOT LIKE "%CANDIDA%" 
AND cs.organism NOT IN ('HAEMOPHILUS INFLUENZAE', 'HAEMOPHILUS PARAINFLUENZAE')
),

-- Finds positive cultures from include list of cultures that were
-- ordered within 24 hours of ER admit
-- Redundent now, but in case we want to change flag for cultures that affect 
-- labels later
has_pos_cult AS (
SELECT DISTINCT fe.pat_enc_csn_id_coded
FROM filtered_examples fe
LEFT JOIN er_admits ea
USING (pat_enc_csn_id_coded)
LEFT JOIN `shc_core.order_proc` op
USING (anon_id)
LEFT JOIN `mining-clinical-decisions.abx.culture_types` culture_types
USING (description)
INNER JOIN `shc_core.culture_sensitivity` cs
USING (order_proc_id_coded)
WHERE culture_types.include = 1 
AND TIMESTAMP_DIFF(op.order_time_jittered_utc, ea.er_admit_time, HOUR) BETWEEN 0 AND 24
AND cs.organism <> "COAG NEGATIVE STAPHYLOCOCCUS"
),

-- finds csns where patient died within 2 weeks of admits
patient_died AS (
SELECT DISTINCT fe.pat_enc_csn_id_coded
FROM filtered_examples fe
LEFT JOIN er_admits ea
USING (pat_enc_csn_id_coded)
LEFT JOIN `shc_core.demographic` demo
USING(anon_id)
WHERE demo.DEATH_DATE_JITTERED IS NOT NULL AND
DATE_DIFF(demo.DEATH_DATE_JITTERED , CAST(ea.er_admit_time AS DATE), DAY) BETWEEN 0 AND 14
)

SELECT fe.*,
CASE WHEN abx.pat_enc_csn_id_coded IS NOT NULL THEN 1 ELSE 0 END abx_not_ordered_or_stopped_within_day_after_cult_result,
CASE WHEN orals.pat_enc_csn_id_coded IS NULL THEN 1 ELSE 0 END not_discharged_with_orals,
CASE WHEN other_cult.pat_enc_csn_id_coded IS NULL THEN 1 ELSE 0 END no_pos_cult_2_weeks_later,
CASE WHEN current_cult.pat_enc_csn_id_coded IS NULL THEN 1 ELSE 0 END no_pos_cult_within_day,
CASE WHEN dead.pat_enc_csn_id_coded IS NULL THEN 1 ELSE 0 END no_death_within_2_weeks
FROM filtered_examples fe
LEFT JOIN started_abx_stopped_within_day_after_cult_results_or_not_started abx
USING (pat_enc_csn_id_coded)
LEFT JOIN discharged_with_orals orals
USING (pat_enc_csn_id_coded)
LEFT JOIN has_pos_cult_ordered_2_weeks other_cult 
USING (pat_enc_csn_id_coded)
LEFT JOIN has_pos_cult current_cult
USING (pat_enc_csn_id_coded)
LEFT JOIN patient_died dead
USING (pat_enc_csn_id_coded) 

