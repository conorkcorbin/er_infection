import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
sns.set(style='white', font_scale=1.0)
import numpy as np
from pulp import *
import os, glob
import pdb
from tqdm import tqdm

from integer_programming import get_clinician_prescribing_patterns

def load_predictions():
    """Helper function that loads predictions from AST classifiers for test set data"""
    abx_options = ["Vancomycin",
                   "Ampicillin",
                   "Cefazolin",
                   "Ceftriaxone",
                   "Cefepime",
                   "Zosyn",
                   "Ciprofloxacin",
                   "Meropenem",
                   "Vancomycin_Meropenem",
                   "Vancomycin_Zosyn",
                   "Vancomycin_Cefepime",
                   "Vancomycin_Ceftriaxone"
                   ]
    base_path="/Users/conorcorbin/repos/er_infection/results/ast_models/testing/{abx}"
    df = pd.DataFrame()
    for i, abx in enumerate(abx_options):
        path = base_path.format(abx=abx)
        f_path = glob.glob(os.path.join(path, '*predictions.csv'))[0]
        if i == 0:
            df = pd.read_csv(f_path)
            df = df[['anon_id', 'pat_enc_csn_id_coded', 'label', 'predictions']]
            df = df.rename(columns={'label' : '%s_label' % abx,
                                    'predictions' : '%s_predictions' % abx})
        else:
            df_preds = pd.read_csv(f_path)
            df_preds = df_preds[['anon_id', 'pat_enc_csn_id_coded', 'label', 'predictions']]
            df_preds = df_preds.rename(columns={'label' : '%s_label' % abx,
                                                'predictions' : '%s_predictions' % abx})
            df = df.merge(df_preds, how='left', on=['anon_id', 'pat_enc_csn_id_coded'])
    
    return df


class AbxDecisionMaker():

    def __init__(self, df_predictions, df_drugs, abx_settings):
        self.df_predictions = df_predictions
        self.df_drugs = df_drugs
        self.df = (df_predictions
            .merge(df_drugs, how='inner', on='pat_enc_csn_id_coded')
        )
        self.abx_settings = abx_settings
        self.abx_options = ["Vancomycin",
                            "Ampicillin",
                            "Cefazolin",
                            "Ceftriaxone",
                            "Cefepime",
                            "Zosyn",
                            "Ciprofloxacin",
                            "Meropenem",
                            "Vancomycin_Meropenem",
                            "Vancomycin_Zosyn",
                            "Vancomycin_Cefepime",
                            "Vancomycin_Ceftriaxone"
                            ]

        self.abx_map = {'Ceftriaxone' : "CEFTRIAXONE",
                        'Vancomycin_Zosyn' : "PIPERACILLIN-TAZOBACTAM VANCOMYCIN",
                        'Zosyn' : "PIPERACILLIN-TAZOBACTAM",
                        'Vancomycin_Ceftriaxone' : "CEFTRIAXONE VANCOMYCIN",
                        'Vancomycin_Cefepime' : "CEFEPIME VANCOMYCIN",
                        'Cefepime' : "CEFEPIME",
                        'Vancomycin' :  "VANCOMYCIN",
                        'Meropenem' : "MEROPENEM",
                        'Vancomycin_Meropenem' : "MEROPENEM VANCOMYCIN",
                        'Cefazolin' : "CEFAZOLIN",
                        'Ciprofloxacin' : "CIPROFLOXACIN",
                        'Ampicillin' : 'AMPICILLIN'
                        }
        self.abx_map_inverse = {self.abx_map[key] : key for key in self.abx_map}

    def set_abx_settings(self, abx_settings):
        self.abx_settings = abx_settings

    def compute_was_covered(self, x, decision_column='med_description'):
        """
        Given med description, find appropriate label column and return whether patient was covered during CSN
        Returns "Not in abx options" if abx regimen isn't in our set of 12 options - useful for filtering later
        """
        if decision_column == 'med_description':
            med_description = x.med_description
        elif decision_column == 'random_med_description':
            med_description = x.random_med_description
        elif decision_column == 'IP_med_description':
            med_description = x.IP_med_description
            
        if med_description == "CEFTRIAXONE":
            return x.Ceftriaxone
        elif med_description == "PIPERACILLIN-TAZOBACTAM VANCOMYCIN":
            return x.Vancomycin_Zosyn
        elif med_description == "PIPERACILLIN-TAZOBACTAM":
            return x.Zosyn
        elif med_description == "CEFTRIAXONE VANCOMYCIN":
            return x.Vancomycin_Ceftriaxone
        elif med_description == "CEFEPIME VANCOMYCIN":
            return x.Vancomycin_Cefepime
        elif med_description == "CEFEPIME":
            return x.Cefepime
        elif med_description == "VANCOMYCIN":
            return x.Vancomycin
        elif med_description == "MEROPENEM":
            return x.Meropenem
        elif med_description == "MEROPENEM VANCOMYCIN":
            return x.Vancomycin_Meropenem
        elif med_description == "CEFAZOLIN":
            return x.Cefazolin
        elif med_description == "CIPROFLOXACIN":
            return x.Ciprofloxacin
        elif med_description == "AMPICILLIN":
            return x.Ampicillin
        else:
            return "Not in abx options"

    def get_coverage_rates(self):
        """
        Create flag for whether clinicians covered the patient during the csn, whether a random assignemnt
        covered patient CSN, and whether optimized assignment covered the patient CSN
        """

        df = (self.df
            .assign(random_med_description=lambda x: np.random.choice(x.IP_med_description,
                                                                      size=len(x.IP_med_description),
                                                                      replace=False))
        )
        df = (df
            .assign(was_covered_dr=df.apply(lambda x: self.compute_was_covered(x), axis=1))
            .assign(was_covered_random=df.apply(lambda x: self.compute_was_covered(x, 
                                                decision_column='random_med_description'),
                                                axis=1))
            .assign(was_covered_IP=df.apply(lambda x: self.compute_was_covered(x, 
                                            decision_column='IP_med_description'),
                                            axis=1))
        )

        clin_covered_rate = df['was_covered_dr'].sum() / len(df)
        random_covered_rate = df['was_covered_random'].sum() / len(df)
        ip_covered_rate = df['was_covered_IP'].sum() / len(df)
        
        return random_covered_rate, clin_covered_rate, ip_covered_rate

    def solve_and_assign(self):

        # Predictions string
        predictions_string = '%s_predictions'
        abx_model = LpProblem("Antibiotics", LpMaximize)

        # Create binary indicators for whether treatment is used
        drug_inds = {}
        for abx in self.abx_options:
            drug_inds[abx] = [LpVariable('%s_%d' % (abx, i), lowBound=0, upBound=1, cat='Binary')
                            for i in range(len(self.df))]

        # Add objective function to model
        per_csn_sum = []
        for i in range(len(self.df)):
            _sum = 0
            for abx in self.abx_options:
                _sum += drug_inds[abx][i] * self.df[predictions_string % abx].values[i]
            per_csn_sum.append(_sum)
            
        abx_model += lpSum(per_csn_sum)

        # Add one selection constraint
        for i in range(len(self.df)):
            selections = []
            for abx in self.abx_options:
                selections.append(drug_inds[abx][i])
            abx_model += lpSum(selections) == 1

        for drug in drug_inds:
            abx_model += lpSum([drug_inds[drug][i] for i in range(len(self.df))]) == self.abx_settings[drug]

        # Solve model
        abx_model.solve()

        # print("Status:", LpStatus[abx_model.status])
        # Save selected antibiotic to df
        abx_decisions = []
        for i in range(len(self.df)):
            abx_decision = None
            for abx in self.abx_options:
                if drug_inds[abx][i].value() == 1:
                    abx_decision = self.abx_map[abx]
            assert abx_decision is not None
            abx_decisions.append(abx_decision)
        self.df['IP_med_description'] = abx_decisions

    
def perform_abx_sweep():
    """
    Simulates antibiotic delivery sweeping through different abx prescribing contraints. 
    For each abx pair, we add the number of times the two were prescribed in total (N) in actual practice
    and then sweep contraints from the extreme where only the first antbiotic prescribed N times to the other
    extreme where the other antibiotic is prescribed N times. We show how coverage rate changes as we change 
    these contraints for both the IP decision and a random decision.
    """
    abx_settings = {"Vancomycin" : 13,
                "Ampicillin" : 0,
                "Cefazolin" : 8,
                "Ceftriaxone" : 367,
                "Cefepime" : 14,
                "Zosyn" : 102,
                "Ciprofloxacin" : 8,
                "Meropenem" : 9,
                "Vancomycin_Meropenem" : 9,
                "Vancomycin_Zosyn" :  113,
                "Vancomycin_Cefepime" : 23,
                "Vancomycin_Ceftriaxone" : 31
                }
    df_predictions = load_predictions()
    df_drugs = get_clinician_prescribing_patterns()
    opt = AbxDecisionMaker(df_predictions, df_drugs, abx_settings)
    opt.solve_and_assign()
    random_covered_rate, clin_covered_rate, ip_covered_rate = opt.get_coverage_rates()

    for abx_to_perturb in abx_settings:
        plt.figure(figsize=(32,24))
        # fig, ax = plt.subplots(3, 4, figsize=(32, 24))
        gs = gridspec.GridSpec(4, 24, wspace=2.0)
        ax1a = plt.subplot(gs[0, 0:6])
        ax1b = plt.subplot(gs[0, 6:12])
        ax1c = plt.subplot(gs[0, 12:18])
        ax1d = plt.subplot(gs[0, 18:24])
        
        ax2a = plt.subplot(gs[1, 3:9])
        ax2b = plt.subplot(gs[1, 9:15])
        ax2c = plt.subplot(gs[1, 15:21])

        ax3a = plt.subplot(gs[2, 0:6])
        ax3b = plt.subplot(gs[2, 6:12])
        ax3c = plt.subplot(gs[2, 12:18])
        ax3d = plt.subplot(gs[2, 18:24])

        axes = [ax1a, ax1b, ax1c, ax1d, ax2a, ax2b, ax2c, ax3a, ax3b, ax3c, ax3d]

        skip = 0
        for ind, abx in enumerate(abx_settings):
            if abx == abx_to_perturb:
                skip = 1
                continue
            else:
                abx_settings_perturbed = {key : abx_settings[key] for key in abx_settings}
                print("Performing %s to %s sweep" % (abx_to_perturb, abx))
                random_rates, ip_rates = [], []
                # Total selections to sweep over
                total_selections = abx_settings[abx] + abx_settings[abx_to_perturb]
                abx_settings_perturbed[abx_to_perturb] = total_selections
                abx_settings_perturbed[abx] = 0
                opt.set_abx_settings(abx_settings_perturbed)
                opt.solve_and_assign()
                r, c, i = opt.get_coverage_rates()
                random_rates.append(r)
                ip_rates.append(i)
                if abx_settings_perturbed == abx_settings:
                    clin_iter = -1 # save point on x axis for clinician performance
                for iter_ in tqdm(range(total_selections)):
                    abx_settings_perturbed[abx_to_perturb] -= 1
                    abx_settings_perturbed[abx] += 1
                    opt.set_abx_settings(abx_settings_perturbed)
                    opt.solve_and_assign()
                    r, c, i = opt.get_coverage_rates()
                    random_rates.append(r)
                    ip_rates.append(i)

                    if abx_settings_perturbed == abx_settings:
                        clin_iter = iter_ # save point on x axis for clinician performance
                    
            axes[ind-skip].plot(range(total_selections+1), random_rates, label='Random Assignment')
            axes[ind-skip].plot(range(total_selections+1), ip_rates, label='Integer Programming')
            axes[ind-skip].plot(clin_iter+1, clin_covered_rate, marker='o', label='Clinician Benchmark')
            # forward = lambda x: total_selections - x
            # backward = lambda x: total_selections - x
            # secax = axes[ind-skip].secondary_xaxis('top', functions=(forward, backward))
            axes[ind-skip].set_xlabel("Num %s Administered" % abx)
            # secaxes[ind-skip].set_xlabel("Num %s Administered" % abx_to_perturb)
            axes[ind-skip].set_ylabel("Coverage Rate")
            axes[ind-skip].set_ylim((0., 1.))
            axes[ind-skip].set_title("%s to %s Sweep" % (abx_to_perturb, abx))
            axes[ind-skip].legend()
        
        fig_name = "%s_sweeps.jpg" % abx_to_perturb
        plt.savefig(fig_name)




def main():
    abx_settings = {"Vancomycin" : 0,
                    "Ampicillin" : 0,
                    "Cefazolin" : 0,
                    "Ceftriaxone" : 0,
                    "Cefepime" : 0,
                    "Zosyn" : 0,
                    "Ciprofloxacin" : 0,
                    "Meropenem" : 0,
                    "Vancomycin_Meropenem" : 697,
                    "Vancomycin_Zosyn" :  0,
                    "Vancomycin_Cefepime" : 0,
                    "Vancomycin_Ceftriaxone" : 0
                    }
    abx_rankings = ['Vancomycin_Meropenem',
                    'Vancomycin_Zosyn',
                    'Vancomycin_Cefepime',
                    'Zosyn',
                    'Vancomycin_Ceftriaxone',
                    'Meropenem',
                    'Cefepime',
                    'Ceftriaxone',
                    'Ciprofloxacin',
                    'Cefazolin',
                    'Ampicillin',
                    'Vancomycin'
                    ]
    
    if os.path.exists('df_predictions.csv'):
        df_predictions = pd.read_csv('df_predictions.csv')
    else:
        df_predictions = load_predictions()
        df_predictions.to_csv('df_predictions.csv', index=None)


    if os.path.exists('df_drugs.csv'):
        df_drugs = pd.read_csv('df_drugs.csv')
    else:
        df_drugs = get_clinician_prescribing_patterns()
        df_drugs.to_csv('df_drugs.csv', index=None)

    deplete_idx = 0
    push_from_idx = 0
    push_to_idx = 1
    random_rates = []
    ip_rates = []
    counter = 0
    while abx_settings['Vancomycin'] != 697: # end of sweep
        opt = AbxDecisionMaker(df_predictions, df_drugs, abx_settings)
        opt.solve_and_assign()
        r, c, i = opt.get_coverage_rates()
        random_rates.append(r)
        ip_rates.append(i)

        if counter == 10:
            print(abx_settings)
            counter = 0

        # Make more narrow spectrum
        abx_settings[abx_rankings[push_from_idx]] -= 1
        abx_settings[abx_rankings[push_to_idx]] += 1
        
        if abx_settings[abx_rankings[deplete_idx]] == 0:
            deplete_idx += 1
            print("Moving deplete index from %s to %s" % (abx_rankings[deplete_idx-1], abx_rankings[deplete_idx]))

        if push_to_idx == len(abx_rankings) - 1:
            push_to_idx = deplete_idx + 1
            push_from_idx = deplete_idx
        else:
            push_from_idx += 1
            push_to_idx += 1

        counter += 1
    
    clin_rates = [c for i in range(len(random_rates))]
    df = pd.DataFrame(data={'clin_rates' : clin_rates,
                            'random_rates' : random_rates,
                            'ip_rates' : ip_rates})
    df.to_csv('summary_sweep.csv')
    fix, ax = plt.subplots(1, 1, figsize=(8,8))
    ax.plot(range(len(random_rates)), random_rates, label='Random Assignment')
    ax.plot(range(len(random_rates)), ip_rates, label='Integer Programming')
    ax.plot(range(len(random_rates)), c, label='Clinician Benchmark')

    fig_name = "summary_sweep.jpg"
    plt.savefig(fig_name)

def test_waterfall():
    abx_settings = {"a" : 2,
                     "b" : 0,
                     "c" : 0
                    }
    abx_rankings = ['a',
                    'b',
                    'c'
                    ]  
    true_settings = [{"a" : 2,
                     "b" : 0,
                     "c" : 0
                    },
                    {"a" : 1,
                     "b" : 1,
                     "c" : 0,
                    },
                    {"a" : 1,
                     "b" : 0,
                     "c" : 1
                    },
                    {"a" : 0,
                     "b" : 1,
                     "c" : 1
                    },
                    {"a" : 0,
                     "b" : 0,
                     "c" : 2
                    }]                 

    settings = []
    deplete_idx = 0
    push_from_idx = 0
    push_to_idx = 1
    import copy
    while abx_settings['c'] != 2: # end of sweep
        settings.append(copy.deepcopy(abx_settings))
        print("Deplete idx:{d} Push From Idx:{f} Push To Idx:{t}".format(d=deplete_idx,
                                                                         f=push_from_idx,
                                                                         t=push_to_idx))
        # Make more narrow spectrum
        abx_settings[abx_rankings[push_from_idx]] -= 1
        abx_settings[abx_rankings[push_to_idx]] += 1

        if abx_settings[abx_rankings[deplete_idx]] == 0:
            deplete_idx += 1

        if push_to_idx == len(abx_rankings)-1:
            push_to_idx = deplete_idx + 1
            push_from_idx = deplete_idx
        else:
            push_from_idx += 1
            push_to_idx += 1

    for i, s in enumerate(settings):
        try:
            assert s == true_settings[i]
        except:
            pdb.set_trace()

if __name__ == '__main__':
    main()