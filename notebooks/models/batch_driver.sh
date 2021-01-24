data_dir=/home/conorcorbin/repos/er_infection/data/
output_dir=/home/conorcorbin/repos/er_infection/results/ast_models/validation/
model_classes=(lightgbm ridge) #lasso random_forest)
labels=(Ampicillin Ciprofloxacin Cefazolin Ceftriaxone Cefepime Zosyn Vancomycin Meropenem 
Vancomycin_Ceftriaxone Vancomycin_Cefepime Vancomycin_Zosyn Vancomycin_Meropenem)
log_dir=/home/conorcorbin/logs/

for mc in ${model_classes[@]}
    do
    for l in ${labels[@]}
        do
        python train_model.py --model_class $mc \
                            --data_dir $data_dir \
                            --label $l \
                            --output_dir ${output_dir}/${mc}/${l} \
                            --val &> ${log_dir}${mc}_${l}_validation.txt & 
        done
    done