#!/bin/bash
yaml_file="run_config_local.yaml" # Replace with your YAML file path

# Check if file exists
if [ ! -f "$yaml_file" ]; then
    echo "Error: File '$yaml_file' not found!"
    exit 1
fi

# Use Python to extract data from YAML
python_command="import yaml; f=open('$yaml_file'); data=yaml.safe_load(f); traits=data.get('traits'); print('\\n'.join(str(i) for i in traits) if isinstance(traits, list) else traits);"

# Use command substitution to capture the output and read into an array
readarray -t output_array < <(python3 -c "$python_command")

# Loop through the array
for trait in "${output_array[@]}"; do
    echo "Trait: $trait"
    sbatch ./bash/run_lm_phenopred.sh "$trait" "$yaml_file"
done

# Consolidate results
# sbatch -p urgent ./bash/run_trainer_test.sh
# sbatch -d "$dependency_string" ./bash/consolidate_results.sh "$yaml_file"