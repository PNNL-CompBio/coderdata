#!/bin/bash

# Define two parallel arrays: one for target classes, another for file names
target_classes=("Sample" "Transcriptomics" "Proteomics" "Mutations" "Copy Number" "Experiments" "Drug")
files=("/tmp/GDSCv2_samples.csv" "/tmp/GDSCv2_transcriptomics.csv" "/tmp/GDSCv2_proteomics.csv" "/tmp/GDSCv2_mutations.csv" "/tmp/GDSCv2_copy_number.csv" "/tmp/GDSCv2_experiments.csv" "/tmp/GDSCv2_drugs.tsv")

# Initialize a flag to track validation status
validation_failed=0

# Get the length of the arrays
array_length=${#target_classes[@]}

# Loop through the arrays
for (( i=0; i<${array_length}; i++ )); do
  target_class=${target_classes[$i]}
  file=${files[$i]}
  echo "Validating $target_class in file $file..."

  # Run the validation command
  linkml-validate --schema schema/coderdata.yaml --target-class "$target_class" "$file"

  # Capture the exit status
  status=$?

  # Check the exit status of the command
  if [ $status -ne 0 ]; then
    echo "Validation failed for $target_class in file $file."
    validation_failed=1
  else
    echo "Validation succeeded for $target_class in file $file."
  fi
done

# Check if any validations failed
if [ $validation_failed -ne 0 ]; then
  echo "One or more validations failed. Exiting with error."
  exit 1
else
  echo "All validations succeeded."
fi

echo "Validation process completed."

