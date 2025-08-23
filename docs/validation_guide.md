# ConTra Pipeline Validation with C. elegans Test Dataset

## Quick Start Guide

### 1. Create Test Dataset
```bash
cd /path/to/ConTra
python3 code/create_celegans_test_dataset.py
```

### 2. Run ConTra Analysis
```bash
# Backup original data (optional)
cp -r data/cleaned_datasets data/cleaned_datasets_backup

# Use test dataset
cp data/test_datasets/celegans_modENCODE/*.csv data/cleaned_datasets/

# Run analysis
cd code
python3 subset_context_dependent_analysis.py
```

### 3. Validate Results
```bash
# Replace with your actual output directory
python3 code/validate_celegans_results.py output/subset_context_dependent_analysis_YYYYMMDD_HHMMSS
```

### 4. Restore Original Data
```bash
rm -rf data/cleaned_datasets
mv data/cleaned_datasets_backup data/cleaned_datasets
```

## Expected Validation Results

When running the validation, you should see:
- **Sensitivity**: ≥ 0.8 (80% of known interactions detected)
- **Specificity**: ≥ 0.8 (low false positive rate)
- **Known interactions detected**:
  - let-7 → lin-41, hbl-1 (developmental timing)
  - lin-4 → lin-14, lin-28 (larval development)

## Troubleshooting

### If validation fails:
1. **Low sensitivity**: ConTra may need parameter tuning for correlation thresholds
2. **Low specificity**: Consider increasing statistical stringency
3. **No interactions detected**: Check that test data was copied correctly

### If analysis crashes:
1. Ensure all dependencies are installed: `pip install -r code/requirements.txt`
2. Check available memory (test dataset is small, should work on any system)
3. Verify file formats match ConTra expectations

## Files Generated

### Test Dataset Creation:
- `data/test_datasets/celegans_modENCODE/`: Complete test dataset
- `validation_metadata.json`: Known interactions and validation criteria
- `README.md`: Dataset documentation

### ConTra Analysis:
- `output/subset_context_dependent_analysis_*/`: Analysis results
- `tables/*.csv`: Interaction tables
- `plots/*.png`: Visualization plots
- `reports/*.{md,html}`: Analysis reports

### Validation:
- `validation_results.png`: Performance metrics plots
- `validation_report_*.md`: Detailed validation report

## Integration into Development Workflow

This test dataset can be used for:

1. **Continuous Integration**: Automated testing of ConTra pipeline
2. **Method Development**: Rapid testing of algorithm changes
3. **Parameter Optimization**: Tuning correlation thresholds and statistical criteria
4. **Performance Benchmarking**: Comparing different versions or configurations
5. **Documentation**: Demonstrating expected pipeline behavior

## Customization

To adapt for other model organisms:

1. **Modify `create_celegans_test_dataset.py`**:
   - Update gene names and known interactions
   - Adjust sample structure and conditions
   - Include organism-specific regulatory patterns

2. **Update `validate_celegans_results.py`**:
   - Define new validation criteria
   - Specify expected regulatory relationships
   - Adjust statistical thresholds

3. **Create new metadata**:
   - Document known regulatory networks
   - Define validation success criteria
   - Include literature references

## Literature References

- **C. elegans timing**: Slack & Ruvkun (1997) Cell 88:635-645
- **let-7 miRNA**: Reinhart et al. (2000) Nature 403:901-906
- **lin-4 cascade**: Lee et al. (1993) Cell 75:843-854
- **modENCODE**: modENCODE Consortium (2010) Science 330:1775-1787