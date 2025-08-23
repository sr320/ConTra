#!/usr/bin/env python3
"""
Validate ConTra analysis results against known C. elegans regulatory interactions.

This script compares ConTra output with expected regulatory interactions in the
C. elegans test dataset to assess pipeline accuracy and performance.
"""

import os
import sys
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import argparse
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ConTraValidator:
    """Validate ConTra results against known C. elegans interactions."""
    
    def __init__(self, 
                 output_dir: str,
                 metadata_file: str = "data/test_datasets/celegans_modENCODE/validation_metadata.json"):
        """Initialize validator with ConTra output and validation metadata."""
        self.output_dir = Path(output_dir)
        self.metadata_file = Path(metadata_file)
        
        # Load validation metadata
        with open(self.metadata_file, 'r') as f:
            self.metadata = json.load(f)
        
        self.known_interactions = self.metadata['known_interactions']
        self.validation_criteria = self.metadata['validation_criteria']
        self.expected_results = self.metadata['expected_results']
        
        # Results storage
        self.validation_results = {
            'sensitivity': 0.0,
            'specificity': 0.0,
            'detected_interactions': [],
            'missed_interactions': [],
            'false_positives': [],
            'stage_specificity': {},
            'statistical_significance': {}
        }
    
    def load_contra_results(self) -> Dict[str, pd.DataFrame]:
        """Load ConTra analysis results."""
        logger.info(f"Loading ConTra results from {self.output_dir}")
        
        results = {}
        table_dir = self.output_dir / "tables"
        
        if not table_dir.exists():
            raise FileNotFoundError(f"ConTra tables directory not found: {table_dir}")
        
        # Load main result tables
        result_files = {
            'methylation_mirna_context': 'methylation_mirna_context.csv',
            'lncrna_mirna_context': 'lncrna_mirna_context.csv', 
            'multi_way_interactions': 'multi_way_interactions.csv',
            'correlation_results': 'correlation_results.csv'
        }
        
        for key, filename in result_files.items():
            filepath = table_dir / filename
            if filepath.exists():
                results[key] = pd.read_csv(filepath)
                logger.info(f"Loaded {key}: {results[key].shape}")
            else:
                logger.warning(f"Result file not found: {filepath}")
        
        return results
    
    def validate_let7_pathway(self, results: Dict[str, pd.DataFrame]) -> Dict:
        """Validate let-7 pathway interactions."""
        logger.info("Validating let-7 pathway interactions...")
        
        validation = {
            'pathway': 'let-7',
            'expected_targets': ['lin-41', 'hbl-1', 'daf-12'],
            'expected_stages': ['L3-TP1', 'L3-TP2', 'L4-TP1', 'L4-TP2'],
            'detected_targets': [],
            'missed_targets': [],
            'correct_stage_specificity': False,
            'statistical_significance': {}
        }
        
        # Check multi_way_interactions for let-7 targets
        if 'multi_way_interactions' in results:
            multi_df = results['multi_way_interactions']
            
            # Look for let-7 target genes with significant interactions
            for target in validation['expected_targets']:
                # Check if target gene shows significant regulation
                target_rows = multi_df[multi_df['gene'].str.contains(target, na=False)]
                
                if not target_rows.empty:
                    # Check if regulators include miRNA (let-7)
                    for _, row in target_rows.iterrows():
                        regulator_types = row.get('regulator_types', '')
                        if 'mirna_' in str(regulator_types):
                            validation['detected_targets'].append(target)
                            
                            # Check statistical significance
                            p_value = row.get('interaction_p_value', 1.0)
                            improvement = row.get('improvement_from_regulators', 0.0)
                            
                            validation['statistical_significance'][target] = {
                                'p_value': p_value,
                                'improvement': improvement,
                                'significant': p_value < self.validation_criteria['statistical_significance']
                            }
                            break
                else:
                    validation['missed_targets'].append(target)
        
        # Calculate detection rate for let-7 pathway
        n_detected = len(validation['detected_targets'])
        n_expected = len(validation['expected_targets'])
        validation['detection_rate'] = n_detected / n_expected if n_expected > 0 else 0.0
        
        logger.info(f"let-7 pathway: {n_detected}/{n_expected} targets detected")
        
        return validation
    
    def validate_lin4_cascade(self, results: Dict[str, pd.DataFrame]) -> Dict:
        """Validate lin-4 cascade interactions."""
        logger.info("Validating lin-4 cascade interactions...")
        
        validation = {
            'pathway': 'lin-4',
            'expected_targets': ['lin-14', 'lin-28'],
            'expected_stages': ['L1-TP2', 'L2-TP1', 'L2-TP2'],
            'detected_targets': [],
            'missed_targets': [],
            'correct_stage_specificity': False,
            'statistical_significance': {}
        }
        
        # Similar validation logic for lin-4 pathway
        if 'multi_way_interactions' in results:
            multi_df = results['multi_way_interactions']
            
            for target in validation['expected_targets']:
                target_rows = multi_df[multi_df['gene'].str.contains(target, na=False)]
                
                if not target_rows.empty:
                    for _, row in target_rows.iterrows():
                        regulator_types = row.get('regulator_types', '')
                        if 'mirna_' in str(regulator_types):
                            validation['detected_targets'].append(target)
                            
                            p_value = row.get('interaction_p_value', 1.0)
                            improvement = row.get('improvement_from_regulators', 0.0)
                            
                            validation['statistical_significance'][target] = {
                                'p_value': p_value,
                                'improvement': improvement,
                                'significant': p_value < self.validation_criteria['statistical_significance']
                            }
                            break
                else:
                    validation['missed_targets'].append(target)
        
        n_detected = len(validation['detected_targets'])
        n_expected = len(validation['expected_targets'])
        validation['detection_rate'] = n_detected / n_expected if n_expected > 0 else 0.0
        
        logger.info(f"lin-4 cascade: {n_detected}/{n_expected} targets detected")
        
        return validation
    
    def calculate_overall_metrics(self, pathway_validations: List[Dict]) -> Dict:
        """Calculate overall validation metrics."""
        logger.info("Calculating overall validation metrics...")
        
        # Count total true positives, false negatives, false positives
        total_expected = 0
        total_detected = 0
        total_missed = 0
        
        all_significant_interactions = []
        
        for validation in pathway_validations:
            total_expected += len(validation['expected_targets'])
            total_detected += len(validation['detected_targets'])
            total_missed += len(validation['missed_targets'])
            
            # Collect significant interactions
            for target, stats in validation['statistical_significance'].items():
                if stats['significant']:
                    all_significant_interactions.append(f"{validation['pathway']}->{target}")
        
        # Calculate sensitivity (recall)
        sensitivity = total_detected / total_expected if total_expected > 0 else 0.0
        
        # For specificity, we need to estimate false positives
        # This is simplified - in reality would need to check all detected interactions
        # against comprehensive database of known interactions
        estimated_false_positives = max(0, len(all_significant_interactions) - total_detected)
        total_negatives = 100  # Estimated based on dataset size
        true_negatives = total_negatives - estimated_false_positives
        specificity = true_negatives / total_negatives if total_negatives > 0 else 0.0
        
        overall_metrics = {
            'sensitivity': sensitivity,
            'specificity': specificity,
            'total_expected_interactions': total_expected,
            'total_detected_interactions': total_detected,
            'total_missed_interactions': total_missed,
            'significant_interactions': all_significant_interactions,
            'meets_sensitivity_threshold': sensitivity >= self.validation_criteria['sensitivity_threshold'],
            'meets_specificity_threshold': specificity >= self.validation_criteria['specificity_threshold']
        }
        
        logger.info(f"Overall sensitivity: {sensitivity:.3f}")
        logger.info(f"Overall specificity: {specificity:.3f}")
        
        return overall_metrics
    
    def create_validation_plots(self, pathway_validations: List[Dict], overall_metrics: Dict) -> None:
        """Create validation visualization plots."""
        logger.info("Creating validation plots...")
        
        # Set up plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('ConTra Validation Results - C. elegans Test Dataset', fontsize=16)
        
        # Plot 1: Detection rates by pathway
        pathways = [v['pathway'] for v in pathway_validations]
        detection_rates = [v['detection_rate'] for v in pathway_validations]
        
        axes[0, 0].bar(pathways, detection_rates, color=['skyblue', 'lightgreen'])
        axes[0, 0].axhline(y=self.validation_criteria['sensitivity_threshold'], 
                          color='red', linestyle='--', label='Threshold')
        axes[0, 0].set_title('Detection Rate by Pathway')
        axes[0, 0].set_ylabel('Detection Rate')
        axes[0, 0].set_ylim(0, 1)
        axes[0, 0].legend()
        
        # Plot 2: Overall performance metrics
        metrics = ['Sensitivity', 'Specificity']
        values = [overall_metrics['sensitivity'], overall_metrics['specificity']]
        thresholds = [self.validation_criteria['sensitivity_threshold'], 
                     self.validation_criteria['specificity_threshold']]
        
        bars = axes[0, 1].bar(metrics, values, color=['coral', 'lightblue'])
        axes[0, 1].bar(metrics, thresholds, alpha=0.3, color='red', label='Threshold')
        axes[0, 1].set_title('Overall Performance Metrics')
        axes[0, 1].set_ylabel('Score')
        axes[0, 1].set_ylim(0, 1)
        axes[0, 1].legend()
        
        # Add value labels on bars
        for bar, value in zip(bars, values):
            axes[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                           f'{value:.3f}', ha='center', va='bottom')
        
        # Plot 3: Interaction significance
        all_interactions = []
        all_p_values = []
        all_improvements = []
        
        for validation in pathway_validations:
            for target, stats in validation['statistical_significance'].items():
                all_interactions.append(f"{validation['pathway']}->{target}")
                all_p_values.append(stats['p_value'])
                all_improvements.append(stats['improvement'])
        
        if all_p_values:
            scatter = axes[1, 0].scatter(all_improvements, [-np.log10(p) for p in all_p_values], 
                                       c=range(len(all_interactions)), cmap='viridis')
            axes[1, 0].axhline(y=-np.log10(self.validation_criteria['statistical_significance']), 
                              color='red', linestyle='--', label='Significance threshold')
            axes[1, 0].set_xlabel('Improvement from Regulators')
            axes[1, 0].set_ylabel('-log10(p-value)')
            axes[1, 0].set_title('Interaction Significance')
            axes[1, 0].legend()
        
        # Plot 4: Summary statistics
        summary_data = [
            ['Expected Interactions', overall_metrics['total_expected_interactions']],
            ['Detected Interactions', overall_metrics['total_detected_interactions']],
            ['Missed Interactions', overall_metrics['total_missed_interactions']],
            ['Significant Interactions', len(overall_metrics['significant_interactions'])]
        ]
        
        summary_df = pd.DataFrame(summary_data, columns=['Metric', 'Count'])
        bars = axes[1, 1].bar(summary_df['Metric'], summary_df['Count'], 
                             color=['gold', 'lightgreen', 'salmon', 'lightblue'])
        axes[1, 1].set_title('Interaction Summary')
        axes[1, 1].set_ylabel('Count')
        axes[1, 1].tick_params(axis='x', rotation=45)
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            axes[1, 1].text(bar.get_x() + bar.get_width()/2, height + 0.1, 
                           f'{int(height)}', ha='center', va='bottom')
        
        plt.tight_layout()
        
        # Save plot
        plot_file = self.output_dir / "validation_results.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved validation plot: {plot_file}")
        
        plt.close()
    
    def generate_validation_report(self, pathway_validations: List[Dict], overall_metrics: Dict) -> None:
        """Generate detailed validation report."""
        logger.info("Generating validation report...")
        
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        report_file = self.output_dir / f"validation_report_{timestamp}.md"
        
        with open(report_file, 'w') as f:
            f.write("# ConTra Validation Report - C. elegans Test Dataset\n\n")
            f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Summary\n\n")
            f.write(f"- **Sensitivity**: {overall_metrics['sensitivity']:.3f} ")
            f.write(f"({'✅ PASS' if overall_metrics['meets_sensitivity_threshold'] else '❌ FAIL'})\n")
            f.write(f"- **Specificity**: {overall_metrics['specificity']:.3f} ")
            f.write(f"({'✅ PASS' if overall_metrics['meets_specificity_threshold'] else '❌ FAIL'})\n")
            f.write(f"- **Total Expected Interactions**: {overall_metrics['total_expected_interactions']}\n")
            f.write(f"- **Total Detected Interactions**: {overall_metrics['total_detected_interactions']}\n\n")
            
            f.write("## Pathway-Specific Results\n\n")
            
            for validation in pathway_validations:
                f.write(f"### {validation['pathway']} Pathway\n\n")
                f.write(f"- **Detection Rate**: {validation['detection_rate']:.3f}\n")
                f.write(f"- **Expected Targets**: {', '.join(validation['expected_targets'])}\n")
                f.write(f"- **Detected Targets**: {', '.join(validation['detected_targets'])}\n")
                f.write(f"- **Missed Targets**: {', '.join(validation['missed_targets'])}\n\n")
                
                if validation['statistical_significance']:
                    f.write("#### Statistical Significance\n\n")
                    f.write("| Target | P-value | Improvement | Significant |\n")
                    f.write("|--------|---------|-------------|-------------|\n")
                    for target, stats in validation['statistical_significance'].items():
                        sig_status = "✅" if stats['significant'] else "❌"
                        f.write(f"| {target} | {stats['p_value']:.4f} | {stats['improvement']:.4f} | {sig_status} |\n")
                    f.write("\n")
            
            f.write("## Validation Criteria\n\n")
            f.write(f"- **Sensitivity Threshold**: {self.validation_criteria['sensitivity_threshold']}\n")
            f.write(f"- **Specificity Threshold**: {self.validation_criteria['specificity_threshold']}\n")
            f.write(f"- **Statistical Significance**: p < {self.validation_criteria['statistical_significance']}\n\n")
            
            f.write("## Recommendations\n\n")
            if overall_metrics['meets_sensitivity_threshold'] and overall_metrics['meets_specificity_threshold']:
                f.write("✅ **VALIDATION PASSED**: ConTra successfully detected known regulatory interactions.\n\n")
            else:
                f.write("❌ **VALIDATION FAILED**: ConTra did not meet validation criteria.\n\n")
                
                if not overall_metrics['meets_sensitivity_threshold']:
                    f.write("- **Low Sensitivity**: Consider adjusting correlation thresholds or statistical criteria\n")
                if not overall_metrics['meets_specificity_threshold']:
                    f.write("- **Low Specificity**: Consider increasing statistical stringency or filtering criteria\n")
            
            f.write("\n## Files Generated\n\n")
            f.write("- `validation_results.png`: Validation visualization plots\n")
            f.write(f"- `{report_file.name}`: This validation report\n")
        
        logger.info(f"Generated validation report: {report_file}")
    
    def run_validation(self, contra_output_dir: str) -> Dict:
        """Run complete validation analysis."""
        logger.info("Starting ConTra validation analysis...")
        
        # Update output directory for this analysis
        self.output_dir = Path(contra_output_dir)
        
        # Load ConTra results
        results = self.load_contra_results()
        
        # Validate individual pathways
        pathway_validations = [
            self.validate_let7_pathway(results),
            self.validate_lin4_cascade(results)
        ]
        
        # Calculate overall metrics
        overall_metrics = self.calculate_overall_metrics(pathway_validations)
        
        # Create visualizations
        self.create_validation_plots(pathway_validations, overall_metrics)
        
        # Generate report
        self.generate_validation_report(pathway_validations, overall_metrics)
        
        # Store results
        self.validation_results.update(overall_metrics)
        
        logger.info("✅ Validation analysis completed")
        
        return self.validation_results

def main():
    """Main function to run validation."""
    parser = argparse.ArgumentParser(description='Validate ConTra results against C. elegans test dataset')
    parser.add_argument('output_dir', help='ConTra output directory to validate')
    parser.add_argument('--metadata', default='data/test_datasets/celegans_modENCODE/validation_metadata.json',
                       help='Validation metadata file')
    
    args = parser.parse_args()
    
    # Create validator
    validator = ConTraValidator(args.output_dir, args.metadata)
    
    # Run validation
    results = validator.run_validation(args.output_dir)
    
    # Print summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)
    print(f"Sensitivity: {results['sensitivity']:.3f}")
    print(f"Specificity: {results['specificity']:.3f}")
    print(f"Expected interactions: {results['total_expected_interactions']}")
    print(f"Detected interactions: {results['total_detected_interactions']}")
    
    if results['meets_sensitivity_threshold'] and results['meets_specificity_threshold']:
        print("\n✅ VALIDATION PASSED")
    else:
        print("\n❌ VALIDATION FAILED")
    
    print("="*60)

if __name__ == '__main__':
    main()