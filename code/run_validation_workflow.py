#!/usr/bin/env python3
"""
Complete ConTra validation workflow using C. elegans test dataset.

This script automates the entire process of:
1. Creating the test dataset
2. Running ConTra analysis  
3. Validating results against known interactions
4. Restoring original data

Usage: python3 run_validation_workflow.py [--keep-test-data]
"""

import os
import sys
import subprocess
import shutil
import argparse
import logging
from pathlib import Path
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ValidationWorkflow:
    """Complete ConTra validation workflow."""
    
    def __init__(self, keep_test_data: bool = False):
        """Initialize workflow."""
        self.keep_test_data = keep_test_data
        self.base_dir = Path(__file__).parent.parent
        self.data_backup_dir = None
        self.output_dir = None
        
    def backup_original_data(self) -> bool:
        """Backup original cleaned datasets."""
        logger.info("Backing up original data...")
        
        original_data = self.base_dir / "data" / "cleaned_datasets"
        self.data_backup_dir = self.base_dir / "data" / "cleaned_datasets_backup"
        
        if not original_data.exists():
            logger.error(f"Original data directory not found: {original_data}")
            return False
            
        if self.data_backup_dir.exists():
            shutil.rmtree(self.data_backup_dir)
            
        shutil.copytree(original_data, self.data_backup_dir)
        logger.info(f"Original data backed up to: {self.data_backup_dir}")
        return True
    
    def create_test_dataset(self) -> bool:
        """Create C. elegans test dataset."""
        logger.info("Creating C. elegans test dataset...")
        
        script_path = self.base_dir / "code" / "create_celegans_test_dataset.py"
        
        try:
            result = subprocess.run([
                sys.executable, str(script_path)
            ], cwd=self.base_dir, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"Test dataset creation failed: {result.stderr}")
                return False
                
            logger.info("Test dataset created successfully")
            return True
            
        except Exception as e:
            logger.error(f"Error creating test dataset: {e}")
            return False
    
    def install_test_data(self) -> bool:
        """Install test data for analysis."""
        logger.info("Installing test data...")
        
        test_data_dir = self.base_dir / "data" / "test_datasets" / "celegans_modENCODE"
        cleaned_data_dir = self.base_dir / "data" / "cleaned_datasets"
        
        if not test_data_dir.exists():
            logger.error(f"Test data directory not found: {test_data_dir}")
            return False
        
        # Copy CSV files
        csv_files = list(test_data_dir.glob("*.csv"))
        
        for csv_file in csv_files:
            dest_file = cleaned_data_dir / csv_file.name
            shutil.copy2(csv_file, dest_file)
            logger.info(f"Copied {csv_file.name}")
        
        logger.info(f"Installed {len(csv_files)} test data files")
        return True
    
    def run_contra_analysis(self) -> bool:
        """Run ConTra analysis on test dataset."""
        logger.info("Running ConTra analysis...")
        
        script_path = self.base_dir / "code" / "subset_context_dependent_analysis.py"
        
        try:
            result = subprocess.run([
                sys.executable, str(script_path)
            ], cwd=self.base_dir / "code", capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"ConTra analysis failed: {result.stderr}")
                return False
            
            # Extract output directory from stdout
            # Look for the pattern in the first few lines
            lines = result.stdout.split('\n')[:10]  # Check first 10 lines only
            
            for line in lines:
                if 'Output directory:' in line:
                    self.output_dir = line.split('Output directory:')[1].strip()
                    # Convert to relative path if needed
                    if self.output_dir.startswith('/home/runner/work/ConTra/ConTra/'):
                        self.output_dir = self.output_dir.replace('/home/runner/work/ConTra/ConTra/', '')
                    logger.info(f"ConTra analysis completed: {self.output_dir}")
                    return True
            
            # Fallback: find the most recent output directory
            output_base = self.base_dir / "output"
            if output_base.exists():
                output_dirs = [d for d in output_base.iterdir() 
                              if d.is_dir() and d.name.startswith('subset_context_dependent_analysis_')]
                if output_dirs:
                    # Get the most recent one
                    latest_dir = max(output_dirs, key=lambda x: x.stat().st_mtime)
                    self.output_dir = str(latest_dir.relative_to(self.base_dir))
                    logger.info(f"Found output directory: {self.output_dir}")
                    return True
            
            logger.error("Could not determine output directory")
            return False
                
        except Exception as e:
            logger.error(f"Error running ConTra analysis: {e}")
            return False
    
    def run_validation(self) -> bool:
        """Run validation analysis."""
        logger.info("Running validation analysis...")
        
        if not self.output_dir:
            logger.error("No output directory specified for validation")
            return False
        
        script_path = self.base_dir / "code" / "validate_celegans_results.py"
        
        try:
            result = subprocess.run([
                sys.executable, str(script_path), self.output_dir
            ], cwd=self.base_dir, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"Validation failed: {result.stderr}")
                return False
            
            # Print validation summary
            if "VALIDATION SUMMARY" in result.stdout:
                summary_start = result.stdout.find("VALIDATION SUMMARY")
                summary_end = result.stdout.find("="*60, summary_start + 20)
                summary = result.stdout[summary_start:summary_end + 60]
                print("\\n" + summary)
            
            logger.info("Validation completed successfully")
            return True
            
        except Exception as e:
            logger.error(f"Error running validation: {e}")
            return False
    
    def restore_original_data(self) -> bool:
        """Restore original data."""
        if self.keep_test_data:
            logger.info("Keeping test data (--keep-test-data specified)")
            return True
            
        logger.info("Restoring original data...")
        
        if not self.data_backup_dir or not self.data_backup_dir.exists():
            logger.error("No backup data found to restore")
            return False
        
        cleaned_data_dir = self.base_dir / "data" / "cleaned_datasets"
        
        if cleaned_data_dir.exists():
            shutil.rmtree(cleaned_data_dir)
        
        shutil.move(str(self.data_backup_dir), str(cleaned_data_dir))
        logger.info("Original data restored")
        return True
    
    def cleanup_backup(self) -> None:
        """Clean up backup directory if it still exists."""
        if self.data_backup_dir and self.data_backup_dir.exists():
            shutil.rmtree(self.data_backup_dir)
            logger.info("Backup directory cleaned up")
    
    def run_complete_workflow(self) -> bool:
        """Run the complete validation workflow."""
        logger.info("Starting complete ConTra validation workflow...")
        
        try:
            # Step 1: Backup original data
            if not self.backup_original_data():
                return False
            
            # Step 2: Create test dataset
            if not self.create_test_dataset():
                return False
            
            # Step 3: Install test data
            if not self.install_test_data():
                return False
            
            # Step 4: Run ConTra analysis
            if not self.run_contra_analysis():
                return False
            
            # Step 5: Run validation
            if not self.run_validation():
                return False
            
            logger.info("✅ Complete validation workflow succeeded!")
            
            if self.output_dir:
                logger.info(f"📊 Results available in: {self.output_dir}")
                logger.info(f"📈 Validation plots: {self.output_dir}/validation_results.png")
                logger.info(f"📄 Validation report: {self.output_dir}/validation_report_*.md")
            
            return True
            
        finally:
            # Always try to restore original data
            self.restore_original_data()
            self.cleanup_backup()

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Run complete ConTra validation workflow')
    parser.add_argument('--keep-test-data', action='store_true',
                       help='Keep test data instead of restoring original data')
    
    args = parser.parse_args()
    
    workflow = ValidationWorkflow(keep_test_data=args.keep_test_data)
    
    success = workflow.run_complete_workflow()
    
    if success:
        print("\\n" + "="*60)
        print("🎉 VALIDATION WORKFLOW COMPLETED SUCCESSFULLY!")
        print("="*60)
        print("ConTra pipeline validation passed using C. elegans test dataset.")
        print("The pipeline successfully detected known regulatory interactions.")
        sys.exit(0)
    else:
        print("\\n" + "="*60)
        print("❌ VALIDATION WORKFLOW FAILED")
        print("="*60)
        print("Check the log messages above for specific error details.")
        sys.exit(1)

if __name__ == '__main__':
    main()