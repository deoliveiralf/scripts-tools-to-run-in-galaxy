#!/usr/bin/env python

import argparse
import csv
import numpy as np
from collections import defaultdict

class SimpleLinearRegression:
    """A simple implementation of linear regression to replace sklearn's LinearRegression"""
    
    def __init__(self):
        self.slope = None
        self.intercept = None
    
    def fit(self, X, y):
        """Fit linear regression model"""
        # Convert to 1D arrays if needed
        X = X.reshape(-1)
        y = y.reshape(-1)
        
        # Calculate the mean of X and y
        x_mean = np.mean(X)
        y_mean = np.mean(y)
        
        # Calculate the slope
        numerator = np.sum((X - x_mean) * (y - y_mean))
        denominator = np.sum((X - x_mean) ** 2)
        
        if denominator == 0:
            self.slope = 0
        else:
            self.slope = numerator / denominator
        
        # Calculate the intercept
        self.intercept = y_mean - self.slope * x_mean
        
        return self
    
    def predict(self, X):
        """Predict using the linear model"""
        # Make sure X is in the right shape
        X = np.array(X).reshape(-1)
        
        # Make predictions
        y_pred = self.slope * X + self.intercept
        
        # Return in the same format as sklearn (2D array)
        return y_pred.reshape(-1, 1)

class PolyamineCalculator:
    """
    A class to calculate polyamine concentrations in both µg/g FW and nmol/g FW
    using the formula:
    Polyamine Content (μg/g FW) = [(Adjusted Peak Area - Intercept) ÷ Slope] × 20 × 40 × Extraction Volume (mL) × Molecular Weight (g/mol) ÷ 1000 ÷ Fresh Weight (g)
    Where:
    * Adjusted Peak Area = (Sample Peak Area ÷ Internal Standard Peak Area) × Reference Value
    """
    
    # Default molecular weights for common polyamines
    DEFAULT_MOLECULAR_WEIGHTS = {
        'Putrescine': 88.15,
        'Spermidine': 145.25,
        'Spermine': 202.34
    }
    
    def __init__(self, peak_area_source, internal_standard_source, initial_fw, volume, 
                 molecular_weight_source=None, preset_molecular_weight=True, 
                 standard_curve_source=None, group_mapping_source=None,
                 injection_volume=20, std_curve_conc=40, dilution_factor=40,
                 default_slope=45548.76, default_intercept=97847.61):
        """
        Initialize the PolyamineCalculator with input data
        
        Parameters:
        - peak_area_source: Path to CSV file with peak area data
        - internal_standard_source: Path to CSV file with internal standard data
        - initial_fw: Initial fresh weight in grams
        - volume: Extraction volume in milliliters (mL)
        - molecular_weight_source: Path to CSV file with custom molecular weights
        - preset_molecular_weight: Whether to use preset molecular weights
        - standard_curve_source: Path to CSV file with standard curve data
        - group_mapping_source: Path to CSV file with sample to group mapping
        - injection_volume: Injection volume in µL (default: 20)
        - std_curve_conc: Standard curve concentration in µL (default: 40)
        - dilution_factor: Dilution factor (default: 40)
        - default_slope: Default slope value for standard curve (default: 45,548.76)
        - default_intercept: Default intercept value for standard curve (default: 97,847.61)
        """
        self.initial_fw = initial_fw
        self.volume = volume / 1000.0 if volume > 100 else volume  # Convert µL to mL if needed
        self.injection_volume = injection_volume
        self.std_curve_conc = std_curve_conc
        self.dilution_factor = dilution_factor
        self.default_slope = default_slope
        self.default_intercept = default_intercept
        self.peak_areas = self._load_peak_areas(peak_area_source)
        self.internal_standards = self._load_internal_standards(internal_standard_source)
        
        # Load molecular weights
        if preset_molecular_weight:
            self.molecular_weights = self.DEFAULT_MOLECULAR_WEIGHTS
        elif molecular_weight_source:
            self.molecular_weights = self._load_molecular_weights(molecular_weight_source)
        else:
            self.molecular_weights = self.DEFAULT_MOLECULAR_WEIGHTS
        
        # Load standard curves if provided
        self.standard_curves = {}
        if standard_curve_source:
            self.standard_curves = self._load_standard_curves(standard_curve_source)
        
        # Load sample groups if provided
        self.sample_groups = {}
        if group_mapping_source:
            self.sample_groups = self._load_group_mapping(group_mapping_source)
    
    def _load_peak_areas(self, file_path):
        """Load peak areas from CSV file"""
        peak_areas = {}
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    sample_name = row.pop(reader.fieldnames[0])
                    peak_areas[sample_name] = {k: float(v) if v else 0 for k, v in row.items()}
            return peak_areas
        except Exception as e:
            print(f"Error loading peak areas: {e}")
            return {}
    
    def _load_internal_standards(self, file_path):
        """Load internal standards from CSV file"""
        internal_standards = {}
        try:
            with open(file_path, 'r') as f:
                reader = csv.reader(f)
                header = next(reader)
                for row in reader:
                    if row:  # Skip empty rows
                        sample_name = row[0]
                        # Assuming internal standard is in the second column
                        internal_standards[sample_name] = float(row[1]) if len(row) > 1 and row[1] else 0
            return internal_standards
        except Exception as e:
            print(f"Error loading internal standards: {e}")
            return {}
    
    def _load_molecular_weights(self, file_path):
        """Load molecular weights from CSV file"""
        molecular_weights = {}
        try:
            with open(file_path, 'r') as f:
                reader = csv.reader(f)
                for row in reader:
                    if len(row) >= 2:
                        try:
                            molecular_weights[row[0]] = float(row[1])
                        except ValueError:
                            pass  # Skip non-numeric values
            return molecular_weights
        except Exception as e:
            print(f"Error loading molecular weights: {e}")
            return {}
    
    def _load_standard_curves(self, file_path):
        """Load standard curves from CSV file"""
        standard_curves = {}
        try:
            with open(file_path, 'r') as f:
                reader = csv.reader(f)
                header = next(reader)
                
                # Find column indices for polyamine name, concentrations, and peak areas
                polyamine_idx = 0  # Default to first column
                conc_start_idx = None
                peak_area_start_idx = None
                
                # Try to find columns by name
                for i, col_name in enumerate(header):
                    col_lower = col_name.lower()
                    if 'polyamine' in col_lower or 'compound' in col_lower or 'name' in col_lower:
                        polyamine_idx = i
                    elif 'concentration' in col_lower and conc_start_idx is None:
                        conc_start_idx = i
                    elif 'peak area' in col_lower and peak_area_start_idx is None:
                        peak_area_start_idx = i
                
                # If couldn't find concentration or peak area columns, make assumptions
                if conc_start_idx is None:
                    conc_start_idx = 1  # Assume concentrations start from second column
                if peak_area_start_idx is None:
                    peak_area_start_idx = conc_start_idx + 1  # Assume peak areas follow concentrations
                
                # Process each row to extract standard curve data
                for row in reader:
                    if len(row) <= polyamine_idx:
                        continue
                    
                    polyamine = row[polyamine_idx]
                    
                    # Extract concentration and peak area values
                    concentrations = []
                    peak_areas = []
                    
                    # Try to get concentrations
                    for i in range(conc_start_idx, peak_area_start_idx):
                        if i < len(row) and row[i]:
                            try:
                                concentrations.append(float(row[i]))
                            except ValueError:
                                pass
                    
                    # Try to get peak areas
                    for i in range(peak_area_start_idx, len(row)):
                        if row[i]:
                            try:
                                peak_areas.append(float(row[i]))
                            except ValueError:
                                pass
                    
                    # If we have matching data points, save them
                    if len(concentrations) > 1 and len(concentrations) == len(peak_areas):
                        # Fit linear regression to get slope and intercept
                        model = SimpleLinearRegression()
                        # Note: We're now fitting peak areas vs concentrations differently
                        # Peak areas on X-axis, concentrations on Y-axis
                        X = np.array(peak_areas).reshape(-1, 1)
                        y = np.array(concentrations).reshape(-1, 1)
                        model.fit(X, y)
                        
                        standard_curves[polyamine] = {
                            'slope': model.slope,
                            'intercept': model.intercept
                        }
            
            return standard_curves
        except Exception as e:
            print(f"Error loading standard curves: {e}")
            return {}
    
    def _load_group_mapping(self, file_path):
        """Load sample to group mapping from CSV file"""
        group_mapping = {}
        try:
            with open(file_path, 'r') as f:
                reader = csv.reader(f)
                header = next(reader, None)  # Skip header row if present
                for row in reader:
                    if len(row) >= 2:
                        sample_name = row[0]
                        group_name = row[1]
                        group_mapping[sample_name] = group_name
            return group_mapping
        except Exception as e:
            print(f"Error loading group mapping: {e}")
            return {}
    
    def calculate_polyamine(self, specific_polyamine=None):
        """
        Calculate polyamine concentrations using the formula:
        Polyamine Content (μg/g FW) = [(Adjusted Peak Area - Intercept) ÷ Slope] × Injection Volume × Dilution Factor × Extraction Volume (mL) × Molecular Weight (g/mol) ÷ 1000 ÷ Fresh Weight (g)
        
        Parameters:
        - specific_polyamine: If provided, only calculate for this polyamine
        
        Returns:
        - Tuple of dictionaries (results_ug, results_nmol)
        """
        # Find the maximum internal standard value as the reference value
        internal_std_values = list(self.internal_standards.values())
        reference_value = max(internal_std_values) if internal_std_values else 1.0
        
        # Process each polyamine
        results_ug = {}
        results_nmol = {}
        
        # Get list of polyamines to process
        polyamines = []
        if specific_polyamine:
            polyamines = [specific_polyamine]
        else:
            # Get all polyamines from peak areas
            if self.peak_areas:
                first_sample = next(iter(self.peak_areas.values()))
                polyamines = list(first_sample.keys())
        
        for polyamine in polyamines:
            # Skip if polyamine not in molecular weights
            if polyamine not in self.molecular_weights:
                print(f"Warning: No molecular weight found for {polyamine}")
                continue
            
            molecular_weight = self.molecular_weights[polyamine]
            
            # Get standard curve data for this polyamine
            curve_data = self.standard_curves.get(polyamine, None)
            
            # Use default slope and intercept if standard curve not available
            slope = curve_data.get('slope', self.default_slope) if curve_data else self.default_slope
            intercept = curve_data.get('intercept', self.default_intercept) if curve_data else self.default_intercept
            
            # Calculate for each sample
            ug_values = []
            nmol_values = []
            sample_results_ug = {}
            sample_results_nmol = {}
            
            for sample, peak_area_dict in self.peak_areas.items():
                if polyamine not in peak_area_dict:
                    continue
                
                peak_area = peak_area_dict[polyamine]
                internal_standard_value = self.internal_standards.get(sample, 1.0)
                
                # Skip calculation if internal standard is zero to avoid division by zero
                if internal_standard_value == 0:
                    continue
                
                # Calculate using the formula:
                # Polyamine Content (μg/g FW) = [(Adjusted Peak Area - Intercept) ÷ Slope] × Injection Volume × Dilution Factor × Extraction Volume (mL) × Molecular Weight (g/mol) ÷ 1000 ÷ Fresh Weight (g)
                
                # 1. Calculate the adjusted peak area
                adjusted_peak_area = (peak_area / internal_standard_value) * reference_value
                
                # 2. Calculate concentration using (Adjusted Peak Area - Intercept) / Slope
                concentration = (adjusted_peak_area - intercept) / slope
		# Add this line to prevent negative concentrations
                concentration = max(0, concentration)
                # 3. Calculate the injection factor (injection volume / std curve concentration)
                injection_factor = self.injection_volume / self.std_curve_conc
                
                # 4. Calculate polyamine content in μg/g FW
                ug_per_g_fw = concentration * injection_factor * self.dilution_factor * self.volume * molecular_weight / 1000 / self.initial_fw
                
                # 5. For nmol/g FW, convert using molecular weight
                nmol_per_g_fw = ug_per_g_fw / molecular_weight * 1000
                
                ug_values.append(ug_per_g_fw)
                nmol_values.append(nmol_per_g_fw)
                
                # Store individual sample results
                sample_results_ug[sample] = ug_per_g_fw
                sample_results_nmol[sample] = nmol_per_g_fw
            
            # Calculate statistics
            if ug_values:
                avg_ug = np.mean(ug_values)
                std_dev_ug = np.std(ug_values) if len(ug_values) > 1 else 0
                cv_percent_ug = (std_dev_ug / avg_ug * 100) if avg_ug != 0 else 0
                
                avg_nmol = np.mean(nmol_values)
                std_dev_nmol = np.std(nmol_values) if len(nmol_values) > 1 else 0
                cv_percent_nmol = (std_dev_nmol / avg_nmol * 100) if avg_nmol != 0 else 0
                
                results_ug[polyamine] = {
                    'average': avg_ug,
                    'std_dev': std_dev_ug,
                    'cv_percent': cv_percent_ug,
                    'values': ug_values,
                    'sample_values': sample_results_ug
                }
                
                results_nmol[polyamine] = {
                    'average': avg_nmol,
                    'std_dev': std_dev_nmol,
                    'cv_percent': cv_percent_nmol,
                    'values': nmol_values,
                    'sample_values': sample_results_nmol
                }
        
        return results_ug, results_nmol
    
    def calculate_group_statistics(self, results):
        """
        Calculate group-level statistics for polyamine results
        
        Parameters:
        - results: Dictionary of polyamine results (either µg/g FW or nmol/g FW)
        
        Returns:
        - Dictionary of group statistics
        """
        # If no sample groups or no results, return empty
        if not self.sample_groups or not results:
            return {}
        
        # Initialize group data structures
        group_data = defaultdict(lambda: defaultdict(list))
        
        # Iterate through each polyamine
        for polyamine, data in results.items():
            if not data or 'sample_values' not in data:
                continue
                
            # Process each sample's value
            for sample, value in data['sample_values'].items():
                # Get group for this sample
                group = self.sample_groups.get(sample)
                if not group:
                    continue
                
                # Add value to group data
                group_data[group][polyamine].append(value)
        
        # Calculate group-level statistics
        group_statistics = {}
        for group, polyamine_data in group_data.items():
            group_statistics[group] = {}
            
            for polyamine, values in polyamine_data.items():
                # Skip if no valid values
                if not values:
                    continue
                
                # Calculate statistics
                mean = np.mean(values)
                std_dev = np.std(values) if len(values) > 1 else 0
                cv_percent = (std_dev / mean * 100) if mean != 0 else 0
                
                group_statistics[group][polyamine] = {
                    'mean': mean,
                    'std_dev': std_dev,
                    'cv_percent': cv_percent,
                    'sample_count': len(values),
                    'sample_values': values  # Store individual values for reference
                }
        
        return group_statistics

def main():
    parser = argparse.ArgumentParser(description='Polyamine Concentration Calculator')
    parser.add_argument('--mol_weight', help='Custom molecular weight file')
    parser.add_argument('--preset_mol_weight', action='store_true', help='Use preset molecular weights')
    parser.add_argument('--standard_curve', help='Standard curve file')
    parser.add_argument('--peak_areas', required=True, help='Peak areas file')
    parser.add_argument('--internal_standard', required=True, help='Internal standard file')
    parser.add_argument('--initial_fw', type=float, required=True, help='Initial fresh weight (g)')
    parser.add_argument('--volume', type=float, required=True, help='Extraction volume (mL or µL)')
    parser.add_argument('--injection_volume', type=float, default=20, help='Injection volume (µL)')
    parser.add_argument('--std_curve_conc', type=float, default=40, help='Standard curve concentration (µL)')
    parser.add_argument('--dilution_factor', type=float, default=40, help='Dilution factor')
    parser.add_argument('--default_slope', type=float, default=45548.76, help='Default slope value')
    parser.add_argument('--default_intercept', type=float, default=97847.61, help='Default intercept value')
    parser.add_argument('--output_ug', required=True, help='Output file for µg/g FW results')
    parser.add_argument('--output_nmol', required=True, help='Output file for nmol/g FW results')
    parser.add_argument('--output_samples_ug', help='Output file for individual sample µg/g FW results')
    parser.add_argument('--output_samples_nmol', help='Output file for individual sample nmol/g FW results')
    parser.add_argument('--polyamine', help='Specific polyamine to calculate')
    parser.add_argument('--group_mapping', help='Group mapping file')
    parser.add_argument('--group_output_ug', help='Group statistics output file for µg/g FW')
    parser.add_argument('--group_output_nmol', help='Group statistics output file for nmol/g FW')
    parser.add_argument('--group_samples_output_ug', help='Group individual samples output file for µg/g FW')
    parser.add_argument('--group_samples_output_nmol', help='Group individual samples output file for nmol/g FW')
    
    args = parser.parse_args()
    
    # Initialize calculator
    calculator = PolyamineCalculator(
        molecular_weight_source=args.mol_weight,
        preset_molecular_weight=args.preset_mol_weight,
        standard_curve_source=args.standard_curve,
        peak_area_source=args.peak_areas,
        internal_standard_source=args.internal_standard,
        initial_fw=args.initial_fw,
        volume=args.volume,
        injection_volume=args.injection_volume,
        std_curve_conc=args.std_curve_conc,
        dilution_factor=args.dilution_factor,
        default_slope=args.default_slope,
        default_intercept=args.default_intercept,
        group_mapping_source=args.group_mapping
    )
    
    # Calculate polyamine concentrations
    results_ug, results_nmol = calculator.calculate_polyamine(args.polyamine)
    
    # Write µg/g FW results
    if args.output_ug:
        with open(args.output_ug, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Polyamine', 'Average (µg/g FW)', 'Std Dev', 'CV%'])
            for polyamine, data in results_ug.items():
                writer.writerow([
                    polyamine, 
                    data['average'], 
                    data['std_dev'], 
                    data['cv_percent']
                ])
    
    # Write nmol/g FW results
    if args.output_nmol:
        with open(args.output_nmol, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Polyamine', 'Average (nmol/g FW)', 'Std Dev', 'CV%'])
            for polyamine, data in results_nmol.items():
                writer.writerow([
                    polyamine, 
                    data['average'], 
                    data['std_dev'], 
                    data['cv_percent']
                ])
    
    # Write individual sample µg/g FW results
    if args.output_samples_ug:
        with open(args.output_samples_ug, 'w', newline='') as f:
            # Get all sample names and polyamines
            all_samples = set()
            all_polyamines = set()
            for polyamine, data in results_ug.items():
                all_polyamines.add(polyamine)
                all_samples.update(data['sample_values'].keys())
            
            # Create header row with all polyamines
            writer = csv.writer(f)
            header = ['Sample']
            header.extend(sorted(all_polyamines))
            writer.writerow(header)
            
            # Write data for each sample
            for sample in sorted(all_samples):
                row = [sample]
                for polyamine in sorted(all_polyamines):
                    value = results_ug.get(polyamine, {}).get('sample_values', {}).get(sample, '')
                    row.append(value)
                writer.writerow(row)
    
    # Write individual sample nmol/g FW results
    if args.output_samples_nmol:
        with open(args.output_samples_nmol, 'w', newline='') as f:
            # Get all sample names and polyamines
            all_samples = set()
            all_polyamines = set()
            for polyamine, data in results_nmol.items():
                all_polyamines.add(polyamine)
                all_samples.update(data['sample_values'].keys())
            
            # Create header row with all polyamines
            writer = csv.writer(f)
            header = ['Sample']
            header.extend(sorted(all_polyamines))
            writer.writerow(header)
            
            # Write data for each sample
            for sample in sorted(all_samples):
                row = [sample]
                for polyamine in sorted(all_polyamines):
                    value = results_nmol.get(polyamine, {}).get('sample_values', {}).get(sample, '')
                    row.append(value)
                writer.writerow(row)
    
    # Calculate and write group statistics if group mapping is provided
    if args.group_mapping:
        # Calculate group statistics for µg/g FW
        group_stats_ug = calculator.calculate_group_statistics(results_ug)
        
        # Write group summary statistics for µg/g FW
        if args.group_output_ug:
            with open(args.group_output_ug, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Group', 'Polyamine', 'Mean (µg/g FW)', 'Std Dev', 'CV%', 'Sample Count'])
                for group, polyamine_data in group_stats_ug.items():
                    for polyamine, stats in polyamine_data.items():
                        writer.writerow([
                            group, 
                            polyamine, 
                            stats['mean'], 
                            stats['std_dev'], 
                            stats['cv_percent'], 
                            stats['sample_count']
                        ])
        
        # Write individual sample values by group for µg/g FW
        if args.group_samples_output_ug:
            with open(args.group_samples_output_ug, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Group', 'Polyamine', 'Sample Number', 'Value (µg/g FW)'])
                for group, polyamine_data in group_stats_ug.items():
                    for polyamine, stats in polyamine_data.items():
                        if 'sample_values' in stats:
                            for i, value in enumerate(stats['sample_values']):
                                writer.writerow([group, polyamine, i+1, value])
        
        # Calculate group statistics for nmol/g FW
        group_stats_nmol = calculator.calculate_group_statistics(results_nmol)
        
        # Write group summary statistics for nmol/g FW
        if args.group_output_nmol:
            with open(args.group_output_nmol, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Group', 'Polyamine', 'Mean (nmol/g FW)', 'Std Dev', 'CV%', 'Sample Count'])
                for group, polyamine_data in group_stats_nmol.items():
                    for polyamine, stats in polyamine_data.items():
                        writer.writerow([
                            group, 
                            polyamine, 
                            stats['mean'], 
                            stats['std_dev'], 
                            stats['cv_percent'], 
                            stats['sample_count']
                        ])
        
        # Write individual sample values by group for nmol/g FW
        if args.group_samples_output_nmol:
            with open(args.group_samples_output_nmol, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Group', 'Polyamine', 'Sample Number', 'Value (nmol/g FW)'])
                for group, polyamine_data in group_stats_nmol.items():
                    for polyamine, stats in polyamine_data.items():
                        if 'sample_values' in stats:
                            for i, value in enumerate(stats['sample_values']):
                                writer.writerow([group, polyamine, i+1, value])

if __name__ == '__main__':
    main()
