import argparse
import csv
import os
import numpy as np
from collections import defaultdict

class AminoAcidCalculator:
    def __init__(self, 
                 molecular_weight_source=None,
                 preset_molecular_weight=None, 
                 standard_curve_source=None, 
                 peak_area_source=None,
                 initial_fw=None, 
                 volume=None,
                 group_mapping_source=None):
        """
        Initialize the Amino Acid Calculator with flexible input sources
        
        Parameters:
        - molecular_weight_source: File path for custom molecular weights
        - preset_molecular_weight: Flag to use preset molecular weights
        - standard_curve_source: File path or equation for standard curve
        - peak_area_source: File path containing peak areas (csv, tsv)
        - initial_fw: Initial fresh weight in grams
        - volume: Volume in milliliters
        - group_mapping_source: File path containing sample to group mappings
        """
        # Updated preset molecular weights based on provided file
        self.preset_molecular_weights = {
            'Aspartic Acid': 133.1,
            'Glutamic Acid': 147.1,
            'Asparagine': 132.1,
            'Serine': 105.1,
            'Glutamine': 146.1,
            'Histidine': 191.7,
            'Glycine': 75.1,
            'Arginine': 210.7,
            'Citruline': 175.2,
            'Threonine': 119.1,
            'Alanine': 89.1,
            'GABA': 103.1,
            'Tyrosine': 181.2,
            'Tryptophan': 204.2,
            'Methionine': 149.2,
            'Valine': 117.2,
            'Phenylalanine': 165.2,
            'Isoleucine': 131.2,
            'Leucine': 131.2,
            'Ornithine': 168.6,
            'Lysine': 182.7
        }
        
        # Load molecular weight from source or use preset
        if molecular_weight_source:
            self.molecular_weight = self._load_molecular_weight(molecular_weight_source)
        elif preset_molecular_weight:
            self.molecular_weight = self.preset_molecular_weights
        else:
            self.molecular_weight = self.preset_molecular_weights
            
        # Load standard curves
        self.standard_curves = self._load_standard_curves(standard_curve_source)
        
        # Load peak areas (now keeping all replicates)
        self.peak_areas, self.sample_names = self._load_peak_areas_with_replicates(peak_area_source)
        
        # Load group mappings if provided
        self.sample_groups = self._load_sample_groups(group_mapping_source)
        
        self.initial_fw = initial_fw
        self.volume = volume
    
    def _load_file_content(self, source, delimiter=None):
        """
        Load file content from CSV or TSV with automatic delimiter detection
        """
        try:
            if delimiter is None:
                # Auto-detect delimiter based on file extension
                if source.endswith('.tsv') or source.endswith('.tabular') or source.endswith('.txt'):
                    delimiter = '\t'
                else:
                    delimiter = ','
            
            with open(source, 'r') as f:
                reader = csv.reader(f, delimiter=delimiter)
                return list(reader)
            
        except Exception as e:
            print(f"Error loading file {source}: {e}")
            return None
    
    def _load_sample_groups(self, source):
        """
        Load sample to group mappings from a file
        
        Expected format:
        sample_name,group_name
        sample1,control
        sample2,control
        sample3,treatment1
        """
        if source is None:
            return None
        
        # Load file content
        rows = self._load_file_content(source)
        
        if not rows or len(rows) < 2:
            return None
        
        # Parse sample to group mappings
        sample_groups = {}
        for row in rows[1:]:  # Skip header row
            if len(row) >= 2:
                sample_name = row[0].strip()
                group_name = row[1].strip()
                sample_groups[sample_name] = group_name
        
        # Validate that all samples have group assignments
        if self.sample_names:
            for sample in self.sample_names:
                if sample not in sample_groups:
                    print(f"Warning: Sample '{sample}' has no group assignment")
        
        return sample_groups
    
    def _load_molecular_weight(self, source):
        """
        Load molecular weight from a file
        """
        if source is None:
            return self.preset_molecular_weights
        
        # Load file content with auto-detected delimiter
        rows = self._load_file_content(source)
        
        if not rows:
            return self.preset_molecular_weights
        
        # Parse molecular weights
        mol_weights = {}
        for row in rows:
            if len(row) >= 2:
                try:
                    amino_acid = row[0].strip()
                    weight = float(row[1].strip())
                    mol_weights[amino_acid] = weight
                except (ValueError, IndexError):
                    continue
        
        return mol_weights if mol_weights else self.preset_molecular_weights
    
    def _calculate_standard_curve(self, x_values, y_values):
        """
        Calculate linear regression for standard curve using numpy
        
        Returns slope, intercept, and r_squared
        """
        if len(x_values) < 2 or len(y_values) < 2:
            return None
        
        # Convert to numpy arrays
        x = np.array(x_values, dtype=float)
        y = np.array(y_values, dtype=float)
        
        # Calculate linear regression
        # Basic formula for linear regression:
        # slope = sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
        # intercept = mean(y) - slope * mean(x)
        n = len(x)
        mean_x = np.mean(x)
        mean_y = np.mean(y)
        
        # Calculate slope
        numerator = np.sum((x - mean_x) * (y - mean_y))
        denominator = np.sum((x - mean_x) ** 2)
        
        if denominator == 0:
            return None
            
        slope = numerator / denominator
        
        # Calculate intercept
        intercept = mean_y - slope * mean_x
        
        # Calculate r_squared
        # r_squared = 1 - sum((y - predicted_y)^2) / sum((y - mean(y))^2)
        y_pred = slope * x + intercept
        ss_total = np.sum((y - mean_y) ** 2)
        ss_residual = np.sum((y - y_pred) ** 2)
        
        if ss_total == 0:
            r_squared = 0
        else:
            r_squared = 1 - (ss_residual / ss_total)
        
        return {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_squared
        }
    
    def _load_standard_curves(self, source):
        """
        Load standard curves for all amino acids from file
        """
        # Default is a simple multiplier of 1.0
        default_curve = {'slope': 1.0, 'intercept': 0.0, 'r_squared': 1.0}
        
        if source is None:
            return {'default': default_curve}
        
        # Load file content
        rows = self._load_file_content(source)
        
        if not rows or len(rows) < 2:
            return {'default': default_curve}
        
        # First row should contain headers (amino acid names)
        headers = rows[0]
        
        # First column should contain concentration values
        concentrations = []
        for row in rows[1:]:
            try:
                concentrations.append(float(row[0]))
            except (ValueError, IndexError):
                continue
        
        # Calculate standard curves for each amino acid
        standard_curves = {}
        for col in range(1, len(headers)):
            if col >= len(headers):
                continue
            
            amino_acid = headers[col].strip()
            
            # Extract peak areas for this amino acid
            peak_areas = []
            for row in rows[1:]:
                if col < len(row):
                    try:
                        peak_areas.append(float(row[col]))
                    except (ValueError, IndexError):
                        peak_areas.append(None)
                else:
                    peak_areas.append(None)
            
            # Filter out None values
            valid_data = [(c, p) for c, p in zip(concentrations, peak_areas) if c is not None and p is not None]
            
            if valid_data:
                x_values = [point[0] for point in valid_data]  # concentrations
                y_values = [point[1] for point in valid_data]  # peak areas
                
                # Calculate standard curve
                curve = self._calculate_standard_curve(x_values, y_values)
                
                if curve:
                    standard_curves[amino_acid] = curve
        
        # If no valid curves were found, use default
        if not standard_curves:
            return {'default': default_curve}
        
        return standard_curves
    
    def _load_peak_areas_with_replicates(self, source):
        """
        Load peak areas from file, preserving replicates
        Returns a dictionary with amino acids as keys, and lists of peak areas as values,
        and a list of sample names
        """
        if source is None:
            return None, None
        
        # Load file content
        rows = self._load_file_content(source)
        
        if not rows:
            return None, None
        
        # Extract headers (first row)
        headers = rows[0]
        
        # Extract sample names (first column, excluding the header)
        sample_names = []
        for row in rows[1:]:
            if row and len(row) > 0:
                sample_names.append(row[0])
        
        # Initialize peak areas dictionary
        # Format: {'Amino Acid': {'raw_values': [val1, val2, ...], 'average': avg_value}}
        peak_areas = {}
        
        # Process each amino acid column
        for col in range(1, len(headers)):
            amino_acid = headers[col].strip()
            raw_values = []
            
            for row in rows[1:]:
                if col < len(row):
                    try:
                        value = float(row[col])
                        raw_values.append(value)
                    except (ValueError, TypeError):
                        raw_values.append(None)  # Add None for invalid values
                else:
                    raw_values.append(None)  # Add None for missing values
            
            if raw_values and any(v is not None for v in raw_values):
                valid_values = [v for v in raw_values if v is not None]
                avg_value = sum(valid_values) / len(valid_values) if valid_values else None
                
                peak_areas[amino_acid] = {
                    'raw_values': raw_values,
                    'average': avg_value
                }
        
        return peak_areas, sample_names
    
    def calculate_amino_acid(self, amino_acid_name=None):
        """
        Calculate amino acid concentration in µg/g FW and nmol/g FW
        
        Parameters:
        - amino_acid_name: Optional specific amino acid to calculate
        
        Returns:
        - Tuple with two dictionaries: one for µg/g FW and one for nmol/g FW,
          each containing amino acids and their calculated concentrations
        """
        # Validate inputs
        if self.peak_areas is None or self.initial_fw is None or self.volume is None:
            raise ValueError("Missing required input parameters")
        
        # If no specific amino acid, calculate for all
        if amino_acid_name is None or amino_acid_name == "":
            results_ug = {}
            results_nmol = {}
            for aa, values in self.peak_areas.items():
                try:
                    result = self._calculate_single_amino_acid_both_units(aa, values['average'])
                    if result is None:
                        continue
                        
                    avg_result_ug, avg_result_nmol = result
                    
                    rep_results_ug = []
                    rep_results_nmol = []
                    
                    # Calculate for each replicate
                    for rep_value in values['raw_values']:
                        if rep_value is not None:
                            rep_result = self._calculate_single_amino_acid_both_units(aa, rep_value)
                            if rep_result is not None:
                                rep_results_ug.append(rep_result[0])
                                rep_results_nmol.append(rep_result[1])
                            else:
                                rep_results_ug.append(None)
                                rep_results_nmol.append(None)
                        else:
                            rep_results_ug.append(None)
                            rep_results_nmol.append(None)
                    
                    # Calculate standard deviation for µg/g FW
                    valid_rep_results_ug = [r for r in rep_results_ug if r is not None]
                    std_dev_ug = np.std(valid_rep_results_ug) if len(valid_rep_results_ug) > 1 else None
                    cv_percent_ug = None
                    if std_dev_ug is not None and avg_result_ug is not None and avg_result_ug != 0:
                        cv_percent_ug = (std_dev_ug / avg_result_ug) * 100
                    
                    # Calculate standard deviation for nmol/g FW
                    valid_rep_results_nmol = [r for r in rep_results_nmol if r is not None]
                    std_dev_nmol = np.std(valid_rep_results_nmol) if len(valid_rep_results_nmol) > 1 else None
                    cv_percent_nmol = None
                    if std_dev_nmol is not None and avg_result_nmol is not None and avg_result_nmol != 0:
                        cv_percent_nmol = (std_dev_nmol / avg_result_nmol) * 100
                    
                    results_ug[aa] = {
                        'average': avg_result_ug,
                        'std_dev': std_dev_ug,
                        'cv_percent': cv_percent_ug,
                        'replicates': rep_results_ug
                    }
                    
                    results_nmol[aa] = {
                        'average': avg_result_nmol,
                        'std_dev': std_dev_nmol,
                        'cv_percent': cv_percent_nmol,
                        'replicates': rep_results_nmol
                    }
                except ValueError as e:
                    print(f"Warning: {e}")
                    continue
            return results_ug, results_nmol
        
        # Calculate for specific amino acid
        if amino_acid_name not in self.peak_areas:
            raise ValueError(f"No peak area found for {amino_acid_name}")
        
        values = self.peak_areas[amino_acid_name]
        result = self._calculate_single_amino_acid_both_units(amino_acid_name, values['average'])
        if result is None:
            return {amino_acid_name: None}, {amino_acid_name: None}
            
        avg_result_ug, avg_result_nmol = result
        
        rep_results_ug = []
        rep_results_nmol = []
        
        # Calculate for each replicate
        for rep_value in values['raw_values']:
            if rep_value is not None:
                rep_result = self._calculate_single_amino_acid_both_units(amino_acid_name, rep_value)
                if rep_result is not None:
                    rep_results_ug.append(rep_result[0])
                    rep_results_nmol.append(rep_result[1])
                else:
                    rep_results_ug.append(None)
                    rep_results_nmol.append(None)
            else:
                rep_results_ug.append(None)
                rep_results_nmol.append(None)
        
        # Calculate standard deviation for µg/g FW
        valid_rep_results_ug = [r for r in rep_results_ug if r is not None]
        std_dev_ug = np.std(valid_rep_results_ug) if len(valid_rep_results_ug) > 1 else None
        cv_percent_ug = None
        if std_dev_ug is not None and avg_result_ug is not None and avg_result_ug != 0:
            cv_percent_ug = (std_dev_ug / avg_result_ug) * 100
        
        # Calculate standard deviation for nmol/g FW
        valid_rep_results_nmol = [r for r in rep_results_nmol if r is not None]
        std_dev_nmol = np.std(valid_rep_results_nmol) if len(valid_rep_results_nmol) > 1 else None
        cv_percent_nmol = None
        if std_dev_nmol is not None and avg_result_nmol is not None and avg_result_nmol != 0:
            cv_percent_nmol = (std_dev_nmol / avg_result_nmol) * 100
        
        return {
            amino_acid_name: {
                'average': avg_result_ug,
                'std_dev': std_dev_ug,
                'cv_percent': cv_percent_ug,
                'replicates': rep_results_ug
            }
        }, {
            amino_acid_name: {
                'average': avg_result_nmol,
                'std_dev': std_dev_nmol,
                'cv_percent': cv_percent_nmol,
                'replicates': rep_results_nmol
            }
        }
    
    def _calculate_single_amino_acid_both_units(self, amino_acid_name, peak_area):
        """
        Calculate concentration for a single amino acid in both µg/g FW and nmol/g FW
        
        Returns:
        - Tuple with two values: (µg/g FW, nmol/g FW)
        """
        if peak_area is None:
            return None
            
        # Get molecular weight for this amino acid
        mol_weight = self.molecular_weight.get(amino_acid_name)
        if mol_weight is None:
            raise ValueError(f"Molecular weight not found for {amino_acid_name}")
        
        # Get standard curve for this amino acid or use default
        curve = self.standard_curves.get(amino_acid_name)
        if curve is None:
            curve = self.standard_curves.get('default', {'slope': 1.0, 'intercept': 0.0})
        
        # Calculate concentration using standard curve: (y - b) / m
        if curve['slope'] == 0:
            raise ValueError(f"Invalid standard curve for {amino_acid_name}: slope is zero")
        
        # Calculate nmol/ml using the standard curve
        nmol_per_ml = (peak_area - curve['intercept']) / curve['slope']
        
        # Convert to µg/ml: nmol/ml * MW / 1000
        ug_per_ml = nmol_per_ml * mol_weight / 1000
        
        # Convert to µg/g FW: (µg/ml * volume) / fresh weight
        ug_per_g_fw = (ug_per_ml * self.volume) / self.initial_fw
        
        # Convert to nmol/g FW: (nmol/ml * volume) / fresh weight
        nmol_per_g_fw = (nmol_per_ml * self.volume) / self.initial_fw
        
        return ug_per_g_fw, nmol_per_g_fw
    
    def calculate_group_statistics(self, results):
        """
        Calculate statistics for each group of samples
        
        Parameters:
        - results: Dictionary with amino acids and their calculated concentrations
        
        Returns:
        - Dictionary with group statistics for each amino acid
        """
        if not self.sample_groups or not self.sample_names:
            return None
        
        # Group samples by their assigned groups
        group_samples = defaultdict(list)
        for i, sample in enumerate(self.sample_names):
            if sample in self.sample_groups:
                group_samples[self.sample_groups[sample]].append(i)
        
        # Calculate statistics for each group
        group_stats = {}
        
        for aa, data in results.items():
            group_stats[aa] = {}
            
            for group, sample_indices in group_samples.items():
                # Get replicate values for this group
                group_values = []
                for idx in sample_indices:
                    if idx < len(data['replicates']) and data['replicates'][idx] is not None:
                        group_values.append(data['replicates'][idx])
                
                if group_values:
                    # Calculate group statistics
                    group_mean = np.mean(group_values)
                    group_std = np.std(group_values) if len(group_values) > 1 else None
                    group_cv = (group_std / group_mean * 100) if group_std is not None and group_mean != 0 else None
                    
                    group_stats[aa][group] = {
                        'mean': group_mean,
                        'std_dev': group_std,
                        'cv_percent': group_cv,
                        'n': len(group_values),
                        'values': group_values
                    }
        
        return group_stats

def save_to_csv(results, sample_names, output_file, unit_label="µg/g FW"):
    """
    Save results to a CSV file, including average, standard deviation (absolute and percent), and replicate values
    
    Parameters:
    - results: Dictionary with amino acids and their calculated concentrations
    - sample_names: List of sample names for the replicates
    - output_file: Path to save the output CSV file
    - unit_label: Unit label to use in the CSV headers (e.g., "µg/g FW" or "nmol/g FW")
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header row with sample names
        header = ['Amino Acid', f'Average ({unit_label})', f'Std Dev ({unit_label})', 'CV (%)']
        for sample in sample_names:
            header.append(f"{sample} ({unit_label})")
        writer.writerow(header)
        
        # Write data rows
        for aa, data in sorted(results.items()):
            row = [
                aa, 
                round(data['average'], 4) if data['average'] is not None else 'N/A',
                round(data['std_dev'], 4) if data['std_dev'] is not None else 'N/A',
                round(data['cv_percent'], 2) if data['cv_percent'] is not None else 'N/A'
            ]
            
            # Add replicate values
            for rep_value in data['replicates']:
                row.append(round(rep_value, 4) if rep_value is not None else 'N/A')
            
            writer.writerow(row)

def save_group_statistics_to_csv(group_stats, output_file, unit_label="µg/g FW"):
    """
    Save group statistics to a CSV file
    
    Parameters:
    - group_stats: Dictionary with group statistics for each amino acid
    - output_file: Path to save the output CSV file
    - unit_label: Unit label to use in the CSV headers (e.g., "µg/g FW" or "nmol/g FW")
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Get all groups
        all_groups = set()
        for aa, groups in group_stats.items():
            all_groups.update(groups.keys())
        all_groups = sorted(all_groups)
        
        # Write header row
        header = ['Amino Acid']
        for group in all_groups:
            header.extend([
                f"{group} - Mean ({unit_label})",
                f"{group} - Std Dev ({unit_label})",
                f"{group} - CV (%)",
                f"{group} - n"
            ])
        writer.writerow(header)
        
        # Write data rows
        for aa, groups in sorted(group_stats.items()):
            row = [aa]
            
            for group in all_groups:
                if group in groups:
                    data = groups[group]
                    row.extend([
                        round(data['mean'], 4) if data['mean'] is not None else 'N/A',
                        round(data['std_dev'], 4) if data['std_dev'] is not None else 'N/A',
                        round(data['cv_percent'], 2) if data['cv_percent'] is not None else 'N/A',
                        data['n']
                    ])
                else:
                    row.extend(['N/A', 'N/A', 'N/A', 'N/A'])
            
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description='Enhanced Amino Acid Calculation with Sample Grouping and Multiple Unit Output')
    
    # Molecular weight source options
    mol_weight_group = parser.add_mutually_exclusive_group()
    mol_weight_group.add_argument('--mol_weight', help='Custom molecular weight file (csv, tsv, txt)', default=None)
    mol_weight_group.add_argument('--preset_mol_weight', action='store_true', help='Use preset molecular weights')
    
    # Other parameters
    parser.add_argument('--standard_curve', help='Standard curve source file (csv, tsv, txt)', default=None)
    parser.add_argument('--peak_areas', help='Peak areas source file (csv, tsv, txt)', required=True)
    parser.add_argument('--initial_fw', type=float, help='Initial Fresh Weight (g)', required=True)
    parser.add_argument('--volume', type=float, help='Volume (ml)', required=True)
    parser.add_argument('--output_ug', help='Output CSV file for µg/g FW results', default='amino_acid_results_ug.csv')
    parser.add_argument('--output_nmol', help='Output CSV file for nmol/g FW results', default='amino_acid_results_nmol.csv')
    parser.add_argument('--amino_acid', help='Specific amino acid to calculate', default=None)
    
    # Sample grouping parameter
    parser.add_argument('--group_mapping', help='Sample to group mapping file (csv, tsv, txt)', default=None)
    parser.add_argument('--group_output_ug', help='Group statistics output CSV file for µg/g FW', default=None)
    parser.add_argument('--group_output_nmol', help='Group statistics output CSV file for nmol/g FW', default=None)
    
    args = parser.parse_args()
    
    try:
        # Initialize calculator
        calculator = AminoAcidCalculator(
            molecular_weight_source=args.mol_weight,
            preset_molecular_weight=args.preset_mol_weight,
            standard_curve_source=args.standard_curve,
            peak_area_source=args.peak_areas,
            initial_fw=args.initial_fw,
            volume=args.volume,
            group_mapping_source=args.group_mapping
        )
        
        # Calculate results in both units
        results_ug, results_nmol = calculator.calculate_amino_acid(args.amino_acid)
        
        # Save results to CSV files
        save_to_csv(results_ug, calculator.sample_names, args.output_ug, "µg/g FW")
        save_to_csv(results_nmol, calculator.sample_names, args.output_nmol, "nmol/g FW")
        print(f"Results saved to {args.output_ug} and {args.output_nmol}")
        
        # If group mapping was provided, calculate and save group statistics
        if calculator.sample_groups:
            # Calculate and save group statistics for µg/g FW
            if args.group_output_ug:
                group_stats_ug = calculator.calculate_group_statistics(results_ug)
                if group_stats_ug:
                    save_group_statistics_to_csv(group_stats_ug, args.group_output_ug, "µg/g FW")
                    print(f"Group statistics (µg/g FW) saved to {args.group_output_ug}")
            
            # Calculate and save group statistics for nmol/g FW
            if args.group_output_nmol:
                group_stats_nmol = calculator.calculate_group_statistics(results_nmol)
                if group_stats_nmol:
                    save_group_statistics_to_csv(group_stats_nmol, args.group_output_nmol, "nmol/g FW")
                    print(f"Group statistics (nmol/g FW) saved to {args.group_output_nmol}")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
