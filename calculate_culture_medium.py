import sys
from datetime import datetime
import math

def calculate_culture_medium(desired_volume_ml, medium_type="MSG", num_flasks=None):
    # Common reagents for all media
    inositol_concentration = 0.1  # g/l
    sucrose_concentration = 30  # g/l
    
    # Medium-specific concentrations
    if medium_type == "MSG":
        macrosalts_concentration = 50  # ml/l
        microsalts_concentration = 10  # ml/l
        f_solution_concentration = 10  # ml/l
        vitamins_concentration = 1  # ml/l
        l_glutamine_concentration = 50  # ml/l
        gelrite_concentration = 3  # g/l
        sorbitol_concentration = 0  # g/l
        aba_concentration = 0  # ml/l
        hydrolyzed_casein_concentration = 0  # g/l
        container_type = "Petri Dishes"
        container_volume = 25  # ml per container (changed from 10 to 25)
    elif medium_type == "MSGm":
        macrosalts_concentration = 50  # ml/l
        microsalts_concentration = 10  # ml/l
        f_solution_concentration = 10  # ml/l
        vitamins_concentration = 1  # ml/l
        l_glutamine_concentration = 50  # ml/l
        gelrite_concentration = 10  # g/l
        sorbitol_concentration = 30  # g/l
        aba_concentration = 12  # ml/l
        hydrolyzed_casein_concentration = 0  # g/l
        container_type = "Petri Dishes"
        container_volume = 25  # ml per container (changed from 10 to 25)
    elif medium_type == "MSG-Liquid":
        macrosalts_concentration = 50  # ml/l
        microsalts_concentration = 10  # ml/l
        f_solution_concentration = 10  # ml/l
        vitamins_concentration = 1  # ml/l
        l_glutamine_concentration = 50  # ml/l  # Total concentration
        # But we'll add only 2.5ml per flask of 50ml (which is 50ml/l)
        gelrite_concentration = 0  # g/l
        sorbitol_concentration = 0  # g/l
        aba_concentration = 0  # ml/l
        hydrolyzed_casein_concentration = 0  # g/l
        container_type = "Erlenmeyer Flasks"
        container_volume = 50  # ml per container
        
        # For MSG-Liquid, if num_flasks is specified, adjust the desired volume
        if num_flasks is not None:
            desired_volume_ml = num_flasks * container_volume
    elif medium_type == "BM":
        macrosalts_concentration = 50  # ml/l
        microsalts_concentration = 20  # ml/l
        f_solution_concentration = 10  # ml/l
        vitamins_concentration = 2  # ml/l
        l_glutamine_concentration = 0  # ml/l  # BM doesn't need L-glutamine after sterilization
        gelrite_concentration = 3  # g/l
        sorbitol_concentration = 0  # g/l
        aba_concentration = 0  # ml/l
        hydrolyzed_casein_concentration = 0.5  # g/l
        container_type = "tubes"
        container_volume = 10  # ml per container
    else:
        raise ValueError(f"Unknown medium type: {medium_type}")

    # Convert desired volume from ml to liters
    desired_volume_l = desired_volume_ml / 1000

    # Calculate the required amounts of each stock solution and reagent
    macrosalts_volume = macrosalts_concentration * desired_volume_l
    microsalts_volume = microsalts_concentration * desired_volume_l
    f_solution_volume = f_solution_concentration * desired_volume_l
    vitamins_volume = vitamins_concentration * desired_volume_l
    inositol_mass = inositol_concentration * desired_volume_l
    sucrose_mass = sucrose_concentration * desired_volume_l
    hydrolyzed_casein_mass = hydrolyzed_casein_concentration * desired_volume_l
    
    # Calculate additional values for the steps
    initial_water_volume = desired_volume_ml / 2
    
    # Special handling for MSG-Liquid
    if medium_type == "MSG-Liquid":
        # For liquid medium, we add 2.5ml glutamine per 50ml flask
        num_flasks = math.ceil(desired_volume_ml / container_volume)
        l_glutamine_volume = 2.5 * num_flasks  # 2.5ml per flask
        pre_additions_volume = desired_volume_ml - l_glutamine_volume
        flask_solution_volume = container_volume - 2.5  # 47.5ml per flask
    elif medium_type == "BM":
        # For BM, we fill to final volume
        pre_additions_volume = desired_volume_ml
        l_glutamine_volume = 0
        gelrite_mass = gelrite_concentration * desired_volume_l
    else:
        l_glutamine_volume = l_glutamine_concentration * desired_volume_l
        gelrite_mass = gelrite_concentration * desired_volume_l
        sorbitol_mass = sorbitol_concentration * desired_volume_l
        aba_volume = aba_concentration * desired_volume_l
        pre_additions_volume = desired_volume_ml - l_glutamine_volume - (aba_volume if medium_type == "MSGm" else 0)

    # Determine if we need to divide the solution
    if medium_type == "MSG-Liquid":
        num_parts = num_flasks
        volume_per_part = flask_solution_volume
        glutamine_per_part = 2.5  # fixed 2.5ml per flask
    elif medium_type == "BM":
        # For BM, we don't divide the solution
        num_parts = 1
        volume_per_part = desired_volume_ml
        gelrite_per_part = gelrite_mass
    elif pre_additions_volume > 500 and medium_type != "MSG-Liquid":
        # FIX: Use pre_additions_volume instead of desired_volume_ml for determining division
        num_parts = math.ceil(pre_additions_volume / 500)
        volume_per_part = pre_additions_volume / num_parts
        gelrite_per_part = gelrite_mass / num_parts
        glutamine_per_part = l_glutamine_volume / num_parts
        aba_per_part = aba_volume / num_parts if medium_type == "MSGm" else 0
    else:
        num_parts = 1
        volume_per_part = pre_additions_volume
        gelrite_per_part = gelrite_mass if medium_type != "MSG-Liquid" else 0
        glutamine_per_part = l_glutamine_volume
        aba_per_part = aba_volume if medium_type == "MSGm" else 0

    # Title and header with larger text using ASCII art
    if medium_type == "MSG":
        medium_name = "MSG (Becwar et al. 1989)"
    elif medium_type == "MSGm":
        medium_name = "MSGm (Maturation)"
    elif medium_type == "BM":
        medium_name = "BM (Gupta and Pullman 1991)"
    else:
        medium_name = "MSG-Liquid (Becwar et al. 1989)"
        
    results = [
        "================================================================",
        f" {medium_name} CULTURE MEDIUM PREPARATION PROTOCOL ",
        "================================================================",
        "",
        f"TO PREPARE {desired_volume_ml:.1f} ml OF {medium_type} CULTURE MEDIUM:",
        "----------------------------------------------------------------",
        "",
        "□ #1. STEP 1",
        "   □ Destilled water initial volume: {:.2f} ml".format(initial_water_volume),
        "   □ Macrosalts stock solution: {:.2f} ml".format(macrosalts_volume),
        "   □ Microsalts stock solution: {:.2f} ml".format(microsalts_volume),
        "   □ F stock solution: {:.2f} ml".format(f_solution_volume),
        "   □ Vitamins stock solution: {:.2f} ml".format(vitamins_volume),
        "   □ Inositol: {:.2f} g".format(inositol_mass),
        "   □ Sucrose: {:.2f} g".format(sucrose_mass),
    ]

    # Add medium-specific components
    if medium_type == "MSGm":
        results.append("   □ Sorbitol: {:.2f} g".format(sorbitol_mass))
    if medium_type == "BM":
        results.append("   □ Hydrolyzed casein: {:.2f} g".format(hydrolyzed_casein_mass))
    
    results.extend([
        "",
        "□ #2. STEP 2",
        f"   □ Fill with water up to: {pre_additions_volume:.2f} ml",
        "",
        "□ #3. STEP 3",
        "   □ Adjust the pH to 5.8",
        "",
    ])
    
    if medium_type == "MSG-Liquid":
        results.extend([
            "□ #4. STEP 4",
            f"   □ Divide the solution in {num_parts} {container_type} with {volume_per_part:.2f} ml each",
            "",
            "□ #5. STEP 5",
            "   □ No Gelrite needed for liquid medium",
            "",
        ])
    elif medium_type == "BM":
        results.extend([
            "□ #4. STEP 4",
            "   □ Keep as a single volume",
            "",
            "□ #5. STEP 5",
            f"   □ Add Gelrite: {gelrite_per_part:.2f} g",
            "",
        ])
    elif num_parts > 1:
        results.extend([
            "□ #4. STEP 4",
            f"   □ Divide the volume in {num_parts} parts of {volume_per_part:.2f} ml each",
            "",
            "□ #5. STEP 5",
            f"   □ Add Gelrite: {gelrite_per_part:.2f} g in each part",
            "",
        ])
    else:
        results.extend([
            "□ #4. STEP 4",
            "   □ Keep as a single volume",
            "",
            "□ #5. STEP 5",
            f"   □ Add Gelrite: {gelrite_per_part:.2f} g",
            "",
        ])
    
    results.extend([
        "□ #6. STEP 6",
        "   □ Autoclave the solutions",
        "",
    ])
    
    # Add Step 7 only for MSG and MSGm and MSG-Liquid
    if medium_type != "BM":
        results.extend([
            "□ #7. STEP 7",
            "   □ After sterilization, add:"
        ])
        
        if medium_type == "MSG-Liquid":
            results.append(f"   □ L-glutamine stock solution: {glutamine_per_part:.2f} ml in each {container_type}")
        elif num_parts > 1:
            results.append(f"   □ L-glutamine stock solution: {glutamine_per_part:.2f} ml in each part")
            if medium_type == "MSGm":
                results.append(f"   □ ABA stock solution: {aba_per_part:.2f} ml in each part")
        else:
            results.append(f"   □ L-glutamine stock solution: {l_glutamine_volume:.2f} ml")
            if medium_type == "MSGm":
                results.append(f"   □ ABA stock solution: {aba_volume:.2f} ml")
        
        results.append("")
    
    results.extend([
        "------------------------------------------------------------",
    ])
    
    if medium_type == "MSG-Liquid":
        results.append(f"Final volume will be {num_parts} {container_type} with {container_volume} ml each")
    else:
        results.append(f"Distribute the solution in {container_type}, with {container_volume} ml in each")
    
    results.extend([
        "",
        "",
        "------------------------------------------------------------",
        "Notes:",
        "- You can check off (☑) each box (□) as you complete the step",
        "- Check off individual reagents as you add them",
        "- Store unused medium at 4°C",
        "- Medium shelf life: 2 weeks after preparation",
        "============================================================",
        "  Developed by Leandro F. Oliveira (2025) ",
        "============================================================",
        "",
        "CITATIONS:",
    ])
    
    if medium_type == "MSG" or medium_type == "MSG-Liquid":
        results.extend([
            "Becwar MR, Noland TL, Wyckoff JL (1989) Maturation, germination, and conversion",
            "of Norway spruce (Picea abies L.) somatic embryos to plants. In Vitro Cell Dev",
            "Biol Plant 26:575–580. https://doi.org/10.1007/BF02623571.",
        ])
    elif medium_type == "MSGm":
        results.extend([
            "Modified from Becwar MR, Noland TL, Wyckoff JL (1989) for maturation purposes.",
            "MSGm includes additional components for embryo maturation.",
        ])
    elif medium_type == "BM":
        results.extend([
            "Gupta PK, Pullman GS (1991) Method for reproducing coniferous plants by somatic",
            "embryogenesis using abscisic acid and osmotic potential variation. US Patent",
            "5,036,007.",
        ])
    
    results.extend([
        "",
        "de Oliveira, LF (2025) Tool Galaxy to calculate and generate a desired quantity",
        f"of {medium_type} culture medium solution (Galaxy Version 1.0.0).",
        "============================================================",
    ])
    
    return results

def save_to_txt(results, output_path):
    with open(output_path, 'w') as file:
        file.write(f"Date and Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        for line in results:
            file.write(line + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python calculate_culture_medium.py <desired_volume_ml> <output_file_path> [medium_type] [num_flasks]")
        print("medium_type options: MSG (default), MSGm, MSG-Liquid, BM")
        print("num_flasks: Optional, number of flasks for MSG-Liquid (each flask will be 50ml)")
        sys.exit(1)
    
    desired_volume_ml = float(sys.argv[1])
    output_file_path = sys.argv[2]
    
    medium_type = "MSG"
    if len(sys.argv) >= 4:
        medium_type = sys.argv[3]

    num_flasks = None
    if len(sys.argv) >= 5 and medium_type == "MSG-Liquid":
        num_flasks = int(sys.argv[4])
    
    results = calculate_culture_medium(desired_volume_ml, medium_type, num_flasks)
    save_to_txt(results, output_file_path)
