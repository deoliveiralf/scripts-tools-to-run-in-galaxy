#!/usr/bin/env python3

import argparse
import csv
import os
import sys
import json

class Expense:
    def __init__(self, description, amount, payer, payee=None, for_both=False):
        self.description = description
        self.amount = float(amount)
        self.payer = payer
        self.payee = payee
        self.for_both = for_both

    def __repr__(self):
        return f"{self.description}: {self.amount} paid by {self.payer}" + (f" for {self.payee}" if self.payee else "") + (f" (shared)" if self.for_both else "")

class HouseholdExpenses:
    def __init__(self, csv_path=None):
        self.csv_path = csv_path
        self.expenses = self.load_expenses() if csv_path and os.path.exists(csv_path) else []
        self.bruno_paid = 0
        self.leandro_paid = 0
        self.bruno_owes = 0
        self.leandro_owes = 0
        self.calculate_balances()

    def load_expenses(self):
        try:
            expenses = []
            with open(self.csv_path, 'r', newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    for_both = False
                    if 'For Both' in row and row['For Both'].lower() in ['yes', 'true', '1']:
                        for_both = True
                    
                    payee = row.get('Payee', '').strip()
                    if not payee:
                        payee = None
                        
                    expenses.append(Expense(
                        row['Description'], 
                        float(row['Amount']), 
                        row['Payer'],
                        payee,
                        for_both
                    ))
            return expenses
        except Exception as e:
            sys.stderr.write(f"Error loading CSV: {e}\n")
            return []

    def save_expenses(self, output_path):
        try:
            with open(output_path, 'w', newline='') as csvfile:
                fieldnames = ['Description', 'Amount', 'Payer', 'Payee', 'For Both']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                for exp in self.expenses:
                    writer.writerow({
                        'Description': exp.description,
                        'Amount': exp.amount,
                        'Payer': exp.payer,
                        'Payee': exp.payee if exp.payee else '',
                        'For Both': 'Yes' if exp.for_both else 'No'
                    })
            sys.stderr.write(f"Expenses saved to {output_path}\n")
        except Exception as e:
            sys.stderr.write(f"Error saving CSV: {e}\n")

    def add_expense(self, description, amount, payer, payee=None, for_both=False):
        if payer not in ["Bruno", "Leandro"]:
            sys.stderr.write("Payer must be either Bruno or Leandro\n")
            return
        if payee and payee not in ["Bruno", "Leandro"]:
            sys.stderr.write("Payee must be either Bruno or Leandro\n")
            return
        if payer == payee:
            sys.stderr.write("Payer and payee cannot be the same person\n")
            return
            
        self.expenses.append(Expense(description, amount, payer, payee, for_both))
        self.calculate_balances()  # Recalculate all balances

    def calculate_balances(self):
        # Reset all counters
        self.bruno_paid = 0
        self.leandro_paid = 0
        self.bruno_owes = 0
        self.leandro_owes = 0
        
        for exp in self.expenses:
            if exp.for_both:
                # Case 1 & 3: Expense is for both (split 50/50)
                if exp.payer == "Bruno":
                    self.bruno_paid += exp.amount
                    self.leandro_owes += exp.amount / 2
                else:  # Leandro paid
                    self.leandro_paid += exp.amount
                    self.bruno_owes += exp.amount / 2
            elif exp.payee:
                # Case 2: One person paid for the other
                if exp.payer == "Bruno" and exp.payee == "Leandro":
                    self.bruno_paid += exp.amount
                    self.leandro_owes += exp.amount
                elif exp.payer == "Leandro" and exp.payee == "Bruno":
                    self.leandro_paid += exp.amount
                    self.bruno_owes += exp.amount
            else:
                # Personal expense, just track who paid
                if exp.payer == "Bruno":
                    self.bruno_paid += exp.amount
                else:
                    self.leandro_paid += exp.amount

    def calculate_settlement(self):
        result = "--- Summary ---\n"
        result += f"Bruno paid in total: ${self.bruno_paid:.2f}\n"
        result += f"Leandro paid in total: ${self.leandro_paid:.2f}\n\n"
        result += f"Bruno owes Leandro: ${self.bruno_owes:.2f}\n"
        result += f"Leandro owes Bruno: ${self.leandro_owes:.2f}\n\n"
        
        # Case 4: Final settlement calculation
        if self.leandro_owes > self.bruno_owes:
            net_amount = self.leandro_owes - self.bruno_owes
            result += f"Final settlement: Leandro owes Bruno ${net_amount:.2f}"
        elif self.bruno_owes > self.leandro_owes:
            net_amount = self.bruno_owes - self.leandro_owes
            result += f"Final settlement: Bruno owes Leandro ${net_amount:.2f}"
        else:
            result += "Final settlement: Everyone is even, no payment needed."
        
        return result

    def get_settlement_json(self):
        """Return settlement data as JSON for Galaxy to process"""
        net_owed = 0
        who_pays = ""
        who_receives = ""
        
        if self.leandro_owes > self.bruno_owes:
            net_amount = self.leandro_owes - self.bruno_owes
            who_pays = "Leandro"
            who_receives = "Bruno"
            net_owed = net_amount
        elif self.bruno_owes > self.leandro_owes:
            net_amount = self.bruno_owes - self.leandro_owes
            who_pays = "Bruno"
            who_receives = "Leandro"
            net_owed = net_amount
        
        return {
            "bruno_paid_total": self.bruno_paid,
            "leandro_paid_total": self.leandro_paid,
            "bruno_owes": self.bruno_owes,
            "leandro_owes": self.leandro_owes,
            "who_pays": who_pays,
            "who_receives": who_receives,
            "net_amount": net_owed
        }

    def display_expenses(self):
        if not self.expenses:
            return "No expenses recorded yet."
        
        # Create a formatted table header
        header = f"{'Description':<30} {'Amount':>10} {'Payer':<10} {'Payee':<10} {'Shared':<8}"
        separator = "-" * len(header)
        
        result = [header, separator]
        
        # Add each expense as a row
        for exp in self.expenses:
            payee_display = exp.payee if exp.payee else "-"
            for_both_display = "Yes" if exp.for_both else "No"
            amount_display = f"${exp.amount:.2f}"
            
            row = f"{exp.description:<30} {amount_display:>10} {exp.payer:<10} {payee_display:<10} {for_both_display:<8}"
            result.append(row)
        
        return "\n".join(result)

    def __repr__(self):
        return self.display_expenses() + "\n\n" + self.calculate_settlement()

def parse_args():
    parser = argparse.ArgumentParser(description='Household Expenses Calculator for Galaxy')
    
    # Input file option
    parser.add_argument('--input', '-i', help='Input CSV file with existing expenses')
    
    # Output files
    parser.add_argument('--output-csv', '-o', help='Output CSV file to save expenses')
    parser.add_argument('--output-summary', '-s', help='Output text file for summary')
    parser.add_argument('--output-json', '-j', help='Output JSON file for settlement details')
    
    # Option to hide settlement JSON from history
    parser.add_argument('--hide-settlement-json', action='store_true', 
                        help='Hide settlement JSON from history output')
    
    # Add new expense parameters - modified to handle multiple expenses
    parser.add_argument('--add-expense', action='append_const', const=True, default=[], 
                        help='Add a new expense - use once for each expense')
    parser.add_argument('--description', action='append', default=[], help='Expense description')
    parser.add_argument('--amount', type=float, action='append', default=[], help='Expense amount')
    parser.add_argument('--payer', action='append', default=[], help='Who paid for the expense')
    parser.add_argument('--expense-type', type=int, action='append', default=[], 
                        help='Expense type: 1=Shared, 2=For other person, 3=Personal')
    parser.add_argument('--process-immediately', action='store_true', 
                        help='Process each expense immediately without waiting for all expenses')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Initialize expenses tracker with input file if provided
    household_expenses = HouseholdExpenses(args.input)
    
    # Add new expenses if requested
    if args.add_expense:
        # Get the number of expenses to add (should match the length of description, amount, etc.)
        num_expenses = len(args.description)
        
        if num_expenses == 0:
            sys.stderr.write("Error: No expense details provided\n")
            sys.exit(1)
            
        # Validate that all lists have the same length
        if not all(len(lst) == num_expenses for lst in [args.description, args.amount, args.payer, args.expense_type]):
            sys.stderr.write("Error: Mismatch in the number of expense parameters\n")
            sys.exit(1)
        
        # Process each expense
        for i in range(num_expenses):
            description = args.description[i]
            amount = args.amount[i]
            payer = args.payer[i]
            expense_type = args.expense_type[i]
            
            for_both = False
            payee = None
            
            if expense_type == 1:  # Shared
                for_both = True
            elif expense_type == 2:  # For other person
                payee = "Leandro" if payer == "Bruno" else "Bruno"
            
            household_expenses.add_expense(description, amount, payer, payee, for_both)
    
    # Save expenses to CSV if output path provided
    if args.output_csv:
        household_expenses.save_expenses(args.output_csv)
    
    # Generate and save summary if requested
    if args.output_summary:
        with open(args.output_summary, 'w') as f:
            f.write(household_expenses.display_expenses())
            f.write("\n\n")
            f.write(household_expenses.calculate_settlement())
    
    # Generate and save JSON if requested
    if args.output_json:
        with open(args.output_json, 'w') as f:
            json.dump(household_expenses.get_settlement_json(), f, indent=2)
    
    # If no output was specified, print to stdout
    if not any([args.output_csv, args.output_summary, args.output_json]):
        print(household_expenses.display_expenses())
        print("\n")
        print(household_expenses.calculate_settlement())

if __name__ == "__main__":
    main()
