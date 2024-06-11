import subprocess
import argparse

'''
Hello I am here.
'''

def main():
    parser = argparse.ArgumentParser(description="Epinano Data Processing and Plotting Pipeline")
    # Assuming initial raw files are required for processing
    parser.add_argument('--file1', type=str, required=True, help="Path to the first raw data CSV file")
    parser.add_argument('--file2', type=str, required=True, help="Path to the second raw data CSV file")
    parser.add_argument('--file3', type=str, default="output.csv", help="Filename for the output of processed data")
    parser.add_argument('--processed_output', type=str, default="processed_data.csv", help="Filename for the output of processed data")
    parser.add_argument('--label1', default='Replicate 1', help='Label for the first processed data')
    parser.add_argument('--label2', default='Replicate 2', help='Label for the second processed data')
    parser.add_argument('--cap', action='store_true', help='Apply cap to the plots')
    args = parser.parse_args()

    # Step 1: Process the raw data files with process_epinano_replicates.py
    print("Processing data...")
    subprocess.run([
        'python', 'process_epinano_replicates.py', 
        '--mode', 'combine',
        '--file1', args.file1, 
        '--file2', args.file2, 
        '--output', args.processed_output
    ], check=True)

    # Assuming the processing step generates two files named like the following (adjust as needed):
    processed_file1 = args.processed_output.replace('.csv', '_1.csv')
    processed_file2 = args.processed_output.replace('.csv', '_2.csv')

    # Step 2: Generate plots with plot_epinano_replicates.py
    # Have to fix plotting schemes
    print("Generating plots...")
    plot_command = [
        'python', 'plot_epinano_replicates.py', 
        '--replicate1', processed_file1, 
        '--replicate2', processed_file2, 
        '--label1', args.label1, 
        '--label2', args.label2
    ]
    if args.cap:
        plot_command.append('--cap')
    subprocess.run(plot_command, check=True)

    print("Pipeline complete.")

if __name__ == "__main__":
    main()