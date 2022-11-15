import os
import time

from iris_validation import generate_report


INPUT_DIR = './example_input/'
OUTPUT_DIR = './example_report/'

PDB_ID = '2a0x'
ROOT_PATH = str(os.path.join(INPUT_DIR, PDB_ID)) + '{suffix}'


if __name__ == '__main__':
    t0 = time.time()

    generate_report(latest_model_path=ROOT_PATH.format(suffix='_final.pdb'),
                    latest_reflections_path=ROOT_PATH.format(suffix='_final.mtz'),
                    latest_sequence_path=ROOT_PATH.format(suffix='.fasta'),
                    latest_distpred_path=ROOT_PATH.format(suffix='.npz'),
                    previous_model_path=ROOT_PATH.format(suffix='_0cyc.pdb'),
                    previous_reflections_path=ROOT_PATH.format(suffix='_0cyc.mtz'),
                    previous_sequence_path=ROOT_PATH.format(suffix='.fasta'),
                    previous_distpred_path=ROOT_PATH.format(suffix='.npz'),
                    run_covariance=False,
                    run_molprobity=False,
                    multiprocessing=True,
                    output_dir=OUTPUT_DIR)

    t1 = time.time()
    print('Time taken:', round(t1-t0, 2), 's')
    print()
