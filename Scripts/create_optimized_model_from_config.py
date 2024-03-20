
import os
os.environ['OMP_NUM_THREADS'] = '16' 
# this must be set before pykeen in imported so that it can take effect before things get started

from pykeen.pipeline import pipeline_from_path
import torch
import argparse

torch.manual_seed(42) 

def get_args():
    parser = argparse.ArgumentParser(description='Create optimized model from config')
    parser.add_argument('--config', type=str, help='Path to config file')
    parser.add_argument('--out', type=str, help='Path to save model')
    return parser.parse_args()

args = get_args()

pipeline_result = pipeline_from_path(args.config,device='cuda')

pipeline_result.save_to_directory(args.out)