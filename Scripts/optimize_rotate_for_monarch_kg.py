import sys
import os
os.environ['OMP_NUM_THREADS'] = '16' 
# this must be set before pykeen in imported so that it can take effect before things get started

from pykeen.pipeline import pipeline
from pykeen.datasets.base import PathDataset
from pykeen.pipeline import pipeline
from pykeen.hpo import hpo_pipeline
from optuna.samplers import RandomSampler, QMCSampler
import torch

job_id = ''
if len(sys.argv) > 1:
    job_id = sys.argv[1]

torch.manual_seed(42)

TEST_PATH =  '/scratch/Shares/layer/workspace/michael_sandbox/LinkPrediction/ELs_for_Rotate/Monarch_KG/test.txt'
TRAIN_PATH = '/scratch/Shares/layer/workspace/michael_sandbox/LinkPrediction/ELs_for_Rotate/Monarch_KG/train.txt'
VALID_PATH = '/scratch/Shares/layer/workspace/michael_sandbox/LinkPrediction/ELs_for_Rotate/Monarch_KG/valid.txt'

study_name = "rotate_monarch_kg_{}".format(job_id)  # Unique identifier of the study.
storage_name = "sqlite:///{}.db".format(study_name) # this for a sqlite database in the current working directory

hpo_pipeline_result = hpo_pipeline(
    # n_trials=30,
    n_trials=4,
    # sampler=RandomSampler, # this was the original one I used, commenting out to try something new
    sampler=QMCSampler,

    model_kwargs={
            "random_seed": 42,
        },

    model_kwargs_ranges=dict(
            embedding_dim=dict(
                type=int,
                low=64,
                high=512,
                q=16,
            ),
        ),

    training_kwargs_ranges=dict(
            num_epochs=dict(
                type=int,
                low=100,
                high=1000,
                q=100,
            ),
        ),
    
    optimizer_kwargs_ranges=dict(
            lr=dict(
                type=float,
                low=0.001,
                high=0.1,
                log=True,
            ),
        ),
    
    negative_sampler_kwargs_ranges=dict(
            num_negs_per_pos=dict(
                type=int,
                low=1,
                high=100,
                q=10,
            ),
        ),
    training=TRAIN_PATH,
    testing=TEST_PATH,
    validation=VALID_PATH,
    model='rotate',
    device='cuda',
    storage=storage_name,
    load_if_exists=True,
    study_name=study_name,
    stopper='early',
    stopper_kwargs={'frequency':10, 'patience':2, 'relative_delta':0.002},
    loss="NSSA",
)

hpo_pipeline_result.save_to_directory('PyKeenOut/{}'.format(study_name))
