# -*- coding: utf-8 -*-

"""
This module originally came from https://github.com/AstraZeneca/biomedical-kg-topological-imbalance/
but has been modified to work with pykeen version 1.10.1 and my specific use cases.
"""

import pandas as pd
from pykeen.datasets.base import Dataset
from pykeen.models import Model
from pykeen import predict

def create_type(row) -> str:    
    if row["in_training"] is True and row["in_testing"] is False:
        return "train"
    elif row["in_training"] is False and row["in_testing"] is True:
        return "test"
    elif row["in_training"] is False and row["in_testing"] is False:
        return "novel"
    else:
        return "unknown"

def annotate_predicted_df(
    df: pd.DataFrame,
    degs: dict,
    position: str,
) -> pd.DataFrame:
    """Annotate a pykeen predictions dataframe."""

    df["entity_type"] = df[position].str.split(":", expand=True)[0]
    df["triple_type"] = df.apply(lambda row: create_type(row), axis=1)
    df["deg"] = [degs[e] for e in list(df[position].values)]

    return df

def make_string_triple(q_entity, q_relation, tail, data) -> str:
     return str([q_entity, q_relation, data.training.entity_id_to_label[tail]])

def get_predictions_tail(
    q_entity: str,
    q_relation: str,
    data: Dataset,
    model: Model,
    degs: dict,
    train_triples: set,
    test_triples: set,
    valid_triples: set
) -> pd.DataFrame:
    """Make a prediction using a a partial triple (missing tail)."""

    pred_df = predict.predict_target(
        model=model,
        head=q_entity,
        relation=q_relation,
        triples_factory=data
    ).df

    pred_df['in_training'] = [ make_string_triple(q_entity, q_relation, x, data) in train_triples or make_string_triple(q_entity, q_relation, x, data)[::-1] in train_triples for x in pred_df['tail_id']]
    pred_df['in_testing'] = [ make_string_triple(q_entity, q_relation, x, data) in test_triples or make_string_triple(q_entity, q_relation, x, data)[::-1] in test_triples for x in pred_df['tail_id']]

    pred_df['tail_label'] = [data.training.entity_id_to_label[x] for x in pred_df['tail_id']]

    pred_df = annotate_predicted_df(pred_df, degs, "tail_label")

    return pred_df


def get_predictions_head(
    q_entity: str,
    q_relation: str,
    data: Dataset,
    model: Model,
    degs: dict,
) -> pd.DataFrame:
    """Make a prediction using a a partial triple (missing head)."""

    pred_df = model.get_head_prediction_df(
        q_relation,
        q_entity,
        triples_factory=data.training,
        testing=data.testing.mapped_triples,
    )
    pred_df = annotate_predicted_df(pred_df, degs, "head_label")

    return pred_df
