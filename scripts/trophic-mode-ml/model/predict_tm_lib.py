#!/usr/bin/env python3

"""
The functions from the code below are from the script predict-tm.py, which is part of the trophic-mode-ml pipeline.
I reorganized the code into a class and functions to make it easier to use in a Jupyter notebook.
Author: Ben Lambert

Load data, extract feature set, train model with feature set, make predictions.

usage: ./predict-tm.py -d [Pfam expression profiles] -t [MMETSP training data] -l [Training labels] -o [Output path and file name for results]
        -f [Feature set] -use-rf [make predictions with random forest model]

"""

# predict_tm_lib.py

import pandas as pd
import numpy as np
import warnings
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier

# Define a warning just to be a bit more specific.
class DataFormatWarning(Warning):
    pass

class PfamLoader():

    # Initialize the class with the path to profiles you want to classify, the path to the MMETSP training dataset,
    # and the list of features (common, all selected features).
    def __init__(self, data, train_data, train_labels, features):
        self.pfam_path = data
        self.train_path = train_data
        self.feat_path = features
        self.labels_path = train_labels

    def load_pfam_profiles(self):
        d = pd.read_csv(self.pfam_path)

        # A common part of pfam annotation leaves you with pfam ids that have decimal + number. We want to
        # get rid of those!
        _cols = [col for col in d.columns if '.' in col]
        if len(_cols) != 0:
            d.columns = d.columns.str.split(".").str[0]

        # Run simple checks on the data content. We want more than 800 non-zero pfams.
        # Very simple check.
        if d.shape[1] < 800:
            warnings.warn("Pfam profiles have less than the suggested 800 non-zero columns.", DataFormatWarning)

        # After the simple check passes, its possible that some rows still have less than 800 non-zero elements.
        counts = list(d.astype(bool).sum(axis=1))
        if min(counts) < 800:
            warnings.warn("Pfam profiles have less than the suggested 800 non-zero columns.", DataFormatWarning)

        self.pfam_data = d

    def get_subset_profiles(self):
        feats = pd.read_csv(self.feat_path)
        train = pd.read_csv(self.train_path)
        l = pd.read_csv(self.labels_path)

        feats = feats.Pfam.str.split(".").str[0]

        _missing = list(set(feats)-set(self.pfam_data.columns))

        for m in _missing:
            self.pfam_data[f'{m}'] = 0

        self.pfam_data = self.pfam_data[feats]
        self.train = train[feats]

        idx = l.index[l['Trophic mode'] == 'Un']
        self.targets = l.drop(idx)
        self.train = self.train.drop(idx)

        return self.train, self.pfam_data, self.targets['Trophic mode']


class TMPredictor():

    # Initialize class with subset profiles obtained with PfamLoader().
    def __init__(self, train_data, train_labels, pfam_profiles, rf):
        self.train_data = train_data
        self.labels = train_labels
        self.profiles = pfam_profiles
        self.use_rf = rf
        self.le = LabelEncoder()

    def create_model(self):
        if self.use_rf:
            model = RandomForestClassifier(n_estimators=100, max_depth=1000)
        else:
            model = XGBClassifier(gamma=0.0, learning_rate=0.5, max_depth=3, n_estimators=10, reg_lambda=0.)

        self.model = model

    def predict(self):
        scaler = MinMaxScaler()
        X = scaler.fit_transform(self.train_data)
        
        # Convert the labels into integer form
        self.labels = self.le.fit_transform(self.labels)

        self.model.fit(X, self.labels)

        scale_profiles = scaler.transform(self.profiles)

        # The model will now output integer predictions, so convert these back to the original labels using inverse_transform
        predictions = self.le.inverse_transform(self.model.predict(scale_profiles))

        return predictions

def run_prediction(data, train, feats, labels, out, rf):
    loader = PfamLoader(data=data, train_data=train, train_labels=labels, features=feats)
    loader.load_pfam_profiles()
    train_data, profiles, targets = loader.get_subset_profiles()
    model = TMPredictor(train_data=train_data, train_labels=targets, pfam_profiles=profiles, rf=rf)
    model.create_model()
    predictions = model.predict()
    p_df = pd.DataFrame(data={'preds': predictions})
    p_df.to_csv(f'{out}', index=False)
    return predictions