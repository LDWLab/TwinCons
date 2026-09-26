#!/usr/bin/env python3
'''Tests for the tools/ scripts. Run from the repository root:
python -m unittest discover -s tools/tests -t tools'''
import csv
import json
import os
import pickle
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from sklearn import svm
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from sklearn.metrics import log_loss
from sklearn.model_selection import KFold

TOOLS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_DIR = os.path.dirname(TOOLS_DIR)
ALIGNMENT_DIR = os.path.join(REPO_DIR, 'tests', 'input_test_data', 'alns')
sys.path.insert(0, TOOLS_DIR)

from twcCalibrate import calibrate
from twcModels import SVCWithProbabilities
from twcSVMtest import fitted_base_classifier, read_features
from twcSVMtrain import train_classifier


def two_classes(n=80, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.random((n, 2))
    y = (X[:, 0] + X[:, 1] + rng.normal(0, 0.2, n) > 1).astype(int)
    return X, y


class TestSVCWithProbabilities(unittest.TestCase):
    def setUp(self):
        self.X, self.y = two_classes()
        self.model = SVCWithProbabilities(C=1, gamma='auto').fit(self.X, self.y)

    def test_decision_function_is_the_plain_svc(self):
        plain = svm.SVC(C=1, gamma='auto').fit(self.X, self.y)
        np.testing.assert_array_equal(self.model.decision_function(self.X), plain.decision_function(self.X))
        np.testing.assert_array_equal(self.model.predict(self.X), plain.predict(self.X))

    def test_probabilities_follow_the_decision_function(self):
        probabilities = self.model.predict_proba(self.X)
        np.testing.assert_allclose(probabilities.sum(axis=1), 1)
        self.assertGreater(np.corrcoef(probabilities[:, 1], self.model.decision_function(self.X))[0, 1], 0.9)

    def test_small_classes(self):
        X, y = two_classes(40, seed=3)
        y = y.copy()
        y[:] = 0
        y[:3] = 1
        self.assertEqual(SVCWithProbabilities(gamma='auto').fit(X, y).predict_proba(X).shape, (40, 2))
        y[:3] = [1, 0, 0]
        self.assertFalse(hasattr(SVCWithProbabilities(gamma='auto').fit(X, y), 'predict_proba'))

    def test_survives_pickling(self):
        restored = pickle.loads(pickle.dumps(self.model))
        np.testing.assert_array_equal(restored.predict_proba(self.X), self.model.predict_proba(self.X))


class TestTrainClassifier(unittest.TestCase):
    def test_sample_weight_forms(self):
        X, y = two_classes()
        for weights in (None, [], list(np.ones(len(y))), np.ones(len(y))):
            classifier = train_classifier(X, y, 1, 'auto', 'rbf', sample_weight=weights)
            self.assertIsInstance(classifier, svm.SVC)

    def test_probabilities_on_request(self):
        X, y = two_classes()
        self.assertFalse(hasattr(train_classifier(X, y, 1, 'auto', 'rbf'), 'predict_proba'))
        self.assertTrue(hasattr(train_classifier(X, y, 1, 'auto', 'rbf', probabilities=True), 'predict_proba'))


class TestCalibrate(unittest.TestCase):
    def test_keeps_the_calibration_with_lower_log_loss(self):
        X, y = two_classes(150, seed=1)
        clf = svm.SVC(gamma='auto').fit(X[:50], y[:50])
        X_valid, y_valid, X_test, y_test = X[50:100], y[50:100], X[100:], y[100:]
        chosen = calibrate(clf, X_test, y_test, X_valid, y_valid, np.ones(len(y_test)))
        losses = {}
        for method in ('isotonic', 'sigmoid'):
            calibrated = CalibratedClassifierCV(FrozenEstimator(clf), method=method, cv=KFold(n_splits=2)).fit(X_valid, y_valid)
            losses[method] = log_loss(y_test, calibrated.predict_proba(X_test))
        self.assertEqual(chosen.method, min(losses, key=losses.get))

    def test_small_validation_set(self):
        X, y = two_classes(60, seed=2)
        clf = svm.SVC(gamma='auto').fit(X[:40], y[:40])
        X_valid = np.array([[0.1, 0.1], [0.2, 0.1], [0.1, 0.3], [0.9, 0.8], [0.8, 0.9], [0.2, 0.2]])
        y_valid = np.array([0, 0, 0, 1, 1, 0])
        chosen = calibrate(clf, X[40:], y[40:], X_valid, y_valid, np.ones(20))
        self.assertEqual(chosen.predict_proba(X[40:]).shape, (20, 2))

    def test_distances_come_from_the_fitted_classifier(self):
        X, y = two_classes()
        clf = svm.SVC(gamma='auto').fit(X, y)
        calibrated = CalibratedClassifierCV(FrozenEstimator(clf), cv=KFold(n_splits=2)).fit(X, y)
        np.testing.assert_array_equal(fitted_base_classifier(calibrated).decision_function(X), clf.decision_function(X))


class TestReadFeatures(unittest.TestCase):
    def test_classifiers_without_minimums(self):
        with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as fh:
            json.dump([{'maxX': 2.0, 'maxY': 3.0}, ['args']], fh)
        self.addCleanup(os.remove, fh.name)
        self.assertEqual(read_features(fh.name), (['args'], [2.0, 3.0, 0.0, 0.0]))


class TestPipeline(unittest.TestCase):
    '''Segments -> SVM training -> scoring of new segments, as run from the command line.'''
    def test_segments_train_and_test(self):
        work = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, work)
        os.makedirs(os.path.join(work, 'alns'))
        # Training labels come from file name prefixes: A_/B_ positive, C_/D_ negative.
        for name, label in (('uL02ab_txid_tagged.fas', 'A_'), ('bS01-RNAP7Ca.fa', 'C_')):
            shutil.copy(os.path.join(ALIGNMENT_DIR, name), os.path.join(work, 'alns', label + name))
        env = dict(os.environ, MPLBACKEND='Agg')
        def run(script, *args):
            subprocess.run([sys.executable, os.path.join(TOOLS_DIR, script), *args], cwd=work, env=env,
                           check=True, capture_output=True)
        run('twcCalculateSegments.py', '-a', 'alns/', 'segments', '-c', '-co', 'cg', 'gt_0.9', 'mx_blosum62')
        run('twcSVMtrain.py', 'segments.csv', 'model.pkl', '-ts', '1', '-twca', 'gt_0.9', 'cg', 'mx_blosum62')
        run('twcSVMtest.py', 'segments.csv', 'scores', 'model.pkl', '-tqa', '-ts', '1')
        with open(os.path.join(work, 'scores')) as fh:
            rows = list(csv.DictReader(fh))
        self.assertGreater(len(rows), 0)
        probabilities = np.array([float(row['Probability']) for row in rows])
        self.assertTrue(np.isfinite(probabilities).all())
        self.assertTrue(((probabilities >= 0) & (probabilities <= 1)).all())


if __name__ == '__main__':
    unittest.main()
