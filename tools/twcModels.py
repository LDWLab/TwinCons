'''Classifiers saved by the tools/ scripts. Kept in their own module so pickled models load from any script.'''
import numpy as np
from sklearn import svm
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold
from sklearn.utils.metaestimators import available_if


class SVCWithProbabilities(svm.SVC):
    '''
    SVC whose predict_proba comes from sigmoid (Platt) calibration of cross-validated decision
    values, as SVC(probability=True) did before scikit-learn deprecated it in 1.9. predict and
    decision_function are those of the SVC fitted on all the data. Up to five folds are used, as
    many as the rarest class allows; with a single example of some class, predict_proba is
    unavailable.
    '''
    def fit(self, X, y, sample_weight=None):
        super().fit(X, y, sample_weight=sample_weight)
        folds = min(5, np.unique(y, return_counts=True)[1].min())
        self.probability_model_ = None
        if folds >= 2:
            self.probability_model_ = CalibratedClassifierCV(svm.SVC(**self.get_params()), method='sigmoid',
                                                             ensemble=False, cv=StratifiedKFold(folds))
            self.probability_model_.fit(X, y, sample_weight=sample_weight)
        return self

    @available_if(lambda self: getattr(self, 'probability_model_', None) is not None)
    def predict_proba(self, X):
        return self.probability_model_.predict_proba(X)
