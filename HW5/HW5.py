'''
IMPORT PART
'''
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer as Imputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import RidgeCV
from sklearn.feature_selection import SelectKBest, mutual_info_regression
import molecular_descriptors  # Importing the molecular descriptors module
import numpy as np
import pandas as pd
import urllib.request
import json

'''
UTILITY FUNCTIONS
'''

def test_train_split(file, test_size_value=0.25, train_size_value=None):
    descriptors, target = descriptor_target_split(file)
    x_train, x_test, y_train, y_test = train_test_split(descriptors, target, test_size=test_size_value,
                                                        train_size=train_size_value)
    train = descriptor_target_join(x_train, y_train)
    train.reset_index(inplace=True)
    train.drop('index', axis=1, inplace=True)
    test = descriptor_target_join(x_test, y_test)
    test.reset_index(inplace=True)
    test.drop('index', axis=1, inplace=True)
    return train, test

def descriptor_target_split(file):
    target = file.loc[:, file.columns == 'Target']
    descriptors = file.loc[:, file.columns != 'Target']
    return descriptors, target

def descriptor_target_join(descriptors, target):
    descriptors['Target'] = target['Target']
    file = descriptors
    return file

'''
DESCRIPTORS PART
'''
def fit_Ridge(X_train, X_test, y_train, y_test):
    pipeline = Pipeline([
        ('impute', Imputer(strategy='median')),
        ('scaling', StandardScaler()),
        ('feature_selection', SelectKBest(score_func=mutual_info_regression)),
        ('model', RidgeCV())
    ])

    param_grid = {'feature_selection__k': [5, 10, 20, 40]}
    grid = GridSearchCV(pipeline, param_grid, cv=10)
    grid.fit(X_train, y_train)
    y_pred = grid.predict(X_test)

    metric = [
        grid.score(X_test, y_test),
        metrics.explained_variance_score(y_test, y_pred),
        metrics.mean_absolute_error(y_test, y_pred),
        metrics.mean_squared_error(y_test, y_pred),
        metrics.median_absolute_error(y_test, y_pred),
        metrics.r2_score(y_test, y_pred)
    ]

    return grid, y_pred, metric

def desc_calc(smiles_list, mode) -> pd.DataFrame:
    return molecular_descriptors.getAllDescriptors(smiles_list, mode)

def sar_model_evaluation(descriptors: pd.DataFrame, target: pd.Series):
    y = target
    X = descriptors
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    model1, y_pred1, metrics1 = fit_Ridge(X_train, X_test, y_train, y_test)
    return model1, y_pred1, metrics1

def sar_model_train(descriptors_train: pd.DataFrame, target_train: pd.Series, indices):
    y_train = target_train
    X_train = descriptors_train.iloc[:, indices]  # Keeping only necessary descriptors according to ANOVA evaluation
    
    a = Imputer(strategy='median')
    b = StandardScaler()
    clf = RidgeCV()
    model = Pipeline([('impute', a), ('scaling', b), ('model', clf)])  # Without ANOVA now

    model.fit(X_train, y_train)
    return model
    
def sar_model_predict(model, descriptors_pred: pd.DataFrame, indices):
    X_pred = descriptors_pred.iloc[:, indices]
    return model.predict(X_pred)

'''
PUBCHEM PART
'''

def pubchem_parsing(url):
    req = urllib.request.Request(url)
    res = urllib.request.urlopen(req).read()
    return json.loads(res.decode())

def get_similar_cids(compound_smiles, threshold=95, maxentries=10):
    pubchem_pug_rest_api_link = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        "fastsimilarity_2d/smiles/%(smiles)s/cids/JSON?Threshold=%(threshold)s&MaxRecords=%(maxentries)s" % {
            "smiles": compound_smiles, "threshold": threshold, "maxentries": maxentries}
    )
    similar_cids = pubchem_parsing(pubchem_pug_rest_api_link)['IdentifierList']['CID']
    return similar_cids

def get_xlogp(compound_cid):
    pubchem_pug_rest_api_link = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        "compound/cid/%s/property/XLogP/JSON" % compound_cid
    )
    try:
        xlogp = pubchem_parsing(pubchem_pug_rest_api_link)['PropertyTable']['Properties'][0]['XLogP']
        return xlogp
    except KeyError:
        return None

'''
MAIN PART
'''

if __name__ == "__main__":
    pd.set_option('use_inf_as_na', True)

    train_data = pd.read_csv('logpfull.csv')  # or 'logp100.csv' for faster computations
    pred_data = pd.read_csv('logp_inputs.csv')
    cpds = pred_data['SMILES'].tolist()

    print("Calculating descriptors for training data...")
    train_descriptors = desc_calc(train_data['SMILES'], mode='train')
    print("Calculating descriptors for prediction data...")
    pred_descriptors = desc_calc(pred_data['SMILES'], mode='predict')

    print("Evaluating regression model parameters...")
    model, y_pred, metrics_values = sar_model_evaluation(train_descriptors, train_data['LogP'])
    print('Best parameters are:', model.best_params_)
    cols = model.best_estimator_.named_steps['feature_selection'].get_support(indices=True)

    print("Training the model with the best parameters...")
    final_model = sar_model_train(train_descriptors, train_data['LogP'], cols)

    for cpd in cpds:
        cpd_descriptors = pred_descriptors[pred_descriptors['SMILES'] == cpd]
        pred = sar_model_predict(final_model, cpd_descriptors, cols)
        print(f"Predicted LogP value for compound {cpd}:", pred)

        result = []

        print("Searching for similar compounds...")
        similarity = get_similar_cids(cpd)

        print("Filtering logP...")
        for cid in similarity:
            xlogp = get_xlogp(cid)
            if xlogp is not None:
                if pred * 0.9 <= xlogp <= pred * 1.1:
                    result.append((cid, xlogp))

        print(f"Request for compound {cpd} completed. Found CIDs with XLogP in the range of {pred} ± 10%: {result}")
