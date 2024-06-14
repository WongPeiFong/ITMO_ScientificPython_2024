'''
IMPORT PART
'''

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import RidgeCV
from sklearn.feature_selection import SelectKBest, mutual_info_regression
import numpy as np
import pandas as pd
import urllib.request
import json
import molecular_descriptors
from regression import fit_Ridge

'''
DESCRIPTORS PART
'''

def desc_calc(data: pd.DataFrame, mode: str) -> pd.DataFrame:
    return molecular_descriptors.getAllDescriptors(data, mode)

def sar_model_evaluation(descriptors: pd.DataFrame):
    y = descriptors['Target']
    X = descriptors.drop(columns=['Target'])
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    model1, y_pred1, metrics1 = fit_Ridge(X_train, X_test, y_train, y_test)
    return model1, y_pred1, metrics1

def sar_model_train(descriptors_train: pd.DataFrame, indices):
    y_train = descriptors_train['Target']
    X_train = descriptors_train.drop(columns=['Target'])
    X_train = X_train.iloc[:, indices]  # Keeping only necessary descriptors according to ANOVA evaluation
    
    a = SimpleImputer(strategy='median')
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

    # loading data
    train_data = pd.read_csv('logpfull.csv')
    pred_data = pd.read_csv('logp_inputs.csv')
    cpds = [row for row in pred_data.loc[:, 'SMILES']]

    # calculating descriptors
    print("Calculating descriptors for training data...")
    train_descriptors = desc_calc(train_data, mode='train')
    print("Calculating descriptors for prediction data...")
    pred_descriptors = desc_calc(pred_data, mode='predict')

    # finding best estimator
    print("Evaluating regression model parameters...")
    model, y_pred, metrics_values = sar_model_evaluation(train_descriptors)
    print('Best parameters are:', model.best_params_)
    cols = model.best_estimator_.named_steps['feature_selection'].get_support(indices=True)  # these are indices from ANOVA

    # train the best estimator and predict values
    print("Training the model with the best parameters...")
    final_model = sar_model_train(train_descriptors, cols)

    for cpd in cpds:      
        cpd_descriptors = pred_descriptors[pred_descriptors['SMILES'] == cpd]
        pred = sar_model_predict(final_model, cpd_descriptors.drop(columns=['SMILES']), cols)
        print(f"Predicted LogP value for compound {cpd}:", pred)

        result = []

        print("Searching for similar compounds...")
        similarity = get_similar_cids(cpd)

        print("Filtering logP...")
        for cid in similarity:
            xlogp = get_xlogp(cid)
            if xlogp:
                if xlogp <= pred * 1.1 and xlogp >= pred * 0.9:
                    result.append((cid, xlogp))

        print(f"Request for compound {cpd} completed. I found the following CIDs in PubChem with XLogP in the range of {pred} ± 10%: {result}")
