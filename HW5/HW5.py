'''
IMPORT PART
'''
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer as Imputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest
from sklearn.linear_model import LassoCV
from sklearn.feature_selection import mutual_info_regression
import numpy as np
import pandas as pd
import urllib.request
import json
import molecular_descriptors
'''
DESCRIPTORS PART
'''
def fit_Lasso(X_train, X_test, y_train, y_test):

    a = Imputer(missing_values=np.nan, strategy='median')
    b = StandardScaler()
    c = SelectKBest(score_func=mutual_info_regression)
    clf = LassoCV(cv=10)
    model = Pipeline([('impute', a), ('scaling', b), ('anova', c), ('rf', clf)])
    parameters = {'anova__k': [5, 10, 20, 40]}
    grid = GridSearchCV(model, parameters)
    grid.fit(X_train, y_train)
    y_pred = grid.predict(X_test)
    metric = [grid.score(X_test, y_test),
               metrics.explained_variance_score(y_test, y_pred),
               metrics.mean_absolute_error(y_test, y_pred),
               metrics.mean_squared_error(y_test, y_pred),
               metrics.median_absolute_error(y_test, y_pred),
               metrics.r2_score(y_test, y_pred)]
    return grid, y_pred, metric

def desc_calc(data, mode='train') -> pd.DataFrame:
    return molecular_descriptors.getAllDescriptors(data, mode)

def sar_model_evaluation(descriptors: pd.DataFrame):

    y = descriptors['Target']
    X = descriptors.drop(columns=['Target'])
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    model1, y_pred1, metrics1 = fit_Lasso(X_train, X_test, y_train, y_test)
    return model1, y_pred1, metrics1

def sar_model_train(descriptors_train: pd.DataFrame, indices):

    y_train = descriptors_train['Target']
    X_train = descriptors_train.drop(columns=['Target'])
    X_train = X_train[X_train.columns[indices]]
    a = Imputer(missing_values=np.nan, strategy='median')
    b = StandardScaler() 
    clf = LassoCV()
    model = Pipeline([('impute', a), ('scaling', b), ('rf', clf)])
    model.fit(X_train, y_train)
    return model
    
def sar_model_predict(model, descriptors_pred, indices):

    X_pred = descriptors_pred
    X_pred = X_pred[X_pred.columns[indices]]
    return model.predict(X_pred)
'''
PUBCHEM PART
'''
def pubchem_parsing(url):

    req = urllib.request.Request(url)
    res = urllib.request.urlopen(req).read()
    fin = json.loads(res.decode())
    return fin

def get_similar_cids(compound_smiles, threshold=95, maxentries=10):

    pubchem_pug_rest_api_link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
    pubchem_pug_rest_api_link+= "fastsimilarity_2d/smiles/%(smiles)s/cids/JSON?Threshold=%(threshold)s&MaxRecords=%(maxentries)s" % {
        "smiles": compound_smiles, "threshold": threshold, "maxentries": maxentries}
    similar_cids = pubchem_parsing(pubchem_pug_rest_api_link)['IdentifierList']['CID']
    return similar_cids

def get_xlogp(compound_cid):

    pubchem_pug_rest_api_link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    pubchem_pug_rest_api_link += "compound/cid/%s/property/XLogP/JSON" % compound_cid
    try:
        xlogp = pubchem_parsing(pubchem_pug_rest_api_link)['PropertyTable']['Properties'][0]['XLogP']
        return xlogp
    except KeyError:
        return None
'''
MAIN PART
'''
if __name__ == "__main__":

    pd.set_option.use_inf_as_na = True
    train_data = pd.read_csv('logp_100.csv')
    pred_data = pd.read_csv('logp_inputs.csv')
    cpds = [row for row in pred_data.loc[:, 'SMILES']]
    print("Calculating descriptors for training data...")
    train_descriptors = desc_calc(train_data)
    print("Calculating descriptors for prediction data...")
    pred_descriptors = desc_calc(pred_data, mode='notTrain')
    print("Evaluating regression model parameters...")
    model = sar_model_evaluation(train_descriptors)
    print('Best parameters are:', model[0].best_params_)
    cols = model[0].best_estimator_.named_steps['anova'].get_support(indices=True)
    print(cols)
    bestParams = model[0].best_params_
    print(bestParams)
    print("Training the model with the best parameters...")
    final_model = sar_model_train(train_descriptors, cols)
    for cpd in cpds:
        cpd_descriptors = pred_descriptors[pred_descriptors['SMILES']==cpd]
        pred = sar_model_predict(final_model, cpd_descriptors, cols)
        print(f"Predicted LogP value for compound {cpd}:", pred)
        result = []
        print("Searching for similar compunds...")
        similarity = get_similar_cids(compound_smiles=cpd)
        print("Filtering logP...")
        for cid in similarity:
            xlogp = get_xlogp(cid)
            if xlogp:
                if xlogp <= pred*1.1 and xlogp >=pred*0.9:
                    result.append((cid, xlogp))
        print(f"Request for compound {cpd} completed. found the following CIDs in PubChem with XLogP in the range of {pred}+- 10%: {result}")
