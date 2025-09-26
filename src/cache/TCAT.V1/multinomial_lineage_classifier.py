#Available local variables are `score_dir` containing the path to the directory containing 
#the scores yaml file, and `usage_for_score` containing the usage DataFrame to be used for
#computing the score. 

#The results are stored in `res` which is expected to be a pandas.Series with a result
#calculated for each cell. 

#Creating functions in this script can lead to confusion with namespace so we recommend
#keeping these scripts short and avoiding defining functions

import os, pickle
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import pandas as pd

model_fn = os.path.join(score_dir, 'TCAT.V1.MultiLogRegression.pkl') 

with open(model_fn, 'rb') as f:
    # Load the object from the pickle file
    model = pickle.load(f)

scaler_tcat = StandardScaler()
scaled_usage = pd.DataFrame(scaler_tcat.fit_transform(usage_for_score),
                                   index = usage_for_score.index,
                                   columns = usage_for_score.columns)
        
res = pd.Series(model.predict(scaled_usage), index = scaled_usage.index)