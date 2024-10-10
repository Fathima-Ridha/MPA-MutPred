from sklearn.ensemble import GradientBoostingRegressor
import joblib
from sklearn.metrics import mean_absolute_error
from scipy.stats import pearsonr, spearmanr
import pandas as pd
import numpy as np
import random

seed_value = 42
np.random.seed(seed_value)
random.seed(seed_value)

# Load your dataset
dataset_path = '/home/MPA-MutPred/train.csv'
df = pd.read_csv(dataset_path)


X = df.drop(columns=['ΔΔG (kcal/mol)'])
y = df['ΔΔG (kcal/mol)']

# GradientBoostingRegressor model with tuning conditions
model = GradientBoostingRegressor(n_estimators=100, max_depth=3, learning_rate=0.1, subsample=1.0, random_state=seed_value)
model.fit(X, y)

predictions = model.predict(X)

# Calculate Mean Absolute Error (MAE)
mae = mean_absolute_error(y, predictions)
#print(f'Mean Absolute Error (MAE): {mae}')

# Calculate Pearson correlation and its p-value
pearson_corr, p_value_pearson = pearsonr(y, predictions)


# Save the trained model to a file
joblib.dump(model, 'trained_gb_model.pkl')



