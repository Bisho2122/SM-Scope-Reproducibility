import numpy as np
from skimage.transform import AffineTransform
import pandas as pd
from pathlib import Path


landmarks_path = "../../../Data/Figure 1/Landmarks/"
out_path = "../../../Data/Figure 1/"


def weighted_rmse(reference_landmarks, transformed_landmarks, weights):
    """
    Calculate the weighted Root Mean Square Error (RMSE) between reference and transformed landmarks.

    Parameters:
    - reference_landmarks: NumPy array of shape (N, 2) with reference landmark coordinates.
    - transformed_landmarks: NumPy array of shape (N, 2) with transformed landmark coordinates.
    - weights: NumPy array of shape (N,) with weights for each landmark.

    Returns:
    - Weighted RMSE value.
    """
    # Ensure the weights are a NumPy array
    weights = np.array(weights)

    # Calculate squared differences
    squared_diff = (reference_landmarks - transformed_landmarks) ** 2

    # Calculate weighted squared differences
    weighted_squared_diff = squared_diff * weights[:, np.newaxis]

    # Calculate the mean of the weighted squared differences
    mean_weighted_squared_diff = np.mean(weighted_squared_diff)

    # Return the square root of the mean weighted squared differences
    return np.sqrt(mean_weighted_squared_diff)
  
def transform_landmarks(landmarks, matrix):
    homogeneous_landmarks = np.hstack([landmarks, np.ones((landmarks.shape[0], 1))])
    transformed_landmarks = homogeneous_landmarks @ matrix.T
    return transformed_landmarks[:, :2]

reference_csv_path = landmarks_path + "Reference_landmarks.csv"
reference_landmarks_df = pd.read_csv(reference_csv_path)
reference_landmarks = reference_landmarks_df[['axis-0', 'axis-1']].to_numpy()

HPF_landmarks = np.array([0,1,2,6,7,8,9])
reference_landmarks = reference_landmarks[HPF_landmarks]


RMSE_res = []
for i in range(66,75):
    moving_csv_path = landmarks_path + "Section_" + str(i) + "_pts.csv"
    moving_landmarks_df = pd.read_csv(moving_csv_path)
    moving_landmarks = moving_landmarks_df[['axis-0', 'axis-1']].to_numpy()

    moving_landmarks =  moving_landmarks[HPF_landmarks]

    affine_transform = AffineTransform()
    affine_transform.estimate(moving_landmarks, reference_landmarks)
    transformation_matrix = affine_transform.params


    #Calculate RMSE
    transformed_landmarks = transform_landmarks(moving_landmarks, transformation_matrix)
    rmse_weights = np.ones(reference_landmarks.shape[0])

    rmse_value_no_extra = weighted_rmse(reference_landmarks, transformed_landmarks, rmse_weights)

    df_dict = {'Section' : 'Section_' + str(i)}
    df_dict["RMSE"] = rmse_value_no_extra

    RMSE_res.append(df_dict)

RMSE_res_df = pd.DataFrame(RMSE_res).drop_duplicates()
RMSE_res_df = RMSE_res_df.sort_values(by='RMSE', ascending=True)
RMSE_res_df.to_csv(out_path + "Atlas_sections_RMSE_results.csv", index=False)
