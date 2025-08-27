import papermill as pm
import os

os.makedirs('temp_qm_scans/', exist_ok=True)

# Define the range of conf_id values
conf_ids = range(0,1)

# Define multiple idx sets
idx_sets = [
    (1,2,3,4),
    (2,3,4,5),
]
# Loop over each idx set and conf_id value
for idx1, idx2, idx3, idx4 in idx_sets:
    for conf_id in conf_ids:
        pm.execute_notebook(
            'qm_scans_forward.ipynb',
            f'temp_qm_scans/temp_qm_scans_forward_{conf_id}.ipynb',
            parameters=dict(idx1=idx1, idx2=idx2, idx3=idx3, idx4=idx4, conf_id=conf_id)
        )
