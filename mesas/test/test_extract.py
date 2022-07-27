
import matplotlib.pyplot as plt

from mesas.sas.model import Model
import os
os.chdir('/Users/esthersida/Documents/Code/mesas')

eg = "hyporheic"
# Create the model
model = Model(data_df=f'./mesas/examples/{eg}/data.csv',
                  config=f'./mesas/examples/{eg}/config.json')
    
# Run the model
model.run()

# Extract results
data_df = model.data_df
flux = model.fluxorder[0]
# make plots
fig = plt.figure()
print(model.solorder)
model.get_mT(2)

