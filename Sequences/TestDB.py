import pandas as pd
SeqDB = pd.read_csv(r'C:\ConditionDependentPBEE\GroundmotionSelection\mainshock_aftershock_file_database.csv')

for i, row in SeqDB.iterrows():
    seq_fn = row['horizontal_1_filename']
    dt = row['dt_sequence_horizontal1']
    npt = row['npt_sequence_horizontal1']
    print("seuqence is: ", seq_fn)
    print("with dt: ", dt)
    print("and total points number: ", npt)