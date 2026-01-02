import pandas as pd
import numpy as np


# 0.0から1.0までを5つの要素で等分割
list_numpy_linspace = np.linspace(0.0, 1.0, num=5).tolist()
print(list_numpy_linspace)

# 0.0から1.0までを0.1刻みで生成 (1.0は含まれない)
list_numpy_arange = np.arange(0.0, 1.0, 0.1).tolist()
print(list_numpy_arange)



time_list = list_pure_python = [i / 10.0 for i in range(11)]

df = pd.DataFrame({
    # NG: 'time': list(float(range(0, 100)))*1/100.0
    'time': list_numpy_arange
}, index = list(range(0, len(list_numpy_arange))))
time_start = 0.01
time_end = 0.70
# NG: cond = (df["time"] => time_start) & (df["time"] <= time_end)
cond = (df["time"] >= time_start) & (df["time"] <= time_end)
new_df = df.loc[cond]

print("new_df: ")
print(new_df)


print("END.")




