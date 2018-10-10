import sys
import datetime
import numpy as np
import pandas as pd

def compare_two_dfs(ref, test):
    ref.set_index('Index', inplace=True)
    test.set_index('Index', inplace=True)
    df_all = pd.concat([ref, test],
        axis='columns', keys=['v1', 'v2'])
    df_final = df_all.swaplevel(axis='columns')[ref.columns[1:]]
    def styleIT(data):
        attr = 'background-color: {}'.format('yellow')
        other = data.xs('v1', axis='columns', level=-1)
        return(pd.DataFrame(np.where(data.ne(other, level=0), attr, ''), index=data.index, columns=data.columns))
    df_final = df_final.style.apply(styleIT, axis=None)
    return df_final

if __name__ == '__main__':
    ref_frame = pd.read_table('MaRSFrame.csv', sep=',', usecols=[0,1,2,3,4,5,8,10])
    ref_frame['Index'] = ref_frame['Sample'] + ':' + ref_frame['Variant']
    new_frame = pd.read_table(sys.argv[1], sep=',', usecols=[0,1,2,3,4,5,8,10])
    new_frame['Index'] = ref_frame['Sample'] + ':' + ref_frame['Variant']
    time = datetime.datetime.now()
    comp_frame = 'MaRSvsNew{0}{1}{2}{3}{4}{5}.xlsx'.format(time.year, time.month, time.day, time.hour, time.minute, time.second)
    res = compare_two_dfs(ref_frame, new_frame)
    res.to_excel(comp_frame)
