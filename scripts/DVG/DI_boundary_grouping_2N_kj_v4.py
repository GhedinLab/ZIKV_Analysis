import argparse
import pandas as pd


parser = argparse.ArgumentParser(description='Group 2N DI boundary records.')
parser.add_argument('--input_file', required=True, help='input file')
parser.add_argument('--output_file', required=True, help='output file')
args = parser.parse_args()


class Group:
    def __init__(self, row):
        self.rows = [row]

    def get_max(self, key, rows=None):
        if rows is None:
            rows = self.rows
        return max([row[key] for row in rows])

    def get_min(self, key, rows=None):
        if rows is None:
            rows = self.rows
        return min([row[key] for row in rows])

    def add(self, new_row):
        new_rows = self.rows + [new_row]

        #if (self.get_max('boundary5p', new_rows) - self.get_min('boundary5p', new_rows)) + (self.get_max('midM_5p', new_rows) - self.get_min('midM_5p', new_rows)) > 10:
        #    return False

        #if (self.get_max('midM_3p', new_rows) - self.get_min('midM_3p', new_rows)) + (self.get_max('boundary3p', new_rows) - self.get_min('boundary3p', new_rows)) > 10:
        #    return False

        if (self.get_max('gap_start1', new_rows) - self.get_min('gap_start1', new_rows)) + (self.get_max('gap_end1', new_rows) - self.get_min('gap_end1', new_rows)) > 10:
            return False

        if (self.get_max('gap_start2', new_rows) - self.get_min('gap_start2', new_rows)) + (self.get_max('gap_end2', new_rows) - self.get_min('gap_end2', new_rows)) > 10:
            return False

        self.rows = new_rows
        return True

    @property
    def df(self):
        return pd.DataFrame(self.rows)


if __name__ == '__main__':
    data = pd.read_csv(args.input_file, keep_default_na=False)
    print(data)
    data.sort_values('freq', inplace=True, ascending=False)
    data.reset_index(drop=True, inplace=True)
    segments = list(set(list(data['segment'])))

    masterDF=pd.DataFrame()

    for seg in segments:
        print(seg)
        groups = []
        df_seg = data[data.segment==seg]
        print(df_seg.shape)

        for i, row in data.iterrows():
            grouped = False
            for group in groups:
                if group.add(row):
                    grouped = True
                    break
            if not grouped:
                groups.append(Group(row))

        dfs = []
        for i, group in enumerate(groups):
            df = group.df
            df['group'] = i + 1

            min_start1 = df['gap_start1'].min()
            max_start1 = df['gap_start1'].max()
            min_end1 = df['gap_end1'].min()
            max_end1 = df['gap_end1'].max()

            min_start2 = df['gap_start2'].min()
            max_start2 = df['gap_start2'].max()
            min_end2 = df['gap_end2'].min()
            max_end2 = df['gap_end2'].max()

            #max_row=df[df.freq==df.freq.max()]
            #newgap1 = max_row['gap1'].iloc[0]
            #newgap2 = max_row['gap2'].iloc[0]

            df['NewGap1'] = list(df[df.freq==df.freq.max()]['gap1'])[0]
            df['NewGap2'] = list(df[df.freq==df.freq.max()]['gap2'])[0]


            df['NewStart1'] = list(df[df.freq==df.freq.max()]['gap_start1'])[0]
            df['NewEnd1'] = list(df[df.freq==df.freq.max()]['gap_end1'])[0]

            df['NewStart2'] = list(df[df.freq==df.freq.max()]['gap_start2'])[0]
            df['NewEnd2'] = list(df[df.freq==df.freq.max()]['gap_end2'])[0]

            df['GapSize1'] = list(df[df.freq==df.freq.max()]['gap_size1'])[0]
            df['GapSize2'] = list(df[df.freq==df.freq.max()]['gap_size2'])[0]

            df['gap1_boundaries'] = '{0}-{1}_{2}-{3}'.format(min_start1,max_start1,min_end1,max_end1)
            df['gap2_boundaries'] = '{0}-{1}_{2}-{3}'.format(min_start2,max_start2,min_end2,max_end2)

            #df['NewGap1'] = newgap1
            #df['NewGap2'] = newgap2
            dfs.append(df)

        segmaster = pd.concat(dfs, ignore_index=True)
        masterDF = masterDF.append(segmaster,ignore_index=True)

    masterDF.to_csv(args.output_file, index=False)
