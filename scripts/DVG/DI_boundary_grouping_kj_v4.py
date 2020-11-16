# edited to work with KJ data
import argparse
import pandas as pd


parser = argparse.ArgumentParser(description='Group 1N DI boundary records.')
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

        #if (self.get_max('boundary5p', new_rows) - self.get_min('boundary5p', new_rows)) + (self.get_max('boundary3p', new_rows) - self.get_min('boundary3p', new_rows)) > 10:
        if (self.get_max('gap_start', new_rows) - self.get_min('gap_start', new_rows)) + (self.get_max('gap_end', new_rows) - self.get_min('gap_end', new_rows)) > 10:
            return False

        self.rows = new_rows
        return True

    @property
    def df(self):
        return pd.DataFrame(self.rows)


if __name__ == '__main__':

    data = pd.read_csv(args.input_file,keep_default_na = False)
    data.sort_values('freq', inplace=True, ascending=False)
    data.reset_index(drop=True, inplace=True)
    segments=(list(set(list(data['segment']))))
    #segments=['HA']

    masterDF=pd.DataFrame()
    for seg in segments:
        groups = []
        print(seg)
        df_seg = data[data.segment==seg]
        print(df_seg.shape)
        for i, row in df_seg.iterrows():
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
            min_start = df['gap_start'].min()
            max_start = df['gap_start'].max()
            min_end = df['gap_end'].min()
            max_end = df['gap_end'].max()

            df['NewGap'] = list(df[df.freq==df.freq.max()]['gap'])[0]
            df['NewStart'] = list(df[df.freq==df.freq.max()]['gap_start'])[0]
            df['NewEnd'] = list(df[df.freq==df.freq.max()]['gap_end'])[0]
            df['GapSize'] = list(df[df.freq==df.freq.max()]['gap_size'])[0]
            df['EstimatedLength'] = list(df[df.freq==df.freq.max()]['estimated_length'])[0]
            df['ORF_Flag'] = list(df[df.freq==df.freq.max()]['ORF_flag'])[0]


            df['boundaries'] = '{0}-{1}_{2}-{3}'.format(min_start,max_start,min_end,max_end)
            dfs.append(df)

        segmaster = pd.concat(dfs, ignore_index=True)
        masterDF = masterDF.append(segmaster,ignore_index=True)

    masterDF.to_csv(args.output_file, index=False)
