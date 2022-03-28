import os

from iris_validation._defs import CONTINUOUS_METRICS


DATA_DIR_PATH = os.path.join(os.path.dirname(__file__), 'data')
RESOLUTION_BINS_PATH = os.path.join(DATA_DIR_PATH, 'resolution_bins.csv')
PERCENTILES_DATA_PATH = os.path.join(DATA_DIR_PATH, 'percentiles_data.csv')
RESOLUTION_BIN_NAMES = ('<10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '>90', 'All')


class PercentileCalculator():
    def __init__(self, resolution=None):
        self.resolution = resolution
        self.percentile_data = None
        self.resolution_bins = None
        self.bin_name = None
        self._load_data()

    def _load_data(self):
        self.percentile_data = { }
        with open(PERCENTILES_DATA_PATH, 'r', encoding='utf8') as infile:
            for i, line in enumerate(infile.readlines()):
                splitline = line.strip().split(',')
                if i == 0:
                    metric_names = splitline[2:]
                    for metric_name in metric_names:
                        self.percentile_data[metric_name] = { }
                        for bin_name in RESOLUTION_BIN_NAMES:
                            self.percentile_data[metric_name][bin_name] = { }
                else:
                    bin_name = splitline[0]
                    percentile = int(splitline[1])
                    metric_values = [ float(x) for x in splitline[2:] ]
                    for metric_name, metric_value in zip(metric_names, metric_values):
                        self.percentile_data[metric_name][bin_name][percentile] = metric_value

        self.resolution_bins = { }
        with open(RESOLUTION_BINS_PATH, 'r', encoding='utf8') as infile:
            infile.readline() # Skip header line
            for line in infile.readlines():
                splitline = line.split(',')
                percentile = int(splitline[0])
                threshold = float(splitline[1])
                self.resolution_bins[percentile] = threshold

        if self.resolution is None:
            return 'All'
        bin_id = 9
        for i, percentile in enumerate(sorted(self.resolution_bins.keys())):
            percentile_resolution = self.resolution_bins[percentile]
            if self.resolution < percentile_resolution:
                bin_id = i
                break
        self.bin_name = RESOLUTION_BIN_NAMES[bin_id]

    def get_percentile(self, metric_id, metric_value, normalise_polarity=True):
        if None in (metric_id, metric_value):
            return None
        metric_name = CONTINUOUS_METRICS[metric_id]['long_name']
        metric_polarity = CONTINUOUS_METRICS[metric_id]['polarity']
        determined_percentile = 100
        for percentile, percentile_value in self.percentile_data[metric_name][self.bin_name].items():
            if metric_value < percentile_value:
                determined_percentile = percentile
                break
        if normalise_polarity and metric_polarity == -1:
            return 101 - determined_percentile
        return determined_percentile
