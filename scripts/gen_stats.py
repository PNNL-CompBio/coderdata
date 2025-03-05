import coderdata as cd
import yaml

class DatasetStatistics:
    def __init__(self, dataset_type):
        self.data = cd.load(dataset_type)

    def count_unique(self, attribute, unique_field):
        if getattr(self.data, attribute) is not None:
            dataset = getattr(self.data, attribute)
            if unique_field in dataset.columns:
                return len(dataset[unique_field].unique())
        return 0

    def count_unique_genes(self):
        gene_ids = set()
        for data_type in ['transcriptomics', 'proteomics', 'mutations', 'copy_number']:
            if getattr(self.data, data_type) is not None:
                dataset = getattr(self.data, data_type)
                if 'entrez_id' in dataset.columns:
                    gene_ids.update(dataset.entrez_id.unique().tolist())
        return len(gene_ids)

def calculate_stats_for_datasets(dataset_types):
    stats = {}
    for dataset_type in dataset_types:
        dataset_stats = DatasetStatistics(dataset_type)
        stats[dataset_type] = {
            'samples': dataset_stats.count_unique('samples', 'improve_sample_id'),
            'genes': dataset_stats.count_unique_genes(),
        }
        # Count drugs if the dataset has a drugs data type
        stats[dataset_type]['drugs'] = dataset_stats.count_unique('drugs', 'improve_drug_id')
        # Count cancer types or cell lines based on dataset type
        if dataset_type in ['cptac', 'hcmi', 'beataml', 'mpnst']:
            stats[dataset_type]['cancer_types'] = dataset_stats.count_unique('samples', 'cancer_type')
        if dataset_type == 'broad_sanger':
            stats[dataset_type]['cell_lines'] = dataset_stats.count_unique('samples', 'cancer_type')

    with open('stats.yml', 'w') as file:
        yaml.dump(stats, file)

# Dataset types
dataset_types = cd.list_datasets(raw=True).keys()
calculate_stats_for_datasets(dataset_types)
