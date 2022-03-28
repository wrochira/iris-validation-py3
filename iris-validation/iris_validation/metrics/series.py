import clipper

from iris_validation import utils
from iris_validation.metrics.model import MetricsModel
from iris_validation.metrics.reflections import ReflectionsHandler


class MetricsModelSeries():
    def __init__(self, metrics_models):
        self.metrics_models = metrics_models
        self.chain_sets = None
        self.chain_alignments = None

        for metrics_model in self.metrics_models:
            metrics_model.parent_series = self

    def align_models(self):
        if len(self.metrics_models) == 0:
            return
        if len(self.metrics_models) > 2:
            raise Exception('Iris currently only supports alignment for two model iterations')

        # Check for and remove chains with no amino acid residues
        bad_chain_ids = set()
        for model in self.metrics_models:
            for chain in model:
                if chain.length == 0:
                    bad_chain_ids.add(chain.chain_id)
        if len(bad_chain_ids) > 0:
            print('WARNING: at least one chain contains no amino acid residues. Ignoring chains: ' + ', '.join(sorted(bad_chain_ids)))
            for model in self.metrics_models:
                for chain_id in bad_chain_ids:
                    model.remove_chain(chain_id)
        if 0 in [ model.chain_count for model in self.metrics_models ]:
            raise Exception('One or more models had no valid chains')

        # Align chains
        chain_id_sets = [ set(chain.chain_id for chain in model) for model in self.metrics_models ]
        common_chain_ids = set.intersection(*chain_id_sets)
        lost_chain_ids = set()
        for model, chain_id_set in zip(self.metrics_models, chain_id_sets):
            model_lost_chain_ids = chain_id_set - common_chain_ids
            lost_chain_ids.update(model_lost_chain_ids)
            if len(model_lost_chain_ids) > 0:
                for chain_id in model_lost_chain_ids:
                    model.remove_chain(chain_id)
        if len(lost_chain_ids) > 0:
            print(f'WARNING: Some chains are not present or valid across all model versions ({sorted(lost_chain_ids)}). These chains will not be represented in the validation report.')

        # Chain sets
        self.chain_sets = { }
        for chain_id in sorted(common_chain_ids):
            self.chain_sets[chain_id] = [ ]
            for model in self.metrics_models:
                matching_chain = [ chain for chain in model if chain.chain_id == chain_id ][0]
                self.chain_sets[chain_id].append(matching_chain)

        # Align residues
        self.chain_alignments = { }
        for chain_id, chain_set in self.chain_sets.items():
            sequences = [ utils.code_three_to_one([ residue.code for residue in chain ]) for chain in chain_set ]
            if len(sequences) == 1:
                self.chain_alignments[chain_id] = (sequences[0], )
                continue
            alignment_pair = utils.needleman_wunsch(sequences[-2], sequences[-1])
            self.chain_alignments[chain_id] = alignment_pair

    def get_raw_data(self):
        if self.chain_alignments is None:
            self.align_models()

        num_versions = len(self.metrics_models)
        has_molprobity = self.metrics_models[0].molprobity_data is not None

        raw_data = [ ]
        for chain_id, chain_set in self.chain_sets.items():
            alignment_strings = self.chain_alignments[chain_id]
            aligned_length = len(alignment_strings[0])
            chain_data = { 'chain_id'           : chain_id,
                           'num_versions'       : num_versions,
                           'has_molprobity'     : has_molprobity,
                           'aligned_length'     : aligned_length,
                           'residue_seqnos'     : [ ],
                           'residue_codes'      : [ ],
                           'residue_validities' : [ ],
                           'discrete_values'    : [ ],
                           'continuous_values'  : [ ],
                           'percentile_values'  : [ ] }

            for alignment_string, chain in zip(alignment_strings, chain_set):
                residue_seqnos = [ ]
                residue_codes = [ ]
                residue_validities = [ ]
                discrete_values = [ ]
                continuous_values = [ ]
                percentile_values = [ ]

                residue_id = -1
                for alignment_char in alignment_string:
                    if alignment_char == '-':
                        residue_seqnos.append(None)
                        residue_codes.append(None)
                        residue_validities.append(False)
                        discrete_values.append(tuple(None for _ in range(3)))
                        continuous_values.append(tuple(None for _ in range(6)))
                        percentile_values.append(tuple(None for _ in range(6)))
                        continue

                    residue_id += 1
                    residue = chain.residues[residue_id]
                    residue_seqnos.append(residue.sequence_number)
                    residue_codes.append(residue.code)
                    residue_validities.append(True)

                    residue_discrete_values = (residue.discrete_indicators['rotamer'],
                                               residue.discrete_indicators['ramachandran'],
                                               residue.discrete_indicators['clash'])
                    residue_continuous_values = (residue.avg_b_factor,
                                                 residue.max_b_factor,
                                                 residue.std_b_factor,
                                                 residue.fit_score,
                                                 residue.mainchain_fit_score,
                                                 residue.sidechain_fit_score)
                    residue_percentile_values = (residue.avg_b_factor_percentile,
                                                 residue.max_b_factor_percentile,
                                                 residue.std_b_factor_percentile,
                                                 residue.fit_score_percentile,
                                                 residue.mainchain_fit_score_percentile,
                                                 residue.sidechain_fit_score_percentile)

                    residue_continuous_values = tuple(round(x, 3) if isinstance(x, float) else x for x in residue_continuous_values)
                    discrete_values.append(residue_discrete_values)
                    continuous_values.append(residue_continuous_values)
                    percentile_values.append(residue_percentile_values)

                discrete_values = list(zip(*discrete_values))
                continuous_values = list(zip(*continuous_values))
                percentile_values = list(zip(*percentile_values))
                chain_data['residue_seqnos'].append(residue_seqnos)
                chain_data['residue_codes'].append(residue_codes)
                chain_data['residue_validities'].append(residue_validities)
                chain_data['discrete_values'].append(discrete_values)
                chain_data['continuous_values'].append(continuous_values)
                chain_data['percentile_values'].append(percentile_values)

            chain_data['discrete_values'] = list(zip(*chain_data['discrete_values']))
            chain_data['continuous_values'] = list(zip(*chain_data['continuous_values']))
            chain_data['percentile_values'] = list(zip(*chain_data['percentile_values']))
            raw_data.append(chain_data)

        return raw_data
