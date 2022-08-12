from iris_validation.metrics.chain import MetricsChain
from iris_validation.metrics.rotamer import RotamerCalculator
from iris_validation.metrics.percentiles import PercentileCalculator


class MetricsModel():
    def __init__(self, mmol_model, covariance_data=None, molprobity_data=None, reflections_data=None):
        self.minimol_model = mmol_model
        self.covariance_data = covariance_data
        self.molprobity_data = molprobity_data
        self.reflections_data = reflections_data

        self._index = -1
        self.minimol_chains = list(mmol_model.model())
        self.chain_count = len(self.minimol_chains)

        self.resolution, self.density_scores = None, None
        if reflections_data is not None:
            self.resolution, self.density_scores = reflections_data
        self.percentile_calculator = PercentileCalculator(self.resolution)
        self.rotamer_calculator = RotamerCalculator()

        self.chains = [ ]
        for mmol_chain in mmol_model:
            chain_id = str(mmol_chain.id().trim())
            chain_covariance_data = None if covariance_data is None else covariance_data[chain_id]
            chain_molprobity_data = None if molprobity_data is None else molprobity_data[chain_id]
            chain_density_scores = None if self.density_scores is None else self.density_scores[chain_id]
            chain = MetricsChain(mmol_chain, self, chain_covariance_data, chain_molprobity_data, chain_density_scores)
            chain.remove_non_aa_residues()
            self.chains.append(chain)

    def __iter__(self):
        return self

    def __next__(self):
        if self._index < len(self.chains)-1:
            self._index += 1
            return self.chains[self._index]
        self._index = -1
        raise StopIteration

    def get_chain(self, chain_id):
        return next(chain for chain in self.chains if chain.chain_id == chain_id)

    def remove_chain(self, chain_id):
        matching_chains = [ chain for chain in self.chains if chain.chain_id == chain_id ]
        if len(matching_chains) == 0:
            print('Error removing chain, no chains matching that ID were found.')
        else:
            for chain in matching_chains:
                self.chains.remove(chain)
                self.chain_count -= 1

    def b_factor_lists(self):
        all_bfs, aa_bfs, mc_bfs, sc_bfs, non_aa_bfs, water_bfs, ligand_bfs, ion_bfs = [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ]
        for chain in self.chains:
            all_bfs_c, aa_bfs_c, mc_bfs_c, sc_bfs_c, non_aa_bfs_c, water_bfs_c, ligand_bfs_c, ion_bfs_c = chain.b_factor_lists()
            for model_li, chain_li in zip((all_bfs, aa_bfs, mc_bfs, sc_bfs, non_aa_bfs, water_bfs, ligand_bfs, ion_bfs),
                                          (all_bfs_c, aa_bfs_c, mc_bfs_c, sc_bfs_c, non_aa_bfs_c, water_bfs_c, ligand_bfs_c, ion_bfs_c)):
                model_li += chain_li
        return all_bfs, aa_bfs, mc_bfs, sc_bfs, non_aa_bfs, water_bfs, ligand_bfs, ion_bfs
