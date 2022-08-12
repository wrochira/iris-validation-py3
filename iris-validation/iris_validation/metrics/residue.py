from math import isnan

import clipper

from iris_validation import utils
from iris_validation._defs import RAMACHANDRAN_THRESHOLDS


class MetricsResidue():
    def __init__(self, mmol_residue, index_in_chain=None, previous_residue=None, next_residue=None, parent_chain=None, covariance_data=None, molprobity_data=None, density_scores=None):
        self.minimol_residue = mmol_residue
        self.initialised_with_context = index_in_chain is not None
        self.index_in_chain = index_in_chain
        self.previous_residue = previous_residue
        self.next_residue = next_residue
        self.parent_chain = parent_chain
        self.covariance_data = covariance_data
        self.molprobity_data = molprobity_data
        self.density_scores = density_scores

        self.atoms = list(mmol_residue)
        self.sequence_number = int(mmol_residue.seqnum())
        self.code = mmol_residue.type().trim()
        self.code_type = utils.code_type(mmol_residue)
        self.backbone_atoms = utils.get_backbone_atoms(mmol_residue)
        self.backbone_atoms_are_correct = None not in self.backbone_atoms
        self.backbone_geometry_is_correct = utils.check_backbone_geometry(mmol_residue) if self.backbone_atoms_are_correct else None
        self.is_aa = utils.check_is_aa(mmol_residue)
        self.is_water = str(mmol_residue.type()).strip() == 'HOH'
        self.is_consecutive_aa = None

        # B-factors
        self.max_b_factor, self.avg_b_factor, self.std_b_factor, self.mc_b_factor, self.sc_b_factor = utils.analyse_b_factors(mmol_residue, self.is_aa, self.backbone_atoms)

        # Backbone torsion angles
        self.phi = clipper.MMonomer.protein_ramachandran_phi(self.previous_residue, mmol_residue) if self.previous_residue else None
        self.psi = clipper.MMonomer.protein_ramachandran_psi(mmol_residue, self.next_residue) if self.next_residue else None
        if self.phi is not None and isnan(self.phi):
            self.phi = None
        if self.psi is not None and isnan(self.psi):
            self.psi = None

        # Side chain torsion angles
        self.chis = utils.calculate_chis(mmol_residue) if self.is_aa else None
        self.is_sidechain_complete = self.chis is not None and None not in self.chis

        # Ramachandran
        self.ramachandran_score = utils.calculate_ramachandran_score(mmol_residue, self.code, self.phi, self.psi)
        self.ramachandran_flags = (None, None, None)
        if self.ramachandran_score is not None:
            if RAMACHANDRAN_THRESHOLDS[0] <= self.ramachandran_score:
                self.ramachandran_flags = (True, False, False)
            elif RAMACHANDRAN_THRESHOLDS[1] <= self.ramachandran_score < RAMACHANDRAN_THRESHOLDS[0]:
                self.ramachandran_flags = (False, True, False)
            elif self.ramachandran_score < RAMACHANDRAN_THRESHOLDS[1]:
                self.ramachandran_flags = (False, False, True)
        self.ramachandran_favoured, self.ramachandran_allowed, self.ramachandran_outlier = self.ramachandran_flags

        # Rotamer
        rotamer_calculator = self.parent_chain.parent_model.rotamer_calculator
        self.rotamer_score = None
        self.rotamer_flags = (None, None, None)
        if self.is_sidechain_complete:
            self.rotamer_score = rotamer_calculator.get_cv_score(self.code, self.chis)
            rotamer_clf_id = rotamer_calculator.get_classification(self.code, self.chis)
            if rotamer_clf_id == 3:
                self.rotamer_flags = (True, False, False)
            elif rotamer_clf_id == 2:
                self.rotamer_flags = (False, True, False)
            elif rotamer_clf_id in (0, 1):
                self.rotamer_flags = (False, False, True)
        self.rotamer_favoured, self.rotamer_allowed, self.rotamer_outlier = self.rotamer_flags

        # MolProbity data
        self.discrete_indicators = self.molprobity_data
        if self.molprobity_data is None:
            ramachandran_indicator = 0 if self.ramachandran_outlier else \
                                     1 if self.ramachandran_allowed else \
                                     2 if self.ramachandran_favoured else None
            rotamer_indicator = 0 if self.rotamer_outlier else \
                                1 if self.rotamer_allowed else \
                                2 if self.rotamer_favoured else None
            self.discrete_indicators = { 'clash' : None,
                                         'c-beta' : None,
                                         'omega' : None,
                                         'ramachandran' : ramachandran_indicator,
                                         'rotamer' : rotamer_indicator }

        # Covariance data
        self.covariance_score, self.cmo_string = None, None
        if self.covariance_data is not None:
            self.covariance_score, self.cmo_string = self.covariance_data
        self.discrete_indicators['cmo'] = self.cmo_string

        # Density fit scores
        self.fit_score, self.mainchain_fit_score, self.sidechain_fit_score = None, None, None
        if self.density_scores is not None:
            self.fit_score, self.mainchain_fit_score, self.sidechain_fit_score = self.density_scores

        # Percentiles
        percentile_calculator = self.parent_chain.parent_model.percentile_calculator
        self.avg_b_factor_percentile = percentile_calculator.get_percentile(0, self.avg_b_factor)
        self.max_b_factor_percentile = percentile_calculator.get_percentile(1, self.max_b_factor)
        self.std_b_factor_percentile = percentile_calculator.get_percentile(2, self.std_b_factor)
        self.fit_score_percentile = percentile_calculator.get_percentile(3, self.fit_score)
        self.mainchain_fit_score_percentile = percentile_calculator.get_percentile(4, self.mainchain_fit_score)
        self.sidechain_fit_score_percentile = percentile_calculator.get_percentile(5, self.sidechain_fit_score)
        self.covariance_score_percentile = percentile_calculator.get_percentile(6, self.covariance_score)
