from math import log

import clipper
from scipy.stats import norm

from iris_validation import utils


class ReflectionsHandler():
    def __init__(self, f_reflections=None, xmap=None, minimol=None):
        self.f_reflections = f_reflections
        self.xmap = xmap
        self.minimol = minimol

        self.hkl = clipper.HKL_info()

        self.grid = None
        self.spacegroup = None
        self.cell = None
        self.resolution = None
        self.resolution_limit = None

        if f_reflections is None:
            if xmap is None:
                raise ValueError('Either a reflections file path or an xmap object must be passed as an argument')
            try:
                self.grid = xmap.grid
            except AttributeError:
                self.grid = None
            self.spacegroup = xmap.spacegroup
            self.cell = xmap.cell
        else:
            extension = f_reflections.split('.')[-1].lower()
            if extension != 'mtz':
                if extension == 'cif':
                    raise ValueError('mmCIF format is not currently supported for reflections data.')
                else:
                    raise ValueError(f'Reflections file has unrecognised extension: {extension}')
            self._load_hkl_data()
            self._calculate_structure_factors()
            self._generate_xmap()
            self._calculate_map_stats()

    def _load_hkl_data(self):
        mtzin = clipper.CCP4MTZfile()
        mtzin.open_read(self.f_reflections)
        mtzin.import_hkl_info(self.hkl)

        mtz_labels_and_types = [ tuple(str(line).strip().split(' ')) for line in mtzin.column_labels() ]
        mtz_column_labels, _ = zip(*mtz_labels_and_types)
        mtz_column_label_suffixes = set([ label.split('/')[-1] for label in mtz_column_labels ])
        # Need a better way to choose the right headers, but I'm not familiar enough with reflections data to know what labels are common
        import_complete = False
        for suffix_pair in ( ('F', 'SIGF'),
                             ('FP', 'SIGFP'),
                             ('FP_ALL', 'SIGFP_ALL') ):
            if len(mtz_column_label_suffixes & set(suffix_pair)) == 2:
                try:
                    self.f_sigf = clipper.HKL_data_F_sigF_float(self.hkl)
                    mtzin.import_hkl_data(self.f_sigf, '/*/*/[' + ','.join(suffix_pair) + ']')
                    import_complete = True
                    break
                except Exception as exception:
                    print('Failed to import HKL data from reflections file')
                    raise exception
        if not import_complete:
            raise ValueError('Reflections file does not contain the required columns')
        mtzin.close_read()

        spacegroup = self.hkl.spacegroup()
        cell = self.hkl.cell()
        resolution = self.hkl.resolution()

        self.spacegroup = spacegroup
        self.cell = cell
        self.resolution = resolution
        self.resolution_limit = self.resolution.limit()

    def _calculate_structure_factors(self, bulk_solvent=True):
        #self.crystal = clipper.MTZcrystal()
        #self.f_phi = clipper.HKL_data_F_phi_float(self.hkl, self.crystal)
        self.f_phi = clipper.HKL_data_F_phi_float(self.hkl)
        atoms = self.minimol.atom_list()
        sf_calc = clipper.SFcalc_obs_bulk_float if bulk_solvent else clipper.SFcalc_obs_base_float
        sf_calc(self.f_phi, self.f_sigf, atoms)

    def _generate_xmap(self):
        self.grid = clipper.Grid_sampling(self.spacegroup, self.cell, self.resolution)
        self.xmap = clipper.Xmap_float(self.spacegroup, self.cell, self.grid)
        try:
            self.xmap.fft_from(self.f_phi)
        except AttributeError:
            self.xmap.fft_from_float(self.f_phi)

    def _calculate_map_stats(self):
        map_stats = clipper.Map_stats(self.xmap)
        self.map_mean = map_stats.mean()
        self.map_std = map_stats.std_dev()

    def get_density_at_point(self, xyz):
        cell = self.xmap.cell()
        grid = self.xmap.grid_sampling()
        co = clipper.Coord_orth(*xyz)
        co_cf_cg = co.coord_frac(cell).coord_grid(grid)
        density = self.xmap.get_data(co_cf_cg)
        return density

    def get_density_at_atom(self, mmol_atom):
        co = mmol_atom.coord_orth()
        xyz = (co.x(), co.y(), co.z())
        return self.get_density_at_point(xyz)

    def get_density_scores_at_residue(self, metrics_residue):
        all_atom_scores, mainchain_atom_scores, sidechain_atom_scores = [ ], [ ], [ ]
        for atom_id, atom in enumerate(metrics_residue.minimol_residue):
            is_mainchain = str(atom.name()).strip() in utils.MC_ATOM_NAMES
            element = str(atom.element()).strip()
            atomic_number = utils.ATOMIC_NUMBERS[element]
            density = self.get_density_at_atom(atom)
            atom_score = None
            density_norm = density / atomic_number
            atom_score = -log(norm.cdf((density_norm - self.map_mean) / self.map_std))
            all_atom_scores.append(atom_score)
            if is_mainchain:
                mainchain_atom_scores.append(atom_score)
            else:
                sidechain_atom_scores.append(atom_score)
        all_score, mainchain_score, sidechain_score = None, None, None
        if len(all_atom_scores) > 0:
            all_score = sum(all_atom_scores) / len(all_atom_scores)
        if metrics_residue.is_aa:
            if len(mainchain_atom_scores) > 0:
                mainchain_score = sum(mainchain_atom_scores) / len(mainchain_atom_scores)
            if len(sidechain_atom_scores) > 0:
                sidechain_score = sum(sidechain_atom_scores) / len(sidechain_atom_scores)
        return all_score, mainchain_score, sidechain_score
