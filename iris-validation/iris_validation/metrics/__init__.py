import clipper

from iris_validation.utils import ONE_LETTER_CODES
from iris_validation.metrics.residue import MetricsResidue
from iris_validation.metrics.chain import MetricsChain
from iris_validation.metrics.model import MetricsModel
from iris_validation.metrics.series import MetricsModelSeries
from iris_validation.metrics.reflections import ReflectionsHandler


def metrics_model_series_from_files(model_paths,
                                    reflections_paths=None,
                                    run_molprobity=False):
    try:
        if isinstance(model_paths, str):
            model_paths = [ model_paths ]
        model_paths = tuple(model_paths)
        if None in model_paths:
            raise TypeError
    except TypeError as exception:
        raise ValueError('Argument \'model_paths\' should be an iterable of filenames') from exception

    path_lists = [ model_paths, reflections_paths ]
    for i in range(1, len(path_lists)):
        if path_lists[i] is None:
            path_lists[i] = tuple([ None for _ in model_paths ])
        if len(path_lists[i]) != len(model_paths) or \
           path_lists[i].count(None) not in (0, len(path_lists[i])):
            raise ValueError('Path arguments should be equal-length iterables of filenames')

    metrics_models = [ ]
    for model_args in zip(*path_lists):
        metrics_model = metrics_model_from_files(*model_args, run_molprobity)
        metrics_models.append(metrics_model)

    metrics_model_series = MetricsModelSeries(metrics_models)
    return metrics_model_series


def metrics_model_from_files(model_path,
                             reflections_path=None,
                             run_molprobity=False):
    fpdb = clipper.MMDBfile()
    minimol = clipper.MiniMol()
    try:
        fpdb.read_file(model_path)
        fpdb.import_minimol(minimol)
    except Exception as exception:
        raise Exception('Failed to import model file') from exception
    reflections_handler = None
    if reflections_path is not None:
        reflections_handler = ReflectionsHandler(reflections_path, minimol=minimol)

    molprobity_data = None
    if run_molprobity:
        molprobity_data = get_molprobity_data(model_path, minimol)

    metrics_model = MetricsModel(minimol, reflections_handler, molprobity_data)
    return metrics_model


def get_molprobity_data(model_path, minimol):
    try:
        from mmtbx.command_line import load_model_and_data
        from mmtbx.command_line.molprobity import get_master_phil
        from mmtbx.validation.molprobity import molprobity, molprobity_flags
    except (ImportError, ModuleNotFoundError):
        print('WARNING: Failed to import MolProbity; continuing without MolProbity analyses')
        return

    try:
        cmdline = load_model_and_data(
            args=[ f'pdb.file_name="{model_path}"', 'quiet=True' ],
            master_phil=get_master_phil(),
            require_data=False,
            process_pdb_file=True)
        validation = molprobity(model=cmdline.model)
    except Exception:
        print('WARNING: Failed to run MolProbity; continuing without MolProbity analyses')
        return

    molprobity_data = { }
    molprobity_data['model_wide'] = { }
    molprobity_data['model_wide']['summary'] = { 'cbeta_deviations' : validation.cbetadev.n_outliers,
                                                 'clashscore' : validation.clashscore(),
                                                 'ramachandran_outliers' : validation.rama_outliers(),
                                                 'ramachandran_favoured' : validation.rama_favored(),
                                                 'rms_bonds' : validation.rms_bonds(),
                                                 'rms_angles' : validation.rms_angles(),
                                                 'rotamer_outliers' : validation.rota_outliers(),
                                                 'molprobity_score' : validation.molprobity_score() }

    molprobity_data['model_wide']['details'] = { 'clash' : [ ],
                                                 'c-beta' : [ ],
                                                 'nqh_flips' : [ ],
                                                 'omega' : [ ],
                                                 'ramachandran' : [ ],
                                                 'rotamer' : [ ] }

    molprobity_results = { 'clash' : validation.clashes.results,
                           'c-beta' : validation.cbetadev.results,
                           'nqh_flips' : validation.nqh_flips.results,
                           'omega' : validation.omegalyze.results,
                           'ramachandran' : validation.ramalyze.results,
                           'rotamer' : validation.rotalyze.results }

    for chain in minimol:
        chain_id = str(chain.id()).strip()
        molprobity_data[chain_id] = { }
        for residue in chain:
            seq_no = int(residue.seqnum())
            molprobity_data[chain_id][seq_no] = { category : None for category in molprobity_results }
            molprobity_data[chain_id][seq_no]['clash'] = 2

    for category, results in molprobity_results.items():
        for result in results:
            if category == 'clash':
                for atom in result.atoms_info:
                    chain_id = atom.chain_id.strip()
                    seq_no = int(atom.resseq.strip())
                    if molprobity_data[chain_id][seq_no][category] > 0:
                        molprobity_data[chain_id][seq_no][category] -= 1
                details_line = [ ' '.join(a.id_str().split()) for a in result.atoms_info ] + [ result.overlap ]
                molprobity_data['model_wide']['details'][category].append(details_line)
                continue

            chain_id = result.chain_id.strip()
            seq_no = int(result.resseq.strip())
            if category in ('ramachandran', 'rotamer'):
                if result.score < 0.3:
                    molprobity_data[chain_id][seq_no][category] = 0
                elif result.score < 2.0:
                    molprobity_data[chain_id][seq_no][category] = 1
                else:
                    molprobity_data[chain_id][seq_no][category] = 2
            else:
                if result.outlier:
                    chain_id = result.chain_id.strip()
                    seq_no = int(result.resseq.strip())
                    molprobity_data[chain_id][seq_no][category] = 0

            if result.outlier:
                score = result.deviation if category == 'c-beta' else result.score
                details_line = [ result.chain_id.strip(), result.resid.strip(), result.resname.strip(), score ]
                molprobity_data['model_wide']['details'][category].append(details_line)

    return molprobity_data
