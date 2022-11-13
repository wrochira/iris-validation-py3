from multiprocessing import Process, Queue

import clipper
import gemmi

from iris_validation.utils import ONE_LETTER_CODES
from iris_validation.metrics.residue import MetricsResidue
from iris_validation.metrics.chain import MetricsChain
from iris_validation.metrics.model import MetricsModel
from iris_validation.metrics.series import MetricsModelSeries
from iris_validation.metrics.reflections import ReflectionsHandler


def _get_minimol_from_path(model_path):
    fpdb = clipper.MMDBfile()
    minimol = clipper.MiniMol()
    try:
        fpdb.read_file(model_path)
        fpdb.import_minimol(minimol)
    except Exception as exception:
        raise Exception('Failed to import model file') from exception
    return minimol


def _get_minimol_seq_nums(minimol):
    seq_nums = { }
    for chain in minimol:
        chain_id = str(chain.id()).strip()
        seq_nums[chain_id] = [ ]
        for residue in chain:
            seq_num = int(residue.seqnum())
            seq_nums[chain_id].append(seq_num)
    return seq_nums


def _get_reflections_data(model_path, reflections_path, model_id=None, out_queue=None):
    minimol = _get_minimol_from_path(model_path)
    reflections_handler = ReflectionsHandler(reflections_path, minimol=minimol)
    resolution = reflections_handler.resolution_limit
    density_scores = reflections_handler.calculate_all_density_scores()
    reflections_data = (resolution, density_scores)
    if out_queue is not None:
        out_queue.put(('reflections', model_id, reflections_data))
    return reflections_data


def _get_molprobity_data(model_path, seq_nums, model_id=None, out_queue=None):
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

    for chain_id, chain_seq_nums in seq_nums.items():
        molprobity_data[chain_id] = { }
        for seq_num in chain_seq_nums:
            molprobity_data[chain_id][seq_num] = { category : None for category in molprobity_results }
            molprobity_data[chain_id][seq_num]['clash'] = 2


    for category, results in molprobity_results.items():
        for result in results:
            if category == 'clash':
                for atom in result.atoms_info:
                    chain_id = atom.chain_id.strip()
                    seq_num = int(atom.resseq.strip())
                    if molprobity_data[chain_id][seq_num][category] > 0:
                        molprobity_data[chain_id][seq_num][category] -= 1
                details_line = [ ' '.join(a.id_str().split()) for a in result.atoms_info ] + [ result.overlap ]
                molprobity_data['model_wide']['details'][category].append(details_line)
                continue

            chain_id = result.chain_id.strip()
            seq_num = int(result.resseq.strip())
            if category in ('ramachandran', 'rotamer'):
                if result.score < 0.3:
                    molprobity_data[chain_id][seq_num][category] = 0
                elif result.score < 2.0:
                    molprobity_data[chain_id][seq_num][category] = 1
                else:
                    molprobity_data[chain_id][seq_num][category] = 2
            else:
                if result.outlier:
                    chain_id = result.chain_id.strip()
                    seq_num = int(result.resseq.strip())
                    molprobity_data[chain_id][seq_num][category] = 0

            if result.outlier:
                score = result.deviation if category == 'c-beta' else result.score
                details_line = [ result.chain_id.strip(), result.resid.strip(), result.resname.strip(), score ]
                molprobity_data['model_wide']['details'][category].append(details_line)

    if out_queue is not None:
        out_queue.put(('molprobity', model_id, molprobity_data))

    return molprobity_data


def _get_covariance_data(model_path,
                         sequence_path,
                         distpred_path,
                         seq_nums,
                         distpred_format='rosettanpz',
                         map_align_exe='map_align',
                         dssp_exe='mkdssp',
                         model_id=None,
                         out_queue=None):
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
        from conkit import applications, command_line, io, plot
    except (ImportError, ModuleNotFoundError):
        print('WARNING: Failed to import Biopython; continuing without covariance analyses')
        return

    parser = PDBParser()
    structure = parser.get_structure('structure', model_path)[0]
    dssp = DSSP(structure, model_path, dssp=dssp_exe, acc_array='Wilke')
    model = io.read(model_path, 'pdb' if model_path.endswith('.pdb') else 'mmcif').top
    prediction = io.read(distpred_path, distpred_format).top
    sequence = io.read(sequence_path, 'fasta').top
    figure = plot.ModelValidationFigure(model, prediction, sequence, dssp, map_align_exe=map_align_exe)

    covariance_data = { }
    for chain_id, chain_seq_nums in seq_nums.items():
        covariance_data[chain_id] = { }
        for seq_num in chain_seq_nums:
            # TODO: by chain
            score = figure.smooth_scores[seq_num] if 0 < seq_num < len(figure.smooth_scores) else None
            alignment = 0 if seq_num in figure.alignment.keys() else 1
            covariance_data[chain_id][seq_num] = (score, alignment)

    if out_queue is not None:
        out_queue.put(('covariance', model_id, covariance_data))

    return covariance_data


def metrics_model_series_from_files(model_paths,
                                    reflections_paths=None,
                                    sequence_paths=None,
                                    distpred_paths=None,
                                    run_covariance=False,
                                    run_molprobity=False,
                                    multiprocessing=True):
    try:
        if isinstance(model_paths, str):
            model_paths = [ model_paths ]
        model_paths = tuple(model_paths)
        if None in model_paths:
            raise TypeError
    except TypeError as exception:
        raise ValueError('Argument \'model_paths\' should be an iterable of filenames') from exception

    path_lists = [ model_paths, reflections_paths, sequence_paths, distpred_paths ]
    for i in range(1, len(path_lists)):
        if path_lists[i] is None:
            path_lists[i] = tuple([ None for _ in model_paths ])
        if len(path_lists[i]) != len(model_paths) or \
           path_lists[i].count(None) not in (0, len(path_lists[i])):
            raise ValueError('Path arguments should be equal-length iterables of filenames')

    all_minimol_data = [ ]
    all_covariance_data = [ ]
    all_molprobity_data = [ ]
    all_reflections_data = [ ]
    num_queued = 0
    results_queue = Queue()
    for model_id, file_paths in enumerate(zip(*path_lists)):
        model_path, reflections_path, sequence_path, distpred_path = file_paths
        minimol = _get_minimol_from_path(model_path)
        seq_nums = _get_minimol_seq_nums(minimol)
        covariance_data = None
        molprobity_data = None
        reflections_data = None
        if run_covariance:
            if multiprocessing:
                p = Process(target=_get_covariance_data,
                            args=(model_path, sequence_path, distpred_path, seq_nums),
                            kwargs={ 'model_id': model_id,
                                     'out_queue': results_queue })
                p.start()
                num_queued += 1
            else:
                covariance_data = _get_covariance_data(model_path, sequence_path, distpred_path)
        if run_molprobity:
            if multiprocessing:
                p = Process(target=_get_molprobity_data,
                            args=(model_path, seq_nums),
                            kwargs={ 'model_id': model_id,
                                     'out_queue': results_queue })
                p.start()
                num_queued += 1
            else:
                molprobity_data = _get_molprobity_data(model_path, seq_nums)
        if reflections_path is not None:
            if multiprocessing:
                p = Process(target=_get_reflections_data,
                            args=(model_path, reflections_path),
                            kwargs={ 'model_id': model_id,
                                     'out_queue': results_queue })
                p.start()
                num_queued += 1
            else:
                reflections_data = _get_reflections_data(model_path, reflections_path)

        all_minimol_data.append(minimol)
        all_covariance_data.append(covariance_data)
        all_molprobity_data.append(molprobity_data)
        all_reflections_data.append(reflections_data)

    if multiprocessing:
        for _ in range(num_queued):
            result_type, model_id, result = results_queue.get()
            if result_type == 'covariance':
                all_covariance_data[model_id] = result
            if result_type == 'molprobity':
                all_molprobity_data[model_id] = result
            elif result_type == 'reflections':
                all_reflections_data[model_id] = result

    metrics_models = [ ]
    for model_id, model_data in enumerate(zip(all_minimol_data, all_covariance_data, all_molprobity_data, all_reflections_data)):
        metrics_model = MetricsModel(*model_data)
        metrics_models.append(metrics_model)

    metrics_model_series = MetricsModelSeries(metrics_models)
    return metrics_model_series
