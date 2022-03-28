import os

from iris_validation.graphics import Panel
from iris_validation.metrics import metrics_model_series_from_files


def generate_report(latest_model_path,
                    previous_model_path=None,
                    latest_reflections_path=None,
                    previous_reflections_path=None,
                    run_molprobity=False,
                    wrap_in_html=True,
                    output_dir=None):

    model_paths = (previous_model_path, latest_model_path)
    reflections_paths = (previous_reflections_path, latest_reflections_path)
    model_series = metrics_model_series_from_files(model_paths,
                                                   reflections_paths,
                                                   run_molprobity)
    model_series_data = model_series.get_raw_data()
    panel = Panel(model_series_data)
    panel_string = panel.dwg.tostring()

    if wrap_in_html:
        panel_string = f'<!DOCTYPE html>\n<html lang="en">\n\t<head>\n\t\t<meta charset="utf-8">\n\t\t<title>Iris Report</title>\n\t</head>\n\t<body>\n\t\t{panel_string}\n\t</body>\n</html>'

    if output_dir is None:
        return panel_string

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    extension = 'html' if wrap_in_html else 'svg'
    with open(os.path.join(output_dir, f'report.{extension}'), 'w', encoding='utf8') as outfile:
        outfile.write(panel_string)
