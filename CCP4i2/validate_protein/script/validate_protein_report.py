import shutil

from report.CCP4ReportParser import *


class validate_protein_report(Report):
    TASKNAME = 'validate_protein'

    def __init__(self, xmlnode=None, jobInfo={ }, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)

        self.append('<title>Multimetric Validation Report</title>')

        results = self.addResults()

        if self.xmlnode.select('//Iris') is None:
            self.add_iris_panel()

        if self.xmlnode.select('//MolProbity') is None:
            self.add_molprobity()
        
        if self.xmlnode.select('//B_factors') is None:
            self.add_b_factors()
        
        if self.xmlnode.select('//Ramachandran') is None:
            self.add_ramachandran()


    def add_iris_panel(self):
        reference_div = self.addDiv(style='margin-top:-90px; float:right;')
        reference = '<h4><a style="text-decoration: none;" href="https://doi.org/10.1002/pro.3955">Rochira W, Agirre J (2021)</a></h4>'
        reference_div.append(reference)
        section_div = self.addDiv()
        iris_panel_svg = self.xmlnode.xpath('//Iris/Panel_svg')[0].text
        section_div.append(iris_panel_svg)
        section_div.addDiv(style='clear:both;')
        self.append('<br/><br/>')


    def add_molprobity(self):
        section_div = self.addDiv()
        section_div.append('<p style="font-weight:bold; font-size: 14px;"><b>MolProbity Analyses</b></p>')

        fold = section_div.addFold(label='Summary', initiallyOpen=True)
        table_summary = fold.addTable(title='MolProbity Summary', select='//MolProbity/Summary', transpose=True)
        table_summary.addData(title='Ramachandran outliers', select='Ramachandran_outliers')
        table_summary.addData(title='Ramachandran favoured', select='Ramachandran_favoured')
        table_summary.addData(title='Rotamer outliers', select='Rotamer_outliers')
        table_summary.addData(title='CBeta deviations', select='CBeta_deviations')
        table_summary.addData(title='Clashscore', select='Clashscore')
        table_summary.addData(title='RMS bonds', select='RMS_bonds')
        table_summary.addData(title='RMS angles', select='RMS_angles')
        table_summary.addData(title='MolProbity score', select='MolProbity_score')
        section_div.append('<br/>')

        # Show tables conditionally
        n_rmchn = len(self.xmlnode.xpath('//MolProbity/Ramachandran_outliers/Outlier'))
        n_omega = len(self.xmlnode.xpath('//MolProbity/Nonplanar_omegas/Outlier'))
        n_rotam = len(self.xmlnode.xpath('//MolProbity/Rotamer_outliers/Outlier'))
        n_cbeta = len(self.xmlnode.xpath('//MolProbity/CBeta_outliers/Outlier'))
        n_clash = len(self.xmlnode.xpath('//MolProbity/Clashes/Outlier'))
        n_flips = len(self.xmlnode.xpath('//MolProbity/Side_chain_flips/Outlier'))

        fold = section_div.addFold(label='Detailed breakdown of outliers', initiallyOpen=False)

        rmchn_div = fold.addDiv()
        if n_rmchn > 0:
            rmchn_div.append('<p><b>Ramachandran outliers: </b></p>')
            table_outliers = rmchn_div.addTable(title='List of outliers', transpose=False)
            table_outliers.addData(title='Chain', data=self.xmlnode.xpath('//MolProbity/Ramachandran_outliers/Outlier/@chain'))
            table_outliers.addData(title='Seqnum', data=self.xmlnode.xpath('//MolProbity/Ramachandran_outliers/Outlier/@seqnum'))
            table_outliers.addData(title='Name', data=self.xmlnode.xpath('//MolProbity/Ramachandran_outliers/Outlier/@name'))
            table_outliers.addData(title='Score', data=self.xmlnode.xpath('//MolProbity/Ramachandran_outliers/Outlier/@score'))

        omega_div = fold.addDiv()
        if n_omega > 0:
            omega_div.append('<p><b>Non-planar omega (peptide bond) angles: </b></p>')
            table_outliers = omega_div.addTable(title='List of outliers', transpose=False)
            table_outliers.addData(title='Chain', data=self.xmlnode.xpath('//MolProbity/Nonplanar_omegas/Outlier/@chain'))
            table_outliers.addData(title='Seqnum', data=self.xmlnode.xpath('//MolProbity/Nonplanar_omegas/Outlier/@seqnum'))
            table_outliers.addData(title='Name', data=self.xmlnode.xpath('//MolProbity/Nonplanar_omegas/Outlier/@name'))

        rotam_div = fold.addDiv()
        if n_rotam > 0:
            rotam_div.append('<p><b>Rotamer outliers: </b></p>')
            table_outliers = rotam_div.addTable(title='List of outliers', transpose=False)
            table_outliers.addData(title='Chain', data=self.xmlnode.xpath('//MolProbity/Rotamer_outliers/Outlier/@chain'))
            table_outliers.addData(title='Seqnum', data=self.xmlnode.xpath('//MolProbity/Rotamer_outliers/Outlier/@seqnum'))
            table_outliers.addData(title='Name', data=self.xmlnode.xpath('//MolProbity/Rotamer_outliers/Outlier/@name'))
            table_outliers.addData(title='Score', data=self.xmlnode.xpath('//MolProbity/Rotamer_outliers/Outlier/@score'))

        cbeta_div = fold.addDiv()
        if n_cbeta > 0:
            cbeta_div.append('<p><b>C<sup><i>&#946;</i></sup> outliers: </b></p>')
            table_outliers = cbeta_div.addTable(title='List of outliers', transpose=False)
            table_outliers.addData(title='Chain', data=self.xmlnode.xpath('//MolProbity/CBeta_outliers/Outlier/@chain'))
            table_outliers.addData(title='Seqnum', data=self.xmlnode.xpath('//MolProbity/CBeta_outliers/Outlier/@seqnum'))
            table_outliers.addData(title='Name', data=self.xmlnode.xpath('//MolProbity/CBeta_outliers/Outlier/@name'))
            table_outliers.addData(title='Deviation', data=self.xmlnode.xpath('//MolProbity/CBeta_outliers/Outlier/@score'))

        flips_div = fold.addDiv()
        if n_flips > 0:
            flips_div.append('<p><b>Suggested side-chain flips: </b></p>')
            table_outliers = flips_div.addTable(title='List of flips', transpose=False)
            table_outliers.addData(title='Chain', data=self.xmlnode.xpath('//MolProbity/Side_chain_flips/Outlier/@chain'))
            table_outliers.addData(title='Seqnum', data=self.xmlnode.xpath('//MolProbity/Side_chain_flips/Outlier/@seqnum'))
            table_outliers.addData(title='Name', data=self.xmlnode.xpath('//MolProbity/Side_chain_flips/Outlier/@name'))

        clash_div = fold.addDiv()
        if n_clash > 0:
            clash_div.append('<p><b>Atomic clashes: </b></p>')
            table_outliers = clash_div.addTable(title='List of clashes', transpose=False)
            table_outliers.addData(title='First atom', data=self.xmlnode.xpath('//MolProbity/Clashes/Outlier/@first_atom'))
            table_outliers.addData(title='Second atom', data=self.xmlnode.xpath('//MolProbity/Clashes/Outlier/@second_atom'))
            table_outliers.addData(title='Overlap', data=self.xmlnode.xpath('//MolProbity/Clashes/Outlier/@overlap'))

        section_div.addDiv(style='clear:both;')
        self.append('<br/><br/>')


    def add_b_factors(self):
        b_factor_averages = { }
        for list_name in ('all', 'amino_acids', 'main_chains', 'side_chains', 'non_amino_acids', 'waters', 'ligands', 'ions'):
            b_factor_averages[list_name] = { }
            chain_ids = self.xmlnode.xpath('//B_factors/' + list_name + '/@chain')
            means = self.xmlnode.xpath('//B_factors/' + list_name + '/@mean')
            stds = self.xmlnode.xpath('//B_factors/' + list_name + '/@std')
            ns = self.xmlnode.xpath('//B_factors/' + list_name + '/@n')
            for chain_id, mean, std, n in zip(chain_ids, means, stds, ns):
                chain_id = int(chain_id) if chain_id != 'All' else 'All'
                mean = round(float(mean), 2) if mean not in ('None', 'nan') else None
                std = round(float(std), 2) if std not in ('None', 'nan') else None
                n = int(n) if n not in ('None', 'nan') else None
                b_factor_averages[list_name][chain_id] = (mean, std, n)

        section_div = self.addDiv()
        section_div.append('<style> th { white-space: nowrap; } </style>')
        section_div.append('<p style="font-weight:bold; font-size: 14px;"><b>B-factor Analyses</b></p>')
        #table_div = section_div.addDiv(style='float:left; margin-left: 30px; width:25%;')
        fold = section_div.addFold(label='Whole model', initiallyOpen=True)
        table = fold.addTable(title='B-factor analysis: whole model', transpose=False)
        table.addData(title='', data=('Mean', 'Stdev.', 'N'))
        table.addData(title='All monomers', data=b_factor_averages['all']['All'])
        table.addData(title='All amino acids', data=b_factor_averages['amino_acids']['All'])
        table.addData(title='Main chains', data=b_factor_averages['main_chains']['All'])
        table.addData(title='Side chains', data=b_factor_averages['side_chains']['All'])
        table.addData(title='All non-amino acids', data=b_factor_averages['non_amino_acids']['All'])
        table.addData(title='Waters', data=b_factor_averages['waters']['All'])
        table.addData(title='Ligands', data=b_factor_averages['ligands']['All'])
        table.addData(title='Ions', data=b_factor_averages['ions']['All'])
        section_div.append('<br/>')

        chain_count = int(self.xmlnode.xpath('//Model_info/Chain_count')[0].text)
        for chain_id in range(chain_count):
            fold = section_div.addFold(label='Chain ' + str(chain_id+1), initiallyOpen=False)
            table = fold.addTable(title='B-factor analysis: chain ' + str(chain_id+1), transpose=False)
            table.addData(title='', data=('Mean', 'Stdev.', 'N'))
            table.addData(title='All monomers', data=b_factor_averages['all'][chain_id])
            table.addData(title='All amino acids', data=b_factor_averages['amino_acids'][chain_id])
            table.addData(title='Main chains', data=b_factor_averages['main_chains'][chain_id])
            table.addData(title='Side chains', data=b_factor_averages['side_chains'][chain_id])
            table.addData(title='All non-amino acids', data=b_factor_averages['non_amino_acids'][chain_id])
            table.addData(title='Waters', data=b_factor_averages['waters'][chain_id])
            table.addData(title='Ligands', data=b_factor_averages['ligands'][chain_id])
            table.addData(title='Ions', data=b_factor_averages['ions'][chain_id])
            section_div.append('<br/>')

        section_div.addDiv(style='clear:both;')
        self.append('<br/><br/>')


    def add_ramachandran(self):
        section_div = self.addDiv()
        section_div.append('<p style="font-weight:bold; font-size:14px;"><b>Ramachandran Analyses</b></p>')

        background_PRO = 'img/rama_pro.png'
        background_GLY = 'img/rama_gly.png'
        background_RST = 'img/rama_rst.png'

        graph_div = self.addDiv(style='float:left; margin-left:50px;')
        rama_graph = graph_div.addFlotGraph(title='Ramachandran Plot', select='//Ramachandran', style='height:500px; width:500px; border:0px; float:left; padding:10px; padding-left:15px;')

        rama_graph.addData(title='PRO_favoured_phi', select='Favoured/Residue[@type=\'PRO\']/Phi')
        rama_graph.addData(title='PRO_favoured_psi', select='Favoured/Residue[@type=\'PRO\']/Psi')
        rama_graph.addData(title='PRO_Allowed_phi' , select='Allowed/Residue[@type=\'PRO\']/Phi')
        rama_graph.addData(title='PRO_Allowed_psi' , select='Allowed/Residue[@type=\'PRO\']/Psi')
        rama_graph.addData(title='PRO_Outliers_phi', select='Outliers/Residue[@type=\'PRO\']/Phi')
        rama_graph.addData(title='PRO_Outliers_psi', select='Outliers/Residue[@type=\'PRO\']/Psi')
        rama_graph.addData(title='GLY_favoured_phi', select='Favoured/Residue[@type=\'GLY\']/Phi')
        rama_graph.addData(title='GLY_favoured_psi', select='Favoured/Residue[@type=\'GLY\']/Psi')
        rama_graph.addData(title='GLY_Allowed_phi' , select='Allowed/Residue[@type=\'GLY\']/Phi')
        rama_graph.addData(title='GLY_Allowed_psi' , select='Allowed/Residue[@type=\'GLY\']/Psi')
        rama_graph.addData(title='GLY_Outliers_phi', select='Outliers/Residue[@type=\'GLY\']/Phi')
        rama_graph.addData(title='GLY_Outliers_psi', select='Outliers/Residue[@type=\'GLY\']/Psi')
        rama_graph.addData(title='RST_favoured_phi', select='Favoured/Residue[@type!=\'PRO\'][@type!=\'GLY\']/Phi')
        rama_graph.addData(title='RST_favoured_psi', select='Favoured/Residue[@type!=\'PRO\'][@type!=\'GLY\']/Psi')
        rama_graph.addData(title='RST_Allowed_phi' , select='Allowed/Residue[@type!=\'PRO\'][@type!=\'GLY\']/Phi')
        rama_graph.addData(title='RST_Allowed_psi' , select='Allowed/Residue[@type!=\'PRO\'][@type!=\'GLY\']/Psi')
        rama_graph.addData(title='RST_Outliers_phi', select='Outliers/Residue[@type!=\'PRO\'][@type!=\'GLY\']/Phi')
        rama_graph.addData(title='RST_Outliers_psi', select='Outliers/Residue[@type!=\'PRO\'][@type!=\'GLY\']/Psi')

        # Add rest graph
        p = rama_graph.addPlotObject()
        p.append('description', 'This graph shows the dihedral Phi and Psi angles for all residues, coloured according to Ramachandran\'s criterion. Source: Richardsons\' Top 500 structures.')
        p.append ('background').text = background_RST
        p.append('title', 'Ramachandran plot [Non-Pro/Gly]')
        p.append('plottype', 'xy')
        p.append('showlegend', 'false')
        p.append('xintegral', 'true')
        p.append('yintegral', 'true')
        p.append('xlabel', 'Phi')
        p.append('ylabel', 'Psi')
        p.append('xrange', min=-180.0, max=180.0)
        p.append('yrange', min=-180.0, max=180.0)
        l = p.append('plotline', xcol=13, ycol=14)
        l.append('colour', 'green')
        l.append('linestyle', '.')
        l.append('symbolsize', '1')
        l.append('symbol', '.')
        l = p.append('plotline', xcol=15, ycol=16)
        l.append('colour', 'orange')
        l.append('linestyle', '.')
        l.append('symbolsize', '3')
        l.append('symbol', 'o')
        l = p.append('plotline', xcol=17, ycol=18)
        l.append('colour', 'red')
        l.append('linestyle', '.')
        l.append('symbolsize', '10')
        l.append('symbol', 'x')

        # Add Proline graph
        p = rama_graph.addPlotObject()
        p.append('description', 'This graph shows the dihedral Phi and Psi angles for all residues, coloured according to Ramachandran\'s criterion. Source: Richardsons\' Top 500 structures.')
        p.append('background').text = background_PRO
        p.append('title', 'Ramachandran plot [Proline]')
        p.append('plottype', 'xy')
        p.append('showlegend', 'false')
        p.append('xintegral', 'true')
        p.append('yintegral', 'true')
        p.append('xlabel', 'Phi')
        p.append('ylabel', 'Psi')
        p.append('xrange', min=-180.0, max=180.0)
        p.append('yrange', min=-180.0, max=180.0)
        l = p.append('plotline', xcol=1, ycol=2)
        l.append('colour', 'green')
        l.append('linestyle', '.')
        l.append('symbolsize', '1')
        l.append('symbol', '.')
        l = p.append('plotline', xcol=3, ycol=4)
        l.append('colour', 'orange')
        l.append('linestyle', '.')
        l.append('symbolsize', '3')
        l.append('symbol', 'o')
        l = p.append('plotline', xcol=5, ycol=6)
        l.append('colour', 'red')
        l.append('linestyle', '.')
        l.append('symbolsize', '10')
        l.append('symbol', 'x')

        # Add Glycine graph
        p = rama_graph.addPlotObject()
        p.append('description', 'This graph shows the dihedral Phi and Psi angles for all residues, coloured according to Ramachandran\'s criterion. Source: Richardsons\' Top 500 structures.')
        p.append('background').text = background_GLY
        p.append('title', 'Ramachandran plot [Glycine]')
        p.append('plottype', 'xy')
        p.append('showlegend', 'false')
        p.append('xintegral', 'true')
        p.append('yintegral', 'true')
        p.append('xlabel', 'Phi')
        p.append('ylabel', 'Psi')
        p.append('xrange', min=-180.0, max=180.0)
        p.append('yrange', min=-180.0, max=180.0)
        l = p.append('plotline', xcol=7, ycol=8)
        l.append('colour', 'green')
        l.append('linestyle', '.')
        l.append('symbolsize', '1')
        l.append('symbol', '.')
        l = p.append('plotline', xcol=9, ycol=10)
        l.append('colour', 'orange')
        l.append('linestyle', '.')
        l.append('symbolsize', '3')
        l.append('symbol', 'o')
        l = p.append('plotline', xcol=11, ycol=12)
        l.append('colour', 'red')
        l.append('linestyle', '.')
        l.append('symbolsize', '10')
        l.append('symbol', 'x')

        n_residues = int(self.xmlnode.xpath('//Ramachandran/Totals/Residues')[0].text)
        n_favoured = int(self.xmlnode.xpath('//Ramachandran/Totals/Favoured')[0].text)
        n_allowed  = int(self.xmlnode.xpath('//Ramachandran/Totals/Allowed')[0].text)
        n_outliers = int(self.xmlnode.xpath('//Ramachandran/Totals/Outliers')[0].text)
        n_na = int(self.xmlnode.xpath('//Ramachandran/Totals/NA')[0].text)

        tab_div = self.addDiv(style='float:left; margin-left:50px;')

        if int(n_outliers) > 0 :
            tab_div.append('<p ><b>Outliers: </b></p>')
            table_outliers = tab_div.addTable(title='List of outliers', transpose=False, downloadable=True)
            table_outliers.addData(title='Chain', data=self.xmlnode.xpath('//Ramachandran/Outliers/Residue/@chain'))
            table_outliers.addData(title='Name',  data=self.xmlnode.xpath('//Ramachandran/Outliers/Residue/@type'))
            table_outliers.addData(title='Residue',    data=self.xmlnode.xpath('//Ramachandran/Outliers/Residue/@seqnum'))

        text ='<b>%s</b> residues have been analysed.<br/><br/>' % n_residues
        text +='In <b><font color=\'green\'>favoured</font></b> regions: <b>%s (%0.2f%%)</b><br/>' % (n_favoured, 100 * n_favoured / n_residues)
        text +='In <b><font color=\'orange\'>allowed</font></b> regions: <b>%s (%0.2f%%)</b><br/>' % (n_allowed, 100 * n_allowed / n_residues)
        text +='In <b><font color=\'red\'>high-energy</font></b> backbone conformations (outliers): <b>%s (%0.2f%%)</b>' % (n_outliers, 100 * n_outliers / n_residues)
        tab_div.append(text)

        section_div.addDiv(style='clear:both;')
        self.append('<br/><br/><br/><br/>')


if __name__ == '__main__':
    validate_protein_report(xmlFile=sys.argv[1], jobId=sys.argv[2])
