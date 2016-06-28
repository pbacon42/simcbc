from glue.ligolw import ligolw
from glue.ligolw import ilwd
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils

def create_empty_row(obj):
    """Create an empty sngl_inspiral row where the
    columns have default values of 0.0 for a float, 0 for an int, '' for
    a string. The ilwd columns have a default where the index is 0.
    """

    # check if sim_inspiral or sngl_inspiral
    row = lsctables.SimInspiral()
    cols = lsctables.SimInspiralTable.validcolumns

    # populate columns with default values
    for entry in cols.keys():
        if cols[entry] in ['real_4','real_8']:
            setattr(row,entry,0.)
        elif cols[entry] == 'int_4s':
            setattr(row,entry,0)
        elif cols[entry] == 'lstring':
            setattr(row,entry,'')
        elif entry == 'process_id':
            row.process_id = ilwd.ilwdchar("sim_inspiral:process_id:0")
        elif entry == 'simulation_id':
            row.simulation_id = ilwd.ilwdchar("sim_inspiral:simulation_id:0")
        else:
            raise ValueError("Column %s not recognized." %(entry) )

    return row

if __name__ == "__main__":

    class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
        pass
    
    lsctables.use_in(LIGOLWContentHandler)

    # Create new SimInspiral table
    table = lsctables.New(lsctables.SimInspiralTable)

    # Insert dummy event in table
    sim = create_empty_row(lsctables.SimInspiral)

    sim.h_end_time = 0.0
    sim.h_end_time_ns = 0.0
    sim.mass1 = 1.0
    sim.mass2 = 1.0
    sim.mchirp = 1.0
    sim.eta = 1.0
    sim.eff_dist_h = 1.0
    sim.coa_phase = 1.0
    sim.spin1z = 0.0
    sim.spin2z = 0.0

    table.append(sim)

    # Create, fill and write XML file
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW()).appendChild(table)
    #xmldoc.write("test.xml")
    
    #xmldoc.childNodes[-1].appendChild(sngl_inspiral_table_curr)
    ligolw_utils.write_filename(xmldoc, "demo.xml", verbose = True)
